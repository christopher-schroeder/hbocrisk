"""High level interface: carrier probabilities and cancer risks for a pedigree.

Risk definitions follow the CanRisk web services ("Cancer Risk Calculations"):

* remaining lifetime risk: from the proband's current age to ages a+1..a+5, a+10 and
  ages divisible by 5 up to 80;
* lifetime risk: proband's age set to 20, risk to age 80 (unaffected probands only);
* 10-year risk: proband's age set to 40, risk to age 50 (unaffected probands only);
* baseline risks: the same calculations for a woman of the same age and year of birth
  (and age at first breast cancer, if affected) from the general population, ignoring
  family history, genetic tests, risk factors and PRS.

Following the reformulated algorithm of BOADICEA v6 (supplement, eq. s.2), the risk is
computed from the proband's posterior genotype distribution given all pedigree data,
combined with the proband's genotype-specific penetrance at future ages.
"""
from copy import deepcopy
from dataclasses import dataclass, field

import numpy as np

from . import params as P
from .genetics import hpm
from .incidence import load_population
from .likelihood import PedigreeLikelihood, genetic_test_factor, phenotype_loglik, _cum
from .models import (BOADICEA, OVARIAN, FEMALE, MALE, cohort_incidences, proband_incidence,
                     model_risk_factors)
from .pedigree import Pedigree, PedigreeError
from .riskfactors import observe, combine, unused_risk_factors

MODELS = {"BC": BOADICEA, "BOADICEA": BOADICEA, "OC": OVARIAN, "OVARIAN": OVARIAN}
CANCER_NAMES = {"bc": "breast", "oc": "ovarian", "prc": "prostate", "pac": "pancreatic"}


@dataclass
class Settings:
    model: str = "BC"
    population: str = "UK"                 # cancer incidence rates
    mutation_frequencies: object = "UK"    # name or dict gene -> allele frequency
    sensitivities: dict = None              # gene -> mutation search sensitivity
    max_rf_categories: int = 1500


@dataclass
class Result:
    model: str
    famid: str
    proband: str
    carrier_probabilities: dict = None
    remaining_lifetime: dict = None        # cancer -> list of (age, risk)
    lifetime: dict = None
    ten_year: dict = None
    baseline_remaining_lifetime: dict = None
    baseline_lifetime: dict = None
    baseline_ten_year: dict = None
    warnings: list = field(default_factory=list)
    settings: dict = None

    def to_dict(self):
        def conv(d):
            if d is None:
                return None
            return {k: [{"age": a, "risk": float(r)} for a, r in v] for k, v in d.items()}
        return {
            "model": self.model, "famid": self.famid, "proband": self.proband,
            "carrier_probabilities": (None if self.carrier_probabilities is None else
                                      {k: float(v) for k, v in self.carrier_probabilities.items()}),
            "remaining_lifetime": conv(self.remaining_lifetime),
            "lifetime": conv(self.lifetime), "ten_year": conv(self.ten_year),
            "baseline_remaining_lifetime": conv(self.baseline_remaining_lifetime),
            "baseline_lifetime": conv(self.baseline_lifetime),
            "baseline_ten_year": conv(self.baseline_ten_year),
            "warnings": list(self.warnings), "settings": self.settings,
        }


class _Context:
    def __init__(self, spec, settings, pedigree):
        self.spec = spec
        self.population = settings.population
        self.table = load_population(settings.population)
        mf = settings.mutation_frequencies
        if isinstance(mf, dict):
            pop = dict(spec.allele_frequencies["UK"])
            pop.update(mf)
            founder = pop
        else:
            name = mf
            if pedigree.is_ashkenazi() and name == "UK":
                name = "Ashkenazi"
            if name not in spec.allele_frequencies:
                raise ValueError(f"unknown mutation frequency set '{name}'")
            founder = spec.allele_frequencies[name]
            # Ashkenazi frequencies are used for founders only; the population
            # incidences are partitioned with the general population frequencies.
            pop = spec.allele_frequencies["UK" if name == "Ashkenazi" else name]
        self.founder_freqs = {g: founder[g] for g in spec.genes}
        self.pop_freq_key = tuple((g, float(pop[g])) for g in spec.genes)
        sens = dict(P.TEST_SENSITIVITY)
        sens.update(settings.sensitivities or {})
        self.sensitivities = sens
        self.max_rf = settings.max_rf_categories
        self.locus = spec.locus

    def cohort(self, yob):
        return self.table.cohort_index(yob)

    def incidences(self, yob):
        return cohort_incidences(self.spec, self.population, self.cohort(yob), self.pop_freq_key)


def _usable(person):
    return bool(person.age) and bool(person.yob)


def _person_incs(ctx, person, proband_primary=None):
    ci = ctx.incidences(person.yob)
    sex = person.sex
    incs = {}
    for (s, c), arr in ci.inc.items():
        if s != sex:
            continue
        if proband_primary is not None and (s, c) == ctx.spec.primary:
            arr = proband_primary
        incs[c] = arr
    return incs, ci


def _proband_primary(ctx, pedigree, proband, current_age):
    """Primary-cancer incidences for the proband given her risk factors (or None)."""
    spec = ctx.spec
    if proband.sex != spec.primary[0]:
        return None, []
    defs = model_risk_factors(spec)
    rfs = observe(defs, pedigree.risk_factors, current_age)
    if not rfs:
        return None, []
    probs, logrr, idx = combine(rfs, ctx.max_rf)
    inc = proband_incidence(spec, ctx.population, ctx.cohort(proband.yob), ctx.pop_freq_key,
                            probs, logrr, idx)
    return inc, rfs


def _prs_factor(ctx, pedigree, proband):
    key = "BC" if ctx.spec.name == "BC" else "OC"
    prs = pedigree.prs.get(key)
    if prs is None or proband.sex != ctx.spec.primary[0]:
        return None
    alpha, z = prs
    return hpm().prs_likelihood(z, alpha)


def proband_posterior(pedigree, ctx, warnings=None):
    """Posterior over (genotype, polygenotype) for the proband and the proband's
    incidences used for risk projection."""
    proband = pedigree.proband
    locus = ctx.locus
    geno_prior = locus.genotype_prior(ctx.founder_freqs)
    primary, rfs = _proband_primary(ctx, pedigree, proband, proband.age)
    pen = {}
    proband_incs = None
    for p in pedigree.people:
        is_proband = p.pid == proband.pid
        table = np.ones((locus.ngeno, hpm().size))
        has_data = False
        if _usable(p):
            incs, ci = _person_incs(ctx, p, primary if is_proband else None)
            ll = phenotype_loglik(p, incs, ci.path if p.sex == FEMALE else None)
            ll = ll - ll.max()
            table = np.exp(ll)[locus.geno_class, :]
            has_data = True
            if is_proband:
                proband_incs = incs

        tf = genetic_test_factor(p, locus, ctx.sensitivities)
        if not np.all(tf == 1):
            table = table * tf[:, None]
            has_data = True
        if is_proband:
            prs = _prs_factor(ctx, pedigree, proband)
            if prs is not None:
                table = table * prs[None, :]
            has_data = True
        if has_data:
            pen[p.pid] = table
    lik = PedigreeLikelihood(pedigree, locus, geno_prior, pen)
    post, loglik = lik.proband_posterior()
    return post, proband_incs, loglik, rfs


def _risk_curve(post, inc, age, ages, locus, factor=1.0):
    """Risk of developing cancer between `age` and each age in `ages`."""
    w = post.T @ locus.class_matrix          # (K, C)
    w = w.T                                    # (C, K)
    cum = _cum(factor * inc)
    out = []
    for a in ages:
        a2 = min(a, P.MAX_AGE)
        r = 1.0 - np.exp(-(cum[..., a2] - cum[..., age]))
        out.append((a, float((w * r).sum())))
    return out


def remaining_ages(age):
    ages = []
    a = age
    while a < P.MAX_AGE:
        a += 1
        if a <= age + 5 or a % 5 == 0 or a == age + 10:
            ages.append(a)
    return ages


def _risks_for(pedigree, ctx, ages, warnings=None, posterior=None):
    """Return dict cancer -> [(age, risk)] for the proband at her current age."""
    proband = pedigree.proband
    post, incs, _, _ = posterior or proband_posterior(pedigree, ctx, warnings)
    out = {}
    age = proband.age
    bc1 = proband.cancer_age("bc1")
    primary = ctx.spec.primary[1]
    for c, inc in incs.items():
        if c != primary:
            continue
        if c == "bc":
            if proband.bc2 is not None:
                continue
            factor = 0.5 if bc1 is not None else 1.0
            name = "contralateral breast" if bc1 is not None else "breast"
            out[name] = _risk_curve(post, inc, age, ages, ctx.locus, factor)
        elif c == "oc" and proband.oc is None:
            out["ovarian"] = _risk_curve(post, inc, age, ages, ctx.locus)
    return out


def _carrier_probs(post, ctx):
    cls = post.sum(axis=1) @ ctx.locus.class_matrix
    out = {"no_pathogenic_variant": cls[0]}
    for i, g in enumerate(ctx.spec.genes, start=1):
        out[g] = cls[i]
    return out


def _baseline_pedigree(pedigree, age=None, keep_bc1=True):
    p = deepcopy(pedigree.proband)
    p.father = p.mother = "0"
    p.mztwin = "0"
    p.tests = {}
    p.pathology = {}
    bc1 = p.bc1 if keep_bc1 else None
    p.bc1, p.bc2, p.oc, p.prc, p.pac = bc1, None, None, None, None
    if age is not None:
        p.age = age
    return Pedigree(pedigree.famid, [p])


def _with_proband_age(pedigree, age):
    ped = deepcopy(pedigree)
    ped.proband.age = age
    return ped


def calculate(pedigree, settings=None, calcs=("carrier_probs", "remaining_lifetime",
                                              "lifetime", "ten_year")):
    """Run the model for a pedigree. Returns a :class:`Result`."""
    settings = settings or Settings()
    spec = MODELS[settings.model.upper()]
    pedigree.validate()
    ctx = _Context(spec, settings, pedigree)
    proband = pedigree.proband
    res = Result(model=spec.name, famid=pedigree.famid, proband=proband.pid,
                 settings={"model": spec.name, "population": settings.population,
                           "mutation_frequencies": dict(ctx.founder_freqs),
                           "sensitivities": {g: ctx.sensitivities[g] for g in spec.genes}})
    warnings = res.warnings
    warnings.extend(unused_risk_factors(spec, pedigree.risk_factors))
    for p in pedigree.people:
        if p.age and p.age >= P.MAX_AGE and p.affected():
            ages = [c for c in ("bc1", "bc2", "oc", "prc", "pac")
                    if p.cancer_age(c) is not None and p.cancer_age(c) >= P.MAX_AGE]
            if ages:
                warnings.append(f"{p.pid}: cancer diagnoses at age 80 or over are not "
                                f"included (individual treated as unaffected at 80)")

    incomplete = [p.pid for p in pedigree.people if not _usable(p)]
    if incomplete:
        warnings.append("year of birth and age at last follow up must be specified in order "
                        f"for {', '.join(incomplete)} to be included in a calculation")
    positive = any(t.tested and t.result == "P" for g, t in proband.tests.items()
                   if g in spec.genes)
    current = proband_posterior(pedigree, ctx, warnings)
    if "carrier_probs" in calcs and not positive:
        res.carrier_probabilities = _carrier_probs(current[0], ctx)

    risks_viable = (proband.sex == spec.primary[0] and not proband.dead
                    and proband.age < P.MAX_AGE and proband.bc2 is None
                    and proband.oc is None and proband.pac is None and proband.prc is None)
    affected = proband.affected()
    if risks_viable and "remaining_lifetime" in calcs:
        ages = remaining_ages(proband.age)
        if ages:
            res.remaining_lifetime = _risks_for(pedigree, ctx, ages, posterior=current)
            base = _baseline_pedigree(pedigree)
            res.baseline_remaining_lifetime = _risks_for(base, _Context(spec, settings, base), ages)
    if risks_viable and not affected:
        for name, (a0, a1) in (("lifetime", (20, 80)), ("ten_year", (40, 50))):
            if name not in calcs:
                continue
            ped = _with_proband_age(pedigree, a0)
            setattr(res, name, _risks_for(ped, ctx, [a1]))
            base = _baseline_pedigree(pedigree, age=a0)
            setattr(res, "baseline_" + name,
                    _risks_for(base, _Context(spec, settings, base), [a1]))
    elif not risks_viable:
        warnings.append("cancer risks are not calculated for this target (male for this "
                        "model, deceased, aged 80 or over, or affected by a cancer other "
                        "than a single breast cancer)")
    return res
