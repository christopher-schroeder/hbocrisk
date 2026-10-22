"""Model specifications and derivation of genotype-specific incidences.

For each cancer, the incidence for an individual of major-gene risk class mu,
polygenotype x and risk-factor category r at age t is

    lambda(t) = h_mu(t) * exp(beta_PG(t) * x + g_r(t))

where beta_PG(t) is the residual polygenic standard deviation and g_r(t) the log-RR of
the risk-factor category. The class baselines h_mu(t) are determined age by age such
that (i) the average incidence of PV carriers in gene mu, over the polygene and risk
factors among those still unaffected, equals RR_mu(t) times the population incidence,
and (ii) the average incidence over all classes equals the population incidence
[ANT2001, LEE2016, LEE2019].
"""
from dataclasses import dataclass, field
from functools import lru_cache
import math

import numpy as np

from . import params as P
from .genetics import hpm, major_locus
from .incidence import load_population

FEMALE, MALE = "F", "M"


def pv_log_rr_variance(q, rr):
    """Variance of the log-RR due to PVs with allele frequency q and RR(t) [LEE2022]."""
    c = q * (2 - q)
    nc = (1 - q) ** 2
    return np.log((nc + c * rr ** 2) / (nc + c * rr) ** 2)


@dataclass(frozen=True, eq=False)
class ModelSpec:
    """Model specification. Hashing is by identity (specs are cache keys), so a modified
    copy (dataclasses.replace) never reuses incidences cached for the original."""
    name: str
    genes: tuple
    allele_frequencies: dict = field(hash=False, compare=False)
    rr: dict = field(hash=False, compare=False)          # (sex, cancer) -> {gene: RR(t)}
    series: dict = field(hash=False, compare=False)      # (sex, cancer) -> incidence series
    primary: tuple = (FEMALE, "bc")                     # cancer with RFs / PRS
    risk_factors: dict = field(hash=False, compare=False, default=None)
    polygenic_variance: dict = field(hash=False, compare=False, default=None)

    @property
    def locus(self):
        return major_locus(self.genes)

    def cancers(self, sex):
        return [c for (s, c) in self.series if s == sex]

    def class_rr(self, sex, cancer):
        rrs = self.rr[(sex, cancer)]
        out = np.ones((len(self.genes) + 1, P.MAX_AGE))
        for i, g in enumerate(self.genes, start=1):
            out[i] = rrs.get(g, P.ONES)
        return out


_SERIES = {(FEMALE, "bc"): "breast_female", (MALE, "bc"): "breast_male",
           (FEMALE, "oc"): "ovary", (MALE, "prc"): "prostate",
           (FEMALE, "pac"): "pancreas_female", (MALE, "pac"): "pancreas_male"}
_RR = {(FEMALE, "bc"): P.BC_RR, (MALE, "bc"): P.MBC_RR, (FEMALE, "oc"): P.OC_RR,
       (MALE, "prc"): P.PRC_RR, (FEMALE, "pac"): P.PAC_RR, (MALE, "pac"): P.PAC_RR}


def _bc_polygenic_variance():
    q = P.BC_ALLELE_FREQUENCIES["UK"]
    female = P.female_bc_polygenic_variance()
    # keep the total variance of the 2008 model: replace the BRCA2 contribution
    female = female + pv_log_rr_variance(q["BRCA2"], P.BC_RR_BRCA2_V6) \
        - pv_log_rr_variance(q["BRCA2"], P.BC_RR["BRCA2"])
    for g in ("PALB2", "CHEK2", "ATM", "BARD1", "RAD51C", "RAD51D"):
        female = female - pv_log_rr_variance(q[g], P.BC_RR[g])
    male = P.MALE_BC_POLYGENIC_VARIANCE - pv_log_rr_variance(q["PALB2"], P.MBC_RR["PALB2"])
    return {(FEMALE, "bc"): np.clip(female, 1e-6, None),
            (MALE, "bc"): np.clip(male, 1e-6, None)}


def _oc_polygenic_variance():
    q = P.OC_ALLELE_FREQUENCIES["UK"]
    v = np.full(P.MAX_AGE, P.EOC_POLYGENIC_VARIANCE)
    for g in ("RAD51D", "RAD51C", "BRIP1", "PALB2"):
        v = v - pv_log_rr_variance(q[g], P.OC_RR[g])
    return {(FEMALE, "oc"): np.clip(v, 1e-6, None)}


BOADICEA = ModelSpec(
    name="BC", genes=tuple(P.BC_GENES), allele_frequencies=P.BC_ALLELE_FREQUENCIES,
    rr=_RR, series=_SERIES, primary=(FEMALE, "bc"),
    risk_factors=None, polygenic_variance=_bc_polygenic_variance())

OVARIAN = ModelSpec(
    name="OC", genes=tuple(P.OC_GENES), allele_frequencies=P.OC_ALLELE_FREQUENCIES,
    rr={k: v for k, v in _RR.items()}, series=_SERIES, primary=(FEMALE, "oc"),
    risk_factors=None, polygenic_variance=_oc_polygenic_variance())


def model_risk_factors(spec):
    from .riskfactors import BC_RISK_FACTORS, OC_RISK_FACTORS
    return BC_RISK_FACTORS if spec.name == "BC" else OC_RISK_FACTORS


# --------------------------------------------------------------------------------------
# Calibration
# --------------------------------------------------------------------------------------


def calibrate(lam_pop, rr_class, class_prior, beta_pg, rf_probs=None, rf_logrr=None):
    """Derive class baselines so that incidences agree with the population.

    lam_pop: (T,) population incidence; rr_class: (C, T) RR relative to population for
    each class (row 0 ignored); class_prior: (C,) class frequencies at birth;
    beta_pg: (T,) residual polygenic SD (0 for no polygene); rf_probs (R,),
    rf_logrr (R, T) risk-factor categories.

    Returns h (C, T), density (C, T) = expected cases per class (unaffected mass times
    class average incidence) used for the tumour subtype constraints.
    """
    H = hpm()
    x = H.values
    K = len(x)
    T = len(lam_pop)
    C = len(class_prior)
    if rf_probs is None:
        rf_probs = np.ones(1)
        rf_logrr = np.zeros((1, T))
    R = len(rf_probs)
    w0 = class_prior[:, None, None] * H.prior[None, :, None] * rf_probs[None, None, :]
    logS = np.zeros((C, K, R))
    h = np.zeros((C, T))
    dens = np.zeros((C, T))
    for t in range(T):
        W = w0 * np.exp(logS)
        M = W.sum(axis=(1, 2))                                   # class masses
        rel = np.exp(beta_pg[t] * x[:, None] + rf_logrr[:, t][None, :])   # (K, R)
        E = (W * rel[None]).sum(axis=(1, 2))
        target = rr_class[:, t] * lam_pop[t]
        Mtot = M.sum()
        carriers = (M[1:] * target[1:]).sum()
        target[0] = max(lam_pop[t] * Mtot - carriers, 1e-3 * lam_pop[t] * M[0]) / M[0]
        with np.errstate(divide="ignore", invalid="ignore"):
            h[:, t] = np.where(E > 0, target * M / E, 0.0)
        dens[:, t] = M * target
        logS -= h[:, t][:, None, None] * rel[None]
    return h, dens


def _pathology_proportions(spec, dens):
    """Class-specific ER-, TN|ER- and CK proportions; non-carrier values constrained so
    that population averages match the general population distributions [MAV2010]."""
    C = len(spec.genes) + 1
    T = P.MAX_AGE
    p = np.zeros((C, T))
    q = np.zeros((C, T))
    ck = np.zeros((C, T, 3))
    for i, g in enumerate(spec.genes, start=1):
        p[i] = P.ER_NEG_PROPORTION.get(g, P.ER_NEG_PROPORTION["population"])
        q[i] = P.TN_GIVEN_ER_NEG_PROPORTION.get(g, P.TN_GIVEN_ER_NEG_PROPORTION["population"])
        ck[i] = P.CK_PROPORTIONS.get(g, P.CK_PROPORTIONS["population"])
    D = np.maximum(dens, 1e-300)
    pp = P.ER_NEG_PROPORTION["population"]
    qp = P.TN_GIVEN_ER_NEG_PROPORTION["population"]
    ckp = np.array(P.CK_PROPORTIONS["population"])
    p[0] = np.clip((pp * D.sum(0) - (D[1:] * p[1:]).sum(0)) / D[0], 1e-3, 0.999)
    Dp = D * p
    q[0] = np.clip((qp * Dp.sum(0) - (Dp[1:] * q[1:]).sum(0)) / Dp[0], 1e-3, 0.999)
    Dpq = Dp * q
    for j in range(3):
        ck[0, :, j] = (ckp[j] * Dpq.sum(0) - (Dpq[1:] * ck[1:, :, j]).sum(0)) / Dpq[0]
    ck[0] = np.clip(ck[0], 1e-3, None)
    ck[0] /= ck[0].sum(axis=1, keepdims=True)
    return p, q, ck


@dataclass
class CohortIncidences:
    """Genotype-specific incidences for one model/population/birth cohort."""
    inc: dict          # (sex, cancer) -> (C, K, T) incidences
    beta_pg: dict      # (sex, cancer) -> (T,)
    path: tuple        # (p_erneg (C,T), p_tn (C,T), ck (C,T,3))
    class_prior: np.ndarray


def _beta_pg(spec, key):
    v = (spec.polygenic_variance or {}).get(key)
    return np.zeros(P.MAX_AGE) if v is None else np.sqrt(v)


@lru_cache(maxsize=256)
def cohort_incidences(spec, population, cohort_idx, freq_key):
    """Incidences for all cancers of a model for a population birth cohort.

    freq_key: tuple of (gene, allele frequency) of the *population* frequencies used
    in the constraint.
    """
    table = load_population(population)
    locus = spec.locus
    class_prior = locus.class_prior(dict(freq_key))
    x = hpm().values
    inc = {}
    betas = {}
    dens_bc = None
    for key, series in spec.series.items():
        lam = table.rates(series, cohort_idx)
        beta = _beta_pg(spec, key)
        h, dens = calibrate(lam, spec.class_rr(*key), class_prior, beta)
        inc[key] = h[:, None, :] * np.exp(beta[None, None, :] * x[None, :, None])
        betas[key] = beta
        if key == (FEMALE, "bc"):
            dens_bc = dens
    path = _pathology_proportions(spec, dens_bc)
    return CohortIncidences(inc, betas, path, class_prior)


def proband_incidence(spec, population, cohort_idx, freq_key, rf_probs, rf_logrr, rf_index):
    """Incidences of the model's primary cancer for the proband's risk factor category."""
    table = load_population(population)
    key = spec.primary
    lam = table.rates(spec.series[key], cohort_idx)
    beta = _beta_pg(spec, key)
    class_prior = spec.locus.class_prior(dict(freq_key))
    h, _ = calibrate(lam, spec.class_rr(*key), class_prior, beta, rf_probs, rf_logrr)
    x = hpm().values
    return h[:, None, :] * np.exp(beta[None, None, :] * x[None, :, None]
                                  + rf_logrr[rf_index][None, None, :])
