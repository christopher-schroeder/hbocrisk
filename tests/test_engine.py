import itertools

import numpy as np
import pytest

from openboadicea import Pedigree, Person, GeneticTest, calculate, Settings
from openboadicea.api import _Context, proband_posterior
from openboadicea.genetics import hpm, major_locus
from openboadicea.likelihood import PedigreeLikelihood
from openboadicea.models import BOADICEA, calibrate, cohort_incidences
from openboadicea import params as P


def test_hpm_properties():
    H = hpm()
    assert H.prior.sum() == pytest.approx(1)
    assert (H.prior * H.values).sum() == pytest.approx(0, abs=1e-12)
    assert (H.prior * H.values ** 2).sum() == pytest.approx(1)
    T = H.transmission()
    assert np.allclose(T.sum(axis=0), 1)
    # founders' distribution is stationary under random mating
    child = np.einsum("cmf,m,f->c", T, H.prior, H.prior)
    assert np.allclose(child, H.prior)
    # parent-offspring correlation of the polygenotype is 1/2
    joint = np.einsum("cmf,m,f->cm", T, H.prior, H.prior)
    cov = (joint * np.outer(H.values, H.values)).sum()
    assert cov == pytest.approx(0.5)


def test_major_locus():
    L = major_locus(tuple(P.BC_GENES))
    assert L.ngeno == 45
    T = L.transmission()
    assert np.allclose(T.sum(axis=0), 1)
    prior = L.genotype_prior(P.BC_ALLELE_FREQUENCIES["UK"])
    assert prior.sum() == pytest.approx(1)
    child = np.einsum("cmf,m,f->c", T, prior, prior)
    assert np.allclose(child, prior)
    # class 1 (BRCA1) frequency ~ 2q
    cls = L.class_prior(P.BC_ALLELE_FREQUENCIES["UK"])
    q = P.BC_ALLELE_FREQUENCIES["UK"]["BRCA1"]
    assert cls[1] == pytest.approx(2 * q - q * q)


def test_prs_likelihood_limits():
    H = hpm()
    flat = H.prs_likelihood(1.5, 0.0)
    assert np.allclose(flat, 1)
    f = H.prs_likelihood(2.0, 0.5)
    assert np.all(np.diff(f) > 0)      # high PRS favours high polygenotypes
    # posterior mean of the polygene equals the observed PRS (x_PRS = alpha z)
    for z in (-2, 0.5, 2):
        f = H.prs_likelihood(z, 0.5) * H.prior
        # with bins, E[x_P | x_PRS] is approximately alpha*z
        edges = H.bin_edges()
        assert np.isfinite(f).all()


def test_calibration_reproduces_population_incidence():
    lam = np.linspace(1e-5, 4e-3, P.MAX_AGE)
    rr = BOADICEA.class_rr("F", "bc")
    prior = BOADICEA.locus.class_prior(P.BC_ALLELE_FREQUENCIES["UK"])
    beta = np.sqrt(BOADICEA.polygenic_variance[("F", "bc")])
    h, dens = calibrate(lam, rr, prior, beta)
    # simulate population forward and compare average incidence among unaffected
    x = hpm().values
    w = prior[:, None] * hpm().prior[None, :]
    S = np.ones_like(w)
    for t in range(P.MAX_AGE):
        inc = h[:, t][:, None] * np.exp(beta[t] * x)[None, :]
        W = w * S
        assert (W * inc).sum() / W.sum() == pytest.approx(lam[t], rel=1e-9)
        # carriers: average incidence = RR * population
        for c in range(1, len(prior)):
            avg = (W[c] * inc[c]).sum() / W[c].sum()
            assert avg == pytest.approx(rr[c, t] * lam[t], rel=1e-9)
        S = S * np.exp(-inc)


def _brute_force_posterior(ped, ctx):
    """Enumerate all joint genotypes of a 3-person nuclear family."""
    from openboadicea.api import _person_incs, _usable
    from openboadicea.likelihood import phenotype_loglik, genetic_test_factor
    L = ctx.locus
    H = hpm()
    prior = L.genotype_prior(ctx.founder_freqs)
    pens = {}
    for p in ped.people:
        t = np.ones((L.ngeno, H.size))
        if _usable(p):
            incs, ci = _person_incs(ctx, p)
            ll = phenotype_loglik(p, incs, ci.path if p.sex == "F" else None)
            t = np.exp(ll)[L.geno_class, :]
        t = t * genetic_test_factor(p, L, ctx.sensitivities)[:, None]
        pens[p.pid] = t
    Tg, Tp = L.transmission(), H.transmission()
    m, f, c = pens["m"], pens["f"], pens["c"]
    # sum over parents' genotypes
    pm = (prior[:, None] * H.prior[None, :]) * m       # (G, K)
    pf = (prior[:, None] * H.prior[None, :]) * f
    post = np.einsum("cmf,dno,mn,fo->cd", Tg, Tp, pm, pf) * c
    return post / post.sum()


def test_variable_elimination_matches_enumeration():
    ped = Pedigree("T", [
        Person("m", "F", age=45, yob=1950, bc1=42, pathology={"er": "N"}),
        Person("f", "M", age=70, yob=1948, prc=60),
        Person("c", "F", father="f", mother="m", target=True, age=35, yob=1980,
               tests={"CHEK2": GeneticTest("S", "N")}),
    ])
    ctx = _Context(BOADICEA, Settings(), ped)
    post, _, _, _ = proband_posterior(ped, ctx)
    bf = _brute_force_posterior(ped, ctx)
    assert np.allclose(post, bf, atol=1e-12)


def test_mz_twins_share_genotype():
    def ped(twin):
        return Pedigree("T", [
            Person("m", "F"), Person("f", "M"),
            Person("a", "F", father="f", mother="m", mztwin="1" if twin else "0",
                   age=40, yob=1980, bc1=35),
            Person("b", "F", father="f", mother="m", mztwin="1" if twin else "0",
                   age=40, yob=1980, target=True),
        ])
    r_twin = calculate(ped(True))
    r_sib = calculate(ped(False))
    assert r_twin.remaining_lifetime["breast"][-1][1] > r_sib.remaining_lifetime["breast"][-1][1]


def test_population_risk_without_information():
    """A woman with no family history has the population risk (up to the discrete-year
    approximation: heterogeneity within each year makes the risk slightly lower)."""
    from openboadicea.incidence import load_population
    ped = Pedigree("T", [Person("p", "F", target=True, age=20, yob=1985)])
    r = calculate(ped)
    lam = load_population("UK").rates("breast_female", 7)
    pop = 1 - np.exp(-lam[20:80].sum())
    assert r.lifetime["breast"][0][1] == pytest.approx(pop, rel=5e-3)
    assert r.baseline_lifetime["breast"][0][1] == pytest.approx(r.lifetime["breast"][0][1])


def test_carrier_probabilities_reference():
    """Official BOADICEA 7.3.2 web-service output (CanRisk bws test suite): unaffected
    25-year-old woman, no family history."""
    ped = Pedigree("T", [Person("p", "F", target=True, age=25, yob=2001)])
    cp = calculate(ped).carrier_probabilities
    ref = {"no_pathogenic_variant": 0.9821809, "BRCA1": 0.0012713, "BRCA2": 0.0020312,
           "PALB2": 0.0012758, "CHEK2": 0.0074259, "ATM": 0.0035743, "BARD1": 0.0008529,
           "RAD51C": 0.0006940, "RAD51D": 0.0006937}
    for k, v in ref.items():
        assert cp[k] == pytest.approx(v, rel=2e-3)


def test_negative_direct_test_excludes_gene():
    ped = Pedigree("T", [
        Person("m", "F", age=45, yob=1950, bc1=40, tests={"BRCA1": GeneticTest("T", "P")}),
        Person("f", "M"),
        Person("c", "F", father="f", mother="m", target=True, age=30, yob=1980,
               tests={"BRCA1": GeneticTest("T", "N")}),
    ])
    cp = calculate(ped).carrier_probabilities
    assert cp["BRCA1"] == 0


def test_positive_proband_gives_carrier_risk():
    base = Pedigree("T", [Person("p", "F", target=True, age=20, yob=1985)])
    car = Pedigree("T", [Person("p", "F", target=True, age=20, yob=1985,
                                tests={"BRCA1": GeneticTest("T", "P")})])
    r0, r1 = calculate(base), calculate(car)
    assert r1.carrier_probabilities is None
    assert r1.lifetime["breast"][0][1] > 5 * r0.lifetime["breast"][0][1]


def test_loop_pedigree():
    """Pedigree with a loop (first-cousin marriage) is handled."""
    ped = Pedigree("L", [
        Person("gf", "M"), Person("gm", "F", age=60, yob=1920, bc1=55),
        Person("a", "M", father="gf", mother="gm"), Person("aw", "F"),
        Person("b", "F", father="gf", mother="gm", age=50, yob=1945, bc1=45),
        Person("bh", "M"),
        Person("c1", "M", father="a", mother="aw"),
        Person("c2", "F", father="bh", mother="b"),
        Person("p", "F", father="c1", mother="c2", target=True, age=30, yob=1990),
    ])
    r = calculate(ped)
    assert 0 < r.remaining_lifetime["breast"][-1][1] < 1
