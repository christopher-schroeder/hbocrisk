from pathlib import Path

import numpy as np
import pytest

from hbocrisk import (Pedigree, Person, GeneticTest, read_canrisk, write_canrisk,
                          calculate, Settings)
from hbocrisk.prs import reference_model, read_prs_file
from hbocrisk.riskfactors import BC_RISK_FACTORS, OC_RISK_FACTORS, observe, combine
from hbocrisk.server import handle_request

CANRISK = """##CanRisk 4.0
##menarche=13
##parity=2
##first_live_birth=28
##oc_use=F:5
##mht_use=N
##BMI=24.5
##alcohol=6
##menopause=51
##birads=c
##height=168
##PRS_BC=alpha=0.501,zscore=1.2
##FamID\tName\tTarget\tIndivID\tFathID\tMothID\tSex\tMZtwin\tDead\tAge\tYob\tBC1\tBC2\tOC\tPRO\tPAN\tAshkn\tBRCA1\tBRCA2\tPALB2\tATM\tCHEK2\tBARD1\tRAD51D\tRAD51C\tBRIP1\tHOXB13\tER:PR:HER2:CK14:CK56
F1\tgm\t0\tgm\t0\t0\tF\t0\t1\t70\t1930\t62\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\tN:N:N:0:0
F1\tgf\t0\tgf\t0\t0\tM\t0\t1\t75\t1928\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0
F1\tmum\t0\tm\tgf\tgm\tF\t0\t0\t58\t1962\t45\t0\t0\t0\t0\t0\tS:N\tS:N\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\tP:0:0:0:0
F1\tdad\t0\tf\t0\t0\tM\t0\t0\t60\t1960\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0
F1\tme\t1\tp\tf\tm\tF\t0\t0\t35\t1990\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0
"""


def test_read_write_roundtrip():
    ped = read_canrisk(CANRISK)[0]
    assert ped.proband.pid == "p"
    assert ped.prs["BC"] == (0.501, 1.2)
    assert ped.by_id["gm"].pathology_code() == 3
    assert ped.by_id["m"].pathology_code() == 0
    assert ped.by_id["m"].tests["BRCA1"].type == "S"
    again = read_canrisk(write_canrisk(ped))[0]
    assert again.to_dict() == ped.to_dict()


def test_dict_roundtrip_and_service():
    ped = read_canrisk(CANRISK)[0]
    d = ped.to_dict()
    out = handle_request({"pedigree": d, "settings": {"model": "BC"}})
    out2 = handle_request({"canrisk": CANRISK})
    assert out[0]["remaining_lifetime"] == out2[0]["remaining_lifetime"]
    assert 0 < out[0]["lifetime"]["breast"][0]["risk"] < 1


def test_risk_factors_parse_and_effect():
    ped = read_canrisk(CANRISK)[0]
    obs = observe(BC_RISK_FACTORS, ped.risk_factors, 35)
    assert {o.name for o in obs} == {"menarche", "parity", "first_live_birth", "oc_use",
                                     "mht_use", "bmi", "alcohol", "menopause", "birads",
                                     "height"}
    probs, logrr, idx = combine(obs)
    assert probs.sum() == pytest.approx(1)
    # dense breasts increase risk relative to no risk factor information
    base = Pedigree("X", [Person("p", "F", target=True, age=40, yob=1985)])
    dense = Pedigree("X", [Person("p", "F", target=True, age=40, yob=1985)],
                     risk_factors={"birads": "d"})
    fatty = Pedigree("X", [Person("p", "F", target=True, age=40, yob=1985)],
                     risk_factors={"birads": "a"})
    r0, r1, r2 = (calculate(x).remaining_lifetime["breast"][-1][1] for x in (base, dense, fatty))
    assert r2 < r0 < r1


def test_average_over_risk_factor_distribution_is_population_risk():
    """Averaging risks over the population distribution of a factor recovers the risk
    without that factor (model internal consistency)."""
    f = BC_RISK_FACTORS["parity"]
    risks = []
    for i, lab in enumerate(["0", "1", "2", "3"]):
        ped = Pedigree("X", [Person("p", "F", target=True, age=20, yob=1985)],
                       risk_factors={"parity": lab})
        risks.append(calculate(ped, calcs=("lifetime",)).lifetime["breast"][0][1])
    avg = float(np.dot(f.probs, risks))
    base = calculate(Pedigree("X", [Person("p", "F", target=True, age=20, yob=1985)]),
                     calcs=("lifetime",)).lifetime["breast"][0][1]
    assert avg == pytest.approx(base, rel=2e-3)


def test_average_over_prs_is_population_risk():
    from scipy.stats import norm
    zs = np.linspace(-6, 6, 241)
    w = norm.pdf(zs)
    w /= w.sum()
    risks = [calculate(Pedigree("X", [Person("p", "F", target=True, age=20, yob=1985)],
                                prs={"BC": (0.5, z)}), calcs=("lifetime",)).lifetime["breast"][0][1]
             for z in zs]
    base = calculate(Pedigree("X", [Person("p", "F", target=True, age=20, yob=1985)]),
                     calcs=("lifetime",)).lifetime["breast"][0][1]
    assert float(np.dot(w, risks)) == pytest.approx(base, rel=2e-3)


def test_prs_reference_models():
    m = reference_model("BCAC 313")
    assert len(m.variants) == 313
    assert m.alpha == 0.501
    # expected raw score under HWE equals the published mean
    p = np.array([v.eaf for v in m.variants])
    b = np.array([v.beta for v in m.variants])
    assert (2 * p * b).sum() == pytest.approx(m.mean, abs=2e-3)
    raw, n = m.raw_score({})
    assert n == 0 and m.zscore(raw) == pytest.approx(0, abs=0.01)
    oc = reference_model("OCAC 36")
    assert oc.disease == "OC" and oc.alpha == 0.223


def test_ovarian_model_risk_factors():
    base = Pedigree("X", [Person("p", "F", target=True, age=30, yob=1990)])
    low = Pedigree("X", [Person("p", "F", target=True, age=30, yob=1990)],
                   risk_factors={"parity": "3", "oc_use": "F:15", "tl": "Y"})
    s = Settings(model="OC")
    r0 = calculate(base, s).remaining_lifetime["ovarian"][-1][1]
    r1 = calculate(low, s).remaining_lifetime["ovarian"][-1][1]
    assert r1 < 0.6 * r0
