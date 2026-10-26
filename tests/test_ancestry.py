"""Ancestry-specific PRS normalisation and UK ethnic groups (BOADICEA v7) [FIC2025E]."""
import numpy as np
import pytest

from hbocrisk import Pedigree, Person, Settings, calculate, read_canrisk
from hbocrisk.ethnicity import ethnic_group
from hbocrisk.params import ETHNICITY_IRR
from hbocrisk.prs import ANCESTRIES, ANCESTRY_PRS, ancestry_file, reference_model


def _data_available():
    try:
        for name in ANCESTRY_PRS:
            ancestry_file(name)
        return True
    except RuntimeError:
        return False


# the ancestry-specific PRS files are downloaded on first use (GPL-3.0, not bundled)
needs_data = pytest.mark.skipif(not _data_available(),
                                reason="ancestry-specific PRS files not available")


@needs_data
def test_ancestry_prs_files():
    m309, m307, m313 = (reference_model(n) for n in ("BCAC 309", "BCAC 307", "BCAC 313"))
    assert len(m309.variants) == 309 and len(m307.variants) == 307
    # the 309/307-SNP scores are subsets of BCAC 313 with identical weights
    w313 = {(v.chrom, v.pos): v.beta for v in m313.variants}
    assert all(w313[(v.chrom, v.pos)] == v.beta for v in m309.variants)
    # Ficorella et al. 2025, Table 2
    table2 = {"European": (-0.335, 0.615, 0.499), "African": (0.292, 0.547, 0.196),
              "East Asian": (0.268, 0.555, 0.326), "South Asian": (-0.160, 0.602, 0.331)}
    for a, (mean, sd, alpha) in table2.items():
        q = m309.ancestry[a]
        assert (q.mean, q.sd, q.alpha) == (mean, sd, alpha)
    assert m309.ancestry["South Asian"].threshold == 0.7
    assert (m309.mean, m309.sd, m309.alpha) == (-0.339, 0.615, 0.499)
    # expected raw score with ancestry-specific frequencies ~ ancestry mean
    for a in ANCESTRIES:
        raw, n = m309.raw_score({}, ancestry=a)
        assert n == 0 and raw == pytest.approx(m309.ancestry[a].mean, abs=0.01)


# reference values from CCGE-BOADICEA/SHARE-PRScalculation vcf2prs (v4.2.12)
@needs_data
@pytest.mark.parametrize("name,raw,props,z,alpha,method", [
    ("BCAC 309", -0.3, (0.4, 0.25, 0.2, 0.15), -0.30988, 0.34210, "mixed"),
    ("BCAC 309", -0.5, (0.79, 0.0, 0.21, 0.0), -0.42257, 0.45176, "mixed"),
    ("BCAC 309", 0.2, (0.85, 0.1, 0.05, 0.0), 0.86992, 0.499, "single"),
    ("BCAC 309", 0.0, (0.2, 0.1, 0.0, 0.7), 0.26578, 0.331, "single"),
    ("BCAC 307", 0.1, (0.5, 0.5, 0.0, 0.0), -0.14589, 0.33070, "mixed"),
    ("OCAC 36", 0.1, (0.5, 0.5, 0.0, 0.0), 1.41592, 0.19267, "mixed"),
    ("OCAC 36", 0.4, (0.0, 1.0, 0.0, 0.0), 3.18504, 0.178, "single"),
])
def test_mixed_ancestry_matches_vcf2prs(name, raw, props, z, alpha, method):
    res = reference_model(name, ancestry=True).normalise(raw, dict(zip(("EUR", "AFR", "EAS", "SAS"), props)))
    assert res["method"] == method
    assert res["zscore"] == pytest.approx(z, abs=1e-5)
    assert res["alpha"] == pytest.approx(alpha, abs=1e-5)


@needs_data
def test_single_and_unknown_ancestry():
    m = reference_model("BCAC 309")
    assert m.normalise(0.292, "African") == {"alpha": 0.196, "zscore": pytest.approx(0),
                                             "ancestry": "African", "method": "single"}
    assert m.normalise(0.0, None)["method"] == "population"
    assert m.normalise(0.0, "unknown")["alpha"] == 0.499
    m313 = reference_model("BCAC 313")
    assert m313.normalise(-0.424, "European")["zscore"] == pytest.approx(0)
    with pytest.raises(ValueError, match="no parameters for African"):
        m313.normalise(0.0, "African")
    with pytest.raises(ValueError, match="add up to 1"):
        m.normalise(0.0, {"EUR": 0.5, "AFR": 0.3})
    with pytest.raises(ValueError, match="unknown ancestry"):
        m.normalise(0.0, "Martian")


def test_ethnic_group_names():
    assert ethnic_group("Asian or Asian British;Chinese") == "East Asian"
    assert ethnic_group("Asian or Asian British;Indian") == "South Asian"
    assert ethnic_group("Black or Black British;African") == "Black"
    assert ethnic_group("White;Irish") == "White"
    assert ethnic_group("Mixed/Multiple ethnic groups;White and Asian") == "Mixed"
    assert ethnic_group("chinese") == "East Asian"
    assert ethnic_group("Unknown") is None and ethnic_group(None) is None
    with pytest.raises(ValueError):
        ethnic_group("Martian")


def test_ethnicity_irr_restandardised():
    # 2016 census proportions (Table S1); 'Other' at the population rate
    props = {"White": 0.851, "Black": 0.035, "South Asian": 0.078, "Mixed": 0.0176,
             "Other": 0.018}
    for series in ("breast_female", "ovary", "prostate"):
        total = sum(p * ETHNICITY_IRR[g].get(series, 1.0) for g, p in props.items())
        assert total == pytest.approx(1.0, abs=2e-3)


def _lifetime(ethnicity, population="UK", pedigree=None):
    ped = pedigree or Pedigree("X", [Person("p", "F", target=True, age=30, yob=1990)])
    res = calculate(ped, Settings(population=population, ethnicity=ethnicity),
                    calcs=("lifetime",))
    return res, res.lifetime["breast"][0][1]


def test_ethnicity_incidences():
    _, pop = _lifetime(None)
    white, w = _lifetime("White")
    black, b = _lifetime("Black")
    other, o = _lifetime("Other")
    # 'Other': population incidences, but the non-European CHEK2 frequency
    assert b < pop < w and o == pytest.approx(pop, rel=1e-3)
    assert other.settings["mutation_frequency_set"] == "UK, non-European"
    assert white.settings["mutation_frequency_set"] == "UK"
    assert black.settings["mutation_frequency_set"] == "UK, non-European"
    assert black.settings["ethnicity"] == "Black"
    # the population baseline is the same ethnic group
    assert black.baseline_lifetime["breast"][0][1] == pytest.approx(b)
    with pytest.raises(ValueError, match="only defined for the UK"):
        _lifetime("Black", population="Germany")


def test_ethnicity_from_pedigree_header():
    text = """##CanRisk 2.0
##ethnicity=Black or Black British;Caribbean
##FamID\tName\tTarget\tIndivID\tFathID\tMothID\tSex\tMZtwin\tDead\tAge\tYob\tBC1\tBC2\tOC\tPRO\tPAN\tAshkn\tBRCA1\tBRCA2\tPALB2\tATM\tCHEK2\tBARD1\tRAD51D\tRAD51C\tBRIP1\tER:PR:HER2:CK14:CK56
F\tp\t1\tp\t0\t0\tF\t0\t0\t30\t1990\t0\t0\t0\t0\t0\t0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0\t0:0:0:0:0
"""
    ped = read_canrisk(text)[0]
    res, risk = _lifetime(None, pedigree=ped)
    assert res.settings["ethnicity"] == "Black"
    assert risk == pytest.approx(_lifetime("Black")[1])
    assert res.baseline_lifetime["breast"][0][1] == pytest.approx(risk)
    # an explicit setting overrides the header
    assert _lifetime("White", pedigree=ped)[1] == pytest.approx(_lifetime("White")[1])


def test_ethnicity_warnings():
    ped = Pedigree("X", [Person("p", "F", target=True, age=30, yob=1990)],
                   risk_factors={"parity": 1, "birads": "c"})
    res = calculate(ped, Settings(ethnicity="South Asian"), calcs=("lifetime",))
    assert any("risk-factor distributions" in w for w in res.warnings)
    assert any("BI-RADS" in w for w in res.warnings)
    res = calculate(ped, Settings(), calcs=("lifetime",))
    assert not any("ethnicity" in w for w in res.warnings)
