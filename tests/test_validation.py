"""Regression tests against published BOADICEA / CanRisk results.

Tolerances reflect differences in the underlying population incidence data (the exact
CanRisk incidence tables are not published) and model version differences.
"""
import pytest

from hbocrisk import Pedigree, Person, GeneticTest, calculate, Settings


def _family(gene=None, mother_bc=None, yob=1985, all_negative=False):
    tests = {}
    if gene:
        tests[gene] = GeneticTest("T", "P")
    if all_negative:
        tests = {g: GeneticTest("S", "N") for g in
                 ["BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM", "BARD1", "RAD51C", "RAD51D"]}
    people = []
    if mother_bc:
        people = [Person("m", "F", age=mother_bc, yob=yob - 25, bc1=mother_bc),
                  Person("f", "M")]
        p = Person("p", "F", father="f", mother="m", target=True, age=20, yob=yob,
                   tests=tests)
    else:
        p = Person("p", "F", target=True, age=20, yob=yob, tests=tests)
    return Pedigree("V", people + [p])


def lifetime(**kw):
    return calculate(_family(**kw), calcs=("lifetime",)).lifetime["breast"][0][1]


@pytest.mark.parametrize("gene,expected", [
    (None, 0.120), ("BARD1", 0.235), ("RAD51C", 0.223), ("RAD51D", 0.208)])
def test_lee2022_table_s1_unknown_family_history(gene, expected):
    # Lee et al. J Med Genet 2022, supplementary table s1 (lifetime risk 20-80, born 1985)
    assert lifetime(gene=gene) == pytest.approx(expected, rel=0.03)


@pytest.mark.parametrize("gene,expected", [
    (None, 0.193), ("BARD1", 0.337), ("RAD51C", 0.322), ("RAD51D", 0.304)])
def test_lee2022_table_s1_mother_affected_50(gene, expected):
    assert lifetime(gene=gene, mother_bc=50) == pytest.approx(expected, rel=0.03)


def test_lee2022_table_s1_negative_panel():
    assert lifetime(all_negative=True) == pytest.approx(0.116, rel=0.03)


def test_ficorella2024_brca2():
    # BRCA2 carrier born 2000: 58% by age 80; 72.5% with a mother affected at 40
    assert lifetime(gene="BRCA2", yob=2000) == pytest.approx(0.58, rel=0.03)
    assert lifetime(gene="BRCA2", yob=2000, mother_bc=40) == pytest.approx(0.725, rel=0.03)


def test_lee2022_palb2_ovarian():
    ped = _family(gene="PALB2")
    r = calculate(ped, Settings(model="OC"), calcs=("lifetime",))
    assert r.lifetime["ovarian"][0][1] == pytest.approx(0.050, rel=0.06)


def test_canrisk_732_population_outputs():
    """Unaffected 25-year-old, no family history (CanRisk bws test suite output of
    'boadicea model 7.3.2'): lifetime 0.1200146, 10-year 0.0171806, risk to 80
    0.1198730."""
    r = calculate(Pedigree("T", [Person("p", "F", target=True, age=25, yob=2001)]))
    assert r.lifetime["breast"][0][1] == pytest.approx(0.1200146, rel=0.02)
    assert r.ten_year["breast"][0][1] == pytest.approx(0.0171806, rel=0.02)
    assert r.remaining_lifetime["breast"][-1][1] == pytest.approx(0.1198730, rel=0.02)
