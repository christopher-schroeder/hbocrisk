"""Compare openboadicea against published BOADICEA/CanRisk results.

Run: python scripts/validate.py
"""
from openboadicea import Pedigree, Person, GeneticTest, calculate, Settings

PANEL = ["BRCA1", "BRCA2", "PALB2", "CHEK2", "ATM", "BARD1", "RAD51C", "RAD51D"]


def family(gene=None, mother_bc=None, yob=1985, age=20, negative=False):
    tests = {gene: GeneticTest("T", "P")} if gene else {}
    if negative:
        tests = {g: GeneticTest("S", "N") for g in PANEL}
    if mother_bc:
        return Pedigree("V", [Person("m", "F", age=mother_bc, yob=yob - 25, bc1=mother_bc),
                              Person("f", "M"),
                              Person("p", "F", father="f", mother="m", target=True,
                                     age=age, yob=yob, tests=tests)])
    return Pedigree("V", [Person("p", "F", target=True, age=age, yob=yob, tests=tests)])


def row(label, ours, ref, source):
    print(f"| {label} | {100 * ours:.4g} | {100 * ref:.4g} | {100 * (ours / ref - 1):+.1f}% | {source} |")


def main():
    print("| scenario | openboadicea (%) | published (%) | rel. diff | source |")
    print("|---|---|---|---|---|")
    r = calculate(Pedigree("T", [Person("p", "F", target=True, age=25, yob=2001)]))
    src = "CanRisk 7.3.2 web-service test output"
    row("population, age 25: lifetime BC (20-80)", r.lifetime["breast"][0][1], 0.1200146, src)
    row("population, age 25: 10-year BC (40-50)", r.ten_year["breast"][0][1], 0.0171806, src)
    row("population, age 25: BC risk to 80", r.remaining_lifetime["breast"][-1][1], 0.1198730, src)
    ref = {"BRCA1": 0.0012713, "BRCA2": 0.0020312, "PALB2": 0.0012758, "CHEK2": 0.0074259,
           "ATM": 0.0035743}
    for g, v in ref.items():
        row(f"population, age 25: P({g} carrier)", r.carrier_probabilities[g], v, src)
    r = calculate(Pedigree("T", [Person("p", "F", target=True, age=64, yob=1962)]))
    row("population, age 64 (yob 1962): lifetime BC", r.lifetime["breast"][0][1], 0.11897, src)
    row("population, age 64 (yob 1962): BC risk to 80", r.remaining_lifetime["breast"][-1][1],
        0.0606968, src)
    ref = {"BRCA1": 0.000408, "BRCA2": 0.0010558, "CHEK2": 0.0068254}
    for g, v in ref.items():
        row(f"population, age 64: P({g} carrier)", r.carrier_probabilities[g], v, src)
    o = calculate(Pedigree("T", [Person("p", "F", target=True, age=64, yob=1962)]),
                  Settings(model="OC"))
    src_o = "CanRisk ovarian model 3.2.2 web-service test output"
    row("population, age 64: lifetime OC (EOC model)", o.lifetime["ovarian"][0][1], 0.0176191, src_o)
    row("population, age 64: OC risk to 80 (EOC model)", o.remaining_lifetime["ovarian"][-1][1],
        0.0098928, src_o)
    src = "Lee et al. JMG 2022, table s1 (born 1985)"
    for gene, ref_u, ref_m in ((None, 0.120, 0.193), ("BARD1", 0.235, 0.337),
                               ("RAD51C", 0.223, 0.322), ("RAD51D", 0.208, 0.304)):
        lab = gene or "untested"
        row(f"{lab}, unknown FH: lifetime BC", calculate(family(gene)).lifetime["breast"][0][1],
            ref_u, src)
        row(f"{lab}, mother BC at 50: lifetime BC",
            calculate(family(gene, mother_bc=50)).lifetime["breast"][0][1], ref_m, src)
    row("negative panel, unknown FH: lifetime BC",
        calculate(family(negative=True)).lifetime["breast"][0][1], 0.116, src)
    row("negative panel, mother BC at 50: lifetime BC",
        calculate(family(negative=True, mother_bc=50)).lifetime["breast"][0][1], 0.179, src)
    src = "Ficorella et al. BJC Rep 2024 (born 2000)"
    row("BRCA2 carrier: BC risk 20-80", calculate(family("BRCA2", yob=2000)).lifetime["breast"][0][1],
        0.58, src)
    row("BRCA2 carrier, mother BC at 40: BC risk 20-80",
        calculate(family("BRCA2", yob=2000, mother_bc=40)).lifetime["breast"][0][1], 0.725, src)
    row("PALB2 carrier: lifetime OC (EOC model)",
        calculate(family("PALB2"), Settings(model="OC")).lifetime["ovarian"][0][1], 0.050,
        "Lee et al. JMG 2022")


if __name__ == "__main__":
    main()
