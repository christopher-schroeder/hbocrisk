"""Command line interface.

    hbocrisk risk family.canrisk [--model BC|OC|both] [--population UK] [--json]
    hbocrisk populations
    hbocrisk fetch-data [--dir DIR]
    hbocrisk serve [--host 127.0.0.1] [--port 8080]
"""
import argparse
import json
import os
import sys

from .api import Settings, calculate
from .incidence import available_populations
from .pedigree import read_canrisk_file
from .prs import ANCESTRY_PRS, ancestry_file, prs_from_vcf, ancestry_name
from .server import serve


def _sens(items):
    out = {}
    for it in items or []:
        g, v = it.split("=", 1)
        out[g.upper()] = float(v)
    return out


def _print_result(r, out):
    print(f"# {r.model} model  family={r.famid}  proband={r.proband}", file=out)
    if r.carrier_probabilities:
        print("Carrier probabilities:", file=out)
        for g, v in r.carrier_probabilities.items():
            print(f"  {g:24s} {100 * v:8.3f} %", file=out)
    for title, attr in (("Remaining lifetime risk", "remaining_lifetime"),
                        ("Lifetime risk (age 20-80)", "lifetime"),
                        ("10-year risk (age 40-50)", "ten_year")):
        d = getattr(r, attr)
        if not d:
            continue
        base = getattr(r, "baseline_" + attr) or {}
        print(f"{title}:", file=out)
        for cancer, rows in d.items():
            b = dict(base.get(cancer, []))
            for age, risk in rows:
                pop = f"   (population {100 * b[age]:6.2f} %)" if age in b else ""
                print(f"  {cancer:22s} to {age:2d}: {100 * risk:6.2f} %{pop}", file=out)
    for w in r.warnings:
        print(f"WARNING: {w}", file=out)


def _ancestry(value):
    """'South Asian' -> 'South Asian'; 'EUR=0.6,AFR=0.4' -> {'EUR': 0.6, 'AFR': 0.4}."""
    if value is None or "=" not in value:
        return value
    out = {}
    for part in value.split(","):
        k, v = part.split("=", 1)
        out[k.strip()] = float(v)
    return out


def ancestry_name_or_none(ancestry):
    if isinstance(ancestry, dict):
        return "mixed"
    return ancestry_name(ancestry)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="hbocrisk", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("risk", help="compute carrier probabilities and cancer risks")
    r.add_argument("pedigree", help="CanRisk (v1-v4) or BOADICEA v4 pedigree file")
    r.add_argument("--model", default="BC", choices=["BC", "OC", "both"])
    r.add_argument("--population", default="UK", help="cancer incidence population")
    r.add_argument("--mutation-frequencies", default="UK",
                   help="UK, 'UK, non-European', Ashkenazi or Iceland")
    r.add_argument("--sensitivity", nargs="*", metavar="GENE=VALUE",
                   help="mutation search sensitivities, e.g. BRCA1=0.9")
    r.add_argument("--prs-vcf", help="VCF with the proband's genotypes/dosages; the PRS "
                   "z-score is computed and added to the pedigree")
    r.add_argument("--prs-name", default=None,
                   help="PRS (e.g. 'BCAC 313', 'OCAC 36') or path to a .prs file; default "
                   "'BCAC 313' for BC and 'OCAC 36' for OC")
    r.add_argument("--prs-sample", default=None, help="sample name in the VCF")
    r.add_argument("--prs-ancestry", default=None,
                   help="genetic ancestry for the PRS normalisation: European, African, "
                   "'East Asian', 'South Asian' or proportions, e.g. EUR=0.6,AFR=0.4; "
                   "non-European ancestry defaults to 'BCAC 309' for BC")
    r.add_argument("--ethnicity", default=None,
                   help="UK ethnic group (White, Black, 'South Asian', 'East Asian', "
                   "Mixed, Other); UK incidences only. Default: ##ethnicity header")
    r.add_argument("--json", action="store_true", help="JSON output")
    sub.add_parser("populations", help="list available incidence populations")
    d = sub.add_parser("fetch-data", help="download the ancestry-specific PRS files "
                       "(GPL-3.0, CCGE-BOADICEA/SHARE-PRScalculation)")
    d.add_argument("--dir", default=None,
                   help="target directory (default $HBOCRISK_DATA or ~/.cache/hbocrisk)")
    s = sub.add_parser("serve", help="run a small JSON HTTP service")
    s.add_argument("--host", default="127.0.0.1")
    s.add_argument("--port", type=int, default=8080)
    args = ap.parse_args(argv)

    if args.cmd == "populations":
        print("\n".join(available_populations()))
        return 0
    if args.cmd == "fetch-data":
        if args.dir:
            os.environ["HBOCRISK_DATA"] = args.dir
        for name in ANCESTRY_PRS:
            print(ancestry_file(name))
        return 0
    if args.cmd == "serve":
        serve(args.host, args.port)
        return 0

    models = ["BC", "OC"] if args.model == "both" else [args.model]
    results = []
    for ped in read_canrisk_file(args.pedigree):
        if args.prs_vcf:
            for m in models:
                ancestry = _ancestry(args.prs_ancestry)
                name = args.prs_name or ("OCAC 36" if m == "OC" else
                                         "BCAC 313" if ancestry_name_or_none(ancestry)
                                         in (None, "European") else "BCAC 309")
                info = prs_from_vcf(args.prs_vcf, name, args.prs_sample, ancestry)
                ped.prs[m] = (info["alpha"], info["zscore"])
                print(f"# {m} PRS {info['name']}: z = {info['zscore']:.3f}, alpha = "
                      f"{info['alpha']:.4f}, {info['n_found']}/{info['n_variants']} "
                      f"variants, ancestry {info['ancestry']} ({info['method']})",
                      file=sys.stderr)
        for m in models:
            st = Settings(model=m, population=args.population,
                          mutation_frequencies=args.mutation_frequencies,
                          sensitivities=_sens(args.sensitivity),
                          ethnicity=args.ethnicity)
            results.append(calculate(ped, st))
    if args.json:
        json.dump([x.to_dict() for x in results], sys.stdout, indent=1)
        print()
    else:
        for x in results:
            _print_result(x, sys.stdout)
            print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
