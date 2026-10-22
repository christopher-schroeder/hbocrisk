"""Command line interface.

    openboadicea risk family.canrisk [--model BC|OC|both] [--population UK] [--json]
    openboadicea populations
    openboadicea serve [--host 127.0.0.1] [--port 8080]
"""
import argparse
import json
import sys

from .api import Settings, calculate
from .incidence import available_populations
from .pedigree import read_canrisk_file
from .prs import prs_from_vcf
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


def main(argv=None):
    ap = argparse.ArgumentParser(prog="openboadicea", description=__doc__,
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
    r.add_argument("--json", action="store_true", help="JSON output")
    sub.add_parser("populations", help="list available incidence populations")
    s = sub.add_parser("serve", help="run a small JSON HTTP service")
    s.add_argument("--host", default="127.0.0.1")
    s.add_argument("--port", type=int, default=8080)
    args = ap.parse_args(argv)

    if args.cmd == "populations":
        print("\n".join(available_populations()))
        return 0
    if args.cmd == "serve":
        serve(args.host, args.port)
        return 0

    models = ["BC", "OC"] if args.model == "both" else [args.model]
    results = []
    for ped in read_canrisk_file(args.pedigree):
        if args.prs_vcf:
            for m in models:
                name = args.prs_name or ("BCAC 313" if m == "BC" else "OCAC 36")
                info = prs_from_vcf(args.prs_vcf, name, args.prs_sample)
                ped.prs[m] = (info["alpha"], info["zscore"])
                print(f"# {m} PRS {info['name']}: z = {info['zscore']:.3f}, alpha = "
                      f"{info['alpha']}, {info['n_found']}/{info['n_variants']} variants",
                      file=sys.stderr)
        for m in models:
            st = Settings(model=m, population=args.population,
                          mutation_frequencies=args.mutation_frequencies,
                          sensitivities=_sens(args.sensitivity))
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
