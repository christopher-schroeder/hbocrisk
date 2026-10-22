"""Build birth-cohort specific, smoothed single-year cancer incidence tables.

Follows the procedure described for BOADICEA v6 (Lee et al., J Med Genet 2022,
supplementary material "Population Incidences") and Lee et al., Br J Cancer 2014:

1. Calendar-year and 5-year-age-group specific incidences are taken from the
   IARC "Cancer Incidence in Five Continents" CI5plus detailed database.
2. Incidences for calendar years before/after the first/last year available are
   the average of the first/last five years of data.
3. For each decadal birth cohort, the age-specific incidence at age ``a`` is the
   average over all birth years ``b`` in the cohort of the incidence applicable to
   someone aged ``a`` in calendar year ``b + a``.
4. The resulting cohort incidences are smoothed over age with LOWESS (local linear
   regression, tricube weights, bandwidth 0.2).

For the UK, BOADICEA uses Office for National Statistics (ONS) registrations for
England from 1992 onwards. ONS registration counts by 5-year age group (parsed with
``scripts/parse_ons_england.py`` from the ONS "Cancer Registration Statistics,
England" reference tables, available here for 2009-2017) replace the CI5plus counts for
those years; CI5plus England person-years are used as denominators.

Usage::

    python scripts/build_incidences.py data_raw/ci5plus/subset.csv \
        src/hbocrisk/data/incidences [data_raw/ons/ons_counts.json]

where ``subset.csv`` is the CI5plus ``data.csv`` restricted to the cancer codes
used here (pancreas 42, breast 67, ovary 86, prostate 101, other female genital
organs 98).
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

# CI5plus cancer codes
CANCERS = {"pancreas": 42, "breast": 67, "ovary": 86, "prostate": 101,
           "other_female_genital": 98}

# Output series: name -> (CI5 cancer code(s), CI5 sex code (1 male, 2 female)).
# Ovarian cancer includes fallopian tube cancers (CanRisk definition; ONS C56-C57),
# taken from CI5plus "Ovary" (C56) plus "Other female genital organs" (C57).
SERIES = {
    "breast_female": (67, 2),
    "breast_male": (67, 1),
    "ovary": ((86, 98), 2),
    "prostate": (101, 1),
    "pancreas_female": (42, 2),
    "pancreas_male": (42, 1),
}

# Populations: name -> list of CI5plus registry ids (cases and person-years summed)
POPULATIONS = {
    "UK": [82600199],  # England
    "Australia": [3600299, 3600399, 3600499, 3600599, 3600699, 3600799, 3600899],
    "Canada": [12400099],
    "USA": [84002010],  # SEER 9 registries, whites (as in BOADICEA)
    "Denmark": [20800099],
    "Estonia": [23300099],
    "France": [25000199, 25000299, 25000399, 25000499, 25000599, 25000699,
               25000899, 25001499, 25001699, 25001799],
    "Iceland": [35200099],
    "Netherlands": [52800099],
    "New-Zealand": [55400099],
    "Norway": [57800099],
    "Slovenia": [70500099],
    "Spain": [72400199, 72400299, 72400399, 72400499, 72400699, 72401099,
              72401199, 72401499],
    # not part of CanRisk, provided as extra options
    "Germany": [27600299, 27601299, 27601699, 27603099],
    "Austria": [4000099],
    "Switzerland": [75600299, 75600599, 75600899, 75601099, 75601199],
}

COHORTS = list(range(1910, 1990, 10))  # 1910-1919 ... 1980-1989 (BOADICEA v6)
MAX_AGE = 85        # ages 0..84 used for smoothing
OUT_AGES = 80       # ages 0..79 stored
BANDWIDTH = 0.2


def lowess(x, y, frac=BANDWIDTH):
    """Locally weighted linear regression (no robustness iterations)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    k = int(np.ceil(frac * n))
    out = np.empty(n)
    for i in range(n):
        d = np.abs(x - x[i])
        h = np.sort(d)[k - 1]
        h = max(h, 1e-12) * 1.0000001
        w = np.clip(1 - (d / h) ** 3, 0, None) ** 3
        sw = w.sum()
        xm = (w * x).sum() / sw
        ym = (w * y).sum() / sw
        sxx = (w * (x - xm) ** 2).sum()
        b = 0.0 if sxx <= 0 else (w * (x - xm) * (y - ym)).sum() / sxx
        out[i] = ym + b * (x[i] - xm)
    return out


def read_subset(path):
    """cases[(id, sex, cancer)][year][agegroup], py[(id, sex)][year][agegroup]"""
    cases = defaultdict(lambda: defaultdict(lambda: np.zeros(18)))
    py = defaultdict(lambda: defaultdict(lambda: np.zeros(18)))
    with open(path) as fh:
        for r in csv.DictReader(fh):
            age = int(r["age"])
            if age > 18:
                continue
            key = (int(r["id_code"]), int(r["sex"]))
            year = int(r["year"])
            cases[key + (int(r["cancer_code"]),)][year][age - 1] += float(r["cases"])
            # person-years are repeated for each cancer: keep one copy
            py[key][year][age - 1] = float(r["py"])
    return cases, py


ONS_SERIES = {"breast_female": ["C50_Females"], "breast_male": ["C50_Males"],
              "ovary": ["C56_Females", "C57_Females"], "prostate": ["C61_Males"],
              "pancreas_female": ["C25_Females"], "pancreas_male": ["C25_Males"]}


def calendar_rates(cases, py, registries, code, sex, override=None):
    """Return (years, rates[year, agegroup]) for a population (registries summed)."""
    years = set()
    for reg in registries:
        years |= set(py.get((reg, sex), {}).keys())
    years = sorted(years)
    if not years:
        return None, None
    c = np.zeros((len(years), 18))
    p = np.zeros((len(years), 18))
    for i, y in enumerate(years):
        for reg in registries:
            pyr = py.get((reg, sex), {}).get(y)
            if pyr is None:
                continue
            p[i] += pyr
            for cc in (code if isinstance(code, tuple) else (code,)):
                c[i] += cases.get((reg, sex, cc), {}).get(y, np.zeros(18))
        if override is not None and y in override:
            c[i] = override[y]
    with np.errstate(invalid="ignore", divide="ignore"):
        rates = np.where(p > 0, c / p, np.nan)
    keep = ~np.all(np.isnan(rates), axis=1)
    return np.array(years)[keep], rates[keep]


def cohort_incidence(years, rates, cohort_start):
    """Average over birth years in the decade, then LOWESS smooth over age."""
    first = np.nanmean(rates[:5], axis=0)
    last = np.nanmean(rates[-5:], axis=0)
    ymin, ymax = years[0], years[-1]
    index = {y: i for i, y in enumerate(years)}
    inc = np.zeros(MAX_AGE)
    for age in range(MAX_AGE):
        grp = min(age // 5, 17)
        vals = []
        for b in range(cohort_start, cohort_start + 10):
            y = b + age
            if y < ymin:
                v = first[grp]
            elif y > ymax:
                v = last[grp]
            elif y in index:
                v = rates[index[y], grp]
            else:  # gap in calendar years: nearest available year
                j = int(np.argmin(np.abs(years - y)))
                v = rates[j, grp]
            if not np.isnan(v):
                vals.append(v)
        inc[age] = np.mean(vals) if vals else 0.0
    smooth = np.clip(lowess(np.arange(MAX_AGE), inc), 0.0, None)
    return smooth[:OUT_AGES]


def load_ons(path):
    """ONS counts -> {series: {year: counts[18]}}"""
    d = json.load(open(path))
    out = {}
    for sname, keys in ONS_SERIES.items():
        out[sname] = {}
        for y, sites in d.items():
            if all(k in sites for k in keys):
                out[sname][int(y)] = np.sum([sites[k] for k in keys], axis=0)
    return out


def main(subset, outdir, ons=None):
    cases, py = read_subset(subset)
    ons = load_ons(ons) if ons else None
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for name, regs in POPULATIONS.items():
        source = "IARC CI5plus detailed (2024 v1.1)"
        if ons and name == "UK":
            source += "; ONS Cancer Registration Statistics England 2009-2017 (counts)"
        table = {"population": name, "source": source,
                 "registries": regs, "cohorts": COHORTS, "ages": [0, OUT_AGES - 1],
                 "units": "cases per person-year", "series": {}, "calendar_years": {}}
        ok = True
        for sname, (code, sex) in SERIES.items():
            override = ons[sname] if (ons and name == "UK") else None
            years, rates = calendar_rates(cases, py, regs, code, sex, override)
            if years is None:
                ok = False
                break
            table["calendar_years"][sname] = [int(years[0]), int(years[-1])]
            table["series"][sname] = [
                [round(float(v), 9) for v in cohort_incidence(years, rates, c)]
                for c in COHORTS]
        if not ok:
            print("skipping", name)
            continue
        fn = outdir / f"{name}.json"
        fn.write_text(json.dumps(table, separators=(",", ":")))
        print("wrote", fn, table["calendar_years"]["breast_female"])


if __name__ == "__main__":
    main(*sys.argv[1:4])
