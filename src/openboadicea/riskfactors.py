"""Questionnaire-based risk factors, mammographic density and height.

Each risk factor is categorical with a population distribution and age-specific
relative risks (RRs) relative to a reference category, as in BOADICEA v5/v6 [LEE2019,
supplementary tables s2.1-s2.10] and the EOC model [LEE2022O, table s2]. Height is
treated as continuous [LEE2022]: the population distribution is discretised into five
bins with binomial B(4, 1/2) masses, and the individual's measured value is added as an
extra category of zero mass whose RR is the RR at the measured height.

Only the risk factors observed for the proband enter the model; unobserved factors are
omitted entirely (the baseline incidence is then the average over their distribution).
"""
from dataclasses import dataclass, field
import math
import re

import numpy as np
from scipy.stats import norm

from .params import AGES, MAX_AGE


@dataclass
class ObservedFactor:
    """A risk factor as used in a calculation: category masses, log-RR per age and the
    index of the proband's category."""
    name: str
    probs: np.ndarray            # (n,)
    logrr: np.ndarray            # (n, MAX_AGE)
    index: int                   # proband's category
    label: str = ""


def _step(rr_young, rr_old, split=20):
    """log-RR rows: RR 1 below `split`, rr_old (list per category) from `split`."""
    out = np.zeros((len(rr_old), MAX_AGE))
    for i, r in enumerate(rr_old):
        out[i, AGES >= split] = math.log(r)
        out[i, AGES < split] = math.log(rr_young)
    return out


@dataclass
class CategoricalFactor:
    name: str
    labels: list
    probs: list
    # function (current_age) -> (n, MAX_AGE) log RR matrix
    logrr_fn: callable = field(repr=False)
    # parser value(str) -> category index (0-based) or None if unobserved
    parse: callable = field(repr=False)
    # extra zero-RR categories with population mass not available for selection
    # (e.g. nulliparous women for age at first live birth)
    implicit_mass: float = 0.0

    def observe(self, value, current_age):
        idx = self.parse(value)
        if idx is None:
            return None
        probs = np.array(self.probs, float)
        logrr = self.logrr_fn(current_age)
        if self.implicit_mass > 0:
            probs = np.append(probs, self.implicit_mass)
            logrr = np.vstack([logrr, np.zeros(MAX_AGE)])
        probs = probs / probs.sum()
        return ObservedFactor(self.name, probs, logrr, idx, self.labels[idx])


@dataclass
class ContinuousNormalFactor:
    """Normally distributed RF with constant log-RR per SD from age `start`."""
    name: str
    mean: float
    sd: float
    logrr_per_sd: float
    nbins: int = 5
    start: int = 20

    def observe(self, value, current_age):
        try:
            x = float(value)
        except (TypeError, ValueError):
            return None
        if not np.isfinite(x) or x <= 0:
            return None
        n = self.nbins
        masses = np.array([math.comb(n - 1, k) for k in range(n)], float) / 2 ** (n - 1)
        edges = norm.ppf(np.concatenate([[0], np.cumsum(masses)]))
        b = self.logrr_per_sd
        # E[exp(b z) | z in bin] for z ~ N(0,1)
        rr = np.exp(b * b / 2) * (norm.cdf(edges[1:] - b) - norm.cdf(edges[:-1] - b)) / masses
        z = (x - self.mean) / self.sd
        rows = np.zeros((n + 1, MAX_AGE))
        rows[:n, AGES >= self.start] = np.log(rr)[:, None]
        rows[n, AGES >= self.start] = b * z
        probs = np.append(masses, 0.0)
        return ObservedFactor(self.name, probs, rows, n, f"{x:g}")


# --------------------------------------------------------------------------------------
# value parsers
# --------------------------------------------------------------------------------------
_NA = {"", "na", "-", "nan", "none", "unknown"}


def _num(value):
    if value is None:
        return None
    s = str(value).strip()
    if s.lower() in _NA:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _range_parser(bounds):
    """bounds: list of (lo_inclusive, hi_exclusive) numeric intervals."""
    def parse(value):
        x = _num(value)
        if x is None:
            return None
        for i, (lo, hi) in enumerate(bounds):
            if lo <= x < hi:
                return i
        return None
    return parse


def _choice_parser(options):
    """options: list of sets of accepted lower-case strings per category."""
    def parse(value):
        if value is None:
            return None
        s = str(value).strip().lower()
        if ":" in s:
            s = s.split(":")[0]
        for i, opts in enumerate(options):
            if s in opts:
                return i
        return None
    return parse


INF = float("inf")

# --------------------------------------------------------------------------------------
# BOADICEA (breast cancer) risk factors  [LEE2019] tables s2.1-s2.10, [LEE2022]
# --------------------------------------------------------------------------------------


def _menopause_logrr(current_age):
    split = [39, 40, 45, 50, 55]
    rr = [0.67, 0.73, 0.86, 1.00, 1.12]
    out = np.zeros((5, MAX_AGE))
    for i in range(5):
        out[i, AGES >= split[i]] = math.log(rr[i])
    return out


def _aflb_logrr(current_age):
    split = [19, 20, 25, 30]
    rr = [1.00, 1.01, 1.11, 1.24]
    out = np.zeros((4, MAX_AGE))
    for i in range(4):
        out[i, AGES >= split[i]] = math.log(rr[i])
    return out


def _bc_oc_use_logrr(current_age):
    rr_lt50 = [1.0, 1.12, 1.33]
    rr_ge50 = [1.0, 1.12, 1.12]
    out = np.zeros((3, MAX_AGE))
    for i in range(3):
        sel = AGES >= current_age
        out[i, sel & (AGES < 50)] = math.log(rr_lt50[i])
        out[i, sel & (AGES >= 50)] = math.log(rr_ge50[i])
    return out


def _bc_mht_logrr(current_age):
    rr = [1.0, 1.57, 2.19]
    out = np.zeros((3, MAX_AGE))
    sel = (AGES >= current_age) & (AGES < current_age + 5)
    for i in range(3):
        out[i, sel] = math.log(rr[i])
    return out


def _age_banded(rr_20_49, rr_ge50):
    def fn(current_age):
        out = np.zeros((len(rr_20_49), MAX_AGE))
        for i in range(len(rr_20_49)):
            out[i, (AGES >= 20) & (AGES < 50)] = math.log(rr_20_49[i])
            out[i, AGES >= 50] = math.log(rr_ge50[i])
        return out
    return fn


def _parse_menarche(value):
    x = _num(value)
    if x is None:
        return None
    if x < 11:
        return 0
    if x > 15:
        return 6
    return int(x) - 10


def _parse_parity_bc(value):
    x = _num(value)
    if x is None:
        return None
    return min(int(x), 3)


def _parse_mht_bc(value):
    s = "" if value is None else str(value).strip().lower()
    if s in {"n", "f", "never/former", "never", "former"}:
        return 0
    if s in {"e", "current e-type"}:
        return 1
    if s in {"c", "current other/unknown type (including combined type)", "current"}:
        return 2
    return None


def _parse_menopause(value):
    s = "" if value is None else str(value).strip().lower()
    if s == "n":          # pre-menopausal: not applicable
        return None
    return _range_parser([(-INF, 40), (40, 45), (45, 50), (50, 55), (55, INF)])(value)


def _parse_birads(value):
    if value is None:
        return None
    s = str(value).strip().lower().replace("bi-rads", "").strip()
    m = {"a": 0, "b": 1, "c": 2, "d": 3, "1": 0, "2": 1, "3": 2, "4": 3}
    return m.get(s)


HEIGHT = ContinuousNormalFactor("height", mean=162.81, sd=6.452, logrr_per_sd=0.10130)

BC_RISK_FACTORS = {
    "menarche": CategoricalFactor(
        "menarche", ["<11", "11", "12", "13", "14", "15", ">15"],
        [0.122, 0.216, 0.268, 0.215, 0.115, 0.044, 0.019],
        lambda c: _step(1.0, [1.19, 1.09, 1.07, 1.00, 0.98, 0.92, 0.82]),
        _parse_menarche),
    "parity": CategoricalFactor(
        "parity", ["0", "1", "2", ">2"], [0.28, 0.16, 0.33, 0.23],
        lambda c: _step(1.0, [1.0, 0.87, 0.81, 0.71]), _parse_parity_bc),
    "first_live_birth": CategoricalFactor(
        "first_live_birth", ["<20", "20-24", "25-29", ">29"], [0.12, 0.23, 0.23, 0.14],
        _aflb_logrr, _range_parser([(-INF, 20), (20, 25), (25, 30), (30, INF)]),
        implicit_mass=0.28),
    "oc_use": CategoricalFactor(
        "oc_use", ["never", "former", "current"], [0.185, 0.684, 0.131],
        _bc_oc_use_logrr,
        _choice_parser([{"n", "never"}, {"f", "former"}, {"c", "current"}])),
    "mht_use": CategoricalFactor(
        "mht_use", ["never/former", "current e-type", "current other/unknown type"],
        [0.913, 0.010962, 0.076038], _bc_mht_logrr, _parse_mht_bc),
    "bmi": CategoricalFactor(
        "bmi", ["<18.5", "18.5-<25", "25-<30", ">=30"], [0.033, 0.357, 0.358, 0.252],
        _age_banded([1.28, 1.0, 0.92, 0.74], [1.0, 1.0, 1.13, 1.25]),
        _range_parser([(-INF, 18.5), (18.5, 25), (25, 30), (30, INF)])),
    "alcohol": CategoricalFactor(
        "alcohol", ["0", "<5", "5-<15", "15-<25", "25-<35", "35-<45", ">=45"],
        [0.106, 0.307, 0.276, 0.133, 0.080, 0.036, 0.061],
        lambda c: _step(1.0, [1.0, 1.01, 1.03, 1.13, 1.21, 1.32, 1.46]),
        lambda v: (None if _num(v) is None else
                   (0 if _num(v) <= 0 else
                    _range_parser([(0, 5), (5, 15), (15, 25), (25, 35), (35, 45),
                                   (45, INF)])(v) + 1))),
    "menopause": CategoricalFactor(
        "menopause", ["<40", "40-44", "45-49", "50-54", ">54"],
        [0.022, 0.070, 0.234, 0.533, 0.14], _menopause_logrr, _parse_menopause),
    "birads": CategoricalFactor(
        "birads", ["BI-RADS a", "BI-RADS b", "BI-RADS c", "BI-RADS d"],
        [0.041, 0.352, 0.470, 0.137],
        _age_banded([0.41, 1.0, 1.75, 3.33], [1.0, 2.04, 2.81, 4.08]), _parse_birads),
    "height": HEIGHT,
}

# --------------------------------------------------------------------------------------
# EOC model risk factors [LEE2022O] table s2 (height continuous as in [LEE2022])
# --------------------------------------------------------------------------------------


def _from_current(rr):
    def fn(current_age):
        out = np.zeros((len(rr), MAX_AGE))
        for i in range(len(rr)):
            out[i, AGES >= current_age] = math.log(rr[i])
        return out
    return fn


def _parse_oc_years(value):
    if value is None:
        return None
    s = str(value).strip().lower()
    if s in _NA:
        return None
    if s in {"n", "never", "<1", "never or <1"}:
        return 0
    if ":" in s:
        s = s.split(":", 1)[1].strip()
    if s in {"<1", "< 1"}:
        return 0
    x = _num(s)
    if x is None:
        return None
    return _range_parser([(-INF, 1), (1, 5), (5, 10), (10, 15), (15, INF)])(x)


def _parse_yes_no(value):
    return _choice_parser([{"n", "no"}, {"y", "yes"}])(value)


def _parse_parity_oc(value):
    x = _num(value)
    if x is None:
        return None
    return min(int(x), 2)


OC_RISK_FACTORS = {
    "parity": CategoricalFactor(
        "parity", ["0", "1", ">1"], [0.28, 0.16, 0.56],
        lambda c: _step(1.0, [1.0, 0.76, 0.58]), _parse_parity_oc),
    "oc_use": CategoricalFactor(
        "oc_use", ["never or <1", "1-4", "5-9", "10-14", ">=15"],
        [0.710, 0.138, 0.086, 0.046, 0.020],
        _from_current([1.0, 0.78, 0.64, 0.56, 0.42]), _parse_oc_years),
    "mht_use": CategoricalFactor(
        "mht_use", ["never", "ever"], [0.73, 0.27], _from_current([1.0, 1.37]),
        _choice_parser([{"n", "never"}, {"e", "f", "c", "ever"}])),
    "tl": CategoricalFactor(
        "tl", ["no", "yes"], [0.77, 0.23], lambda c: _step(1.0, [1.0, 0.74]),
        _parse_yes_no),
    "endo": CategoricalFactor(
        "endo", ["no", "yes"], [0.90, 0.10], lambda c: _step(1.0, [1.0, 1.53]),
        _parse_yes_no),
    "bmi": CategoricalFactor(
        "bmi", ["<22.5", "22.5-<30", ">=30"], [0.325, 0.529, 0.146],
        lambda c: _step(1.0, [1.0, 1.066, 1.13]),
        _range_parser([(-INF, 22.5), (22.5, 30), (30, INF)])),
    "height": HEIGHT,
}

# synonyms accepted in CanRisk file headers / dicts
SYNONYMS = {
    "menarche_age": "menarche", "age_of_first_live_birth": "first_live_birth",
    "oral_contraception": "oc_use", "mht": "mht_use", "alcohol_intake": "alcohol",
    "age_of_menopause": "menopause", "tubal_ligation": "tl", "endometriosis": "endo",
    "mammographic_density": "birads", "bi-rads": "birads",
}


def normalise_name(name):
    n = re.sub(r"[^a-z0-9_\-]", "", str(name).strip().lower())
    return SYNONYMS.get(n, n)


# header entries that are not risk factors of either model
_IGNORED = {"ethnicity"}
_NOT_IMPLEMENTED = {"volpara", "stratus"}


def unused_risk_factors(spec, values):
    """Warnings for supplied risk factors that are not used by the model."""
    from .models import model_risk_factors
    defs = model_risk_factors(spec)
    other = OC_RISK_FACTORS if defs is BC_RISK_FACTORS else BC_RISK_FACTORS
    out = []
    for raw, value in (values or {}).items():
        name = normalise_name(raw)
        if name in _NOT_IMPLEMENTED:
            out.append(f"risk factor '{raw}' (continuous mammographic density) is not "
                       f"implemented and was ignored; use BI-RADS")
        elif name in defs:
            if str(value).strip().lower() not in _NA and \
                    defs[name].observe(value, 40) is None and name != "menopause":
                out.append(f"risk factor '{raw}={value}' could not be interpreted and "
                           f"was ignored")
        elif name not in other and name not in _IGNORED:
            out.append(f"unknown risk factor '{raw}' ignored")
    return out


def observe(factor_defs, values, current_age):
    """Return list of ObservedFactor for the observed risk factors."""
    out = []
    for raw_name, value in (values or {}).items():
        name = normalise_name(raw_name)
        f = factor_defs.get(name)
        if f is None:
            continue
        o = f.observe(value, current_age)
        if o is not None:
            out.append(o)
    out.sort(key=lambda o: o.name)
    return out


# --------------------------------------------------------------------------------------
# Combination of factors into joint categories (with amalgamation)
# --------------------------------------------------------------------------------------


def combine(factors, max_categories=1500):
    """Combine independent factors into joint categories.

    Returns (probs (m,), logrr (m, MAX_AGE), proband_index). When the number of joint
    categories exceeds `max_categories` they are amalgamated: categories are ordered
    by their mean log-RR over ages 20-79 and adjacent categories merged into bins of
    approximately equal mass, with the RR of a merged category the mass-weighted
    average of its constituents (as in BOADICEA v5, supplement section 2). The
    proband's category is never merged.
    """
    probs = np.ones(1)
    logrr = np.zeros((1, MAX_AGE))
    idx = 0
    for f in factors:
        n = len(f.probs)
        probs = (probs[:, None] * f.probs[None, :]).reshape(-1)
        logrr = (logrr[:, None, :] + f.logrr[None, :, :]).reshape(-1, MAX_AGE)
        idx = idx * n + f.index
        if len(probs) > max_categories:
            probs, logrr, idx = _amalgamate(probs, logrr, idx, max_categories // 2)
    return probs, logrr, idx


def _amalgamate(probs, logrr, idx, nbins):
    m = len(probs)
    keep = np.ones(m, bool)
    keep[idx] = False
    own_p, own_l = probs[idx], logrr[idx]
    others = np.nonzero(keep)[0]
    key = logrr[others, 20:].mean(axis=1)
    order = others[np.argsort(key, kind="stable")]
    p = probs[order]
    cum = np.cumsum(p) / max(p.sum(), 1e-300)
    bins = np.minimum((cum * nbins - 1e-12).astype(int), nbins - 1)
    starts = np.flatnonzero(np.r_[True, bins[1:] != bins[:-1]])
    pb = np.add.reduceat(p, starts)
    rrsum = np.add.reduceat(p[:, None] * np.exp(logrr[order]), starts, axis=0)
    ok = pb > 0
    newp = np.append(pb[ok], own_p)
    newl = np.vstack([np.log(rrsum[ok] / pb[ok, None]), own_l[None, :]])
    return newp, newl, len(newp) - 1
