"""Population cancer incidence tables (birth-cohort and country specific).

Tables are produced by ``scripts/build_incidences.py`` from IARC CI5plus data and
stored as JSON in ``hbocrisk/data/incidences``. Users can supply their own tables
(e.g. national registry data) with :func:`register_population`.
"""
from functools import lru_cache
from importlib import resources
import json

import numpy as np

from .params import MAX_AGE, ETHNICITY_IRR

SERIES = ("breast_female", "breast_male", "ovary", "prostate", "pancreas_female",
          "pancreas_male")

_CUSTOM = {}

# CanRisk names -> table names
ALIASES = {"Other": "UK", "New Zealand": "New-Zealand", "United Kingdom": "UK",
           "England": "UK", "The Netherlands": "Netherlands"}


class IncidenceTable:
    def __init__(self, name, cohorts, series, meta=None):
        self.name = name
        self.cohorts = list(cohorts)            # first birth year of each decadal cohort
        self.series = {k: np.asarray(v, float) for k, v in series.items()}
        for k, v in self.series.items():
            if v.shape != (len(self.cohorts), MAX_AGE):
                raise ValueError(f"incidence series {k} must have shape "
                                 f"({len(self.cohorts)}, {MAX_AGE})")
        self.meta = meta or {}

    def cohort_index(self, yob):
        """Decadal cohort for a year of birth; earlier/later births use the first/last."""
        idx = 0
        for i, c in enumerate(self.cohorts):
            if yob >= c:
                idx = i
        return idx

    def rates(self, series, cohort_idx):
        return self.series[series][cohort_idx]


def available_populations():
    names = set(_CUSTOM)
    for f in resources.files("hbocrisk").joinpath("data", "incidences").iterdir():
        if f.name.endswith(".json"):
            names.add(f.name[:-5])
    return sorted(names)


def register_population(name, cohorts, series, meta=None):
    """Register custom incidences. `series` maps each of SERIES to an array of shape
    (n_cohorts, 80) with annual incidence rates (cases per person-year) at ages 0..79."""
    _CUSTOM[name] = IncidenceTable(name, cohorts, series, meta)
    load_population.cache_clear()


def ethnic_population(name, ethnicity):
    """Name of the incidence table for a UK ethnic group (see :mod:`hbocrisk.ethnicity`),
    e.g. 'UK (Black)'. Only the UK incidences can be combined with an ethnicity."""
    if ALIASES.get(name, name) != "UK":
        raise ValueError(f"{name} cancer rates with a UK ethnicity parameter ({ethnicity}) "
                         "are not valid; ethnicity-specific incidences are only defined "
                         "for the UK")
    return f"UK ({ethnicity})"


@lru_cache(maxsize=None)
def load_population(name="UK"):
    if name.startswith("UK (") and name.endswith(")"):
        ethnicity = name[4:-1]
        base = load_population("UK")
        irr = ETHNICITY_IRR[ethnicity]
        series = {k: v * irr.get(k, 1.0) for k, v in base.series.items()}
        return IncidenceTable(name, base.cohorts, series,
                              dict(base.meta, ethnicity=ethnicity, rate_ratios=irr))
    name = ALIASES.get(name, name)
    if name in _CUSTOM:
        return _CUSTOM[name]
    path = resources.files("hbocrisk").joinpath("data", "incidences", f"{name}.json")
    if not path.is_file():
        raise ValueError(f"unknown cancer incidence population '{name}'; available: "
                         f"{available_populations()}")
    d = json.loads(path.read_text())
    return IncidenceTable(name, d["cohorts"], d["series"],
                          {k: v for k, v in d.items() if k != "series"})
