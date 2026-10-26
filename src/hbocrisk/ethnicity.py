"""UK ethnic groups of BOADICEA BC v7 / EOC v3 [FIC2025E].

Self-reported ethnicity selects the population cancer incidences (UK incidences scaled
by ethnicity-specific rate ratios) and the CHEK2 pathogenic-variant frequency. CanRisk
offers this only together with the UK incidences. Genetic ancestry for the PRS is a
separate input (see :mod:`hbocrisk.prs`).

Not implemented: the ethnicity-specific distributions of questionnaire risk factors and
BI-RADS density (Tables S18, S21) and breast tumour subtypes (Tables S4-S6); the
default (v6) distributions are used for every group.
"""
from .params import ETHNICITY_IRR

GROUPS = tuple(ETHNICITY_IRR)

_NAMES = {
    "white": "White", "black": "Black", "south asian": "South Asian",
    "asian": "South Asian", "east asian": "East Asian", "chinese": "East Asian",
    "mixed": "Mixed", "other": "Other",
}
# ONS 2011/2021 census groups (CanRisk "##ethnicity=Group;Background")
_ONS = {
    "white": "White", "mixed/multiple ethnic groups": "Mixed",
    "asian or asian british": "South Asian", "black or black british": "Black",
    "other ethnic group": "Other",
}
_UNSPECIFIED = ("", "unknown", "unspecified", "na", "none", "not specified")


def ethnic_group(value):
    """UK ethnic group (White, Black, South Asian, East Asian, Mixed, Other) or None
    for unspecified ethnicity. Accepts these names, UK Biobank names (asian, chinese)
    and CanRisk ONS strings such as 'Asian or Asian British;Chinese'."""
    if value is None:
        return None
    group, _, background = str(value).partition(";")
    group, background = group.strip().lower(), background.strip().lower()
    if group in _UNSPECIFIED:
        return None
    if group in _ONS:
        if group == "asian or asian british" and background == "chinese":
            return "East Asian"
        return _ONS[group]
    if group in _NAMES:
        return _NAMES[group]
    raise ValueError(f"unknown ethnicity '{value}'; use one of {', '.join(GROUPS)} or "
                     "an ONS group (e.g. 'Black or Black British;African')")
