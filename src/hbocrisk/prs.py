"""Polygenic risk scores: raw score from genotypes/dosages and normalisation.

BOADICEA takes a PRS as a standard-normal z-score together with alpha, the square root
of the proportion of the model's polygenic variance explained by the PRS [LEE2019,
MAV2023]. For the PRS supported in CanRisk the raw score is normalised with fixed,
published population parameters (European ancestry; CanRisk FAQ "What variants are
used in the PRS?", Mavaddat et al. CEBP 2023; Ficorella et al. BJC 2025 for other
ancestries). The variant files in ``data/prs`` are the CanRisk reference files
(effect allele, effect allele frequency, log odds ratio per allele).

Ancestry [FIC2025E]: the raw score is the same for every ancestry, but its mean, SD and
alpha depend on genetic ancestry (European, African, East Asian, South Asian). Files
with ancestry-specific parameters (``*-UKB_all.prs`` from CCGE-BOADICEA/
SHARE-PRScalculation, for BCAC 309, BCAC 307 and OCAC 36) carry them in their header,
together with ancestry-specific effect-allele frequencies. They are distributed there
under GPL-3.0 and therefore not bundled: they are downloaded on first use (pinned
commit, SHA-256 checked) into ``$HBOCRISK_DATA`` (default ``~/.cache/hbocrisk``); see
:func:`ancestry_file` and ``hbocrisk fetch-data``. Individuals whose main
ancestry proportion reaches the threshold (0.8; 0.7 for South Asian) are standardised
with that ancestry's parameters; for mixed genetic ancestry the z-score and alpha are
derived from the ancestry proportions as in vcf2prs (Ficorella et al. 2025, eqs. for
Z, alpha and sigma). Unknown ancestry uses the population parameters.
"""
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
import gzip
import hashlib
import math
import os
import re
import tempfile
import urllib.request

import numpy as np

# name: (bundled file or None, disease, mean, sd, alpha)
REFERENCE_PRS = {
    "BCAC 77": ("BCAC_77_PRS.prs", "BC", -0.892, 0.449, 0.396),
    "BCAC 313": ("BCAC_313_PRS.prs", "BC", -0.424, 0.611, 0.501),
    "BCAC 3820": ("BCAC_3820_PRS.prs", "BC", -0.445, 0.460, 0.516),
    "BRIDGES 306": ("BRIDGES_306_PRS.prs", "BC", -0.422, 0.608, 0.498),
    "EGLH-CEN 301": ("EGLH-CEN_301_PRS.prs", "BC", -0.378, 0.604, 0.496),
    "EGLH-CEN 303": ("EGLH-CEN_303_PRS.prs", "BC", -0.407, 0.606, 0.497),
    "MAINZ 309": ("MAINZ_309_PRS.prs", "BC", -0.351, 0.607, 0.499),
    "PERSPECTIVE 295": ("PERSPECTIVE_295_PRS.prs", "BC", -0.448, 0.599, 0.492),
    "PRISMA 268": ("PRISM_268_PRS.prs", "BC", 0.322, 0.558, 0.451),
    "WISDOM 75": ("WISDOM_75_PRS.prs", "BC", -1.057, 0.567, 0.360),
    "WISDOM 128": ("WISDOM_128_PRS.prs", "BC", -0.180, 0.464, 0.452),
    "OCAC 36": ("OCAC_36_PRS.prs", "OC", -0.259, 0.315, 0.223),
    "OC-EGLH-CEN 34": ("OC_EGLH-CEN_34_PRS.prs", "OC", -0.092, 0.309, 0.219),
    # BCAC 313 without the 4 SNPs not imputed in UK Biobank (309), and additionally
    # without the 2 SNPs correlated with CHEK2 c.1100del (307) [FIC2025E]; with
    # ancestry-specific parameters only (downloaded, see ANCESTRY_PRS).
    "BCAC 309": (None, "BC", -0.339, 0.615, 0.499),
    "BCAC 307": (None, "BC", -0.007, 0.613, 0.497),
}

# PRS with ancestry-specific parameters: (file, sha256) in SHARE-PRScalculation v4.2.12
ANCESTRY_PRS = {
    "BCAC 309": ("BCAC_309_PRS-UKB_all.prs",
                 "fb6ea386ed254d84225d7aa305add5cca7cf204cab625e48ff5c3f33d7cd4ba9"),
    "BCAC 307": ("BCAC_307_PRS-UKB_all.prs",
                 "992dab272a940ed2c9d98314d7b365b6b58be822618176e356aa3a2e5d0be0f4"),
    "OCAC 36": ("OCAC_36_PRS-UKB_all.prs",
                "f2c25085c4475a0cdd335de91813da64372f172753e973917cec6f9e27496bb6"),
}
ANCESTRY_PRS_URL = ("https://raw.githubusercontent.com/CCGE-BOADICEA/SHARE-PRScalculation/"
                    "5edd79feb9f986b6cfb0d4d5616ee3c86c563235/PRSmodels_CanRisk/{}")

ANCESTRIES = ("European", "African", "East Asian", "South Asian")
_ANCESTRY_COLUMNS = ("eaf_eur", "eaf_afr", "eaf_eas", "eaf_sas")
_ANCESTRY_NAMES = {
    "european": "European", "eur": "European", "white": "European",
    "african": "African", "afr": "African", "black": "African",
    "east asian": "East Asian", "eastasian": "East Asian", "eas": "East Asian",
    "south asian": "South Asian", "southasian": "South Asian", "sas": "South Asian",
}
_UNKNOWN_ANCESTRY = ("", "unknown", "pop", "population", "na", "other", "mixed")


def ancestry_name(name):
    """Canonical ancestry name (European, African, East Asian, South Asian) or None for
    unknown/unspecified ancestry. Accepts EUR/AFR/EAS/SAS and similar spellings."""
    if name is None:
        return None
    key = re.sub(r"[\s_-]+", " ", str(name).strip().lower())
    if key in _UNKNOWN_ANCESTRY:
        return None
    if key in _ANCESTRY_NAMES:
        return _ANCESTRY_NAMES[key]
    if key.replace(" ", "") in _ANCESTRY_NAMES:
        return _ANCESTRY_NAMES[key.replace(" ", "")]
    raise ValueError(f"unknown ancestry '{name}'; use one of {', '.join(ANCESTRIES)}, "
                     "'unknown', or a dict of ancestry proportions")


def ancestry_proportions(props):
    """Normalise a mapping ancestry -> proportion to a vector over ANCESTRIES."""
    vec = np.zeros(len(ANCESTRIES))
    for k, v in props.items():
        name = ancestry_name(k)
        if name is None:
            raise ValueError(f"ancestry proportion for unsupported group '{k}'")
        vec[ANCESTRIES.index(name)] += float(v)
    if np.any(vec < 0):
        raise ValueError("negative ancestry proportions")
    if not 0.99 < vec.sum() < 1.01:
        raise ValueError(f"ancestry proportions must add up to 1 (sum {vec.sum():.3f})")
    return vec / vec.sum()


@dataclass
class AncestryParameters:
    threshold: float     # minimum proportion of the main ancestry for single ancestry
    mean: float
    sd: float
    alpha: float


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    effect: str
    eaf: float
    beta: float
    eaf_ancestry: tuple = None   # effect-allele frequencies in ANCESTRIES order


@dataclass
class PRSModel:
    name: str
    variants: list
    alpha: float
    mean: float = None
    sd: float = None
    disease: str = "BC"
    ancestry: dict = None        # ancestry name -> AncestryParameters

    def __post_init__(self):
        p = np.array([v.eaf for v in self.variants])
        b = np.array([v.beta for v in self.variants])
        if self.mean is None:
            self.mean = float((2 * p * b).sum())       # expected raw score under HWE
        if self.sd is None:
            self.sd = float(np.sqrt((2 * p * (1 - p) * b * b).sum()))

    def raw_score(self, dosages, impute_missing=True, ancestry=None):
        """dosages: mapping from variant key to effect-allele dosage (0..2). Keys may be
        (chrom, pos) tuples or 'chrom:pos' / 'chrom:pos:ref:effect' strings. Missing
        variants contribute their expected value 2*eaf*beta if `impute_missing`, with
        the ancestry-specific (or proportion-weighted) frequency if `ancestry` is given
        and the model has ancestry-specific frequencies. Returns (raw_score, n_found)."""
        eafs = self._imputation_eafs(ancestry)
        norm = {}
        for k, v in dosages.items():
            if isinstance(k, tuple):
                norm[(_chrom(k[0]), int(k[1]))] = float(v)
            else:
                parts = str(k).split(":")
                norm[(_chrom(parts[0]), int(parts[1]))] = float(v)
        s, n = 0.0, 0
        for var, eaf in zip(self.variants, eafs):
            d = norm.get((var.chrom, var.pos))
            if d is None or (isinstance(d, float) and math.isnan(d)):
                if impute_missing:
                    s += 2 * eaf * var.beta
                continue
            s += d * var.beta
            n += 1
        return s, n

    def _imputation_eafs(self, ancestry):
        pop = [v.eaf for v in self.variants]
        if ancestry is None or not self.ancestry:
            return pop
        if isinstance(ancestry, dict):
            w = ancestry_proportions(ancestry)
        else:
            name = ancestry_name(ancestry)
            if name is None:
                return pop
            w = np.eye(len(ANCESTRIES))[ANCESTRIES.index(name)]
        return [float(np.dot(v.eaf_ancestry, w)) for v in self.variants]

    def zscore(self, raw):
        return (raw - self.mean) / self.sd

    def normalise(self, raw, ancestry=None):
        """alpha and z-score for a raw score given the genetic ancestry.

        ancestry: None/'unknown' (population parameters), an ancestry name, or a dict of
        ancestry proportions (e.g. {'EUR': 0.6, 'AFR': 0.4}). Returns a dict with alpha,
        zscore, the ancestry used and the method ('population', 'single', 'mixed')."""
        if isinstance(ancestry, dict):
            if not self.ancestry:
                raise ValueError(f"PRS '{self.name}' has no ancestry-specific parameters")
            props = ancestry_proportions(ancestry)
            params = [self.ancestry[a] for a in ANCESTRIES]
            single = [a for a, p, q in zip(ANCESTRIES, props, params) if p >= q.threshold]
            if single:
                ancestry = single[0]
            else:
                return self._mixed(raw, props, params)
        name = ancestry_name(ancestry)
        if name is None:
            return {"alpha": self.alpha, "zscore": self.zscore(raw), "ancestry": "unknown",
                    "method": "population"}
        if not self.ancestry:
            if name == "European":     # CanRisk reference PRS are European parameters
                return {"alpha": self.alpha, "zscore": self.zscore(raw),
                        "ancestry": name, "method": "single"}
            raise ValueError(f"PRS '{self.name}' has no parameters for {name} ancestry; "
                             f"use one with ancestry-specific parameters "
                             f"({', '.join(n for n in REFERENCE_PRS if _has_ancestry(n))})")
        q = self.ancestry[name]
        return {"alpha": q.alpha, "zscore": (raw - q.mean) / q.sd, "ancestry": name,
                "method": "single"}

    def _mixed(self, raw, props, params):
        mean = np.array([q.mean for q in params])
        sd = np.array([q.sd for q in params])
        alpha = np.array([q.alpha for q in params])
        eaf = np.array([v.eaf_ancestry for v in self.variants])
        beta = np.array([v.beta for v in self.variants])
        # variance of the raw score in the admixed individual: within-ancestry variance
        # plus the variance of each SNP's frequency over the ancestry mixture
        var = (props * sd ** 2).sum() + 4 * (((eaf ** 2) @ props - (eaf @ props) ** 2)
                                            * beta ** 2).sum()
        sigma = math.sqrt(var)
        z = (alpha * props * (raw - mean) / sd).sum() / (alpha * props / sd).sum() / sigma
        a = (alpha * props * sd).sum() / sigma
        return {"alpha": float(a), "zscore": float(z),
                "ancestry": {n: float(p) for n, p in zip(ANCESTRIES, props)},
                "method": "mixed"}


def _chrom(c):
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


_ANCESTRY_HEADER = re.compile(r"^([-+0-9.eE,\s]+)#\s*(threshold|mean|sd|alpha)\s*$", re.I)


def read_prs_file(path_or_text, name=None, mean=None, sd=None, disease="BC"):
    """Read a CanRisk .prs reference file ("alpha = x" header line followed by CSV with
    Chromosome,Position,Reference_Allele,Effect_Allele,Effect_Allele_Frequency,
    Log_Odds_Ratio).

    Multi-ancestry files (vcf2prs v4) instead have header lines
    ``Pop,Eur,Afr,EAs,SAs  #threshold`` (and #mean, #sd, #alpha) and the columns
    eaf_eur, eaf_afr, eaf_eas and eaf_sas; their population values are used as the
    default normalisation, and the ancestry-specific ones are kept in ``ancestry``."""
    text = path_or_text
    if "\n" not in str(path_or_text):
        with open(path_or_text) as fh:
            text = fh.read()
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    alpha = None
    multi = {}
    variants = []
    header = None
    for ln in lines:
        m = _ANCESTRY_HEADER.match(ln)
        if m:
            multi[m.group(2).lower()] = [float(x) for x in m.group(1).split(",") if x.strip()]
            continue
        if ln.startswith("#"):
            continue
        m = re.match(r"^\s*alpha\s*=\s*([0-9.eE+-]+)", ln, re.I)
        if m:
            alpha = float(m.group(1))
            continue
        cols = [c.strip() for c in ln.split(",")]
        if header is None and any(c.lower() == "position" for c in cols):
            header = [c.lower() for c in cols]
            continue
        if header is None:
            continue
        d = dict(zip(header, cols))
        eaf_anc = None
        if all(c in d for c in _ANCESTRY_COLUMNS):
            eaf_anc = tuple(float(d[c]) for c in _ANCESTRY_COLUMNS)
        variants.append(Variant(_chrom(d["chromosome"]), int(d["position"]),
                                d.get("reference_allele", ""), d["effect_allele"],
                                float(d["effect_allele_frequency"]),
                                float(d["log_odds_ratio"]), eaf_anc))
    ancestry = None
    if multi:
        if set(multi) != {"threshold", "mean", "sd", "alpha"} or \
                any(len(v) != 1 + len(ANCESTRIES) for v in multi.values()):
            raise ValueError("ancestry header must give threshold, mean, sd and alpha for "
                             "Pop,Eur,Afr,EAs,SAs")
        if any(v.eaf_ancestry is None for v in variants):
            raise ValueError("multi-ancestry PRS file without eaf_eur/eaf_afr/eaf_eas/eaf_sas")
        ancestry = {a: AncestryParameters(*(multi[k][i + 1] for k in
                                            ("threshold", "mean", "sd", "alpha")))
                    for i, a in enumerate(ANCESTRIES)}
        mean = multi["mean"][0] if mean is None else mean
        sd = multi["sd"][0] if sd is None else sd
        alpha = multi["alpha"][0] if alpha is None else alpha
    return PRSModel(name or "custom", variants, alpha, mean, sd, disease, ancestry)


def data_dir():
    """Directory for downloaded data: $HBOCRISK_DATA or ~/.cache/hbocrisk."""
    d = os.environ.get("HBOCRISK_DATA")
    if d:
        return Path(d)
    return Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "hbocrisk"


def ancestry_file(name, download=True):
    """Path of the ancestry-specific file of a reference PRS, downloaded if missing."""
    f, sha256 = ANCESTRY_PRS[name]
    path = data_dir() / f
    if not path.exists():
        if not download:
            raise FileNotFoundError(f"{path} missing; run 'hbocrisk fetch-data'")
        url = ANCESTRY_PRS_URL.format(f)
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = r.read()
        except OSError as e:
            raise RuntimeError(f"cannot download {url} ({e}); download it into "
                               f"{path.parent} or set HBOCRISK_DATA") from e
        if hashlib.sha256(data).hexdigest() != sha256:
            raise RuntimeError(f"checksum mismatch for {url}")
        path.parent.mkdir(parents=True, exist_ok=True)
        # atomic, so that concurrent processes never read a partial file
        fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f + ".")
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        os.replace(tmp, path)
    elif hashlib.sha256(path.read_bytes()).hexdigest() != sha256:
        raise RuntimeError(f"checksum mismatch for {path}; delete it to download again")
    return path


def reference_model(name, ancestry=False):
    """Load one of the PRS supported by CanRisk with its published normalisation.

    With ``ancestry=True`` (always for BCAC 309/307) the file with ancestry-specific
    parameters and allele frequencies is used (downloaded on first use)."""
    f, disease, mean, sd, alpha = REFERENCE_PRS[name]
    if f is None or (ancestry and name in ANCESTRY_PRS):
        text = ancestry_file(name).read_text()
    else:
        text = resources.files("hbocrisk").joinpath("data", "prs", f).read_text()
    model = read_prs_file(text, name=name, mean=mean, sd=sd, disease=disease)
    model.alpha = alpha
    return model


def _has_ancestry(name):
    return name in ANCESTRY_PRS


def ancestry_parameters(name):
    """Ancestry-specific (threshold, mean, sd, alpha) of a reference PRS."""
    return reference_model(name, ancestry=True).ancestry


def vcf_dosages(path, sample=None, model=None):
    """Effect-allele dosages from a VCF (GT or DS fields) for one sample.

    If a PRSModel is given, dosages are oriented to the model's effect alleles."""
    opener = gzip.open if str(path).endswith(".gz") else open
    effects = {}
    if model is not None:
        effects = {(v.chrom, v.pos): v for v in model.variants}
    out = {}
    with opener(path, "rt") as fh:
        samples = None
        for ln in fh:
            if ln.startswith("##"):
                continue
            f = ln.rstrip("\n").split("\t")
            if ln.startswith("#"):
                samples = f[9:]
                col = 9 + (samples.index(sample) if sample is not None else 0)
                continue
            chrom, pos, ref, alt = _chrom(f[0]), int(f[1]), f[3], f[4].split(",")[0]
            fmt = f[8].split(":")
            vals = dict(zip(fmt, f[col].split(":")))
            dose = None
            if "DS" in vals and vals["DS"] not in (".", ""):
                dose = float(vals["DS"])
            elif "GT" in vals:
                al = re.split(r"[/|]", vals["GT"])
                if "." not in al:
                    dose = float(sum(1 for a in al if a != "0"))
            if dose is None:
                continue
            key = (chrom, pos)
            var = effects.get(key)
            if var is not None:
                if var.effect == alt:
                    pass
                elif var.effect == ref:
                    dose = 2.0 - dose
                else:
                    continue                 # allele mismatch
            out[key] = dose
    return out


def prs_from_vcf(path, name="BCAC 313", sample=None, ancestry=None):
    """Return dict(alpha, zscore, raw, n_found, n_variants, ancestry, method) for a
    sample. `ancestry` as in :meth:`PRSModel.normalise`."""
    if name in REFERENCE_PRS:
        needs = isinstance(ancestry, dict) or ancestry_name(ancestry) is not None
        model = reference_model(name, ancestry=needs)
    else:
        model = read_prs_file(name)
    dos = vcf_dosages(path, sample, model)
    raw, n = model.raw_score(dos, ancestry=ancestry)
    norm = model.normalise(raw, ancestry)
    return {"name": model.name, "alpha": norm["alpha"], "zscore": norm["zscore"],
            "raw": raw, "n_found": n, "n_variants": len(model.variants),
            "ancestry": norm["ancestry"], "method": norm["method"]}
