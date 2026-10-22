"""Polygenic risk scores: raw score from genotypes/dosages and normalisation.

BOADICEA takes a PRS as a standard-normal z-score together with alpha, the square root
of the proportion of the model's polygenic variance explained by the PRS [LEE2019,
MAV2023]. For the PRS supported in CanRisk the raw score is normalised with fixed,
published population parameters (European ancestry; CanRisk FAQ "What variants are
used in the PRS?", Mavaddat et al. CEBP 2023; Ficorella et al. BJC 2025 for other
ancestries). The variant files in ``data/prs`` are the CanRisk reference files
(effect allele, effect allele frequency, log odds ratio per allele).
"""
from dataclasses import dataclass
from importlib import resources
import gzip
import math
import re

import numpy as np

# name: (file, disease, mean, sd, alpha)
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
}

# Ancestry-specific parameters (Ficorella et al. BJC 2025, Table 2): mean, sd, alpha.
# 309-SNP and 307-SNP versions of BCAC 313 (UK Biobank), and OCAC 36.
ANCESTRY_PRS = {
    "BCAC 309": {"European": (-0.335, 0.615, 0.499), "African": (0.292, 0.547, 0.196),
                 "East Asian": (0.268, 0.555, 0.326), "South Asian": (-0.160, 0.602, 0.331)},
    "BCAC 307": {"European": (-0.005, 0.613, 0.497), "African": (0.634, 0.547, 0.196),
                 "East Asian": (0.611, 0.555, 0.326), "South Asian": (0.179, 0.601, 0.329)},
    "OCAC 36": {"European": (0.007, 0.322, 0.223), "African": (-0.152, 0.254, 0.178),
                "East Asian": (-0.161, 0.235, 0.171), "South Asian": (-0.180, 0.278, 0.194)},
}


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    effect: str
    eaf: float
    beta: float


@dataclass
class PRSModel:
    name: str
    variants: list
    alpha: float
    mean: float = None
    sd: float = None
    disease: str = "BC"

    def __post_init__(self):
        p = np.array([v.eaf for v in self.variants])
        b = np.array([v.beta for v in self.variants])
        if self.mean is None:
            self.mean = float((2 * p * b).sum())       # expected raw score under HWE
        if self.sd is None:
            self.sd = float(np.sqrt((2 * p * (1 - p) * b * b).sum()))

    def raw_score(self, dosages, impute_missing=True):
        """dosages: mapping from variant key to effect-allele dosage (0..2). Keys may be
        (chrom, pos) tuples or 'chrom:pos' / 'chrom:pos:ref:effect' strings. Missing
        variants contribute their expected value 2*eaf*beta if `impute_missing`.
        Returns (raw_score, n_found)."""
        norm = {}
        for k, v in dosages.items():
            if isinstance(k, tuple):
                norm[(_chrom(k[0]), int(k[1]))] = float(v)
            else:
                parts = str(k).split(":")
                norm[(_chrom(parts[0]), int(parts[1]))] = float(v)
        s, n = 0.0, 0
        for var in self.variants:
            d = norm.get((var.chrom, var.pos))
            if d is None or (isinstance(d, float) and math.isnan(d)):
                if impute_missing:
                    s += 2 * var.eaf * var.beta
                continue
            s += d * var.beta
            n += 1
        return s, n

    def zscore(self, raw):
        return (raw - self.mean) / self.sd


def _chrom(c):
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def read_prs_file(path_or_text, name=None, mean=None, sd=None, disease="BC"):
    """Read a CanRisk .prs reference file ("alpha = x" header line followed by CSV with
    Chromosome,Position,Reference_Allele,Effect_Allele,Effect_Allele_Frequency,
    Log_Odds_Ratio)."""
    text = path_or_text
    if "\n" not in str(path_or_text):
        with open(path_or_text) as fh:
            text = fh.read()
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    alpha = None
    variants = []
    header = None
    for ln in lines:
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
        variants.append(Variant(_chrom(d["chromosome"]), int(d["position"]),
                                d.get("reference_allele", ""), d["effect_allele"],
                                float(d["effect_allele_frequency"]),
                                float(d["log_odds_ratio"])))
    return PRSModel(name or "custom", variants, alpha, mean, sd, disease)


def reference_model(name):
    """Load one of the PRS supported by CanRisk with its published normalisation."""
    f, disease, mean, sd, alpha = REFERENCE_PRS[name]
    text = resources.files("openboadicea").joinpath("data", "prs", f).read_text()
    model = read_prs_file(text, name=name, mean=mean, sd=sd, disease=disease)
    model.alpha = alpha
    return model


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


def prs_from_vcf(path, name="BCAC 313", sample=None):
    """Return dict(alpha, zscore, raw, n_variants_found, n_variants) for a sample."""
    model = reference_model(name) if name in REFERENCE_PRS else read_prs_file(name)
    dos = vcf_dosages(path, sample, model)
    raw, n = model.raw_score(dos)
    return {"name": model.name, "alpha": model.alpha, "zscore": model.zscore(raw),
            "raw": raw, "n_found": n, "n_variants": len(model.variants)}
