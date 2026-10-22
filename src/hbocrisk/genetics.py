"""Major-gene locus and hypergeometric polygenic model (HPM).

Major genes
-----------
Pathogenic variants (PVs) in the model genes are represented as alleles of a single
"major gene" locus: allele 0 is wild type and allele i (1..G) a PV in the i-th gene in
dominance order [LEE2016]. A genotype is an unordered pair of alleles. The risk class of
a genotype is the highest-ranked (lowest index) PV it carries, or 0 for non-carriers.

Polygenic component
-------------------
The polygene is approximated by the hypergeometric polygenic model of Cannings et al.
(1978) / Lange (1997) with 2N alleles (N = 3): the polygenotype R (number of "high
risk" alleles) is Binomial(2N, 1/2) in founders and each parent transmits N of their 2N
alleles sampled without replacement. The standardised polygenotype is
x = (R - N) / sqrt(N/2), with mean 0 and variance 1.
"""
from functools import lru_cache
import math

import numpy as np
from scipy.stats import norm, binom


class MajorLocus:
    def __init__(self, genes):
        self.genes = list(genes)
        self.nalleles = len(genes) + 1
        pairs = [(a, b) for a in range(self.nalleles) for b in range(a, self.nalleles)]
        self.pairs = pairs
        self.ngeno = len(pairs)
        self.index = {p: i for i, p in enumerate(pairs)}
        cls = []
        for a, b in pairs:
            nz = [x for x in (a, b) if x > 0]
            cls.append(min(nz) if nz else 0)
        self.geno_class = np.array(cls)                 # (ngeno,)
        self.nclass = self.nalleles
        # indicator: genotype carries allele of gene g (g = 1..G)
        self.carries = np.zeros((self.nalleles, self.ngeno), bool)
        for i, (a, b) in enumerate(pairs):
            self.carries[a, i] = True
            self.carries[b, i] = True
        self.class_matrix = np.zeros((self.ngeno, self.nclass))
        self.class_matrix[np.arange(self.ngeno), self.geno_class] = 1.0
        self._transmission = None

    def allele_freqs(self, freq_by_gene):
        q = np.array([float(freq_by_gene.get(g, 0.0)) for g in self.genes])
        return np.concatenate([[1.0 - q.sum()], q])

    def genotype_prior(self, freq_by_gene):
        """Hardy-Weinberg genotype frequencies."""
        f = self.allele_freqs(freq_by_gene)
        out = np.empty(self.ngeno)
        for i, (a, b) in enumerate(self.pairs):
            out[i] = f[a] * f[b] * (1 if a == b else 2)
        return out

    def class_prior(self, freq_by_gene):
        return self.genotype_prior(freq_by_gene) @ self.class_matrix

    def transmission(self):
        """T[c, m, f] = P(child genotype c | mother m, father f)."""
        if self._transmission is None:
            n = self.ngeno
            T = np.zeros((n, n, n))
            for m, (m1, m2) in enumerate(self.pairs):
                for f, (f1, f2) in enumerate(self.pairs):
                    for x in (m1, m2):
                        for y in (f1, f2):
                            c = self.index[(min(x, y), max(x, y))]
                            T[c, m, f] += 0.25
            self._transmission = T
        return self._transmission


class HPM:
    def __init__(self, n=3):
        self.n = n
        self.size = 2 * n + 1
        r = np.arange(self.size)
        self.values = (r - n) / math.sqrt(n / 2.0)           # standardised polygenotype
        self.prior = binom.pmf(r, 2 * n, 0.5)
        self._transmission = None

    def transmission(self):
        """T[c, m, f] = P(child polygenotype c | mother m, father f)."""
        if self._transmission is None:
            n2 = 2 * self.n
            s = self.size
            # gamete distribution: number of high alleles transmitted j by parent with R
            g = np.zeros((s, self.n + 1))
            for R in range(s):
                for j in range(self.n + 1):
                    if j <= R and self.n - j <= n2 - R:
                        g[R, j] = math.comb(R, j) * math.comb(n2 - R, self.n - j) / math.comb(n2, self.n)
            T = np.zeros((s, s, s))
            for m in range(s):
                for f in range(s):
                    for jm in range(self.n + 1):
                        for jf in range(self.n + 1):
                            T[jm + jf, m, f] += g[m, jm] * g[f, jf]
            self._transmission = T
        return self._transmission

    def bin_edges(self):
        """Standard normal quantile bins matching the binomial masses [LEE2019] (s1.13)."""
        cdf = np.cumsum(self.prior)
        cdf[-1] = 1.0
        edges = np.concatenate([[-np.inf], norm.ppf(cdf[:-1]), [np.inf]])
        return edges

    def prs_likelihood(self, zscore, alpha):
        """f(x_PRS | x_P in bin k) for each polygenotype bin k [LEE2019] (s1.14).

        x_PRS = alpha * zscore is the observed PRS with variance alpha^2 and the total
        polygene x_P ~ N(0, 1). Returns an array over the 2N+1 bins; only ratios
        matter (constant factors cancel in pedigree likelihood ratios).
        """
        alpha = float(alpha)
        if alpha <= 0:
            return np.ones(self.size)
        alpha = min(alpha, 0.999999)
        x = alpha * float(zscore)
        edges = self.bin_edges()
        s = math.sqrt(1.0 - alpha * alpha)
        # joint density integrated over bin: phi_alpha(x) * P(x_P in bin | x_PRS = x)
        dens = norm.pdf(x, scale=alpha) * (norm.cdf((edges[1:] - x) / s) - norm.cdf((edges[:-1] - x) / s))
        return dens / self.prior


@lru_cache(maxsize=None)
def hpm(n=3):
    return HPM(n)


@lru_cache(maxsize=None)
def major_locus(genes):
    return MajorLocus(genes)
