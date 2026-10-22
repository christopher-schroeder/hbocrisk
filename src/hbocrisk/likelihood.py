"""Penetrance functions and pedigree likelihood (peeling by variable elimination).

Each individual has two hidden genetic variables: the major-gene genotype G (unordered
allele pair at the single "major gene" locus) and the polygenotype R of the
hypergeometric polygenic model. They are transmitted independently. The phenotype
likelihood ("penetrance") depends on (G, R) through the genotype's risk class.

The pedigree likelihood is computed by sum-product variable elimination with a greedy
minimum-size elimination order, which reduces to Elston-Stewart peeling for pedigrees
without loops and also handles loops (at higher cost).
"""
import string

import numpy as np

from .genetics import hpm
from .params import MAX_AGE

_LETTERS = string.ascii_letters


# --------------------------------------------------------------------------------------
# penetrance
# --------------------------------------------------------------------------------------

def _cum(inc):
    """cumulative hazard before age a: cum[..., a] = sum_{u<a} inc[..., u], a=0..80."""
    c = np.zeros(inc.shape[:-1] + (MAX_AGE + 1,))
    c[..., 1:] = np.cumsum(inc, axis=-1)
    return c


def log_event(inc, cum, age):
    """log P(survive to age `age` then event within year `age`)."""
    return -cum[..., age] + np.log(-np.expm1(-inc[..., age]))


def log_surv(cum, age):
    return -cum[..., age]


def censor_age(person):
    """Age at last follow-up censored at 80 (diagnoses at 80 or later are ignored)."""
    return min(int(person.age), MAX_AGE)


def phenotype_loglik(person, incs, pathology=None, is_proband_inc=None):
    """Log-likelihood over (class, polygenotype) of a person's cancer phenotypes.

    incs: dict cancer -> (C, K, T) incidences applicable to this person (by sex).
    pathology: (p_erneg, p_tn, ck) class-specific tumour proportions (females).
    """
    A = censor_age(person)
    ll = 0.0
    for cancer, inc in incs.items():
        cum = _cum(inc)
        if cancer == "bc":
            d1 = person.cancer_age("bc1")
            if d1 is not None and d1 >= MAX_AGE:
                d1 = None
            if d1 is None:
                ll = ll + log_surv(cum, A)
                continue
            ll = ll + log_event(inc, cum, d1)
            if pathology is not None and person.sex == "F":
                ll = ll + np.log(pathology_factor(person.pathology_code(), pathology, d1))[:, None]
            # contralateral breast cancer: half the breast incidence (one breast at risk)
            half = 0.5 * inc
            hc = _cum(half)
            d2 = person.cancer_age("bc2")
            if d2 is not None and d2 >= MAX_AGE:
                d2 = None
            if d2 is not None:
                d2 = max(d2, d1)
                if d2 == d1:
                    ll = ll + np.log(-np.expm1(-half[..., d1]))
                else:
                    ll = ll + (-(hc[..., d2] - hc[..., d1 + 1])
                               + np.log(-np.expm1(-half[..., d2])))
            elif A > d1 + 1:
                ll = ll - (hc[..., A] - hc[..., d1 + 1])
        else:
            col = {"oc": "oc", "prc": "prc", "pac": "pac"}[cancer]
            d = person.cancer_age(col)
            if d is not None and d < MAX_AGE:
                ll = ll + log_event(inc, cum, d)
            else:
                ll = ll + log_surv(cum, A)
    return ll


def pathology_factor(code, path, age):
    """P(observed tumour pathology | class) at age of diagnosis, shape (C,)."""
    p_er, p_tn, ck = path
    p = p_er[:, age]
    q = p_tn[:, age]
    c = ck[:, age, :]
    if code == 0:
        return 1 - p
    if code == 1:
        return p
    if code == 2:
        return p * (1 - q)
    if code == 3:
        return p * q
    if code == 4:
        return p * q * c[:, 2]
    if code == 5:
        return p * q * c[:, 1]
    if code == 6:
        return p * q * c[:, 0]
    return np.ones_like(p)


def genetic_test_factor(person, locus, sensitivity):
    """Likelihood of genetic test results over genotypes, shape (ngeno,).

    Tests are assumed 100% specific. A negative mutation search ('S') in gene g has
    probability 1 - sensitivity_g for genotypes carrying a PV in g; a negative direct
    test ('T') excludes them. A positive test identifies the risk class as the
    highest-ranked gene with a positive result; positive results in lower ranked
    genes are ignored (dominance model, CanRisk FAQ).
    """
    f = np.ones(locus.ngeno)
    positives = []
    for gi, gene in enumerate(locus.genes, start=1):
        t = person.tests.get(gene)
        if t is None or not t.tested:
            continue
        if t.result == "P":
            positives.append(gi)
        elif t.result == "N":
            if t.type == "T":
                f[locus.carries[gi]] = 0.0
            else:
                f[locus.carries[gi]] *= 1.0 - float(sensitivity.get(gene, 1.0))
    if positives:
        f = f * (locus.geno_class == min(positives))
    return f


# --------------------------------------------------------------------------------------
# variable elimination
# --------------------------------------------------------------------------------------

class Factor:
    __slots__ = ("vars", "table", "logscale")

    def __init__(self, vars_, table, logscale=0.0):
        self.vars = tuple(vars_)
        self.table = table
        self.logscale = logscale


def _normalise(f):
    m = f.table.max()
    if m <= 0 or not np.isfinite(m):
        if m <= 0:
            raise ZeroDivisionError("pedigree data are inconsistent (zero likelihood)")
        raise FloatingPointError("non-finite likelihood")
    f.table = f.table / m
    f.logscale += np.log(m)
    return f


def eliminate(factors, sizes, keep):
    """Sum out all variables not in `keep`. Returns (remaining factors, logscale)."""
    store = {i: f for i, f in enumerate(factors)}
    next_id = len(factors)
    index = {}
    for i, f in store.items():
        for u in f.vars:
            index.setdefault(u, set()).add(i)
    to_elim = set(index) - set(keep)
    logscale = 0.0

    def cost(v):
        scope = set()
        for i in index[v]:
            scope.update(store[i].vars)
        scope.discard(v)
        c = 1
        for u in scope:
            c *= sizes[u]
        return c

    while to_elim:
        v = min(to_elim, key=lambda u: (cost(u), u))
        to_elim.discard(v)
        ids = sorted(index.pop(v))
        fs = [store.pop(i) for i in ids]
        for f in fs:
            for u in f.vars:
                if u != v:
                    index[u].difference_update(ids)
        scope = []
        for f in fs:
            for u in f.vars:
                if u != v and u not in scope:
                    scope.append(u)
        letters = {}
        for f in fs:
            for u in f.vars:
                if u not in letters:
                    letters[u] = _LETTERS[len(letters)]
        sub = ",".join("".join(letters[u] for u in f.vars) for f in fs)
        out = "".join(letters[u] for u in scope)
        if len(fs) == 1:
            table = np.einsum(sub + "->" + out, fs[0].table)
        else:
            table = np.einsum(sub + "->" + out, *[f.table for f in fs], optimize="greedy")
        nf = _normalise(Factor(scope, np.asarray(table), sum(f.logscale for f in fs)))
        if not scope:
            logscale += nf.logscale + np.log(float(nf.table))
        else:
            store[next_id] = nf
            for u in scope:
                index[u].add(next_id)
            next_id += 1
    return list(store.values()), logscale


class PedigreeLikelihood:
    """Builds the factor graph of a pedigree for a given model context."""

    def __init__(self, pedigree, locus, geno_prior, penetrance):
        """penetrance: dict pid -> array (ngeno, K) (missing = no phenotype data)."""
        H = hpm()
        self.locus = locus
        self.K = H.size
        people = pedigree.people
        by_id = {p.pid: p for p in people}
        # merge MZ twins into one genetic unit
        unit = {}
        twin_rep = {}
        for p in people:
            if p.mztwin not in ("0", "", None):
                key = (p.mztwin, p.father, p.mother)
                if key in twin_rep:
                    unit[p.pid] = unit[twin_rep[key]]
                    continue
                twin_rep[key] = p.pid
            unit[p.pid] = len(set(unit.values()))
        self.unit = unit
        nunits = len(set(unit.values()))
        rep = {}
        for p in people:
            rep.setdefault(unit[p.pid], p)
        # data per unit
        pen = {}
        for pid, arr in penetrance.items():
            u = unit[pid]
            pen[u] = arr if u not in pen else pen[u] * arr
        proband_unit = unit[pedigree.proband.pid]
        # prune units without data and without (kept) descendants
        children = {u: set() for u in range(nunits)}
        parents = {}
        for u, p in rep.items():
            if not p.is_founder:
                parents[u] = (unit[p.mother], unit[p.father])
                children[unit[p.mother]].add(u)
                children[unit[p.father]].add(u)
        alive = set(range(nunits))
        changed = True
        while changed:
            changed = False
            for u in list(alive):
                if u == proband_unit or u in pen:
                    continue
                if not (children[u] & alive):
                    alive.discard(u)
                    changed = True
        # variables: G_u = 2u, P_u = 2u + 1
        self.sizes = {}
        factors = []
        Tg = locus.transmission()
        Tp = H.transmission()
        for u in alive:
            g, r = 2 * u, 2 * u + 1
            self.sizes[g] = locus.ngeno
            self.sizes[r] = self.K
            if u in parents and parents[u][0] in alive and parents[u][1] in alive:
                m, f = parents[u]
                factors.append(Factor((g, 2 * m, 2 * f), Tg))
                factors.append(Factor((r, 2 * m + 1, 2 * f + 1), Tp))
            else:
                factors.append(Factor((g,), geno_prior))
                factors.append(Factor((r,), H.prior))
            if u in pen:
                factors.append(_normalise(Factor((g, r), pen[u])))
        self.factors = factors
        self.proband_vars = (2 * proband_unit, 2 * proband_unit + 1)

    def proband_posterior(self):
        """Posterior distribution over (genotype, polygenotype) of the proband and the
        log-likelihood of the pedigree."""
        remaining, logscale = eliminate(self.factors, self.sizes, self.proband_vars)
        g, r = self.proband_vars
        post = np.ones((self.sizes[g], self.sizes[r]))
        for f in remaining:
            if f.vars == (g, r):
                post = post * f.table
            elif f.vars == (r, g):
                post = post * f.table.T
            elif f.vars == (g,):
                post = post * f.table[:, None]
            elif f.vars == (r,):
                post = post * f.table[None, :]
            else:
                raise RuntimeError(f"unexpected factor scope {f.vars}")
            logscale += f.logscale
        tot = post.sum()
        if not tot > 0:
            raise ZeroDivisionError("pedigree data are inconsistent (zero likelihood)")
        return post / tot, logscale + np.log(tot)
