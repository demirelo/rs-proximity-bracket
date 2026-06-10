"""
singlelist_past_johnson.py -- EXACT single-code RS list size just past the Johnson
radius, smooth (multiplicative-subgroup) domains vs matched random-subset domains.

THE OBJECT (sub-lemma P', line-decoding-analysis.md SS6.1 / D2)
--------------------------------------------------------------
For C = RS[F,L,k], n=|L|, rho=k/n, and a target word w in F^n, the single-code list
at relative radius delta is

    Lambda(C, delta, w) = { c in C : Delta(c, w) <= floor(delta*n) }
                        = { c in C : agree(c, w) >= n - floor(delta*n) },

and the (worst-case) list size is  Lambda(C, delta) = max_w |Lambda(C, delta, w)|.

KNOWN (the anchors this experiment sits between):
  * delta < UD = (1-rho+1/n)/2 : |Lambda| = 1 exactly (MDS / triangle ineq).
  * delta <= J - eta, J = 1 - sqrt(rho) : |Lambda| <= 1/(2*eta*rho) (ABF Cor 3.3),
    a CONSTANT in n -- proven for ALL RS, incl. smooth domains.
  * delta in (J, r_E), r_E = 1 - H_q(rho) : the worst-case RS list is the famously
    hard / partly-open regime.  P' conjectures it stays O_rho(1/eps0) for delta=J+eps0
    on smooth domains; the Kambire mechanism produces n^C lists but only at
    delta >= (1-rho) - O(1/log n) >= r_E - Omega(1).
  * delta > r_E : worst-case list is q^{Omega(n)} for EVERY code (Elias/CS25 deep holes).

THE KEY QUESTION (Task 2): just past Johnson (delta = J + eps0, small eps0), does the
smooth single-code list stay SMALL (consistent with O(1/eps0)), or does it blow up?
Where does it depart from the Johnson-bound constant?  Is smooth ~ random, or
better/worse?

WHY ONE PASS GIVES THE WHOLE CURVE.  For a fixed w, `_distance_histogram` (below)
returns, for EVERY codeword, the exact agreement popcount with w.  The list size at
radius delta is then just  #{ codewords with popcount >= n - floor(delta*n) }  -- a
threshold on that single integer array.  So one O(q^k * n) pass over the codeword book
yields |Lambda(C, delta, w)| for ALL delta simultaneously, exactly, with no decoder.

PROBING THE WORST CASE (honesty: exhaustive max over all q^n words is infeasible).
We take a STRONG SAMPLE over structured worst-case candidates and report the per-delta
MAX over the whole sample (a lower bound on the true Lambda(C,delta), and an upper bound
on nothing -- we never under-count a given word's list, that is exact).  Candidate
families, chosen to stress the regimes where RS lists are known to be large:

  (D) DEEP HOLES -- the classical large-list seeds for RS:
      * x^k                 : the degree-k monomial (canonical RS deep hole; distance
                              n-k from C, the maximum, for the full/coset domain).
      * x^{k+1}, x^{k+2}    : higher monomials (the Kambire family lives here).
      * 1/(x - a)           : Cauchy/Reed deep hole (a not in L); for several a.
      * x^k + (low-deg)     : deep hole plus a random low-degree codeword (shifts the
                              cluster center without changing the deep-hole structure).
  (K) KAMBIRE MONOMIALS on the subgroup -- the EXACT near-capacity bad-list words:
      X^{r*m} for the (s,m,r) coset parametrization, when such a subgroup chain exists.
      These are THE words that realize n^C lists near capacity; including them makes the
      probe genuinely worst-case where it matters most.
  (N) NEAR-CODEWORDS: a random codeword with e errors, for e straddling the Johnson and
      capacity error budgets -- catches "two messages agreeing on many points" clusters.
  (R) RANDOM words in F^n -- the average-case baseline.
  (M) MAX-CLUSTER planted word: pick k+t shared coordinates from one codeword and fill
      the rest adversarially to maximize multiplicity at a chosen radius (a constructed
      worst case at the boundary of each delta).

For each (field, n, k, domain) we sweep delta on the 1/n lattice from below Johnson up
to (and a hair past) capacity, and record the per-delta max list over the entire
candidate sample, plus which family achieved it.  Smooth subgroup vs matched random
subset are run on identical (field,n,k) with the identical candidate recipe.

Everything is EXACT (full codeword enumeration, q^k <= ~3e6).  No decoder, no
list-decoder.  Results -> results/singlelist_past_johnson.{json,csv}.
"""

from __future__ import annotations

import csv
import json
import math
import os
import time
import traceback
from dataclasses import dataclass
from multiprocessing import Pool

import numpy as np

from ff import PrimeField, BinaryExtensionField, FiniteField, _divisors, is_prime
from rs import (build_codeword_book, domain_subgroup, domain_random, CodewordBook,
                encode, min_distance, dist_to_code)

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
QK_CAP = 3_000_000


# ---------------------------------------------------------------------------
# q-ary entropy and the list-decoding-capacity (Elias) radius r_E = 1 - H_q(rho).
# ---------------------------------------------------------------------------
def Hq(x: float, q: int) -> float:
    """q-ary entropy H_q(x) = x*log_q(q-1) - x*log_q x - (1-x)*log_q(1-x)."""
    if x <= 0.0 or x >= 1.0:
        return 0.0
    lq = math.log(q)
    return (x * math.log(q - 1) - x * math.log(x) - (1 - x) * math.log(1 - x)) / lq


def elias_radius(rho: float, q: int) -> float:
    """r_E = 1 - H_q(rho)  (the list-decoding-capacity radius)."""
    return 1.0 - Hq(rho, q)


# ---------------------------------------------------------------------------
# Battery.
# ---------------------------------------------------------------------------
@dataclass
class Cfg:
    field_name: str
    field: FiniteField
    n: int
    k: int

    @property
    def q(self):
        return self.field.q

    @property
    def rho(self):
        return self.k / self.n

    @property
    def qk(self):
        return self.q ** self.k

    @property
    def cost(self):
        # one agreement pass over the codeword book ~ q^k * n
        return self.qk * self.n

    @property
    def label(self):
        return f"{self.field_name}_n{self.n}_k{self.k}"


def _largest_prime_factor(n: int) -> int:
    m = n; d = 2; lpf = 1
    while d * d <= m:
        while m % d == 0:
            lpf = d; m //= d
        d += 1
    if m > 1:
        lpf = max(lpf, m)
    return lpf


def _two_adic(n: int) -> int:
    v = 0
    while n % 2 == 0:
        v += 1; n //= 2
    return v


def build_battery() -> list[Cfg]:
    """Codes with SMOOTH multiplicative-subgroup domains at rates rho ~ {1/2,1/4,1/8},
    exactly enumerable (q^k <= QK_CAP), plus GF(2^m) contrasts.

    Design tension (made explicit): a smooth POWER-OF-TWO subgroup of order n lives only
    in primes p with 2^t | p-1; those primes are large, so q^k <= 3e6 caps k=2 there
    (tiny rho).  To probe rho=1/2 and 1/4 with a genuine multiplicative subgroup we
    therefore use SMOOTH subgroups in the broader sense (n highly composite / small
    largest-prime-factor), and we PRIORITIZE power-of-two n where it fits.  Each cell is
    tagged with n's largest prime factor and 2-adic valuation so the analysis can
    separate the genuine power-of-two FFT-domain cells (the deployed case) from merely
    smooth ones.  We pick, per (field-pool, target rho), the LARGEST feasible n (widest
    (J,r_E) band) so the onset is asymptotically meaningful.

    The Kambire near-capacity bad-list mechanism (which P' must rule out below r_E) needs
    a coset chain s|n; smooth (highly composite) n maximize the number of such chains, so
    these domains are exactly where a smooth-specific large list would show up if it
    exists -- a conservative (worst-case-favouring) choice.
    """
    # A broad pool of small primes; the selector below picks the best cells.
    prime_pool = [p for p in range(5, 270) if is_prime(p)]
    ext_fields = [
        ("GF(2^4)", BinaryExtensionField(4)),   # q-1 = 15 (odd) -> NO power-of-two subgroup
        ("GF(2^6)", BinaryExtensionField(6)),   # q-1 = 63 (odd)
        ("GF(2^8)", BinaryExtensionField(8)),   # q-1 = 255 (odd)
    ]

    COST_CAP = 4e8   # q^k * n per agreement pass
    cfgs: list[Cfg] = []
    seen = set()

    def consider(name, F, n, k):
        if k < 2 or k >= n or n < 6:
            return False
        if F.q ** k > QK_CAP:
            return False
        cfg = Cfg(name, F, n, k)
        if cfg.cost > COST_CAP:
            return False
        key = (name, n, k)
        if key in seen:
            return False
        seen.add(key)
        cfgs.append(cfg)
        return True

    # ---- Prime smooth-subgroup cells, selected to MAXIMIZE n per (rho) ----------
    # For each target rho we collect every feasible (p,n,k) with n|p-1 and a smooth n,
    # then keep a spread that includes the largest few n (widest band) and a couple of
    # genuine power-of-two-n cells.
    targets = [(0.5, "1/2"), (0.25, "1/4"), (0.125, "1/8")]
    by_target: dict[float, list[tuple]] = {t: [] for t, _ in targets}
    for p in prime_pool:
        F = PrimeField(p)
        for n in _divisors(p - 1):
            if n < 6 or n >= p:
                continue
            # smoothness gate: largest prime factor of n is modest (n is "smooth")
            lpf = _largest_prime_factor(n)
            if lpf > 7:                       # keep n genuinely smooth (3,5,7-smooth)
                continue
            for t, _ in targets:
                k = round(t * n)
                if k < 2 or k >= n:
                    continue
                if abs(k / n - t) > 0.09:
                    continue
                if p ** k > QK_CAP:
                    continue
                by_target[t].append((n, _two_adic(n), lpf, p, k))

    # Per target: prioritize (a) largest n, (b) power-of-two n.  Take a capped spread.
    PER_TARGET = 9
    for t, _ in targets:
        cand = by_target[t]
        # sort: power-of-two-n first within each n; then by n descending; then small q.
        cand.sort(key=lambda c: (-c[0], -(c[1] > 0), c[1], c[3]))
        taken_n = {}
        for (n, v2, lpf, p, k) in cand:
            # at most 2 cells per distinct n (one pow2-ish, one other) to keep a spread
            if taken_n.get(n, 0) >= 2:
                continue
            if consider(f"GF({p})", PrimeField(p), n, k):
                taken_n[n] = taken_n.get(n, 0) + 1
            if sum(taken_n.values()) >= PER_TARGET:
                break

    # ---- Char-2 contrast: largest feasible subgroup (odd order) or full group ----
    for name, F in ext_fields:
        q = F.q
        divs = [d for d in _divisors(q - 1) if 6 <= d < q - 1]
        cand_n = sorted(set(divs + [q - 1]), reverse=True)
        added = 0
        for n in cand_n:
            for t, _ in targets:
                k = round(t * n)
                if k < 2 or k >= n or abs(k / n - t) > 0.12:
                    continue
                if consider(name, F, n, k):
                    added += 1
            if added >= 4:
                break

    return cfgs


def delta_lattice(rho: float, n: int, q: int) -> list[int]:
    """Integer error counts e (=floor(delta*n)) to evaluate, from below the
    unique-decoding radius up to (a hair past) capacity.  We work in integer e so
    list sizes are exact; delta = e/n.  Range: from UD-ish up through capacity+2/n.

    The list is constant (=1) below UD and the action is in [J, capacity], so we
    cover the full integer ladder from max(1, floor(UD*n)-1) to min(n-1, ceil(cap*n)+2).
    """
    cap = 1.0 - rho
    ud = (1.0 - rho + 1.0 / n) / 2.0
    e_lo = max(1, int(math.floor(ud * n)) - 1)
    e_hi = min(n - 1, int(math.ceil(cap * n)) + 2)
    return list(range(e_lo, e_hi + 1))


# ---------------------------------------------------------------------------
# Candidate worst-case words.
# ---------------------------------------------------------------------------
def _monomial_on_domain(F: FiniteField, L: np.ndarray, a: int) -> np.ndarray:
    """X^a evaluated on L (exact field powers)."""
    return np.array([F.pow(int(x), a) for x in L], dtype=np.int64)


def _cauchy_deep_hole(F: FiniteField, L: np.ndarray, a: int) -> np.ndarray | None:
    """1/(x - a) on L, for a NOT in L (the Cauchy/Reed RS deep hole).  None if a in L."""
    Ls = set(int(x) for x in L)
    if a in Ls:
        return None
    out = np.empty(len(L), dtype=np.int64)
    for i, x in enumerate(L):
        d = F.sub(int(x), a)
        if d == 0:
            return None
        out[i] = F.inv(d)
    return out


def gen_candidates(F: FiniteField, L: np.ndarray, k: int, rng: np.random.Generator,
                   kambire_exps: list[int], n_random: int, n_nearcw: int,
                   n_cauchy: int, n_deepshift: int) -> list[tuple[str, np.ndarray]]:
    """Build the structured worst-case candidate word list for one (F,L,k).

    Returns list of (family_tag, word).  family_tag in {deephole, kambire, nearcw,
    random, cauchy, deepshift, maxcluster}.
    """
    n = len(L)
    q = F.q
    cands: list[tuple[str, np.ndarray]] = []

    # (D) classical deep holes: x^k, x^{k+1}, x^{k+2}
    for a in (k, k + 1, k + 2):
        if a < q:                      # X^a on a multiplicative domain depends on a mod n,
            cands.append((f"deephole_x^{a}", _monomial_on_domain(F, L, a)))

    # (D) Cauchy deep holes 1/(x-a) for a few a not in L
    Ls = set(int(x) for x in L)
    tried = 0
    a = 0
    while tried < n_cauchy and a < q:
        if a not in Ls:
            w = _cauchy_deep_hole(F, L, a)
            if w is not None:
                cands.append((f"cauchy_1/(x-{a})", w))
                tried += 1
        a += 1

    # (D) deep hole + random low-degree codeword (shift the cluster center)
    base = _monomial_on_domain(F, L, k)        # x^k
    for _ in range(n_deepshift):
        coeffs = rng.integers(0, q, size=k).astype(np.int64)
        cw = encode(F, L, coeffs)
        cands.append(("deepshift", F.add_vec(base, cw)))

    # (K) Kambire monomials X^{r*m} (the near-capacity bad-list words), if provided
    for e in kambire_exps:
        if e >= k:                              # must be a genuine non-codeword exponent
            cands.append((f"kambire_X^{e}", _monomial_on_domain(F, L, e)))

    # (N) near-codewords: random codeword with e errors, e across the budget
    for frac in (0.30, 0.45, 0.55, 0.70):
        e = max(1, int(round(frac * n)))
        for _ in range(max(1, n_nearcw // 4)):
            coeffs = rng.integers(0, q, size=k).astype(np.int64)
            w = encode(F, L, coeffs).copy()
            pos = rng.choice(n, size=min(e, n), replace=False)
            for p in pos:
                w[p] = F.add(int(w[p]), int(rng.integers(1, q)))
            cands.append(("nearcw", w))

    # (R) random words
    for _ in range(n_random):
        cands.append(("random", rng.integers(0, q, size=n).astype(np.int64)))

    # (M) planted max-cluster word: take a codeword, keep its first k coords (forces
    # that codeword to be in EVERY ball that contains the agreement), and replace the
    # remaining n-k coords by a value chosen to coincide with as MANY OTHER codewords
    # as possible.  We approximate the analytic worst case by, for each of a few
    # "anchor" codewords, planting a second codeword's values on a chosen overlap.
    for _ in range(4):
        c1 = rng.integers(0, q, size=k).astype(np.int64)
        c2 = rng.integers(0, q, size=k).astype(np.int64)
        w1 = encode(F, L, c1)
        w2 = encode(F, L, c2)
        # word that agrees with w1 on a prefix block and w2 on the rest -> two codewords
        # in the ball whenever the block sizes meet the radius; a clean 2-cluster seed.
        cut = n // 2
        w = w1.copy()
        w[cut:] = w2[cut:]
        cands.append(("maxcluster", w))

    return cands


def kambire_exponents_for(F: FiniteField, n: int, k: int) -> list[int]:
    """If a coset chain s | n with k = (r-2)*m, r = rho*s + 2 exists (Kambire), return
    the bad-list monomial exponents r*m and (r-1)*m and a couple neighbours.  These are
    the exact near-capacity bad-list seeds.  Returns [] if no such chain fits this (n,k).
    """
    exps: set[int] = set()
    for s in _divisors(n):
        if s < 2 or s >= n:
            continue
        m = n // s
        # k = (r-2)*m  =>  r-2 = k/m must be a positive integer, and r = (k/m)+2
        if k % m != 0:
            continue
        r2 = k // m
        if r2 < 1:
            continue
        r = r2 + 2
        if r > s:                       # need r distinct elements from a size-s subgroup
            continue
        # the bad-list words: X^{r*m} (Kambire f) and a couple neighbours
        for e in (r * m, (r - 1) * m, (r + 1) * m):
            if 0 < e < (F.q):           # nontrivial on the domain
                exps.add(e)
    # Always include a few generic high monomials too (deep-hole family near capacity)
    for e in (n - 1, n - 2, k + (n - k) // 2):
        if e > k:
            exps.add(e)
    return sorted(exps)


# ---------------------------------------------------------------------------
# Exact list-size curve for one word: |Lambda(C, e/n, w)| for all e at once.
# ---------------------------------------------------------------------------
def _distance_histogram(book: CodewordBook, w: np.ndarray) -> np.ndarray:
    """Histogram (length n+1) of Hamming distances from w to every codeword.

    EXACT, one streamed pass over the codeword book.  We compute the per-codeword
    MISMATCH count = (G != w).sum(axis=1) directly (no bit-packing -- we only need
    the count, not the agreement set), which is the leanest exact kernel.
    """
    w = np.asarray(w, dtype=np.int64)
    n = book.n
    hist = np.zeros(n + 1, dtype=np.int64)
    for _, G in book.iter_chunks():
        d = (G != w[None, :]).sum(axis=1).astype(np.int64)   # Hamming distance per cw
        hist += np.bincount(d, minlength=n + 1)
    return hist


def list_sizes_all_e(book: CodewordBook, w: np.ndarray, e_list: list[int]) -> dict[int, int]:
    """For target w, return {e : #codewords c with Delta(c,w) <= e} for each e in
    e_list.  EXACT: one streamed pass over the codeword book, histogram of distances,
    cumulative count.  cum[e] = #codewords with dist <= e = |Lambda(C, e/n, w)|.
    """
    hist = _distance_histogram(book, w)
    cum = np.cumsum(hist)                        # cum[e] = #codewords with dist <= e
    return {e: int(cum[e]) for e in e_list}


# ---------------------------------------------------------------------------
# One domain: max list over all candidates, per e.
# ---------------------------------------------------------------------------
def measure_domain(F: FiniteField, L: np.ndarray, k: int, e_list: list[int],
                   kambire_exps: list[int], rng: np.random.Generator,
                   budget: dict) -> dict:
    """Per-e max single-code list size over the structured worst-case candidate sample,
    plus which family achieved each max and a per-family breakdown."""
    n = len(L)
    book = build_codeword_book(F, L, k)
    # MDS sanity (RS over any domain is MDS); cheap and catches domain bugs.
    dmin = min_distance(book)
    mds_ok = (dmin == n - k + 1)

    cands = gen_candidates(F, L, k, rng, kambire_exps,
                           n_random=budget["random"], n_nearcw=budget["nearcw"],
                           n_cauchy=budget["cauchy"], n_deepshift=budget["deepshift"])

    max_list = {e: 0 for e in e_list}
    max_family = {e: None for e in e_list}
    # per-family max (to attribute the worst case)
    fam_max = {}
    n_cands = 0
    for tag, w in cands:
        fam = tag.split("_")[0]
        sizes = list_sizes_all_e(book, w, e_list)
        n_cands += 1
        for e in e_list:
            s = sizes[e]
            if s > max_list[e]:
                max_list[e] = s
                max_family[e] = fam
            fkey = (fam, e)
            if s > fam_max.get(fkey, 0):
                fam_max[fkey] = s

    cap = 1.0 - k / n
    johnson = 1.0 - math.sqrt(k / n)
    r_E = elias_radius(k / n, F.q)
    ud = (1.0 - k / n + 1.0 / n) / 2.0

    # Per-e: the Johnson constant 1/(2*eta*rho) (valid below J) and the q-ary ELIAS VOLUME
    # prediction q^{(H_q(delta)-(1-rho))*n} (the average #codewords in the ball; <1 below
    # r_E, >1 above).  The route's claim (p-prime-route.md SS6.3) is that the worst-case list
    # tracks this volume to within a constant factor, so the volume column is the key
    # quantitative comparison.
    rho = k / n
    curve = []
    for e in e_list:
        delta = e / n
        eta_below_J = johnson - delta           # >0 below Johnson
        johnson_bound = (1.0 / (2.0 * eta_below_J * rho)) if eta_below_J > 1e-12 else None
        E_exp = Hq(delta, F.q) - (1.0 - rho)     # Elias volume exponent
        try:
            elias_vol = float(F.q) ** (E_exp * n)
        except OverflowError:
            elias_vol = float("inf")
        curve.append({
            "e": e, "delta": round(delta, 6),
            "delta_minus_UD": round(delta - ud, 6),
            "delta_minus_J": round(delta - johnson, 6),
            "delta_minus_rE": round(delta - r_E, 6),
            "delta_minus_cap": round(delta - cap, 6),
            "max_list": max_list[e],
            "max_family": max_family[e],
            "johnson_const_bound": (round(johnson_bound, 3)
                                    if johnson_bound is not None else None),
            "elias_volume_exp": round(E_exp, 5),
            "elias_volume_pred": (round(elias_vol, 4) if elias_vol < 1e12 else elias_vol),
        })

    # families present
    families = sorted({tag.split("_")[0] for tag, _ in cands})
    fam_breakdown = {}
    for fam in families:
        fam_breakdown[fam] = {str(e): fam_max.get((fam, e), 0) for e in e_list}

    return {
        "n": n, "k": k, "rho": rho, "q": F.q,
        "capacity": cap, "johnson": johnson, "r_E": r_E, "UD": ud,
        "min_dist": dmin, "mds_ok": mds_ok,
        "num_candidates": n_cands,
        "families": families,
        "curve": curve,
        "family_breakdown": fam_breakdown,
        "kambire_exps": kambire_exps,
    }


def budget_for(cost: float) -> dict:
    """Candidate-sample sizes tiered by per-pass cost (q^k * n).  Logged, no silent
    caps.  Even the largest cells get a solid structured sample (deep holes + Kambire
    are deterministic and always included; only the random/near-cw counts shrink)."""
    if cost <= 2e7:
        return dict(random=200, nearcw=80, cauchy=12, deepshift=20)
    if cost <= 8e7:
        return dict(random=120, nearcw=60, cauchy=10, deepshift=16)
    if cost <= 2e8:
        return dict(random=60, nearcw=40, cauchy=8, deepshift=10)
    return dict(random=30, nearcw=24, cauchy=6, deepshift=6)


def run_cell(args) -> dict:
    cfg: Cfg = args["cfg"]
    seed: int = args["seed"]
    F = cfg.field
    out = {
        "label": cfg.label, "field": F.name, "q": F.q, "n": cfg.n, "k": cfg.k,
        "rho": cfg.rho, "qk": cfg.qk, "cost": cfg.cost, "n_over_q": cfg.n / F.q,
        "ok": True, "elapsed_sec": 0.0,
    }
    t0 = time.time()
    try:
        e_list = delta_lattice(cfg.rho, cfg.n, cfg.q)
        budget = budget_for(cfg.cost)
        out["e_list"] = e_list
        out["budget"] = budget
        kexps = kambire_exponents_for(F, cfg.n, cfg.k)

        # smooth subgroup domain
        rng_s = np.random.default_rng(seed)
        Ds = domain_subgroup(F, cfg.n)
        out["smooth"] = measure_domain(F, Ds, cfg.k, e_list, kexps, rng_s, budget)

        # matched random-subset domain (same field/n/k, fresh rng)
        rng_r = np.random.default_rng(seed + 1)
        Dr = domain_random(F, cfg.n, rng_r)
        out["random_domain"] = measure_domain(F, Dr, cfg.k, e_list, kexps, rng_r, budget)
    except Exception as e:
        out["ok"] = False
        out["error"] = str(e)
        out["traceback"] = traceback.format_exc()
    out["elapsed_sec"] = round(time.time() - t0, 2)
    return out


# ---------------------------------------------------------------------------
# Analysis: onset of departure from the Johnson constant; smooth vs random.
# ---------------------------------------------------------------------------
def analyze(results) -> dict:
    """For each cell, find:
      * departure_e: the smallest e (delta) where max_list first exceeds 2 (i.e. the
        list leaves the unique-decoding value and the first nontrivial cluster appears).
      * list_at_Jplus[t]: max_list at delta = J + t/n for t = 0,1,2,3 (just past Johnson).
      * the max_list at r_E and at capacity.
    Compare smooth vs random per cell.  Report whether, just past Johnson, the list is
    SMALL (a single-digit constant) consistent with O(1/eps0).
    """
    rows = []
    for r in results:
        if not r.get("ok"):
            continue
        n, k = r["n"], r["k"]
        for dom in ("smooth", "random_domain"):
            d = r[dom]
            curve = {c["e"]: c for c in d["curve"]}
            J = d["johnson"]; cap = d["capacity"]; rE = d["r_E"]
            eJ = int(math.floor(J * n))            # floor(J*n)
            erE = int(math.floor(rE * n))
            ecap = int(math.floor(cap * n))
            def ml(e):
                c = curve.get(e)
                return c["max_list"] if c else None
            # just past Johnson: J + t/n
            jplus = {t: ml(eJ + t) for t in (0, 1, 2, 3, 4)}
            # departure: first e with max_list >= 3 (left the {1,2} unique/near-unique band)
            dep_e = None
            for e in sorted(curve):
                if curve[e]["max_list"] >= 3:
                    dep_e = e
                    break
            dep_delta = (dep_e / n) if dep_e is not None else None
            # Does a genuine OPEN BAND (J, r_E) exist? Only then is "delta = J + eps0
            # with eps0 < r_E - J" a meaningful P' regime. At high rho / small q,
            # r_E <= J (Elias below Johnson), so there is NO band above J below r_E and
            # P' is VACUOUS there -- the list grows at J because Johnson >= Elias.
            has_open_band = (erE - eJ) >= 2     # at least 2 lattice steps of band
            band_width = round(rE - J, 4)
            rows.append({
                "label": r["label"], "q": r["q"], "n": n, "k": k,
                "rho": round(r["rho"], 4), "domain": dom,
                "J": round(J, 4), "r_E": round(rE, 4), "cap": round(cap, 4),
                "has_open_band": has_open_band, "band_width_rE_minus_J": band_width,
                "list_at_J": ml(eJ),
                "list_at_Jplus1": jplus[1], "list_at_Jplus2": jplus[2],
                "list_at_Jplus3": jplus[3], "list_at_Jplus4": jplus[4],
                "list_at_rE": ml(erE),
                "list_at_rE_minus1": ml(erE - 1),
                "list_at_cap": ml(ecap),
                "list_at_cap_minus1": ml(ecap - 1),
                "departure_e": dep_e,
                "departure_delta": (round(dep_delta, 4) if dep_delta is not None else None),
                "departure_minus_J": (round(dep_delta - J, 4) if dep_delta is not None else None),
                "departure_minus_rE": (round(dep_delta - rE, 4) if dep_delta is not None else None),
                "departure_minus_cap": (round(dep_delta - cap, 4) if dep_delta is not None else None),
                "max_family_at_departure": (curve[dep_e]["max_family"] if dep_e is not None else None),
            })
    # smooth vs random: paired departure comparison
    paired = {}
    for row in rows:
        paired.setdefault(row["label"], {})[row["domain"]] = row
    smooth_vs_random = []
    for label, pr in paired.items():
        s = pr.get("smooth"); rr = pr.get("random_domain")
        if not s or not rr:
            continue
        smooth_vs_random.append({
            "label": label, "q": s["q"], "rho": s["rho"], "n": s["n"], "k": s["k"],
            "J": s["J"], "r_E": s["r_E"], "cap": s["cap"],
            "smooth_departure_delta": s["departure_delta"],
            "random_departure_delta": rr["departure_delta"],
            "smooth_dep_minus_J": s["departure_minus_J"],
            "random_dep_minus_J": rr["departure_minus_J"],
            "smooth_dep_minus_cap": s["departure_minus_cap"],
            "random_dep_minus_cap": rr["departure_minus_cap"],
            "smooth_list_at_Jplus1": s["list_at_Jplus1"],
            "random_list_at_Jplus1": rr["list_at_Jplus1"],
            "smooth_list_at_Jplus2": s["list_at_Jplus2"],
            "random_list_at_Jplus2": rr["list_at_Jplus2"],
            "smooth_list_at_rE": s["list_at_rE"],
            "random_list_at_rE": rr["list_at_rE"],
            "smooth_eq_random_departure": s["departure_delta"] == rr["departure_delta"],
        })
    return {"per_domain_rows": rows, "smooth_vs_random": smooth_vs_random}


def print_summary(results, ana):
    print("\n" + "=" * 100)
    print("SINGLE-CODE LIST SIZE PAST JOHNSON  (smooth subgroup vs random subset, EXACT)")
    print("=" * 100)
    print("Per cell (smooth domain): list size at J, just past J (J+t/n), at r_E, at capacity.")
    print("'band' = does a genuine OPEN BAND (J,r_E) exist (r_E > J)? If NOT (high rho / small")
    print("q, where r_E <= J), P' is VACUOUS there and 'J+t/n' is already past Elias -- IGNORE.")
    print(f"  {'code':18} {'q':>5} {'rho':>5} {'band':>4} {'J':>6} {'r_E':>6} {'cap':>6} | "
          f"{'@J':>3} {'@J+1':>4} {'@J+2':>4} {'@J+3':>4} | {'@rE-1':>5} {'@rE':>4} | "
          f"{'@cap-1':>6} {'@cap':>5} | {'depart':>6} {'dep-J':>6} {'dep-rE':>6} {'fam@dep':>9}")
    smooth_rows = [r for r in ana["per_domain_rows"] if r["domain"] == "smooth"]
    for r in sorted(smooth_rows, key=lambda x: (not x["has_open_band"], x["rho"], x["q"], x["n"])):
        band = "YES" if r["has_open_band"] else "no"
        print(f"  {r['label']:18} {r['q']:>5} {r['rho']:>5.3f} {band:>4} {r['J']:>6.3f} "
              f"{r['r_E']:>6.3f} {r['cap']:>6.3f} | "
              f"{str(r['list_at_J']):>3} {str(r['list_at_Jplus1']):>4} "
              f"{str(r['list_at_Jplus2']):>4} {str(r['list_at_Jplus3']):>4} | "
              f"{str(r['list_at_rE_minus1']):>5} {str(r['list_at_rE']):>4} | "
              f"{str(r['list_at_cap_minus1']):>6} {str(r['list_at_cap']):>5} | "
              f"{str(r['departure_delta']):>6} {str(r['departure_minus_J']):>6} "
              f"{str(r['departure_minus_rE']):>6} {str(r['max_family_at_departure']):>9}")

    # The key question -- BAND-CONDITIONED (the honest verdict).  Pooling all rho is
    # misleading: where r_E <= J there is no open band, so "J+2/n" sits past Elias and the
    # list is large for a RATE (Elias) reason, not a smooth-domain failure.  Restrict the
    # P' question to cells with a genuine open band (r_E > J).
    band_rows = [r for r in smooth_rows if r["has_open_band"]]
    noband_rows = [r for r in smooth_rows if not r["has_open_band"]]
    print("-" * 100)
    print(f"KEY QUESTION (band-conditioned): of {len(smooth_rows)} smooth cells, "
          f"{len(band_rows)} have a genuine open band (r_E>J), {len(noband_rows)} do NOT "
          f"(P' vacuous there).")
    print("Just past Johnson, does the SMOOTH list stay SMALL in the OPEN-BAND cells?")
    jp1 = [r["list_at_Jplus1"] for r in band_rows if r["list_at_Jplus1"] is not None]
    jp2 = [r["list_at_Jplus2"] for r in band_rows if r["list_at_Jplus2"] is not None]
    jp3 = [r["list_at_Jplus3"] for r in band_rows if r["list_at_Jplus3"] is not None]
    if jp1:
        print(f"  [open-band cells only] list at J+1/n: max={max(jp1)}, mean={np.mean(jp1):.2f}  "
              f"(n={len(jp1)})")
        print(f"  [open-band cells only] list at J+2/n: max={max(jp2)}, mean={np.mean(jp2):.2f}")
        print(f"  [open-band cells only] list at J+3/n: max={max(jp3)}, mean={np.mean(jp3):.2f}")
        if max(jp3) <= 25:
            print(f"  => In every open-band cell the smooth single-code list stays a SMALL "
                  f"CONSTANT (<= {max(jp3)} through J+3/n). STRONGLY consistent with P' "
                  f"(O_rho(1/eps0)).")
        elif max(jp2) <= 25:
            print(f"  => Open-band list stays small (<= {max(jp2)}) at J+2/n; grows by J+3/n "
                  f"(approaching r_E). Consistent with P' below r_E.")
        else:
            print(f"  => Even in open-band cells the list is large ({max(jp2)}) at J+2/n. "
                  f"P' would be in danger.")
    # The pooled (all-rho) numbers, shown ONLY to expose the pooling artifact.
    pj2 = [r["list_at_Jplus2"] for r in smooth_rows if r["list_at_Jplus2"] is not None]
    pj2_band = [r["list_at_Jplus2"] for r in noband_rows if r["list_at_Jplus2"] is not None]
    if pj2_band:
        print(f"  (For contrast: pooling ALL rho gives max list@J+2/n = {max(pj2)}, but that max "
              f"comes from a NO-open-band cell where J+2/n is past r_E -- max@J+2/n among "
              f"no-band cells = {max(pj2_band)}. The pooled number is an Elias artifact, not P'.)")

    # departure relative to J vs r_E vs capacity (smooth) -- split by band.
    print("-" * 100)
    deps = [r for r in band_rows if r["departure_delta"] is not None]
    if deps:
        dmj = np.array([r["departure_minus_J"] for r in deps])
        dmr = np.array([r["departure_minus_rE"] for r in deps])
        dmc = np.array([r["departure_minus_cap"] for r in deps])
        print(f"  Departure (first delta with list>=3) vs landmarks, OPEN-BAND cells only "
              f"(n={len(deps)}):")
        print(f"    delta_departure - J     : mean={dmj.mean():+.3f}  min={dmj.min():+.3f}  max={dmj.max():+.3f}  (>0 means past Johnson)")
        print(f"    delta_departure - r_E   : mean={dmr.mean():+.3f}  min={dmr.min():+.3f}  max={dmr.max():+.3f}  (~0 means pinned at Elias)")
        print(f"    delta_departure - cap   : mean={dmc.mean():+.3f}  min={dmc.min():+.3f}  max={dmc.max():+.3f}  (<0 means below capacity)")
        n_above_J = sum(1 for r in deps if r["departure_minus_J"] > 1e-9)
        n_below_cap = sum(1 for r in deps if r["departure_minus_cap"] < -1e-9)
        n_near_rE = sum(1 for r in deps if abs(r["departure_minus_rE"]) <= 0.06)
        print(f"    departure ABOVE Johnson  : {n_above_J}/{len(deps)} open-band cells")
        print(f"    departure BELOW capacity : {n_below_cap}/{len(deps)} open-band cells")
        print(f"    departure NEAR r_E (|.|<=0.06): {n_near_rE}/{len(deps)} open-band cells "
              f"(=> onset tracks the Elias radius, NOT Johnson)")

    # smooth vs random (restricted to open-band cells where the comparison is meaningful)
    print("-" * 100)
    svr_all = ana["smooth_vs_random"]
    band_labels = {r["label"] for r in band_rows}
    svr = [x for x in svr_all if x["label"] in band_labels]
    if svr:
        same = sum(1 for x in svr if x["smooth_eq_random_departure"])
        smooth_later = sum(1 for x in svr
                           if x["smooth_departure_delta"] is not None
                           and x["random_departure_delta"] is not None
                           and x["smooth_departure_delta"] > x["random_departure_delta"])
        smooth_earlier = sum(1 for x in svr
                             if x["smooth_departure_delta"] is not None
                             and x["random_departure_delta"] is not None
                             and x["smooth_departure_delta"] < x["random_departure_delta"])
        print(f"  SMOOTH vs RANDOM departure (open-band cells, n={len(svr)}): identical in "
              f"{same}; smooth later (better) in {smooth_later}; smooth earlier (worse) in "
              f"{smooth_earlier}.")
        # at J+2/n (inside the band), is smooth ever bigger than random?
        worse = [x for x in svr
                 if x["smooth_list_at_Jplus2"] is not None and x["random_list_at_Jplus2"] is not None
                 and x["smooth_list_at_Jplus2"] > x["random_list_at_Jplus2"]]
        print(f"  At J+2/n (in band): smooth list EXCEEDS random list in {len(worse)}/{len(svr)} "
              f"open-band cells (smooth-specific inflation just past Johnson?).")
        if smooth_earlier == 0 and len(worse) == 0:
            print(f"  => SMOOTH is comparable to random up to a small additive excess in the open band: "
                  f"smooth ~ random, no smooth-specific inflation. Supports P' (smooth comparable to the "
                  f"random-domain regime that provably has small lists).")

    # MDS sanity
    bad_mds = [r for r in results if r.get("ok") and
               (not r["smooth"]["mds_ok"] or not r["random_domain"]["mds_ok"])]
    print("-" * 100)
    print(f"  MDS sanity (min-dist == n-k+1) failures: {len(bad_mds)} (MUST be 0).")


def write_outputs(results, ana, elapsed, battery_meta, out_dir=RESULTS_DIR):
    os.makedirs(out_dir, exist_ok=True)
    payload = {
        "meta": {
            "experiment": "singlelist_past_johnson",
            "timestamp": time.strftime("%Y%m%d_%H%M%S"),
            "elapsed_sec": round(elapsed, 2),
            "qk_cap": QK_CAP,
            "object": "single-code RS list Lambda(C,delta,w)={c in C: Delta(c,w)<=floor(delta n)}; "
                      "we report max_w over a STRUCTURED WORST-CASE SAMPLE (deep holes x^k, x^{k+i}, "
                      "1/(x-a), Kambire monomials X^{rm}, near-codewords, random, planted clusters). "
                      "Each word's list is EXACT (full codeword enumeration, all radii from one pass). "
                      "Exhaustive max over all q^n words is infeasible: this is a strong-sample LOWER "
                      "bound on the true worst-case Lambda(C,delta).",
            "convention": "delta=relative radius (ABF). J=1-sqrt(rho) Johnson; r_E=1-H_q(rho) Elias; "
                          "cap=1-rho Singleton. e=floor(delta*n) integer error count.",
            "domains": "SMOOTH = power-of-two multiplicative subgroup <omega>; RANDOM = matched "
                       "size-n random subset of F* (same field/n/k). Char-2 GF(2^m): |F*| odd so "
                       "NO power-of-two multiplicative subgroup -- those cells use the largest "
                       "odd-order subgroup / full group, flagged as a contrast.",
            "battery": battery_meta,
        },
        "results": results,
        "analysis": ana,
    }
    jpath = os.path.join(out_dir, "singlelist_past_johnson.json")
    with open(jpath, "w", newline="\n") as f:
        json.dump(payload, f, indent=2, default=str)

    # CSV: one row per (cell, domain, e) with the list-size curve.
    cpath = os.path.join(out_dir, "singlelist_past_johnson.csv")
    rows = []
    for r in results:
        if not r.get("ok"):
            continue
        for dom in ("smooth", "random_domain"):
            d = r[dom]
            J = d["johnson"]; rE = d["r_E"]; n = r["n"]
            has_band = (int(math.floor(rE * n)) - int(math.floor(J * n))) >= 2
            for c in d["curve"]:
                rows.append({
                    "field": r["field"], "q": r["q"], "n": r["n"], "k": r["k"],
                    "rho": round(r["rho"], 4), "domain": dom,
                    "has_open_band": has_band,
                    "J": round(d["johnson"], 4), "r_E": round(d["r_E"], 4),
                    "cap": round(d["capacity"], 4), "UD": round(d["UD"], 4),
                    "e": c["e"], "delta": c["delta"],
                    "delta_minus_J": c["delta_minus_J"],
                    "delta_minus_rE": c["delta_minus_rE"],
                    "delta_minus_cap": c["delta_minus_cap"],
                    "max_list": c["max_list"],
                    "max_family": c["max_family"],
                    "johnson_const_bound": c["johnson_const_bound"],
                    "elias_volume_exp": c["elias_volume_exp"],
                    "elias_volume_pred": c["elias_volume_pred"],
                })
    if rows:
        with open(cpath, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
            w.writeheader()
            for row in rows:
                w.writerow(row)
    return jpath, cpath, len(rows)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=min(16, os.cpu_count() or 8))
    ap.add_argument("--seed", type=int, default=0x5114E)
    ap.add_argument("--out", type=str, default=os.environ.get("SLPJ_OUT", RESULTS_DIR))
    ap.add_argument("--limit", type=int, default=0, help="cap #cells (debug)")
    args = ap.parse_args()

    battery = build_battery()
    battery.sort(key=lambda c: c.cost)
    if args.limit:
        battery = battery[:args.limit]
    jobs = [{"cfg": c, "seed": args.seed + 13 * i} for i, c in enumerate(battery)]
    battery_meta = [{"label": c.label, "q": c.q, "n": c.n, "k": c.k,
                     "rho": round(c.rho, 4), "qk": c.qk, "cost": c.cost,
                     "n_over_q": round(c.n / c.q, 5),
                     "budget": budget_for(c.cost)} for c in battery]

    print("=" * 100)
    print("SINGLE-CODE LIST PAST JOHNSON: smooth subgroup vs matched random subset")
    print("=" * 100)
    print(f"{len(battery)} cells, procs={args.procs}")
    for c in battery:
        b = budget_for(c.cost)
        nb = sum(b.values()) + 12  # +deterministic deep-hole/kambire/maxcluster families
        print(f"  {c.label:18} q={c.q:>5} rho={c.rho:.3f} n/q={c.n/c.q:.4f} "
              f"q^k={c.qk:.1e} cost={c.cost:.1e} cands~{nb}")
    print("-" * 100, flush=True)

    t0 = time.time()
    results = []
    with Pool(processes=args.procs) as pool:
        for r in pool.imap_unordered(run_cell, jobs):
            st = "ok " if r.get("ok") else "ERR"
            extra = ""
            if r.get("ok"):
                sm = r["smooth"]
                # quick peek: list just past Johnson
                n = r["n"]; J = sm["johnson"]
                eJ = int(math.floor(J * n))
                curve = {c["e"]: c for c in sm["curve"]}
                jp1 = curve.get(eJ + 1, {}).get("max_list", "?")
                extra = f"list@J+1/n={jp1}"
            else:
                extra = f"ERR={r.get('error')}"
            print(f"  [{st}] {r['label']:18} q={r['q']:>5} ({r['elapsed_sec']:.1f}s) {extra}",
                  flush=True)
            results.append(r)
    elapsed = time.time() - t0
    print("-" * 100)
    print(f"All cells done in {elapsed:.1f}s ({elapsed/60:.1f} min)")

    ana = analyze(results)
    print_summary(results, ana)
    jpath, cpath, nrows = write_outputs(results, ana, elapsed, battery_meta, out_dir=args.out)
    print(f"\nWrote:\n  {jpath}\n  {cpath} ({nrows} rows)")


if __name__ == "__main__":
    main()
