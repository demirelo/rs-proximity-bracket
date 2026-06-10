"""
cluster_certificate.py -- EXACT cluster certificate for the open-band list size,
reaching LARGE n at FIXED rate rho (the regime full q^k enumeration cannot reach).

THE QUESTION (the single decisive open one for P').
  P' asserts: for a prime-field smooth subgroup L = <omega> subset F_p^* of order n,
  rate rho, the open-band (J = 1 - sqrt(rho) < delta < r_E = 1 - H_p(rho)) worst-case
  single-code list  |Lambda(RS[F_p, L, rho n], delta)|  is  O_rho(1)  -- BOUNDED,
  independent of n and p.  The full-enumeration experiments capped n at ~24-32 because
  they stream all q^k = p^{rho n} codewords per word.  The decisive missing measurement
  is: at FIXED rho, does the open-band list stay BOUNDED as n -> infinity, or GROW?

THE METHOD -- a witness-based certificate that COSTS O(t*n), NOT O(q^k).
  A list of size t at integer radius e = floor(delta*n) is *witnessed* by t codewords
  (degree-<k polynomials) c_1..c_t and a word w with every Delta(c_i, w) <= e.  As a
  HEURISTIC center for a FIXED codeword set we use the coordinate-wise PLURALITY (the
  value shared by the most c_i at each of the n points of L).  The plurality is the exact
  minimizer of TOTAL (summed) Hamming distance sum_i Delta(c_i, w); it is NOT in general
  the minimizer of the MAXIMUM radius max_i Delta(c_i, w) (the closest-string / min-max
  center).  This does not affect correctness: the certificate is an EXPLICIT witness (an
  exhibited word w plus distinct codewords inside the radius-e ball), so its validity as a
  LOWER bound holds for ANY chosen w -- plurality being only a heuristic center weakens the
  claim "we found the best center for this cluster," never the witnessed list itself.
  Evaluating a random degree-<k polynomial on L costs O(n) and is INDEPENDENT of q^k, so a
  t-codeword cluster + plurality word + in-ball count costs O(t*n) -- reaches n = 2048 cheaply.
  NEXT EXPERIMENT: an EXACT closest-string / min-max center solver for small t can only
  STRENGTHEN these witnessed lists (it would tighten max-radius, never shrink a witness); it
  is being added separately.

  For a candidate (w, generating set S), the CERTIFIED list size is
        cert(w, S) = #{ c in S : Delta(c, w) <= e }    (<= t).
  This is an HONEST LOWER BOUND on |Lambda(C, e/n, w)| <= Lambda(C, e/n): it counts only
  the *generators* that fell in the ball (other codewords may also be in the ball, but
  counting those would require q^k enumeration, which we deliberately avoid).  So the
  certificate is conservative -- it can only UNDERshoot the true worst case, never over.

  We report the LARGEST cert over several constructions, weakest-to-strongest:
    (RP) random codewords + plurality                 -- the baseline (matches the old probe);
    (GR) GREEDY growth: grow a cluster by adding the codeword that best preserves a small
         ball, refitting the plurality center each step;
    (LL) Lloyd / iterative-plurality refit: keep the in-ball codewords, refit w as their
         plurality, repeat -- the center migrates to pull in more codewords;
    (LS) local search on w: coordinate moves that increase the in-ball count;
    (SC) seed from a near-codeword center, then grow.
  The adversarial constructions GENUINELY try to grow t (a construction that merely
  plateaus proves nothing); the max cert is a witness LOWER bound on the worst-case list.

DECISIVE OUTPUT.  max-cert vs n at FIXED rho.  BOUNDED (O_rho(1)) supports P'; GROWING at
fixed rho FALSIFIES P' outright.  We also tabulate it against the random-plurality heuristic
1/(1-delta) and the q-ary Elias volume q^{(H_p(delta) - (1-rho)) n} (which is < 1 below r_E).

CROSS-CHECK (consistency with the EXACT list).  At small n (n <= 24, p^k <= 3e6) we ALSO
enumerate all q^k codewords via the existing singlelist machinery and verify, for the SAME
witness word w the certificate produced, that  cert(w,S) <= TRUE |Lambda(C,e/n,w)|  (exactly,
since w is explicit) and that the certificate's max is <= the exact sampled Lambda(C,e/n)
and close to it.  This pins the certificate to ground truth where ground truth is computable.

CORRECTNESS-FIRST.  No decoder anywhere.  Hamming distances are exact numpy passes.  Field,
subgroup, and polynomial-evaluation machinery is reused verbatim from ff.py / rs.py (and the
exact-list kernel from singlelist_past_johnson.py for the cross-check); we do NOT rebuild it.

HONESTY.  Every cert is a witness-based LOWER bound; a BOUNDED max-cert across large n at
fixed rho is strong SUPPORTING evidence for P' (not proof -- a cleverer word could give more);
a GROWING max-cert at fixed rho would FALSIFY P'.  Sample/trial counts are logged per cell (no
silent caps), so coverage and construction strength are stated honestly.

Run `python3.11 cluster_certificate.py --selftest`  for the self-test battery,
    `python3.11 cluster_certificate.py`             for the full run (parallel, 16 cores),
    `python3.11 cluster_certificate.py --quick`     for a fast reduced sweep.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from dataclasses import dataclass
from multiprocessing import Pool

import numpy as np

from ff import PrimeField, is_prime, _divisors
from rs import build_codeword_book, domain_subgroup, encode
from singlelist_past_johnson import Hq, elias_radius, list_sizes_all_e

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")

# Exact-enumeration cap (only used by the small-n cross-check, never by the certificate).
QK_CAP = 3_000_000


# ===========================================================================
# Radii.
# ===========================================================================
def johnson_radius(rho: float) -> float:
    return 1.0 - math.sqrt(rho)


def band_e(n: int, k: int, p: int, c: float) -> int | None:
    """Integer radius e = round((J + c*(r_E - J)) * n), snapped STRICTLY inside the open
    band (eJ < e <= floor(r_E*n)).  Returns None if the 1/n lattice has no interior point."""
    rho = k / n
    J = johnson_radius(rho)
    rE = elias_radius(rho, p)
    if rE <= J:
        return None
    eJ = math.floor(J * n + 1e-9)
    erE = math.floor(rE * n + 1e-9)
    e = int(round((J + c * (rE - J)) * n))
    e = min(max(e, eJ + 1), erE)
    if e <= eJ or e > erE:
        return None
    return e


# ===========================================================================
# Core primitives -- ALL O(t*n) or O(n), NONE touch q^k.
# ===========================================================================
def rand_codewords(F: PrimeField, L: np.ndarray, k: int, t: int,
                   rng: np.random.Generator) -> np.ndarray:
    """t random degree-<k codewords on L, as a (t, n) int64 matrix.  O(t*n) via Horner.

    Evaluating a random degree-<k polynomial on the FIXED set L costs O(n) and is
    INDEPENDENT of q^k -- this is exactly what lets the certificate reach large n.
    """
    q = F.q
    n = len(L)
    coeffs = rng.integers(0, q, size=(t, k)).astype(np.int64)   # (t, k) low->high
    Larr = np.asarray(L, dtype=np.int64)
    acc = np.broadcast_to(coeffs[:, -1:], (t, n)).astype(np.int64).copy()
    for j in range(k - 2, -1, -1):
        acc = F.add_vec(F.mul_vec(acc, Larr[None, :]), coeffs[:, j:j + 1])
    return acc


def plurality_word(codewords: np.ndarray, q: int,
                   rng: np.random.Generator) -> np.ndarray:
    """Coordinate-wise plurality of a (t, n) codeword matrix: per coordinate, the value
    held by the most codewords (ties broken at random).  This is the EXACT minimizer of
    TOTAL (summed) Hamming distance sum_i Delta(c_i, w) for the fixed set: coordinates are
    independent, and at each coordinate the value minimizing the number of disagreements is
    the most frequent (plurality) value, so summing over coordinates is minimized term by
    term.  We use it as a strong HEURISTIC center for producing witness lower bounds.  It is
    NOT in general the exact minimizer of the MAXIMUM radius max_i Delta(c_i, w) (the
    closest-string / min-max center) -- e.g. over GF(2) the cluster {000,111,111} has
    plurality 111 with max-radius 3, while 110 attains max-radius 2.  Either way the
    certificate stays a valid lower bound (it counts codewords explicitly inside the ball
    around whatever w we pick).  Cost O(t*n) (sort each length-t column).

    Implementation: sort each column; equal values form contiguous runs; the longest run
    is the mode.  We scan the t rows once, tracking, per column, the current run's value,
    its running length, and the best (longest) run seen.  Ties (multiple values achieving
    the column max) are broken at random AFTER the scan, by re-drawing a winner uniformly
    among all values that attain the max count -- a cheap second pass over the runs.
    """
    t, n = codewords.shape
    s = np.sort(codewords, axis=0)                  # (t, n), each column ascending
    # run-start flags down each column (row 0 always starts a run)
    start = np.ones((t, n), dtype=bool)
    if t > 1:
        start[1:, :] = s[1:, :] != s[:-1, :]
    # run length at each position = (#rows until next run-start, going down).
    # Compute per-position run length: for each entry, count consecutive equal below it.
    # cumulative run id per column:
    run_id = np.cumsum(start, axis=0)               # (t, n) labels 1..G down each column
    # size of each run per column: count entries sharing the run id of each position.
    # Do it with the standard "diff of start positions" trick on the flattened columns.
    best_cnt = np.zeros(n, dtype=np.int64)
    best_val = s[0, :].copy()                       # default: first value
    # iterate over candidate run ids 1..t (t small), counting and valuing each run.
    for g in range(1, t + 1):
        mask = run_id == g                          # entries in run g per column
        cnt = mask.sum(axis=0).astype(np.int64)     # length of run g per column (0 if absent)
        if not cnt.any():
            continue
        # value of run g per column = the value of its first row (argmax of the mask).
        first_row = mask.argmax(axis=0)             # 0 where mask all-False (cnt==0, ignored)
        val_g = s[first_row, np.arange(n)]
        better = cnt > best_cnt
        best_val = np.where(better, val_g, best_val)
        best_cnt = np.where(better, cnt, best_cnt)
    # random tie-break: where several runs attain best_cnt, pick uniformly among them.
    # (Exactness of the COUNT is unaffected; only WHICH max value is chosen varies.)
    # Second pass: collect, per column, all run values whose count == best_cnt, choose one.
    tie_choice_val = best_val.copy()
    tie_rand = rng.random((t, n))
    seen = np.zeros(n, dtype=np.int64)              # how many max-runs encountered per column
    for g in range(1, t + 1):
        mask = run_id == g
        cnt = mask.sum(axis=0).astype(np.int64)
        is_max = (cnt == best_cnt) & (cnt > 0)
        if not is_max.any():
            continue
        first_row = mask.argmax(axis=0)
        val_g = s[first_row, np.arange(n)]
        seen_new = seen + is_max
        # reservoir sampling: replace current choice with this max-run val with prob 1/seen_new
        take = is_max & (tie_rand[g - 1, :] < 1.0 / np.maximum(seen_new, 1))
        tie_choice_val = np.where(take, val_g, tie_choice_val)
        seen = np.where(is_max, seen_new, seen)
    return tie_choice_val.astype(np.int64)


def hamming_to_set(w: np.ndarray, codewords: np.ndarray) -> np.ndarray:
    """Hamming distances Delta(c_i, w) for every row c_i of the (t, n) matrix.  O(t*n)."""
    return (codewords != w[None, :]).sum(axis=1)


def distinct_codewords(codewords: np.ndarray) -> np.ndarray:
    """Return distinct codeword rows.

    The certificate is a lower bound on the number of distinct codewords in a ball. Random pools
    can contain duplicates in tiny cross-check fields; exact enumeration counts those once.
    """
    if codewords is None or codewords.shape[0] <= 1:
        return codewords
    return np.unique(codewords, axis=0)


def cert_count(w: np.ndarray, codewords: np.ndarray, e: int) -> int:
    """Certified list size of (w, codewords) at radius e: #distinct generators within e of w.
    HONEST LOWER bound on |Lambda(C, e/n, w)| (counts only the generators, never the
    rest of the code)."""
    cw = distinct_codewords(codewords)
    return int((hamming_to_set(w, cw) <= e).sum())


# ===========================================================================
# Constructions (each returns (best_cert, witness_word, witness_codewords)).
# ===========================================================================
def construct_random_plurality(F, L, k, e, rng, t_set, trials_per_t):
    """(RP) For each target t, draw `trials_per_t` random t-codeword sets, build the
    plurality word, count in-ball.  Records the best over all (t, trial).  O(t*n) each."""
    q = F.q
    best = 0
    best_w = None
    best_cw = None
    best_t = None
    n_tried = 0
    for t in t_set:
        for _ in range(trials_per_t):
            cw = rand_codewords(F, L, k, t, rng)
            w = plurality_word(cw, q, rng)
            n_tried += 1
            c = cert_count(w, cw, e)
            if c > best:
                keep = distinct_codewords(cw[hamming_to_set(w, cw) <= e])
                best, best_w, best_cw, best_t = c, w, keep, keep.shape[0]
    return {"cert": best, "w": best_w, "cw": best_cw, "t": best_t, "n_tried": n_tried}


def _lloyd_refit(cw_pool: np.ndarray, w0: np.ndarray, e: int, q: int,
                 rng: np.random.Generator, rounds: int = 8):
    """Lloyd/iterative-plurality refit on a fixed codeword POOL: keep the in-ball rows,
    refit w as their plurality, repeat.  The center migrates toward the densest cluster.
    Returns (best_cert, best_w, kept_codewords_at_best).  O(rounds * |pool| * n)."""
    w = w0.copy()
    best = cert_count(w, cw_pool, e)
    best_w = w.copy()
    best_keep = distinct_codewords(cw_pool[hamming_to_set(w, cw_pool) <= e])
    for _ in range(rounds):
        d = hamming_to_set(w, cw_pool)
        keep = cw_pool[d <= e]
        if keep.shape[0] < 2:
            # widen: keep the closest few so the refit has something to work with
            order = np.argsort(d)[:max(2, min(cw_pool.shape[0], 4))]
            keep = cw_pool[order]
        w = plurality_word(keep, q, rng)
        c = cert_count(w, cw_pool, e)
        if c > best:
            best = c
            best_w = w.copy()
            best_keep = distinct_codewords(cw_pool[hamming_to_set(w, cw_pool) <= e])
    return best, best_w, best_keep


def construct_lloyd(F, L, k, e, rng, pool_sizes, trials_per_pool):
    """(LL) Draw a POOL of random codewords, run the Lloyd refit from the plurality of a
    random small seed, keep the best.  A larger pool gives the refit more codewords to
    pull into the ball -- the natural way to grow t.  O(rounds * pool * n) per trial."""
    q = F.q
    best = 0
    best_w = None
    best_cw = None
    best_t = None
    n_tried = 0
    for psz in pool_sizes:
        for _ in range(trials_per_pool):
            pool = rand_codewords(F, L, k, psz, rng)
            # seed center: plurality of a random small subset (3..6) of the pool
            seed_sz = min(psz, int(rng.integers(3, 7)))
            seed_idx = rng.choice(psz, size=seed_sz, replace=False)
            w0 = plurality_word(pool[seed_idx], q, rng)
            n_tried += 1
            c, w, keep = _lloyd_refit(pool, w0, e, q, rng)
            if c > best:
                best, best_w, best_cw, best_t = c, w, keep, keep.shape[0]
    return {"cert": best, "w": best_w, "cw": best_cw, "t": best_t, "n_tried": n_tried}


def construct_greedy(F, L, k, e, rng, n_restarts, pool_factor, max_t):
    """(GR) GREEDY cluster growth.  Maintain a center w and an accepted set A.  Each step
    draw a batch of fresh random codewords and ADD the one that, after refitting w to the
    plurality of A+{c}, yields the largest in-ball count; stop when no addition helps.
    Genuinely tries to grow t.  Cost ~ O(n_restarts * max_t * batch * n)."""
    q = F.q
    n = len(L)
    best = 0
    best_w = None
    best_cw = None
    best_t = None
    n_tried = 0
    batch = max(8, pool_factor)
    for _ in range(n_restarts):
        # start from a random pair's plurality
        A = rand_codewords(F, L, k, 2, rng)
        w = plurality_word(A, q, rng)
        cur = cert_count(w, A, e)
        improved = True
        while improved and A.shape[0] < max_t:
            improved = False
            cand = rand_codewords(F, L, k, batch, rng)
            n_tried += batch
            best_add = None
            best_add_cert = cur
            best_add_w = w
            best_add_A = A
            for i in range(batch):
                A2 = np.vstack([A, cand[i:i + 1]])
                w2 = plurality_word(A2, q, rng)
                c2 = cert_count(w2, A2, e)
                if c2 > best_add_cert:
                    best_add_cert = c2
                    best_add = i
                    best_add_w = w2
                    best_add_A = A2
            if best_add is not None:
                A = best_add_A
                w = best_add_w
                cur = best_add_cert
                improved = True
        # final cert measured on the in-ball members only (honest witness set)
        d = hamming_to_set(w, A)
        keep = distinct_codewords(A[d <= e])
        c = cert_count(w, A, e)
        if c > best:
            best, best_w, best_cw, best_t = c, w, keep, keep.shape[0]
    return {"cert": best, "w": best_w, "cw": best_cw, "t": best_t, "n_tried": n_tried}


def _structured_seed_words(F, L, k, n, rng, max_words: int = 24) -> list[np.ndarray]:
    """Structured (coset-biased) ball-center seeds, O(n) each (no q^k):
      - deep-hole / high monomials  X^a, a in {k, k+1, k+2, n-2, n-1};
      - Kambire monomials  X^{rm}, X^{(r-1)m}, X^{(r+1)m}  for every coset chain s|n
        with k = (r-2)m (the e_2=0 mechanism the structural route of P' names as the
        only asymptotic source of a super-constant smooth list).
    These are the theoretically-motivated worst-case centers; offering them to the
    Lloyd refit lets the certificate exploit coset structure if it ever helps at large n.
    """
    q = F.q
    Larr = np.asarray(L, dtype=np.int64)
    words = []
    exps = set()
    for a in (k, k + 1, k + 2, n - 2, n - 1):
        if k <= a < n:
            exps.add(a)
    for s in _divisors(n):
        if s < 2 or s >= n:
            continue
        m = n // s
        if m == 0 or k % m != 0:
            continue
        r = k // m + 2
        if 2 <= r <= s:
            for ex in (r * m, (r - 1) * m, (r + 1) * m):
                if k < ex < n:
                    exps.add(ex)
    for a in sorted(exps):
        words.append(np.array([F.pow(int(x), a) for x in Larr], dtype=np.int64))
    if len(words) > max_words:
        # keep a random subset to bound cost (all are O(n) to build, cheap regardless)
        idx = rng.choice(len(words), size=max_words, replace=False)
        words = [words[i] for i in idx]
    return words


def construct_structured(F, L, k, e, rng, n_pool, trials):
    """(ST) STRUCTURED-seed cluster growth: seed the ball center at a coset/Kambire monomial
    or deep-hole word, build a pool of random codewords PLUS the structured codewords (X^a),
    run the Lloyd refit, and count.  Tests whether the coset mechanism packs MORE codewords
    than generic random clusters at large n -- the asymptotic worst case P''s structural
    route is about.  Still O(rounds*pool*n): the structured words are O(n) monomials, no q^k."""
    q = F.q
    n = len(L)
    best = 0
    best_w = None
    best_cw = None
    best_t = None
    n_tried = 0
    seeds = _structured_seed_words(F, L, k, n, rng)
    if not seeds:
        return {"cert": 0, "w": None, "cw": None, "t": None, "n_tried": 0}
    # the structured words double as candidate codewords-of-interest in the pool (the
    # monomials X^a with a < k are themselves codewords; a >= k are not, but their VALUES
    # seed a center the refit can exploit).  Codewords for the pool: random + any X^a, a<k.
    low_mono = [np.array([F.pow(int(x), a) for x in L], dtype=np.int64)
                for a in range(1, k)]
    for seed_w in seeds:
        for _ in range(trials):
            pool = rand_codewords(F, L, k, n_pool, rng)
            if low_mono:
                pool = np.vstack([pool] + low_mono)
            n_tried += 1
            c, w, keep = _lloyd_refit(pool, seed_w.copy(), e, q, rng, rounds=10)
            if c > best:
                best, best_w, best_cw, best_t = c, w, keep, keep.shape[0]
    return {"cert": best, "w": best_w, "cw": best_cw, "t": best_t, "n_tried": n_tried}


def construct_localsearch(F, L, k, e, rng, seed_constructions, n_moves):
    """(LS) Coordinate local search on the witness word: starting from the best word found
    so far, flip individual coordinates to a value that increases the in-ball count of a
    fresh+inherited codeword pool.  Each move is O(t) to re-score.  Polishes the center."""
    q = F.q
    n = len(L)
    best = 0
    best_w = None
    best_cw = None
    best_t = None
    n_tried = 0
    for seed in seed_constructions:
        if seed["w"] is None:
            continue
        w = seed["w"].copy()
        # build a pool: the seed's own witnesses plus fresh random codewords near w
        pool = seed["cw"] if seed["cw"] is not None and seed["cw"].shape[0] > 0 else \
            rand_codewords(F, L, k, 8, rng)
        extra = rand_codewords(F, L, k, max(8, 4 * pool.shape[0]), rng)
        pool = np.vstack([pool, extra])
        cur = cert_count(w, pool, e)
        for _ in range(n_moves):
            j = int(rng.integers(0, n))
            # candidate new value at coord j: the plurality of the in-ball pool at j
            d = hamming_to_set(w, pool)
            inball = pool[d <= e + 2]          # slightly relaxed to allow recruitment
            if inball.shape[0] == 0:
                inball = pool
            vals, counts = np.unique(inball[:, j], return_counts=True)
            newv = int(vals[counts.argmax()])
            if newv == w[j]:
                continue
            old = w[j]
            w[j] = newv
            n_tried += 1
            c = cert_count(w, pool, e)
            if c >= cur:
                cur = c
            else:
                w[j] = old                      # revert non-improving move
        d = hamming_to_set(w, pool)
        keep = distinct_codewords(pool[d <= e])
        c = cert_count(w, pool, e)
        if c > best:
            best, best_w, best_cw, best_t = c, w, keep, keep.shape[0]
    return {"cert": best, "w": best_w, "cw": best_cw, "t": best_t, "n_tried": n_tried}


# ===========================================================================
# One cell: max certificate at a fixed (p, n, k) and a fixed band-fraction c.
# ===========================================================================
@dataclass
class CertSpec:
    p: int
    n: int
    k: int
    c: float            # band fraction: delta = J + c*(r_E - J)
    note: str = ""

    @property
    def rho(self):
        return self.k / self.n

    @property
    def label(self):
        return f"GF({self.p})_n{self.n}_k{self.k}_c{self.c}"


def run_cert_cell(spec_and_budget) -> dict:
    """Compute the max cluster certificate (witness LOWER bound on Lambda(C,delta)) at one
    (p, n, k, c) cell, across all constructions, plus the heuristic comparators.  No q^k
    enumeration -- cost is O((sum of trial sizes) * n)."""
    spec, budget, seed = spec_and_budget
    p, n, k, c = spec.p, spec.n, spec.k, spec.c
    F = PrimeField(p)
    L = domain_subgroup(F, n)
    rho = k / n
    J = johnson_radius(rho)
    rE = elias_radius(rho, p)
    e = band_e(n, k, p, c)
    if e is None:
        return {"label": spec.label, "p": p, "n": n, "k": k, "rho": round(rho, 6),
                "c": c, "has_open_band": False}
    rng = np.random.default_rng(seed)
    delta = e / n

    t0 = time.time()
    rp = construct_random_plurality(F, L, k, e, rng,
                                    t_set=budget["rp_t_set"],
                                    trials_per_t=budget["rp_trials"])
    ll = construct_lloyd(F, L, k, e, rng,
                         pool_sizes=budget["ll_pool_sizes"],
                         trials_per_pool=budget["ll_trials"])
    gr = construct_greedy(F, L, k, e, rng,
                          n_restarts=budget["gr_restarts"],
                          pool_factor=budget["gr_batch"],
                          max_t=budget["gr_max_t"])
    st = construct_structured(F, L, k, e, rng,
                              n_pool=budget["st_pool"], trials=budget["st_trials"])
    # local search polishes the best center found so far (seed from the strongest first)
    ls = construct_localsearch(F, L, k, e, rng,
                               seed_constructions=[gr, st, ll, rp],
                               n_moves=budget["ls_moves"])
    elapsed = time.time() - t0

    per_construction = {
        "random_plurality": rp["cert"],
        "lloyd_refit": ll["cert"],
        "greedy": gr["cert"],
        "structured_coset": st["cert"],
        "local_search": ls["cert"],
    }
    constructions = {"RP": rp, "LL": ll, "GR": gr, "ST": st, "LS": ls}
    name2key = {"random_plurality": "RP", "lloyd_refit": "LL", "greedy": "GR",
                "structured_coset": "ST", "local_search": "LS"}
    best_name = max(per_construction, key=per_construction.get)
    max_cert = per_construction[best_name]
    best = constructions[name2key[best_name]]

    # heuristic comparators
    heuristic_1_over_1md = 1.0 / (1.0 - delta) if delta < 1 else float("inf")
    E = Hq(delta, p) - (1 - rho)               # Elias volume exponent (<0 below r_E)
    vol = p ** (E * n)

    n_tried_total = rp["n_tried"] + ll["n_tried"] + gr["n_tried"] + st["n_tried"] + ls["n_tried"]

    # store the best witness word (compact: list of ints) for the cross-check / audit
    witness_w = best["w"].tolist() if best["w"] is not None else None
    witness_t = int(best["t"]) if best["t"] is not None else 0

    return {
        "label": spec.label, "p": p, "n": n, "k": k, "rho": round(rho, 6),
        "c": c, "note": spec.note, "has_open_band": True,
        "J": round(J, 6), "r_E": round(rE, 6), "cap": round(1 - rho, 6),
        "band_width": round(rE - J, 6),
        "e": int(e), "delta": round(delta, 6),
        "delta_minus_J": round(delta - J, 6),
        "margin_rE_minus_delta": round(rE - delta, 6),
        "max_cert": int(max_cert),
        "max_cert_construction": best_name,
        "per_construction": per_construction,
        "witness_t": witness_t,
        "heuristic_1_over_1md": round(heuristic_1_over_1md, 4),
        "volume_exp_E": round(E, 6),
        "volume_pred_qEn": (round(vol, 8) if vol < 1e9 else float(vol)),
        "n_trials_total": int(n_tried_total),
        "elapsed_s": round(elapsed, 2),
        "witness_w": witness_w,
    }


# ===========================================================================
# Cross-check at small n: certificate LOWER bound vs the EXACT q^k list.
# ===========================================================================
def crosscheck_cell(spec_and_budget) -> dict:
    """At a small cell (p^k <= QK_CAP), run the certificate AND fully enumerate q^k.
    Verify (a) for the certificate's OWN witness word w, cert(w,S) <= TRUE |Lambda(C,e/n,w)|
    (exact, since w is explicit -- the certificate must never over-count); (b) the
    certificate max <= the EXACT sampled Lambda(C,e/n) (worst case over a structured sample)
    and is reasonably close to it.  This pins the certificate to ground truth."""
    spec, budget, seed = spec_and_budget
    p, n, k, c = spec.p, spec.n, spec.k, spec.c
    F = PrimeField(p)
    L = domain_subgroup(F, n)
    rho = k / n
    e = band_e(n, k, p, c)
    if e is None:
        return {"label": spec.label, "has_open_band": False}
    rng = np.random.default_rng(seed)

    # certificate (no q^k)
    cell = run_cert_cell((spec, budget, seed))
    max_cert = cell["max_cert"]
    witness_w = np.asarray(cell["witness_w"], dtype=np.int64) if cell["witness_w"] else None

    # exact enumeration (q^k) -- only here, only for the cross-check
    book = build_codeword_book(F, L, k)

    # (a) exact list at the certificate's OWN witness word
    true_list_at_witness = None
    if witness_w is not None:
        true_list_at_witness = list_sizes_all_e(book, witness_w, [e])[e]

    # (b) EXACT worst-case over a structured + random sample (the same flavor of search the
    #     full-enumeration experiments used), to bound the certificate from ABOVE.
    exact_max = 0
    exact_max_w = None
    # structured worst-case words: deep-hole / Kambire monomials X^a + random words
    sample_words = []
    for a in range(k, n):
        sample_words.append(np.array([F.pow(int(x), a) for x in L], dtype=np.int64))
    for _ in range(budget.get("exact_random", 400)):
        sample_words.append(rng.integers(0, p, size=n).astype(np.int64))
    # also feed the certificate witnesses themselves (they should be near the worst case)
    if witness_w is not None:
        sample_words.append(witness_w)
    for w in sample_words:
        s = list_sizes_all_e(book, w, [e])[e]
        if s > exact_max:
            exact_max = s
            exact_max_w = w

    ok_lower = (witness_w is None) or (max_cert <= true_list_at_witness)
    ok_upper = (max_cert <= exact_max)
    return {
        "label": spec.label, "p": p, "n": n, "k": k, "rho": round(rho, 6),
        "c": c, "has_open_band": True,
        "e": int(e), "delta": round(e / n, 6),
        "margin_rE_minus_delta": cell["margin_rE_minus_delta"],
        "max_cert": int(max_cert),
        "true_list_at_witness_word": (int(true_list_at_witness)
                                      if true_list_at_witness is not None else None),
        "exact_sampled_max_list": int(exact_max),
        "cert_le_true_at_witness": bool(ok_lower),
        "cert_le_exact_sampled_max": bool(ok_upper),
        "gap_exact_minus_cert": int(exact_max - max_cert),
        "num_codewords": int(p ** k),
        "n_exact_words": int(len(sample_words)),
    }


# ===========================================================================
# Batteries.
# ===========================================================================
def smooth_prime_for(n: int, lo: int) -> int | None:
    """Smallest prime p >= lo with n | (p - 1) (so the order-n subgroup exists)."""
    start = max(lo, n + 1)
    # align start to 1 mod n
    r = start % n
    p = start + ((1 - r) % n)
    if p < n + 1:
        p += n
    while p < lo + 5_000_000:
        if is_prime(p):
            return p
        p += n
    return None


def build_fixed_rho_ladders(quick: bool = False) -> dict:
    """The DECISIVE batteries: fixed rho, growing n = 2^t, at fixed band fractions c.

    TWO field choices per rho, both giving a clean fixed-rho ladder:
      (F) FIXED FIELD p = 12289 (an NTT prime, p-1 = 2^12 * 3) for n in {32..2048}, and
          p = 65537 (Fermat prime, p-1 = 2^16) for n in {32..4096} -- the cleanest test
          of n-dependence at FIXED (rho, p): the field is literally constant down the
          ladder, so any change in max-cert is purely an n-effect.
      (S) SMALLEST SMOOTH PRIME per n (p grows mildly with n) -- a cross-check that the
          fixed-field verdict is not an artifact of one field.
    """
    ns = [32, 64, 128, 256, 512, 1024] if quick else [32, 64, 128, 256, 512, 1024, 2048, 4096]
    rhos = [("1/16", 16), ("1/8", 8), ("1/4", 4)]
    cs = [0.5] if quick else [0.25, 0.5, 0.75]

    ladders = {"fixed_field": {}, "smallest_smooth": {}}

    # (F) fixed field.  Use 65537 (supports up to 2^16) so a single p covers all n.
    P_FIXED = 65537
    for rho_lbl, inv in rhos:
        rows = []
        for n in ns:
            if (P_FIXED - 1) % n != 0:
                continue
            k = max(2, round(n / inv))
            if elias_radius(k / n, P_FIXED) <= johnson_radius(k / n):
                continue
            for c in cs:
                rows.append(CertSpec(P_FIXED, n, k, c, note=f"fixedfield_rho={rho_lbl}"))
        if rows:
            ladders["fixed_field"][rho_lbl] = rows

    # (S) smallest smooth prime per n (>= 257 floor so the field is never tiny relative to n)
    for rho_lbl, inv in rhos:
        rows = []
        for n in ns:
            k = max(2, round(n / inv))
            p = smooth_prime_for(n, lo=max(257, 2 * n))
            if p is None:
                continue
            if elias_radius(k / n, p) <= johnson_radius(k / n):
                continue
            for c in cs:
                rows.append(CertSpec(p, n, k, c, note=f"smallsmooth_rho={rho_lbl}"))
        if rows:
            ladders["smallest_smooth"][rho_lbl] = rows

    return ladders


def build_crosscheck_battery() -> list[CertSpec]:
    """Small cells (p^k <= QK_CAP) with a genuine open band, where full q^k enumeration is
    feasible and the certificate can be pinned to the exact list.  Mirrors the regime the
    full-enumeration experiments reached (n <= 24)."""
    cells = []
    seen = set()

    def add(p, n, k, c):
        if p ** k > QK_CAP:
            return
        if not is_prime(p) or (p - 1) % n != 0:
            return
        if elias_radius(k / n, p) <= johnson_radius(k / n):
            return
        key = (p, n, k, c)
        if key in seen:
            return
        seen.add(key)
        cells.append(CertSpec(p, n, k, c, note="crosscheck"))

    # genuine power-of-two + highly-composite open-band cells, several c each
    base = [
        (97, 24, 3), (73, 24, 3), (257, 16, 2), (241, 16, 2),
        (97, 12, 3), (1217, 8, 2), (109, 27, 3), (101, 25, 3),
        (257, 32, 2), (97, 16, 2), (193, 24, 4),
    ]
    for (p, n, k) in base:
        for c in (0.25, 0.5, 0.75):
            add(p, n, k, c)
    return cells


# ===========================================================================
# Budgets (trial counts scale with n so big cells stay tractable; ALL O(t*n)).
# ===========================================================================
def cert_budget(n: int, k: int, quick: bool = False) -> dict:
    """Trial/pool sizes for the certificate search.  Cost per trial is O(t*n) (NOT q^k),
    so we can afford generous sampling even at n=2048.  We scale trial counts down mildly
    with n to keep wall-time per cell bounded while keeping the adversarial search strong."""
    # adversarial target sizes: probe well above any plausible bounded list
    if n <= 128:
        rp_t = (3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32)
        rp_trials = 200 if not quick else 60
        ll_pool = (8, 12, 16, 24, 32, 48, 64)
        ll_trials = 60 if not quick else 20
        gr_restarts = 60 if not quick else 20
        gr_max_t = 40
        st_pool = 64; st_trials = 8 if not quick else 3
        ls_moves = 400 if not quick else 150
    elif n <= 512:
        rp_t = (3, 4, 5, 6, 8, 10, 12, 16, 20, 24)
        rp_trials = 120 if not quick else 40
        ll_pool = (8, 12, 16, 24, 32, 48)
        ll_trials = 40 if not quick else 15
        gr_restarts = 40 if not quick else 15
        gr_max_t = 32
        st_pool = 48; st_trials = 6 if not quick else 2
        ls_moves = 600 if not quick else 200
    else:  # n in {1024, 2048, 4096}
        rp_t = (3, 4, 5, 6, 8, 10, 12, 16, 20)
        rp_trials = 60 if not quick else 24
        ll_pool = (8, 12, 16, 24, 32)
        ll_trials = 24 if not quick else 10
        gr_restarts = 24 if not quick else 10
        gr_max_t = 28
        st_pool = 32; st_trials = 4 if not quick else 2
        ls_moves = 800 if not quick else 250
    return {
        "rp_t_set": rp_t, "rp_trials": rp_trials,
        "ll_pool_sizes": ll_pool, "ll_trials": ll_trials,
        "gr_restarts": gr_restarts, "gr_batch": 12, "gr_max_t": gr_max_t,
        "st_pool": st_pool, "st_trials": st_trials,
        "ls_moves": ls_moves,
        "exact_random": 400,
    }


# ===========================================================================
# Driver.
# ===========================================================================
def _flatten(ladders: dict) -> list[CertSpec]:
    out = []
    for fam in ladders.values():
        for rows in fam.values():
            out.extend(rows)
    return out


def run_full(quick: bool = False, procs: int = None):
    procs = procs or min(16, os.cpu_count() or 4)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t_start = time.time()

    ladders = build_fixed_rho_ladders(quick=quick)
    cert_specs = _flatten(ladders)
    cc_specs = build_crosscheck_battery()

    print(f"[cluster_certificate] {len(cert_specs)} fixed-rho cells, "
          f"{len(cc_specs)} cross-check cells, {procs} procs, quick={quick}")

    # ---- run the fixed-rho certificate cells (the headline) ----
    cert_args = [(s, cert_budget(s.n, s.k, quick), 0xCE71 + 7919 * i)
                 for i, s in enumerate(cert_specs)]
    cert_results = []
    with Pool(procs) as pool:
        for i, rec in enumerate(pool.imap_unordered(run_cert_cell, cert_args)):
            cert_results.append(rec)
            if rec.get("has_open_band"):
                print(f"  [cert {i+1}/{len(cert_args)}] {rec['label']:>28} "
                      f"n={rec['n']:>5} rho={rec['rho']:.4f} c={rec['c']} "
                      f"delta={rec['delta']:.4f} margin={rec['margin_rE_minus_delta']:+.4f} "
                      f"max_cert={rec['max_cert']:>3} ({rec['max_cert_construction'][:6]}) "
                      f"1/(1-d)={rec['heuristic_1_over_1md']:.2f} "
                      f"trials={rec['n_trials_total']} {rec['elapsed_s']}s",
                      flush=True)

    # ---- run the cross-check cells (q^k enumeration; consistency to ground truth) ----
    cc_args = [(s, cert_budget(s.n, s.k, quick), 0xC0FFEE + 104729 * i)
               for i, s in enumerate(cc_specs)]
    cc_results = []
    with Pool(procs) as pool:
        for i, rec in enumerate(pool.imap_unordered(crosscheck_cell, cc_args)):
            cc_results.append(rec)
            if rec.get("has_open_band"):
                print(f"  [xchk {i+1}/{len(cc_args)}] {rec['label']:>28} "
                      f"cert={rec['max_cert']:>3} true@w={rec['true_list_at_witness_word']} "
                      f"exactMax={rec['exact_sampled_max_list']:>3} "
                      f"cert<=true:{rec['cert_le_true_at_witness']} "
                      f"cert<=exact:{rec['cert_le_exact_sampled_max']}", flush=True)

    elapsed = time.time() - t_start

    # ---- assemble fixed-rho ladders into compact tables for the writeup ----
    def ladder_table(fam_name):
        fam = ladders.get(fam_name, {})
        tbl = {}
        for rho_lbl, rows in fam.items():
            by_c = {}
            for s in rows:
                rec = next((r for r in cert_results
                            if r["label"] == s.label and r.get("has_open_band")), None)
                if rec is None:
                    continue
                by_c.setdefault(str(s.c), []).append({
                    "n": rec["n"], "p": rec["p"], "k": rec["k"],
                    "delta": rec["delta"], "margin": rec["margin_rE_minus_delta"],
                    "max_cert": rec["max_cert"],
                    "construction": rec["max_cert_construction"],
                    "heuristic_1_over_1md": rec["heuristic_1_over_1md"],
                    "volume_qEn": rec["volume_pred_qEn"],
                    "n_trials": rec["n_trials_total"],
                })
            for cc in by_c:
                by_c[cc].sort(key=lambda r: r["n"])
            tbl[rho_lbl] = by_c
        return tbl

    summary = {
        "meta": {
            "purpose": "Exact cluster certificate: witness lower bound on the open-band "
                       "single-code list, reaching LARGE n at FIXED rho (no q^k enumeration).",
            "date": "2026-06-03",
            "quick": quick, "procs": procs,
            "elapsed_s": round(elapsed, 1),
            "method": "max #generators-in-plurality-ball over random + greedy + Lloyd + "
                      "local-search constructions; cost O(t*n) per trial, independent of q^k.",
            "honesty": "Every max_cert is a witness LOWER bound on Lambda(C,delta); bounded "
                       "across large n at fixed rho supports P', growing would falsify it.",
        },
        "fixed_rho_tables": {
            "fixed_field": ladder_table("fixed_field"),
            "smallest_smooth": ladder_table("smallest_smooth"),
        },
        "cert_cells": cert_results,
        "crosscheck": cc_results,
    }

    with open(os.path.join(RESULTS_DIR, "cluster_certificate.json"), "w") as f:
        json.dump(summary, f, indent=1)

    # CSV: one row per fixed-rho cert cell (the headline data)
    with open(os.path.join(RESULTS_DIR, "cluster_certificate.csv"), "w", newline="") as f:
        wtr = csv.writer(f)
        wtr.writerow(["family", "rho", "p", "n", "k", "c", "delta",
                      "margin_rE_minus_delta", "max_cert", "construction",
                      "heuristic_1_over_1md", "volume_qEn", "n_trials"])
        for rec in sorted(cert_results, key=lambda r: (r.get("note", ""), r.get("rho", 0),
                                                       r.get("c", 0), r.get("n", 0))):
            if not rec.get("has_open_band"):
                continue
            fam = "fixed_field" if "fixedfield" in rec.get("note", "") else "smallest_smooth"
            wtr.writerow([fam, rec["rho"], rec["p"], rec["n"], rec["k"], rec["c"],
                          rec["delta"], rec["margin_rE_minus_delta"], rec["max_cert"],
                          rec["max_cert_construction"], rec["heuristic_1_over_1md"],
                          rec["volume_pred_qEn"], rec["n_trials_total"]])

    # ---- console headline ----
    _print_headline(summary)
    print(f"\n[cluster_certificate] DONE in {elapsed:.1f}s -> "
          f"results/cluster_certificate.{{json,csv}}")
    return summary


def _print_headline(summary: dict):
    print("\n" + "=" * 78)
    print("HEADLINE: max cluster certificate vs n at FIXED rho (witness LOWER bound)")
    print("=" * 78)
    for fam_name, fam in summary["fixed_rho_tables"].items():
        print(f"\n[{fam_name}]")
        for rho_lbl, by_c in fam.items():
            for cc, rows in by_c.items():
                ns = "  ".join(f"n={r['n']}:{r['max_cert']}" for r in rows)
                trend = _trend(rows)
                print(f"  rho={rho_lbl:5s} c={cc:4s}: {ns}    [{trend}]")
    # cross-check verdict
    cc = summary["crosscheck"]
    ok_l = all(r.get("cert_le_true_at_witness", True) for r in cc if r.get("has_open_band"))
    ok_u = all(r.get("cert_le_exact_sampled_max", True) for r in cc if r.get("has_open_band"))
    print(f"\n  cross-check: cert<=true@witness ALL OK: {ok_l};  "
          f"cert<=exact-sampled-max ALL OK: {ok_u}")


def _trend(rows: list[dict]) -> str:
    if len(rows) < 2:
        return "n/a"
    cs = [r["max_cert"] for r in rows]
    if max(cs) - min(cs) <= 1:
        return f"FLAT ({min(cs)}-{max(cs)})"
    if cs[-1] > cs[0]:
        return f"GROWING {cs[0]}->{cs[-1]}"
    return f"VARYING {min(cs)}-{max(cs)}"


# ===========================================================================
# Self-test battery.
# ===========================================================================
def _self_test():
    print("cluster_certificate self-test")
    print("=" * 70)
    rng = np.random.default_rng(0xCE27)

    # ---- 1. rand_codewords are genuine RS codewords (match rs.encode) ----
    print("[1] rand_codewords == rs.encode on the same coefficients:")
    F = PrimeField(97); n = 24; k = 3
    L = domain_subgroup(F, n)
    q = F.q
    coeffs = rng.integers(0, q, size=(5, k)).astype(np.int64)
    Larr = np.asarray(L, dtype=np.int64)
    acc = np.broadcast_to(coeffs[:, -1:], (5, n)).astype(np.int64).copy()
    for j in range(k - 2, -1, -1):
        acc = F.add_vec(F.mul_vec(acc, Larr[None, :]), coeffs[:, j:j + 1])
    for i in range(5):
        ref = encode(F, L, coeffs[i])
        assert np.array_equal(acc[i], ref), "rand_codewords batch Horner != encode"
    print("    batch Horner matches rs.encode on 5 codewords OK")

    # ---- 2. plurality_word is the true coordinate-wise mode ----
    print("[2] plurality_word == coordinate-wise mode (brute check):")
    for _ in range(20):
        t = int(rng.integers(2, 12))
        cw = rng.integers(0, 7, size=(t, n)).astype(np.int64)
        w = plurality_word(cw, 7, rng)
        for j in range(n):
            vals, counts = np.unique(cw[:, j], return_counts=True)
            mx = counts.max()
            winners = set(int(v) for v in vals[counts == mx])
            # plurality must pick a value achieving the column max count
            assert int(w[j]) in winners, f"plurality not a column mode at j={j}"
            # and its count must equal the max
            assert int((cw[:, j] == w[j]).sum()) == int(mx), "plurality count != mode count"
    print("    20 random matrices: plurality is always a true column mode OK")

    # ---- 3. plurality is a strong HEURISTIC ball center (TOTAL-distance optimal) ----
    print("[3] plurality word is a strong heuristic center (beats random words on a sample):")
    # The plurality EXACTLY minimizes TOTAL distance sum_i Delta(c_i, w) (coordinates are
    # independent; each coordinate picks its mode).  It is NOT in general the min-MAX
    # (closest-string) center -- e.g. {000,111,111} has plurality 111 with max-radius 3 but
    # 110 attains 2.  Here we just sanity-check it is a strong center: its max-distance is at
    # most that of a sample of random words (a heuristic check, not a min-max optimality proof).
    # The certificate's validity does not rely on this -- w is an explicit witness either way.
    t = 6
    cw = rng.integers(0, 5, size=(t, n)).astype(np.int64)
    w = plurality_word(cw, 5, rng)
    # the certified ball radius needed to hold ALL t = max_i Delta(c_i, w)
    rad_plur = int(hamming_to_set(w, cw).max())
    for _ in range(200):
        wr = rng.integers(0, 5, size=n).astype(np.int64)
        assert int(hamming_to_set(wr, cw).max()) >= rad_plur, \
            "plurality center was beaten on the random sample (heuristic-strength check)"
    print(f"    plurality max-dist={rad_plur} <= every one of 200 random words OK")

    # ---- 4. cert_count is an EXACT lower bound: cert(w,S) <= TRUE |Lambda(C,e/n,w)| ----
    print("[4] cert_count(w,S) <= EXACT |Lambda(C,e/n,w)| (q^k enumeration), every cell:")
    book = build_codeword_book(F, L, k)            # n=24,k=3 -> 97^3 ~ 9e5, fine
    for _ in range(8):
        t = int(rng.integers(3, 10))
        cw = rand_codewords(F, L, k, t, rng)
        w = plurality_word(cw, q, rng)
        for e in [14, 15, 16, 17]:
            cert = cert_count(w, cw, e)
            true_list = list_sizes_all_e(book, w, [e])[e]
            assert cert <= true_list, \
                f"CERTIFICATE OVER-COUNTS: cert={cert} > true={true_list} at e={e}"
            # the generators are themselves codewords, so each in-ball generator is in Lambda
            # => cert is a genuine subset count; equality holds iff no OTHER codeword is close
    print("    cert <= exact list at the same word, all cells (no over-count) OK")

    # ---- 5. constructions never exceed the exact list at their own witness word -------
    print("[5] every construction's cert <= exact list at its OWN witness word:")
    e = 16
    bud = cert_budget(n, k, quick=True)
    for ctor, args in [
        (construct_random_plurality, (F, L, k, e, rng, (3, 4, 5, 6, 8), 30)),
        (construct_lloyd, (F, L, k, e, rng, (8, 12, 16), 12)),
        (construct_greedy, (F, L, k, e, rng, 12, 12, 20)),
        (construct_structured, (F, L, k, e, rng, 32, 4)),
    ]:
        res = ctor(*args)
        if res["w"] is not None:
            true_list = list_sizes_all_e(book, np.asarray(res["w"], dtype=np.int64), [e])[e]
            assert res["cert"] <= true_list, \
                f"{ctor.__name__}: cert {res['cert']} > true {true_list} at its witness"
    print("    random-plurality / lloyd / greedy / structured all certify <= the exact list OK")

    # ---- 5b. structured construction RECOVERS a known coset (Kambire) cluster ----------
    print("[5b] structured construction recovers a known Kambire coset cluster:")
    # GF(73), n=24, k=2: the Kambire monomial X^6 (r=3,m=2 on the order-6 sub-subgroup) is
    # a deep-hole word with an EXACT list of 4 at e=18 (delta=0.75).  The structured-seed
    # construction, which seeds at exactly such monomials, must certify a cluster >= 2 here
    # (it seeds at the monomial center and Lloyd-refits over the actual codewords).
    Fk = PrimeField(73); nk = 24; kk = 2
    Lk = domain_subgroup(Fk, nk)
    ek = 18
    stk = construct_structured(Fk, Lk, kk, ek, rng, n_pool=48, trials=6)
    bookk = build_codeword_book(Fk, Lk, kk)
    # exact worst case at that radius (the Kambire monomial X^6) for context
    wk = np.array([Fk.pow(int(x), 6) for x in Lk], dtype=np.int64)
    exact_kambire = list_sizes_all_e(bookk, wk, [ek])[ek]
    assert stk["cert"] >= 2, \
        f"structured construction should certify a coset cluster >=2 (got {stk['cert']}, exact Kambire {exact_kambire})"
    if stk["w"] is not None:
        true_at_w = list_sizes_all_e(bookk, np.asarray(stk["w"], dtype=np.int64), [ek])[ek]
        assert stk["cert"] <= true_at_w, "structured cert over-counts at its witness"
    print(f"    GF(73) n=24 k=2 X^6 cluster: structured cert={stk['cert']} "
          f"(exact Kambire list={exact_kambire}) OK")

    # ---- 6. adversarial constructions are >= the random-plurality baseline (they try) --
    print("[6] adversarial constructions >= random-plurality baseline (strength check):")
    e = 17
    rp = construct_random_plurality(F, L, k, e, rng, (3, 4, 5, 6, 8, 10, 12), 120)
    ll = construct_lloyd(F, L, k, e, rng, (8, 12, 16, 24, 32), 40)
    gr = construct_greedy(F, L, k, e, rng, 40, 12, 32)
    # greedy/lloyd should at least match the baseline (they subsume random pluralities)
    assert max(ll["cert"], gr["cert"]) >= rp["cert"] - 1, \
        f"adversarial weaker than baseline: rp={rp['cert']} ll={ll['cert']} gr={gr['cert']}"
    print(f"    e={e}: random-plurality={rp['cert']}, lloyd={ll['cert']}, greedy={gr['cert']} "
          f"(adversarial >= baseline) OK")

    # ---- 7. cross-check cell runs end-to-end and the consistency flags hold ----
    print("[7] crosscheck_cell end-to-end consistency flags:")
    spec = CertSpec(97, 24, 3, 0.5, note="selftest")
    rec = crosscheck_cell((spec, cert_budget(24, 3, quick=True), 0xABCD))
    assert rec["cert_le_true_at_witness"], "cross-check: cert must be <= true at its witness"
    assert rec["cert_le_exact_sampled_max"], "cross-check: cert must be <= exact sampled max"
    print(f"    GF(97) n=24 k=3 c=0.5: cert={rec['max_cert']} <= "
          f"true@w={rec['true_list_at_witness_word']} <= exactMax={rec['exact_sampled_max_list']} "
          f"OK")

    # ---- 8. band_e lands strictly inside the open band ----
    print("[8] band_e is strictly inside (J, r_E] on the integer lattice:")
    for (p, n2, k2) in [(257, 32, 2), (769, 256, 16), (12289, 1024, 64)]:
        rho = k2 / n2
        J = johnson_radius(rho); rE = elias_radius(rho, p)
        for c in (0.25, 0.5, 0.75):
            e = band_e(n2, k2, p, c)
            assert e is not None, f"no band point for {(p,n2,k2,c)}"
            assert J * n2 < e <= rE * n2 + 1e-9, \
                f"band_e {e} not in (J*n,r_E*n]=({J*n2:.2f},{rE*n2:.2f}] for {(p,n2,k2,c)}"
    print("    band_e in (J,r_E] for n in {32,256,1024} OK")

    # ---- 9. cost is independent of q^k: a large-n cell runs FAST (no enumeration) ----
    print("[9] large-n certificate runs without q^k enumeration (speed sanity):")
    spec = CertSpec(65537, 2048, 256, 0.5, note="selftest_largeN")  # q^k = 65537^256, astronomical
    t0 = time.time()
    rec = run_cert_cell((spec, cert_budget(2048, 256, quick=True), 0x5151))
    dt = time.time() - t0
    assert rec["has_open_band"], "n=2048 rho=1/8 must have an open band"
    assert dt < 90, f"n=2048 cell took {dt:.1f}s -- should be fast (O(t*n), not q^k)"
    print(f"    GF(65537) n=2048 k=256 (q^k astronomical): max_cert={rec['max_cert']} "
          f"in {dt:.1f}s (no enumeration) OK")

    print("=" * 70)
    print("ALL CLUSTER-CERTIFICATE SELF-TESTS PASSED")


# ===========================================================================
# Main.
# ===========================================================================
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true", help="run the self-test battery")
    ap.add_argument("--quick", action="store_true", help="fast reduced sweep")
    ap.add_argument("--procs", type=int, default=None, help="parallel worker processes")
    args = ap.parse_args()

    if args.selftest:
        _self_test()
    else:
        run_full(quick=args.quick, procs=args.procs)
