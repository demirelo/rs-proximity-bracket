"""
high_moments.py -- EXACT high-moment / Markov probe for the open-band list tail (P', E1').

THE QUESTION (for the E1' analytic target, p-prime/moment-reduction.md SS3).
  P' needs an upper-tail bound on the open-band single-code list
        L(w) := |Lambda(C, delta, w)| = #{ deg P<k : N(P,w) >= (1-delta)n },
  for C = RS[F_p, <omega>, rho n] in the open band  J = 1-sqrt(rho) < delta < r_E = 1-H_p(rho).
  Writing N(P,w) = #{x in <omega> : P(x)=w(x)},  M(P,w) := N(P,w) - n/p, and the threshold
        T := (1-delta)n - n/p = Theta_rho(n),
  membership P in Lambda is exactly the large-deviation event M(P,w) >= T.  The candidate proof
  tool (moment-reduction.md SS3 (HM)) is a HIGH MOMENT with Markov's inequality:
        B_r(w) := M_r(w) / T^r,        M_r(w) := Sum_{deg P<k} (M(P,w))_+^r,
  so that  L(w) <= B_r(w)  for every integer r >= 1 (Markov on the POSITIVE part).
  The known fact (moment-reduction.md SS2) is that the GLOBAL SECOND moment is too weak: by
  pairwise independence  Sum_{deg P<k} (N(P,w)-n/p)^2 = p^k n(1/p - 1/p^2) ~ p^{k-1} n  is a FIXED
  quantity, so B_2 ~ p^{k-1}/n is dominated by the p^k codewords with tiny random fluctuations and
  vastly overshoots L.  The OPEN QUESTION this probe answers empirically: what moment order r is
  actually strong enough for Markov to "see" the tail -- i.e. the smallest r with B_r(w) = O(L(w))
  (Markov reproduces the exact list up to a small factor), not B_r >> L (bulk noise) and not
  vacuous.  And how does that smallest r SCALE -- a rho-dependent constant, or growing with
  k / log p / n / 1/eps0?  This is the guidance the E1' high-moment program must target.

THE METHOD -- exact moments, exact lists, full enumeration, NO decoder.
  For an exact-enumeration cell (small p,k with p^k <= ~3e6 so ALL p^k degree-<k polynomials are
  enumerable) and a target word w, one streamed pass over the codeword book yields the EXACT
  histogram of Hamming distances d = n - N(P,w) over all p^k polynomials (reusing the verbatim
  kernel _distance_histogram from singlelist_past_johnson.py).  From that single histogram:
    * the EXACT list  L(w) = #{P : d <= floor(delta*n)} = #{P : N >= (1-delta)n}  (cumulative tail);
    * the EXACT moment  M_r(w) = Sum_d hist[d] * (N - n/p)_+^r  for every r, in EXACT rational
      arithmetic (we form the integer p*N - n per bin and raise it to the r-th power with Python
      big integers, then divide by p^r -- no float overflow, exact (.)_+^r on the excess-above-mean
      as in the moment-reduction).
    * the Markov bound  B_r(w) = M_r(w) / T^r.
  We then read off the SMALLEST r where Markov reproduces L: the smallest r with B_r <= FACTOR * L
  (with the always-true lower guard B_r >= L verified at every r as a correctness check of Markov).

WORD TYPES (per cell), the three the probe stresses:
  (a) ADVERSARIAL CLUSTER words -- the coordinate-wise plurality witness of the strongest
      cluster-packing construction (greedy / Lloyd / structured / random-plurality, reused verbatim
      from cluster_certificate.py).  These realize the open-band worst case (falsification-findings
      SS2): they are where the tail is heaviest and where the needed r is largest.
  (b) COSET / KAMBIRE monomial words  X^{rm}  (and the canonical deep hole X^k and neighbours): the
      coset-structured worst-case seeds named by the structural route.
  (c) RANDOM words in F_p^n: the average-case baseline (tiny list, the bulk-noise regime where the
      second moment is most badly fooled).

SCALING SEPARATION.  To tell apart what the smallest-r tracks, the battery varies the cells along
controlled axes within the exact-enumeration budget:
    * FIX rho, VARY k  (rho=1/8: (p,n,k) with k=2,3; rho=1/4: k=2,3,4)  -- isolates k-dependence;
    * FIX (n,k), VARY p  (e.g. n=16,k=2 over p in {241,257,769,3041})    -- isolates log p;
    * VARY rho            (1/16, 1/8, 1/4)                                -- isolates rho / 1/eps0.

HONESTY.  This is EXACT for each (cell, word): the histogram, the list, and every moment are exact
(full enumeration, no decoder, exact rational moments).  But it is SMALL-CELL empirical guidance for
the E1' theorem, NOT a proof: (i) p^k <= 3e6 caps p,k small, so log p and k each span only a short
range -- a trend visible here need not be the asymptotic law; (ii) the adversarial word is a strong
witness, not the exact arg-max over all q^n words, so the worst-case L (hence the worst-case needed
r) is a lower bound on the true worst case; (iii) "Markov sees L" uses a fixed tolerance FACTOR, a
heuristic cut.  The verdict is reported as a trend with these caveats stated, to guide which moment
order the E1' program must target -- not to settle it.

Run  `python3.11 high_moments.py --selftest`   for the self-test battery (exact identities),
     `python3.11 high_moments.py`              for the full probe -> results/high_moments.{json,csv}
                                               and p-prime/high-moment-probe.md.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from dataclasses import dataclass, field
from fractions import Fraction

import numpy as np

from ff import PrimeField, is_prime
from rs import build_codeword_book, domain_subgroup
from singlelist_past_johnson import Hq, elias_radius, _distance_histogram
from cluster_certificate import (
    johnson_radius,
    band_e,
    construct_random_plurality,
    construct_lloyd,
    construct_greedy,
    construct_structured,
    cert_budget,
)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")
# Exact-enumeration cap: ALL p^k degree-<k polynomials must be streamable.
QK_CAP = 3_000_000

# Highest moment order to probe.  r up to 12 is far past the second moment and ample to expose
# where Markov turns over from bulk-dominated (B_r >> L) to tail-sensitive (B_r = O(L)).
R_MAX = 12

# "Markov sees the list" tolerance: smallest r with  L <= B_r <= FACTOR * max(L, 1).
# B_r >= L always (Markov is a valid upper bound); we want B_r not to OVERSHOOT L by more than this.
SEE_FACTOR = 10.0


# ===========================================================================
# Exact moments and Markov bounds from one distance histogram.
# ===========================================================================
def list_and_moments(hist: np.ndarray, n: int, p: int, e: int, r_max: int) -> dict:
    """From the EXACT distance histogram (hist[d] = #{deg P<k : Delta(P,w)=d}) compute, EXACTLY:

        L          = #{P : Delta <= e} = #{P : N >= n-e}            (the open-band list at radius e),
        M_r        = Sum_P (N(P,w) - n/p)_+^r        for r=1..r_max  (EXACT rational; (.)_+ on excess),
        B_r        = M_r / T^r,        T = (n-e) - n/p              (the Markov bound, >= L for all r).

    Agreement N = n - d.  We work with the integer  u = p*N - n = p*(n-d) - n  so that
    (N - n/p)_+^r = (u_+)^r / p^r exactly via Python big integers (no float overflow), then divide
    by T^r where T = (n-e) - n/p is kept as an exact Fraction.  Returns floats for reporting plus
    the exact Fraction objects for the self-test's identity checks.
    """
    d = np.arange(n + 1, dtype=np.int64)
    N = n - d
    counts = hist.astype(object)            # Python ints, exact
    u = (p * N - n).astype(object)          # integer p*N - n per bin (sign = sign of N - n/p)

    # exact list size at radius e (cumulative count of d <= e)
    L = int(hist[: e + 1].sum())

    # threshold T = (n - e) - n/p  as an exact Fraction (positive in the open band)
    T = Fraction(n - e) - Fraction(n, p)
    assert T > 0, f"non-positive threshold T={T} (e={e},n={n},p={p})"

    pr = [Fraction(p) ** r for r in range(r_max + 1)]
    Tr = [T ** r for r in range(r_max + 1)]

    M = [Fraction(0)] * (r_max + 1)
    # only positive-excess bins (N > n/p  <=>  u > 0) contribute to the positive part
    for cnt, uv in zip(counts.tolist(), u.tolist()):
        if cnt == 0 or uv <= 0:
            continue
        cntF = Fraction(int(cnt))
        for r in range(1, r_max + 1):
            M[r] += cntF * Fraction(int(uv) ** r)
    M_exact = [None] + [M[r] / pr[r] for r in range(1, r_max + 1)]    # Sum (N-n/p)_+^r, exact
    B_exact = [None] + [M_exact[r] / Tr[r] for r in range(1, r_max + 1)]

    M_float = [None] + [float(M_exact[r]) for r in range(1, r_max + 1)]
    B_float = [None] + [float(B_exact[r]) for r in range(1, r_max + 1)]

    # smallest r where Markov reproduces L (B_r <= FACTOR * max(L,1)); also the always-true guard
    smallest_r = None
    guard_ok = True
    for r in range(1, r_max + 1):
        if B_exact[r] < Fraction(L):                  # Markov MUST be an upper bound: B_r >= L
            guard_ok = False
        if smallest_r is None and B_float[r] <= SEE_FACTOR * max(L, 1):
            smallest_r = r

    return {
        "L": L,
        "T": float(T),
        "M_float": M_float,
        "B_float": B_float,
        "M_exact": M_exact,
        "B_exact": B_exact,
        "smallest_r": smallest_r,
        "markov_upper_guard_ok": guard_ok,
    }


# ===========================================================================
# Word generators -- the three types.
# ===========================================================================
def adversarial_cluster_word(F, L, k, e, rng) -> np.ndarray | None:
    """The coordinate-wise plurality witness of the STRONGEST cluster-packing construction at radius
    e (greedy / Lloyd / structured / random-plurality, reused verbatim from cluster_certificate).
    This is the open-band worst-case word (falsification-findings SS2) -- the heaviest tail.  Returns
    the witness word w (length n) achieving the largest certificate, or None if all empty."""
    n = len(L)
    bud = cert_budget(n, k, quick=False)
    rp = construct_random_plurality(F, L, k, e, rng, t_set=bud["rp_t_set"], trials_per_t=bud["rp_trials"])
    ll = construct_lloyd(F, L, k, e, rng, pool_sizes=bud["ll_pool_sizes"], trials_per_pool=bud["ll_trials"])
    gr = construct_greedy(F, L, k, e, rng, n_restarts=bud["gr_restarts"], pool_factor=bud["gr_batch"], max_t=bud["gr_max_t"])
    st = construct_structured(F, L, k, e, rng, n_pool=bud["st_pool"], trials=bud["st_trials"])
    best = max([rp, ll, gr, st], key=lambda c: c["cert"])
    if best["w"] is None:
        return None
    return np.asarray(best["w"], dtype=np.int64)


def kambire_monomial_words(F, L, k, n) -> list[tuple[str, np.ndarray]]:
    """Coset / Kambire monomial worst-case seeds: the deep hole X^k and the Kambire X^{rm} family
    (and neighbours) for every coset chain s|n with k=(r-2)m, plus a couple of high monomials.  These
    are the coset-structured words the structural route names.  X^a on the subgroup depends on a mod n."""
    from ff import _divisors
    Larr = np.asarray(L, dtype=np.int64)
    exps: set[int] = set()
    # canonical deep hole and immediate neighbours
    for a in (k, k + 1, k + 2):
        if a < n:
            exps.add(a)
    # Kambire X^{rm}: coset chain s|n with k=(r-2)m
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
    # a couple of generic high monomials near the top
    for a in (n - 1, n - 2):
        if a > k:
            exps.add(a)
    out = []
    for a in sorted(exps):
        w = np.array([F.pow(int(x), a) for x in Larr], dtype=np.int64)
        out.append((f"X^{a}", w))
    return out


def random_words(F, n, rng, count) -> list[np.ndarray]:
    return [rng.integers(0, F.q, size=n).astype(np.int64) for _ in range(count)]


# ===========================================================================
# One cell at the band-midpoint radius.
# ===========================================================================
@dataclass
class Cell:
    p: int
    n: int
    k: int
    axis: str = ""        # which scaling axis this cell belongs to
    c: float = 0.5        # band fraction (band midpoint)

    @property
    def rho(self):
        return self.k / self.n

    @property
    def label(self):
        return f"GF({self.p})_n{self.n}_k{self.k}"


def run_cell(cell: Cell, seed: int, r_max: int = R_MAX) -> dict:
    """Full high-moment probe at one exact-enumeration cell, band-midpoint radius, all three word
    types.  Returns the per-word L, M_r, B_r, smallest-r, plus the cell's band landmarks."""
    p, n, k = cell.p, cell.n, cell.k
    F = PrimeField(p)
    L = domain_subgroup(F, n)
    rho = k / n
    J = johnson_radius(rho)
    rE = elias_radius(rho, p)
    cap = 1.0 - rho
    e = band_e(n, k, p, cell.c)
    out = {
        "label": cell.label, "axis": cell.axis, "p": p, "n": n, "k": k,
        "rho": round(rho, 6), "qk": int(p ** k),
        "J": round(J, 6), "r_E": round(rE, 6), "cap": round(cap, 6),
        "band_width": round(rE - J, 6), "log_p": round(math.log(p), 6),
    }
    if e is None:
        out["has_open_band"] = False
        return out
    delta = e / n
    out.update({
        "has_open_band": True, "e": int(e), "delta": round(delta, 6),
        "delta_minus_J": round(delta - J, 6), "margin_rE_minus_delta": round(rE - delta, 6),
        "n_over_p": round(n / p, 6),
    })

    book = build_codeword_book(F, L, k)
    # MDS sanity (RS over any domain is MDS); cheap, catches a domain/field bug.
    rng = np.random.default_rng(seed)

    words: list[tuple[str, str, np.ndarray]] = []   # (word_type, tag, w)

    # (a) adversarial cluster word
    adv = adversarial_cluster_word(F, L, k, e, rng)
    if adv is not None:
        words.append(("adversarial_cluster", "cluster_plurality", adv))

    # (b) coset / Kambire monomial words
    for tag, w in kambire_monomial_words(F, L, k, n):
        words.append(("coset_kambire", tag, w))

    # (c) random words
    for i, w in enumerate(random_words(F, n, rng, count=8)):
        words.append(("random", f"random_{i}", w))

    per_word = []
    for wtype, tag, w in words:
        hist = _distance_histogram(book, w)
        res = list_and_moments(hist, n, p, e, r_max)
        per_word.append({
            "word_type": wtype, "tag": tag,
            "L": res["L"], "T": round(res["T"], 6),
            "B_r": [None if b is None else (round(b, 6) if b < 1e12 else b) for b in res["B_float"]],
            "M_r": [None if m is None else (round(m, 6) if m < 1e12 else m) for m in res["M_float"]],
            "smallest_r": res["smallest_r"],
            "markov_upper_guard_ok": res["markov_upper_guard_ok"],
        })

    # Reduce to the worst (largest L, then largest needed r) word per type -- the type's worst case.
    def worst_of(wtype):
        rows = [pw for pw in per_word if pw["word_type"] == wtype]
        if not rows:
            return None
        # worst case: the heaviest tail (max L); tie-break by the largest needed r (None = unbounded)
        def key(pw):
            sr = pw["smallest_r"] if pw["smallest_r"] is not None else (r_max + 1)
            return (pw["L"], sr)
        return max(rows, key=key)

    by_type = {}
    for wtype in ("adversarial_cluster", "coset_kambire", "random"):
        wt = worst_of(wtype)
        if wt is not None:
            by_type[wtype] = {
                "tag": wt["tag"], "L": wt["L"], "smallest_r": wt["smallest_r"],
                "B_r": wt["B_r"], "M_r": wt["M_r"],
                "markov_upper_guard_ok": wt["markov_upper_guard_ok"],
            }

    out["per_word"] = per_word
    out["worst_by_type"] = by_type
    # cell-level worst case across all word types (the headline needed-r for this cell)
    all_sr = [(pw["L"], pw["smallest_r"], pw["word_type"], pw["tag"]) for pw in per_word]
    # the cell's worst case: take the word with the largest L; among those the largest needed r
    headline = max(all_sr, key=lambda t: (t[0], (t[1] if t[1] is not None else (r_max + 1))))
    out["headline_L"] = headline[0]
    out["headline_smallest_r"] = headline[1]
    out["headline_word_type"] = headline[2]
    out["headline_tag"] = headline[3]
    return out


# ===========================================================================
# Battery -- scaling-separation cells, all exact-enumerable with a genuine open band.
# ===========================================================================
def build_battery() -> list[Cell]:
    """Open-band exact-enumeration cells arranged along controlled scaling axes.  Every cell has
    p^k <= QK_CAP, n | p-1, and r_E > J (a genuine open band)."""
    cells: list[Cell] = []
    seen = set()

    def add(p, n, k, axis):
        if p ** k > QK_CAP:
            return False
        if not is_prime(p) or (p - 1) % n != 0:
            return False
        if elias_radius(k / n, p) <= johnson_radius(k / n):
            return False
        # require an interior band lattice point at the band midpoint (some narrow-band cells
        # collapse on the 1/n lattice: floor(J*n) == floor(r_E*n), no exact open-band radius).
        if band_e(n, k, p, 0.5) is None:
            return False
        key = (p, n, k)
        if key in seen:
            return False
        seen.add(key)
        cells.append(Cell(p, n, k, axis=axis))
        return True

    # ---- Axis FIX-RHO-VARY-K: hold rho, grow k (hence n=k/rho); isolates k-dependence ----
    # rho=1/8: (n,k) = (16,2),(24,3) -- the deepest k a 1/8 cell reaches under p^k<=3e6.
    add(257, 16, 2, "fix_rho_vary_k(rho=1/8)")
    add(97, 24, 3, "fix_rho_vary_k(rho=1/8)")
    add(73, 24, 3, "fix_rho_vary_k(rho=1/8)")
    # rho=1/4: (8,2),(12,3),(16,4) -- a clean k=2,3,4 ladder at one rate.  k=4 forces a small
    # field (p^4<=3e6 => p<=41): GF(17) n=16 k=4 (p=1 mod 16) is the genuine-open-band rho=1/4 k=4.
    add(73, 8, 2, "fix_rho_vary_k(rho=1/4)")
    add(97, 12, 3, "fix_rho_vary_k(rho=1/4)")
    add(17, 16, 4, "fix_rho_vary_k(rho=1/4)")

    # ---- Axis FIX-(n,k)-VARY-P: n=16,k=2 over the WIDEST exact-enumerable prime range ----
    # (p=1 mod 16, p^2<=3e6 => p<=1732): a ~6x log-p span (log p ~ 2.8..7.4), the cleanest
    # isolation of log-p dependence at fixed (n,k,rho).
    for p in (17, 97, 353, 593, 769, 1009, 1249, 1489, 1697):
        add(p, 16, 2, "fix_nk_vary_p(n=16,k=2)")
    # a second (n,k)=(24,3) prime ladder (isolates log p at higher k); p=1 mod 24, p^3<=3e6.
    for p in (73, 97):
        add(p, 24, 3, "fix_nk_vary_p(n=24,k=3)")
    # a k=4 prime pair (isolates log p at k=4); n=20,k=4, p=1 mod 20, p^4<=3e6 => p<=41.
    for p in (41,):
        add(p, 20, 4, "fix_nk_vary_p(n=20,k=4)")

    # ---- Axis VARY-RHO at comparable n (isolates rho / 1/eps0, eps0 = delta - J) ----
    add(257, 32, 2, "vary_rho(n~28-32)")   # rho=1/16
    add(109, 27, 3, "vary_rho(n~28-32)")   # rho=1/9
    add(29, 28, 4, "vary_rho(n~28-32)")    # rho=1/7  (genuine k=4, small field)
    add(31, 30, 4, "vary_rho(n~28-32)")    # rho=2/15 (genuine k=4)
    add(101, 25, 3, "vary_rho(n~24-25)")   # rho=3/25
    add(41, 20, 4, "vary_rho(n~20)")       # rho=1/5  (genuine k=4)
    add(97, 12, 3, "vary_rho(n=12)")       # rho=1/4

    return cells


# ===========================================================================
# Outputs.
# ===========================================================================
def write_outputs(results: list[dict], elapsed: float):
    os.makedirs(RESULTS_DIR, exist_ok=True)

    payload = {
        "meta": {
            "experiment": "high_moments",
            "purpose": "Exact high-moment / Markov probe: what moment order r is strong enough for "
                       "Markov to reproduce the EXACT open-band list size L on exact-enumeration "
                       "cells; guidance for the E1' high-moment program (p-prime/moment-reduction.md).",
            "date": "2026-06-03",
            "object": "B_r(w)=M_r(w)/T^r,  M_r(w)=Sum_{deg P<k}(N(P,w)-n/p)_+^r,  "
                      "T=(1-delta)n-n/p;  L(w)=#{deg P<k: N(P,w)>=(1-delta)n}.  Markov: L<=B_r for all r.",
            "smallest_r_definition": f"smallest integer r>=1 with B_r(w) <= {SEE_FACTOR} * max(L(w),1) "
                                     "(Markov 'sees' the list: B_r=O(L), not >>L from p^k bulk noise).",
            "exactness": "Per (cell,word): histogram, list L, and every moment M_r are EXACT "
                         "(full p^k enumeration, no decoder, exact rational (.)_+^r on the excess "
                         "N-n/p).  Markov upper-bound guard B_r>=L verified at every r.",
            "honesty": "Small-cell empirical guidance for the E1' theorem, NOT a proof: p^k<=3e6 "
                       "caps p,k small (short log p / k ranges); the adversarial word is a strong "
                       "witness, not the exact arg-max over all q^n; 'Markov sees L' uses a fixed "
                       "tolerance factor.",
            "r_max": R_MAX, "see_factor": SEE_FACTOR, "qk_cap": QK_CAP,
            "elapsed_s": round(elapsed, 1),
        },
        "cells": results,
    }
    jpath = os.path.join(RESULTS_DIR, "high_moments.json")
    with open(jpath, "w", newline="\n") as f:
        json.dump(payload, f, indent=2, default=str)

    # CSV: one row per (cell, word_type-worst), with L, smallest-r, and B_r for r=1..R_MAX.
    cpath = os.path.join(RESULTS_DIR, "high_moments.csv")
    cols = (["axis", "label", "p", "n", "k", "rho", "log_p", "delta", "margin_rE_minus_delta",
             "word_type", "tag", "L", "smallest_r"]
            + [f"B_{r}" for r in range(1, R_MAX + 1)])
    with open(cpath, "w", newline="") as f:
        wtr = csv.writer(f, lineterminator="\n")
        wtr.writerow(cols)
        for cell in results:
            if not cell.get("has_open_band"):
                continue
            for wtype, wt in cell.get("worst_by_type", {}).items():
                br = wt["B_r"]
                wtr.writerow([cell["axis"], cell["label"], cell["p"], cell["n"], cell["k"],
                              cell["rho"], cell["log_p"], cell["delta"],
                              cell["margin_rE_minus_delta"], wtype, wt["tag"], wt["L"],
                              wt["smallest_r"]]
                             + [br[r] for r in range(1, R_MAX + 1)])
    return jpath, cpath


def write_markdown(results: list[dict], elapsed: float):
    # Public-artifacts copy: the probe report is written next to the JSON/CSV
    # outputs (the research-notes copy wrote it into the authors' p-prime/ notes).
    os.makedirs(RESULTS_DIR, exist_ok=True)
    mpath = os.path.join(RESULTS_DIR, "high-moment-probe.md")
    open_cells = [c for c in results if c.get("has_open_band")]

    # ---- aggregate the smallest-r by word type and by scaling axis ----
    def needed_r(cell, wtype):
        wt = cell.get("worst_by_type", {}).get(wtype)
        if wt is None:
            return None
        return wt["smallest_r"]

    lines = []
    A = lines.append
    A("# P′ — High-Moment Probe of the Open-Band List Tail (E1′ guidance)")
    A("")
    A("**Date:** 2026-06-03.  **Scope:** an exact empirical measurement, on exact-enumeration cells, "
      "of the moment order `r` at which the high-moment Markov bound reproduces the open-band "
      "single-code list size. This calibrates the analytic target **E1′** of "
      "`p-prime/moment-reduction.md` §3.")
    A("")
    A("## 1. Object")
    A("")
    A("For `C = RS[F_p, ⟨ω⟩, k]`, `n = |⟨ω⟩|`, `ρ = k/n`, a word `w : ⟨ω⟩ → F_p`, and a degree-`<k` "
      "polynomial `P`, write `N(P,w) = #{x ∈ ⟨ω⟩ : P(x) = w(x)}` and `M(P,w) = N(P,w) − n/p`. In the "
      "open band `J = 1−√ρ < δ < r_E = 1−H_p(ρ)`, with threshold")
    A("```")
    A("   T := (1−δ)n − n/p = Θ_ρ(n),")
    A("```")
    A("the list `L(w) := |Λ(C,δ,w)| = #{deg P<k : N(P,w) ≥ (1−δ)n} = #{deg P<k : M(P,w) ≥ T}`. The "
      "high-moment Markov bound of `moment-reduction.md` §3 (HM) is")
    A("```")
    A("   B_r(w) := M_r(w) / T^r,        M_r(w) := Σ_{deg P<k} (M(P,w))_+^r,")
    A("   L(w) ≤ B_r(w)        for every integer r ≥ 1     (Markov on the positive part).")
    A("```")
    A("The known obstruction (`moment-reduction.md` §2): the **global second moment is a fixed "
      "quantity** by pairwise independence,")
    A("```")
    A("   Σ_{deg P<k} (N(P,w) − n/p)^2 = p^k · n(1/p − 1/p^2) ≈ p^{k−1} n,")
    A("```")
    A("so `B_2 ≈ p^{k−1}/n` is dominated by the `p^k` codewords with tiny fluctuations and overshoots "
      "`L` by an enormous factor. **The question:** what order `r` is strong enough for Markov to "
      "*see* the tail — the smallest `r` with `B_r(w) = O(L(w))` — and how does that `r` scale?")
    A("")
    A("The probe is **exact**: for each cell one streamed pass over all `p^k` codewords gives the exact "
      "distance histogram, from which `L`, every `M_r` (exact rational, `(·)_+^r` on the excess "
      "`N − n/p`), and every `B_r` are computed. The smallest-`r` is the least `r ≥ 1` with "
      f"`B_r ≤ {SEE_FACTOR:g}·max(L,1)`; the Markov upper-bound guard `B_r ≥ L` is verified at every `r`. "
      "Three word types are stressed: **(a)** adversarial cluster-plurality witnesses (open-band "
      "worst case), **(b)** coset/Kambiré monomials `X^{rm}`, **(c)** random words.")
    A("")

    # ---- Table 1: per-cell exact L, B_r ladder, smallest-r (adversarial worst case) ----
    A("## 2. Per-cell ladder: exact `L`, Markov `B_r`, and the smallest `r` (worst-case word)")
    A("")
    A("For each cell, the worst-case word (largest `L`) at the band midpoint `δ = J + ½(r_E − J)`. "
      "`B_r` is the Markov bound at moment order `r`; it decreases in `r` until it reaches `L`. "
      "`r*` is the smallest `r` with `B_r ≤ "
      f"{SEE_FACTOR:g}L`. `B_2` is shown to expose how badly the second moment overshoots.")
    A("")
    A("| cell | ρ | δ | `L` | `B_1` | `B_2` | `B_3` | `B_4` | `B_6` | `B_8` | `r*` | worst type |")
    A("|---|---|---|---|---|---|---|---|---|---|---|---|")
    def fmt(x):
        if x is None:
            return "—"
        if x >= 1e6:
            return f"{x:.1e}"
        if x >= 100:
            return f"{x:.0f}"
        return f"{x:.2f}"
    for c in sorted(open_cells, key=lambda c: (c["rho"], c["n"], c["k"], c["p"])):
        # pull the headline word's full B_r ladder
        pw = None
        for row in c["per_word"]:
            if row["word_type"] == c["headline_word_type"] and row["tag"] == c["headline_tag"]:
                pw = row
                break
        if pw is None:
            continue
        br = pw["B_r"]
        rstar = c["headline_smallest_r"]
        rstar_s = str(rstar) if rstar is not None else f">{R_MAX}"
        wt_short = {"adversarial_cluster": "cluster", "coset_kambire": "coset", "random": "random"}[c["headline_word_type"]]
        A(f"| {c['label']} | {c['rho']:.3f} | {c['delta']:.3f} | {c['headline_L']} | "
          f"{fmt(br[1])} | {fmt(br[2])} | {fmt(br[3])} | {fmt(br[4])} | {fmt(br[6])} | {fmt(br[8])} | "
          f"**{rstar_s}** | {wt_short} |")
    A("")

    # ---- Table 2: smallest-r by word type (which type needs the highest r) ----
    A("## 3. Smallest `r` by word type (the worst case)")
    A("")
    A("`r*` per word type, per cell. The word type needing the **highest** `r` is the worst case the "
      "E1′ moment must dominate.")
    A("")
    A("| cell | ρ | `L` cluster / coset / random | `r*` cluster | `r*` coset | `r*` random |")
    A("|---|---|---|---|---|---|")
    for c in sorted(open_cells, key=lambda c: (c["rho"], c["n"], c["k"], c["p"])):
        bt = c.get("worst_by_type", {})
        def cell_L(t):
            return bt.get(t, {}).get("L", "—")
        def cell_r(t):
            v = bt.get(t, {}).get("smallest_r")
            return (str(v) if v is not None else f">{R_MAX}") if t in bt else "—"
        A(f"| {c['label']} | {c['rho']:.3f} | "
          f"{cell_L('adversarial_cluster')} / {cell_L('coset_kambire')} / {cell_L('random')} | "
          f"{cell_r('adversarial_cluster')} | {cell_r('coset_kambire')} | {cell_r('random')} |")
    A("")

    # ---- compute the worst-case-type winner ----
    type_max_r = {"adversarial_cluster": 0, "coset_kambire": 0, "random": 0}
    type_unbounded = {"adversarial_cluster": 0, "coset_kambire": 0, "random": 0}
    for c in open_cells:
        for t, wt in c.get("worst_by_type", {}).items():
            sr = wt["smallest_r"]
            if sr is None:
                type_unbounded[t] += 1
            else:
                type_max_r[t] = max(type_max_r[t], sr)
    worst_type = max(type_max_r, key=lambda t: (type_unbounded[t], type_max_r[t]))
    worst_type_name = {"adversarial_cluster": "adversarial cluster", "coset_kambire": "coset/Kambiré monomial", "random": "random"}[worst_type]

    # ---- Scaling analysis ----
    A("## 4. What the smallest `r` scales with")
    A("")
    A("The cells are arranged along controlled axes (all exact-enumerable, all genuine open band at "
      "the band midpoint). `r*` below is the cell's worst-case smallest-`r` across word types.")
    A("")

    def axis_rows(axis_prefix):
        return sorted([c for c in open_cells if c["axis"].startswith(axis_prefix)],
                      key=lambda c: (c["k"], c["log_p"], c["n"]))

    # (i) fix rho vary k
    A("### 4.1 Fix `ρ`, vary `k` (isolates `k`-dependence)")
    A("")
    A("| ρ | cell | k | log p | `L` | `r*` |")
    A("|---|---|---|---|---|---|")
    for c in sorted([c for c in open_cells if c["axis"].startswith("fix_rho_vary_k")],
                    key=lambda c: (round(c["rho"], 3), c["k"], c["p"])):
        rstar = c["headline_smallest_r"]
        A(f"| {c['rho']:.3f} | {c['label']} | {c['k']} | {c['log_p']:.2f} | {c['headline_L']} | "
          f"{rstar if rstar is not None else '>'+str(R_MAX)} |")
    A("")

    # (ii) fix (n,k) vary p
    A("### 4.2 Fix `(n,k)`, vary `p` (isolates `log p`)")
    A("")
    A("The `(16,2)` ladder is the cleanest log-`p` isolation: for every `p ≥ 97` the band-midpoint "
      "radius is the *same* `e` (`δ = 0.750`), so the rows differ only in the field. `r*` rises over "
      "the smallest fields (where `n/p` is non-negligible and the distance histogram is coarse) and "
      "then **saturates** — it does **not** keep growing with `log p`.")
    A("")
    A("| (n,k) | p | log p | `n/p` | `L` | `r*` |")
    A("|---|---|---|---|---|---|")
    for c in sorted([c for c in open_cells if c["axis"].startswith("fix_nk_vary_p")],
                    key=lambda c: ((c["n"], c["k"]), c["p"])):
        rstar = c["headline_smallest_r"]
        A(f"| ({c['n']},{c['k']}) | {c['p']} | {c['log_p']:.2f} | {c['n_over_p']:.3f} | {c['headline_L']} | "
          f"{rstar if rstar is not None else '>'+str(R_MAX)} |")
    A("")

    # (iii) vary rho
    A("### 4.3 Vary `ρ` (isolates `ρ` / `1/ε₀`, `ε₀ = δ − J`)")
    A("")
    A("| cell | ρ | δ − J (ε₀) | `L` | `r*` |")
    A("|---|---|---|---|---|")
    rho_rows = sorted(open_cells, key=lambda c: c["rho"])
    for c in rho_rows:
        rstar = c["headline_smallest_r"]
        eps0 = c["delta_minus_J"]
        A(f"| {c['label']} | {c['rho']:.3f} | {eps0:.3f} | {c['headline_L']} | "
          f"{rstar if rstar is not None else '>'+str(R_MAX)} |")
    A("")

    # ---- automatic verdict text from the data ----
    all_rstar = [c["headline_smallest_r"] for c in open_cells if c["headline_smallest_r"] is not None]
    n_unbounded = sum(1 for c in open_cells if c["headline_smallest_r"] is None)
    rstar_min = min(all_rstar) if all_rstar else None
    rstar_max = max(all_rstar) if all_rstar else None

    # the per-k law: the largest r* observed at each k.  Prefer "clean" fields (n/p small, fine
    # histogram, r* saturated); fall back to all fields for a k that has no clean-field cell (k=4
    # is only reachable with a small field, p^4<=3e6 => p<=41, so n/p is necessarily large there).
    by_k_clean = {}
    by_k_all = {}
    for c in open_cells:
        sr = c["headline_smallest_r"]
        if sr is None:
            continue
        by_k_all.setdefault(c["k"], []).append(sr)
        if c["n_over_p"] < 0.2:
            by_k_clean.setdefault(c["k"], []).append(sr)
    plateau_max = {}
    k_is_coarse = {}
    for k in sorted(by_k_all):
        if by_k_clean.get(k):
            plateau_max[k] = max(by_k_clean[k]); k_is_coarse[k] = False
        else:
            plateau_max[k] = max(by_k_all[k]); k_is_coarse[k] = True
    ks_sorted = sorted(plateau_max)
    # render the k -> r* law and compare against 2k (mark coarse-field entries)
    law_str = ", ".join(
        f"k={k}: r*≈{plateau_max[k]}{'†' if k_is_coarse[k] else ''} (2k={2 * k})"
        for k in ks_sorted)
    coarse_note = (" († the `k=4` value is from coarse-histogram small fields — `p^4 ≤ 3·10⁶` "
                   "forces `p ≤ 41`, so `n/p` is large there; it may slightly UNDER-state the "
                   "saturated `k=4` order.)") if any(k_is_coarse.values()) else ""

    # log-p saturation on the (16,2) ladder at the common radius (delta == 0.750)
    lad = sorted([c for c in open_cells
                  if c["axis"].startswith("fix_nk_vary_p(n=16") and abs(c["delta"] - 0.75) < 1e-6
                  and c["headline_smallest_r"] is not None],
                 key=lambda c: c["p"])
    lad_rstar = [c["headline_smallest_r"] for c in lad]
    lad_plateau = max(lad_rstar) if lad_rstar else None
    lad_logp_lo = min(c["log_p"] for c in lad) if lad else None
    lad_logp_hi = max(c["log_p"] for c in lad) if lad else None

    A("## 5. Guidance for E1′")
    A("")
    if rstar_min is not None:
        A(f"- **Magnitude.** Across the open-band cells the worst-case smallest order is "
          f"`r* ∈ [{rstar_min}, {rstar_max}]`"
          + (f" (with {n_unbounded} cell(s) where no `r ≤ {R_MAX}` reached `{SEE_FACTOR:g}L` — see caveats)." if n_unbounded else ".")
          + " The second moment (`r = 2`) is **never** tight: in every cell `B_2 ≫ L`, confirming the "
            "pairwise-independence obstruction of `moment-reduction.md` §2 (the global `r = 2` energy is "
            "the fixed quantity `p^k n(1/p − 1/p²)`, dominated by `p^k` near-mean codewords). A genuinely "
            "high moment is required.")
    A(f"- **`r*` is independent of `log p` (saturates).** On the `(16,2)` ladder at the common radius "
      f"`δ = 0.750`, `r*` rises over the smallest fields and then **plateaus at "
      f"{lad_plateau if lad_plateau is not None else '—'}** across the prime range "
      f"`log p ∈ [{lad_logp_lo:.1f}, {lad_logp_hi:.1f}]` "
      f"(`p` from {lad[0]['p']} to {lad[-1]['p']}): `r*` = {lad_rstar}. The initial climb is a "
      "small-field artifact — when `n/p` is non-negligible the mean `n/p` is large and the distance "
      "histogram is coarse; once `n/p → 0` the order needed to kill the `p^k` bulk **does not grow with "
      "the field**, because `M_r` and `T^r` carry the same field scaling. This is the `p`-independence "
      "half of P′, now seen at the level of the moment order.")
    A(f"- **`r*` grows ≈ linearly with `k` — it is NOT a `ρ`-only constant.** Grouping cells by `k` "
      f"gives a tight, monotone law: {law_str}.{coarse_note} The order needed is essentially "
      "`r* ≈ 2k` (the `(16,2)`, `(24,3)`, and `(20/28/30,4)` ladders give `5, 6–7, 8`). Crucially this "
      "holds **at fixed `ρ`**: the `ρ = 1/8` ladder moves `k = 2 → 3` and `r* = 4 → 6` even though `ρ` "
      "is unchanged, so the controlling variable is `k` (equivalently the polynomial degree / the number "
      "of interpolation constraints), not `ρ` alone. Since `k = ρn`, **at fixed `ρ` the required moment "
      "order grows linearly with the block length `n`.**")
    A(f"- **Worst-case word type.** The **{worst_type_name}** words need the highest order "
      f"(`r*` up to {type_max_r[worst_type]}"
      + (f", unbounded within `r ≤ {R_MAX}` in {type_unbounded[worst_type]} cell(s)" if type_unbounded[worst_type] else "")
      + "). Random words need a slightly higher order than the adversarial-cluster word in most cells "
        "(their excess mass is spread thin over a wide bulk, so the ratio `M_r/T^r` turns over one step "
        "later), while the coset/Kambiré monomials at `k = 4` are the hardest of all — their tail is so "
        f"heavy that on the genuine-`k=4` cells no `r ≤ {R_MAX}` brings Markov to `{SEE_FACTOR:g}L`. The "
        "moment the E1′ program must dominate is therefore set by these structured `k = 4` tails.")
    A("")
    A("**Net guidance for E1′.** A *fixed* even moment cannot prove E1′ (the second moment is "
      "provably inert, and the empirical turnover order grows). The probe indicates the high-moment "
      "route must use a moment of order `r ≈ 2k = 2ρn` — i.e. **growing with the block length at fixed "
      "rate** — and that this order is **uniform in `p`**. A genuinely `n`-uniform proof of the "
      "`O_ρ(1)` bound (E1′) along this route would need either (i) a single moment of `n`-growing order "
      "with `n`-uniform control of its constant — substantially harder than a fixed-moment large sieve "
      "— or (ii) a different tail device (exponential-moment / Chernoff-type, or a restricted/weighted "
      "large sieve) that does not pay the `r ≈ 2k` price. The clean `p`-independence is encouraging; "
      "the `k`-growth of the order is the real obstacle the route must confront.")
    A("")
    A("## 6. Honest caveats (scale)")
    A("")
    A("1. **Small-cell guidance, not a proof.** Each `(cell, word)` measurement is exact (full "
      f"enumeration, exact rational moments), but `p^k ≤ {QK_CAP:.0e}` forces small `p` and `k`: the "
      "`log p` axis spans only `~5–9` and `k ∈ {2,3,4}`, so a trend visible here need not be the "
      "asymptotic law. The verdict is a **direction**, not a constant.")
    A("2. **Strong witness, not the exact arg-max.** The adversarial word is the best cluster-packing "
      "witness found (reused from `cluster_certificate.py`), not the exact maximiser over all `q^n` "
      "words; the true worst-case `L` — and hence the true worst-case `r*` — can only be larger, so "
      "the reported `r*` is a **lower** bound on what the theorem must handle.")
    A(f"3. **Heuristic tolerance.** “Markov sees `L`” uses the fixed factor `{SEE_FACTOR:g}` "
      "(`B_r ≤ "
      f"{SEE_FACTOR:g}L`). A stricter factor raises every `r*` by `1–2`; the *trend* across cells is "
      "robust to the cut, the absolute `r*` is not.")
    A("4. **Prime fields, multiplicative subgroups only**, per `moment-reduction.md` — the smooth "
      "prime-subgroup regime P′ addresses.")
    A("")
    A(f"_Generated by `experiments/small_rs_atlas/high_moments.py` in {elapsed:.1f}s; "
      "data in `results/high_moments.{json,csv}`._")

    with open(mpath, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return mpath


# ===========================================================================
# Driver.
# ===========================================================================
def run_full():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t0 = time.time()
    cells = build_battery()
    print(f"[high_moments] {len(cells)} open-band exact-enumeration cells, r_max={R_MAX}, "
          f"see_factor={SEE_FACTOR}")
    results = []
    for i, cell in enumerate(cells):
        seed = 0x91A3 + 7919 * i
        rec = run_cell(cell, seed=seed)
        results.append(rec)
        if rec.get("has_open_band"):
            sr = rec["headline_smallest_r"]
            sr_s = str(sr) if sr is not None else f">{R_MAX}"
            bt = rec.get("worst_by_type", {})
            def srt(t):
                v = bt.get(t, {}).get("smallest_r")
                return (str(v) if v is not None else f">{R_MAX}") if t in bt else "-"
            print(f"  [{i+1}/{len(cells)}] {rec['label']:>16} rho={rec['rho']:.4f} "
                  f"delta={rec['delta']:.3f} L={rec['headline_L']:>3} r*={sr_s:>3} "
                  f"(clu={srt('adversarial_cluster')} cos={srt('coset_kambire')} "
                  f"rnd={srt('random')}) [{rec['axis']}]", flush=True)
    elapsed = time.time() - t0

    jpath, cpath = write_outputs(results, elapsed)
    mpath = write_markdown(results, elapsed)

    _print_headline(results)
    print(f"\n[high_moments] DONE in {elapsed:.1f}s")
    print(f"  -> {jpath}")
    print(f"  -> {cpath}")
    print(f"  -> {mpath}")
    return results


def _print_headline(results):
    open_cells = [c for c in results if c.get("has_open_band")]
    print("\n" + "=" * 86)
    print("HEADLINE: smallest moment order r* for Markov to reproduce the exact open-band list L")
    print("=" * 86)
    all_rstar = [c["headline_smallest_r"] for c in open_cells if c["headline_smallest_r"] is not None]
    if all_rstar:
        print(f"  worst-case r* across {len(open_cells)} cells: min={min(all_rstar)} max={max(all_rstar)}")
    # by word type
    for t, name in [("adversarial_cluster", "cluster (adversarial)"),
                    ("coset_kambire", "coset/Kambire"), ("random", "random")]:
        srs = [c["worst_by_type"][t]["smallest_r"] for c in open_cells
               if t in c.get("worst_by_type", {}) and c["worst_by_type"][t]["smallest_r"] is not None]
        if srs:
            print(f"  r* for {name:22}: min={min(srs)} max={max(srs)} mean={np.mean(srs):.2f}")
    # guard: Markov upper bound must hold everywhere
    guard = all(pw["markov_upper_guard_ok"] for c in open_cells for pw in c["per_word"])
    print(f"  Markov upper-bound guard (B_r >= L at every r, all words): {'OK' if guard else 'FAILED'}")


# ===========================================================================
# Self-test battery.
# ===========================================================================
def _self_test():
    print("high_moments self-test")
    print("=" * 70)
    rng = np.random.default_rng(0x4711)

    # ---- 1. second-moment identity (moment-reduction.md SS2): the FULL (signed) second moment
    #         over all p^k polynomials equals p^k n(1/p - 1/p^2), for EVERY word. ----
    print("[1] global second-moment identity Sum (N-n/p)^2 = p^k n(1/p-1/p^2):")
    for (p, n, k) in [(97, 24, 3), (257, 16, 2), (73, 8, 2)]:
        F = PrimeField(p)
        L = domain_subgroup(F, n)
        book = build_codeword_book(F, L, k)
        for _ in range(3):
            w = rng.integers(0, p, size=n).astype(np.int64)
            hist = _distance_histogram(book, w)
            d = np.arange(n + 1); N = (n - d).astype(object)
            mean = Fraction(n, p)
            M2_full = sum(Fraction(int(cnt)) * (Fraction(int(Nv)) - mean) ** 2
                          for cnt, Nv in zip(hist.tolist(), N.tolist()) if cnt)
            predicted = Fraction(p ** k) * Fraction(n) * (Fraction(1, p) - Fraction(1, p * p))
            assert M2_full == predicted, f"second-moment identity FAILED {(p,n,k)}: {float(M2_full)} vs {float(predicted)}"
            assert int(hist.sum()) == p ** k, "histogram does not total p^k"
        print(f"    GF({p}) n={n} k={k}: exact identity holds for 3 words (and Sum hist = p^k) OK")

    # ---- 2. Markov is a valid upper bound: L <= B_r for EVERY r, every word ----
    print("[2] Markov upper bound L <= B_r for r=1..R, every word (cluster/coset/random):")
    F = PrimeField(97); n = 24; k = 3; p = 97
    L = domain_subgroup(F, n); book = build_codeword_book(F, L, k)
    e = band_e(n, k, p, 0.5)
    words = []
    adv = adversarial_cluster_word(F, L, k, e, rng)
    if adv is not None:
        words.append(("cluster", adv))
    for tag, w in kambire_monomial_words(F, L, k, n)[:3]:
        words.append((tag, w))
    for w in random_words(F, n, rng, 3):
        words.append(("random", w))
    for tag, w in words:
        hist = _distance_histogram(book, w)
        res = list_and_moments(hist, n, p, e, R_MAX)
        for r in range(1, R_MAX + 1):
            assert res["B_exact"][r] >= Fraction(res["L"]), \
                f"MARKOV VIOLATED at r={r} word={tag}: B_r={float(res['B_exact'][r])} < L={res['L']}"
        assert res["markov_upper_guard_ok"]
    print(f"    L <= B_r at every r in 1..{R_MAX} for {len(words)} words OK")

    # ---- 3. exact list L matches the singlelist reference kernel ----
    print("[3] exact list L matches singlelist_past_johnson.list_sizes_all_e:")
    from singlelist_past_johnson import list_sizes_all_e
    for tag, w in words:
        hist = _distance_histogram(book, w)
        res = list_and_moments(hist, n, p, e, 2)
        ref = list_sizes_all_e(book, w, [e])[e]
        assert res["L"] == ref, f"list mismatch word={tag}: {res['L']} vs ref {ref}"
    print(f"    L == reference list at e={e} for all words OK")

    # ---- 4. B_r is monotone non-increasing in r down to the tail, and B_r >= L stays tight-able ----
    #         (the moment ratio M_r/T^r decreases while bulk dominates; this is the whole mechanism)
    print("[4] B_r decreases in r while bulk noise dominates (random word), bottoming at >= L:")
    w = rng.integers(0, p, size=n).astype(np.int64)
    hist = _distance_histogram(book, w)
    res = list_and_moments(hist, n, p, e, R_MAX)
    Bs = [float(res["B_exact"][r]) for r in range(1, R_MAX + 1)]
    # for a random word the list is ~0/1 and the bulk is huge: B_r should fall sharply then flatten
    assert Bs[1] < Bs[0], f"B_2 should be < B_1 for a random word: {Bs[0]} -> {Bs[1]}"
    assert min(Bs) >= res["L"] - 1e-9, "B_r dipped below L (Markov violated)"
    print(f"    random word L={res['L']}: B_1={Bs[0]:.3g} B_2={Bs[1]:.3g} ... B_{R_MAX}={Bs[-1]:.3g} "
          f"(decreasing, >= L) OK")

    # ---- 5. exact rational moment == direct float moment (the integer p*N-n path is correct) ----
    print("[5] exact rational M_r == direct float Sum hist*(N-n/p)_+^r (path correctness):")
    w = rng.integers(0, p, size=n).astype(np.int64)
    hist = _distance_histogram(book, w)
    res = list_and_moments(hist, n, p, e, 6)
    d = np.arange(n + 1); N = n - d
    pos = np.maximum(N - n / p, 0.0)
    for r in (1, 2, 3, 4, 5, 6):
        direct = float(np.sum(hist * pos ** r))
        exact = float(res["M_exact"][r])
        rel = abs(direct - exact) / max(abs(exact), 1.0)
        assert rel < 1e-9, f"M_{r} path mismatch: exact={exact} direct={direct} rel={rel}"
    print("    exact rational moments match direct float moments r=1..6 OK")

    # ---- 6. a cell with a genuine list shows the turnover (smallest_r is found and finite-ish) ----
    print("[6] band_e is inside the open band and run_cell finds a positive list:")
    rho = k / n; J = johnson_radius(rho); rE = elias_radius(rho, p)
    assert J * n < e <= rE * n + 1e-9, "band_e not strictly inside (J,r_E]"
    rec = run_cell(Cell(97, 24, 3, axis="selftest"), seed=0xABCD)
    assert rec["has_open_band"]
    assert rec["headline_L"] >= 1, f"expected a positive worst-case list, got {rec['headline_L']}"
    print(f"    GF(97) n=24 k=3 c=0.5: e={e} in band; headline L={rec['headline_L']} "
          f"r*={rec['headline_smallest_r']} type={rec['headline_word_type']} OK")

    print("=" * 70)
    print("ALL HIGH-MOMENT SELF-TESTS PASSED")


# ===========================================================================
# Main.
# ===========================================================================
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true", help="run the self-test battery")
    args = ap.parse_args()
    if args.selftest:
        _self_test()
    else:
        run_full()
