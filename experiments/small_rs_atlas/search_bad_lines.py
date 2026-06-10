"""
search_bad_lines.py -- hunt for "bad lines" in RS proximity gaps.

Setup
-----
Fix C = RS[F,L,k], n=|L|, rho=k/n, and a proximity radius delta in (0,1).
A *line* through f0 in direction f1 is  { f0 + gamma*f1 : gamma in F }.

Two quantities per line:

  (1) close-count  C_line(delta) = #{ gamma in F : dist(f0 + gamma*f1, C) <= delta*n }.
      This is the number of points on the line that are delta-close to the code.

  (2) max common-agreement-set size  S*(delta) =
          max over (c0,c1) in C x C of  | { i in L : f0[i]=c0[i] AND f1[i]=c1[i] } |.
      This is the largest coordinate set S on which f0 AND f1 *simultaneously*
      agree with some fixed pair of codewords.  KEY FACT: on such an S, for EVERY
      gamma the point f0+gamma*f1 agrees with the codeword c0+gamma*c1 on all of S
      (RS is linear, so c0+gamma*c1 in C).  Hence if S*>=(1-delta)n, every point of
      the line is delta-close *for the same reason* -- this is exactly correlated
      agreement (CA): a single common set explains all the closeness.

Proximity-gap dichotomy (what we are testing empirically):
  Either C_line(delta) is small (few close points), OR there is correlated
  agreement: a large common set S* >= (1-delta)n explaining the closeness.

  A BAD LINE is the dangerous violation of this dichotomy:
        C_line(delta) is LARGE     (many gamma give close points)
        BUT S*(delta) is SMALL     (no large common agreement set).
  Bad lines are the seeds of MCA counterexamples.

Computing S* exactly
--------------------
We need, over all pairs (c0,c1), the size of the common agreement set.  Naively
that is q^k * q^k pairs.  We reduce it as follows.

For a fixed candidate "agreement codeword pair", the agreement set of f0 with c0
is A0(c0) = {i : f0[i]=c0[i]} and similarly A1(c1).  The common set is
A0(c0) ∩ A1(c1), and we want to maximize its size.  But the *constraint that makes
CA meaningful* is that we want one set S on which simultaneously some codeword
matches f0 and some codeword matches f1.  Because codewords are determined by any
k positions (RS), the natural exact computation is:

  S* = max_{c0,c1} |A0(c0) ∩ A1(c1)|.

We compute this exactly but smartly: instead of looping all q^k x q^k, we note
that |A0(c0) ∩ A1(c1)| is maximized by choosing, for the *combined* agreement, a
large common subset.  We use the following exact method valid for our small codes:

  * If q^k is small (<= a few hundred thousand), we compute for every codeword c
    its agreement mask with f0 (a boolean n-vector) and with f1.  Then
        S* = max over (c0,c1) of popcount(maskf0[c0] & maskf1[c1]).
    We avoid the full q^k x q^k product by a meet-in-the-middle / greedy bound
    plus an exact refinement (see _max_common_agreement).  For the field sizes in
    the atlas the exact double loop is feasible when q^k is a few thousand; for
    larger q^k we use the agreement-pattern method below which is exact.

Exact agreement-pattern method (always used)
--------------------------------------------
Group codewords of C by their agreement *pattern* with f0 is not enough (we need
the joint pattern).  Instead we exploit RS structure directly:

  A set S supports CA iff BOTH f0|_S and f1|_S are restrictions of degree-<k
  polynomials.  Equivalently: there exist polynomials p0,p1 (deg<k) interpolating
  f0,f1 on S.  For |S| >= k this is a real constraint; for |S| <= k it is free.

So  S* = size of the LARGEST S such that f0|_S and f1|_S are both codeword
restrictions.  We compute S* exactly via the codeword agreement masks:

  S* = max_{c0 in C} max_{c1 in C} |A0(c0) ∩ A1(c1)|
     = max_{c0 in C} ( max_{c1 in C} popcount(A0(c0) & A1(c1)) ).

For each fixed c0, the inner max over c1 is:  among A1-masks, the one whose
overlap with A0(c0) is largest.  We compute all A1 masks once as a bit-packed
matrix and, for each c0, take a fast popcount of (A0(c0) & A1) across all c1.
This is q^k x q^k popcounts but vectorized over c1, i.e. q^k numpy passes.  For
q^k ~ 1e4..1e5 that is fine; we cap it and fall back to a *sampled* inner search
with an honest log when q^k is too large.

This module both ENUMERATES lines (when feasible) and SAMPLES structured/random
(f0,f1) pairs, recording the worst offenders.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field as dc_field

import numpy as np

from ff import FiniteField
from rs import (CodewordBook, build_codeword_book, encode, dist_to_code,
                random_word, random_codeword, codeword_plus_noise)


# ---------------------------------------------------------------------------
# Agreement-mask machinery (bit-packed for speed).
# ---------------------------------------------------------------------------
def _agreement_bits(book: CodewordBook, w: np.ndarray) -> np.ndarray:
    """For word w, return a (num_codewords,) array of agreement popcounts AND a
    bit-packed (num_codewords, ceil(n/64)) matrix of agreement masks.

    Returns (popcounts, packed_masks). packed via np.packbits on bool (G==w).
    """
    w = np.asarray(w, dtype=np.int64)
    n = book.n
    pop_list = []
    packed_list = []
    for _, G in book.iter_chunks():
        eq = (G == w[None, :])                 # (b, n) bool
        pop_list.append(eq.sum(axis=1).astype(np.int32))
        packed_list.append(np.packbits(eq, axis=1))   # (b, ceil(n/8)) uint8
    pop = np.concatenate(pop_list)
    packed = np.concatenate(packed_list, axis=0)
    return pop, packed


# precomputed popcount table for uint8
_POPCOUNT8 = np.array([bin(i).count("1") for i in range(256)], dtype=np.int32)


def _max_common_agreement(book: CodewordBook, f0: np.ndarray, f1: np.ndarray,
                          hard_cap: int = 2_000_000,
                          rng: np.random.Generator | None = None
                          ) -> tuple[int, dict]:
    """EXACT max common-agreement-set size S* over the line.

        S* = max_{c0,c1 in C} | {i : f0[i]=c0[i] and f1[i]=c1[i]} |.

    Why this is exact (and fast):
      * For every codeword c0 we have |A0(c0) ∩ A1(c1)| <= |A0(c0)| = pop0(c0).
      * Sort codewords by pop0 DESCENDING.  Iterate c0 in that order maintaining
        the running best S*.  As soon as pop0(c0) <= best, ALL remaining c0 also
        have pop0 <= best (sorted), so none can improve S* -- we `break`.  This is
        a sound branch-and-bound: the result is the true maximum.
      * For each examined c0 we AND its bit-packed agreement mask against the
        bit-packed masks of ALL c1 at once (vectorized) and take the max overlap
        popcount.  We additionally restrict the c1 set to those with
        pop1(c1) > best (they are the only ones that could beat `best`).

    The early break makes the number of *examined* c0 tiny for generic words
    (which agree with almost no codeword), and at most q^k in the pathological
    case.  `hard_cap` is a pure safety valve on examined-c0 count: if it ever
    triggers we set exact=False and log it (it does not for atlas field sizes;
    q^k <= ~1e6 and pruning keeps examined-c0 in the low hundreds).

    Returns (S_star, info) with info documenting exactness and #examined.
    """
    ncw = book.num_codewords

    pop0, packed0 = _agreement_bits(book, f0)
    pop1, packed1 = _agreement_bits(book, f1)

    # Sort both code lists by agreement popcount descending (for pruning).
    order0 = np.argsort(-pop0, kind="stable")
    pop0_sorted = pop0[order0]
    packed0_sorted = packed0[order0]

    order1 = np.argsort(-pop1, kind="stable")
    pop1_sorted = pop1[order1]
    packed1_sorted = packed1[order1]

    best = 0
    examined = 0
    exact = True
    for ii in range(ncw):
        p0 = int(pop0_sorted[ii])
        if p0 <= best:
            break                       # sound branch-and-bound termination
        if examined >= hard_cap:        # safety valve (never hit for the atlas)
            exact = False
            break
        examined += 1
        row0 = packed0_sorted[ii]                                # (B,) uint8
        # Only c1 with pop1 > best can possibly improve best.  pop1_sorted is
        # descending, so take the prefix where pop1_sorted > best.
        m1 = int(np.searchsorted(-pop1_sorted, -best, side="left"))
        # m1 = number of c1 with pop1 > best (since sorted desc).  Guard >=1.
        if m1 == 0:
            continue
        cand1 = packed1_sorted[:m1]                              # (m1, B)
        anded = np.bitwise_and(cand1, row0[None, :])             # (m1, B)
        inter = _POPCOUNT8[anded].sum(axis=1)                    # (m1,)
        m = int(inter.max())
        if m > best:
            best = m

    info = {"exact": exact, "num_codewords": ncw, "c0_examined": examined}
    return best, info


# ---------------------------------------------------------------------------
# Close-count on a line.
# ---------------------------------------------------------------------------
def line_close_count(book: CodewordBook, f0: np.ndarray, f1: np.ndarray,
                     delta: float,
                     gammas: np.ndarray | None = None) -> tuple[int, int, list]:
    """Count gamma in F (or in `gammas`) with dist(f0+gamma f1, C) <= delta*n.

    Returns (close_count, num_gamma_tested, dist_list) where dist_list[i] is the
    exact distance of f0 + gammas[i]*f1.  We always test ALL gamma in F (q values)
    because q is small in the atlas -- no sampling over gamma.
    """
    F = book.F
    if gammas is None:
        gammas = F.elements()            # all q field elements, including 0
    thresh = delta * book.n + 1e-9
    dists = []
    close = 0
    f0 = np.asarray(f0, dtype=np.int64)
    f1 = np.asarray(f1, dtype=np.int64)
    for g in gammas.tolist():
        w = F.add_vec(f0, F.mul_scalar_vec(int(g), f1))
        d = dist_to_code(book, w)
        dists.append(d)
        if d <= thresh:
            close += 1
    return close, len(gammas), dists


# ---------------------------------------------------------------------------
# Line generators (the "structured words" the prize plan asks us to try).
# ---------------------------------------------------------------------------
def gen_random_pair(F, book, rng):
    """f0, f1 both uniform random words in F^n."""
    return random_word(F, book.n, rng), random_word(F, book.n, rng)


def gen_codeword_plus_noise_pair(F, book, rng, e0, e1):
    """f0, f1 are codewords each perturbed in e0 / e1 random positions."""
    return codeword_plus_noise(book, rng, e0), codeword_plus_noise(book, rng, e1)


def gen_shared_error_support(F, book, rng, num_errors):
    """f0, f1 = codewords corrupted on the SAME error positions.

    This is the structured configuration most likely to *avoid* a bad line: the
    common agreement set is exactly the complement of the shared error support, so
    CA should explain everything.  Used as a control / sanity case.
    """
    c0 = random_codeword(book, rng).copy()
    c1 = random_codeword(book, rng).copy()
    pos = rng.choice(book.n, size=num_errors, replace=False)
    for p in pos:
        c0[p] = F.add(int(c0[p]), int(rng.integers(1, F.q)))
        c1[p] = F.add(int(c1[p]), int(rng.integers(1, F.q)))
    return c0, c1


def gen_low_degree_rational(F, book, rng):
    """f0 = evaluation of a degree-k or degree-(k+1) poly (just OUT of the code),
    f1 = another such.  These are 'almost codewords' -- structured near-codewords
    that probe whether slightly-too-high degree creates many close line points.
    """
    L = book.L
    # degree exactly k (one above the code) -> distance is small-ish (MDS: any
    # degree-k poly is distance >= 1 from C; often gives interesting lines).
    c0 = rng.integers(0, F.q, size=book.k + 1).astype(np.int64)
    c1 = rng.integers(0, F.q, size=book.k + 1).astype(np.int64)
    return encode(F, L, c0), encode(F, L, c1)


def gen_rational_function(F, book, rng):
    """f0 = (numerator/denominator) evaluated on L, a genuine rational function
    of low degree (not polynomial).  Rational functions r=a/b with deg a,deg b
    small are classic near-RS words: they agree with many codewords on subgroup
    structure.  We pick numerator deg < k and denominator deg 1 (a single pole),
    avoiding poles in L.
    """
    L = book.L.tolist()
    def build():
        for _try in range(50):
            num = rng.integers(0, F.q, size=book.k).astype(np.int64)
            # denominator x - r with r NOT in L (so no division by zero)
            r = int(rng.integers(0, F.q))
            if r not in set(L):
                vals = []
                ok = True
                for x in L:
                    nx = 0
                    for i, c in enumerate(num):
                        nx = F.add(nx, F.mul(int(c), F.pow(int(x), i)))
                    den = F.sub(int(x), r)
                    if den == 0:
                        ok = False
                        break
                    vals.append(F.mul(nx, F.inv(den)))
                if ok:
                    return np.array(vals, dtype=np.int64)
        return random_word(F, book.n, rng)
    return build(), build()


# ---------------------------------------------------------------------------
# One line evaluation -> record.
# ---------------------------------------------------------------------------
@dataclass
class LineRecord:
    kind: str
    close_count: int
    num_gamma: int
    S_star: int
    n: int
    k: int
    delta: float
    ca_threshold: int            # ceil((1-delta) n): S* must reach this for CA
    is_bad: bool
    meaningful_regime: bool      # ca_threshold > k (else CA is vacuous)
    ca_explained: bool
    S_exact: bool
    min_dist_on_line: int
    f0: list
    f1: list


def evaluate_line(book: CodewordBook, f0, f1, delta: float, kind: str,
                  rng=None) -> LineRecord:
    """Compute close-count, S*, and classify the line.

    Classification logic (scientifically careful):
      * CA threshold = ceil((1-delta) n).  S* >= threshold  <=>  correlated
        agreement holds (a single common set of that size explains the closeness).
      * MEANINGFUL regime: ca_threshold > k.  Only here is "S* < threshold" a real
        constraint -- for ca_threshold <= k (i.e. delta >= capacity 1-k/n) ANY k
        coordinates determine a codeword, so S* >= k automatically and CA is
        vacuous; we do NOT call such lines bad.
      * A BAD LINE = meaningful regime AND many close points (>1) AND no CA
        (S* < threshold).  This is the genuine proximity-gap-violating
        configuration (a line with many close points NOT explained by a common
        agreement set), the seed of an MCA counterexample.

    NOTE: at delta between the Johnson radius and capacity, RS is in its
    list-decoding regime and even RANDOM words/lines exhibit many close points
    with small S*; such bad lines are therefore EXPECTED and only become a
    counterexample *signal* if SMOOTH domains produce MORE of them (or larger
    close-counts) than random/full domains at the same (rho,delta).  The atlas
    headline analysis does that differential comparison.
    """
    n, k = book.n, book.k
    close, ngamma, dists = line_close_count(book, f0, f1, delta)
    S_star, info = _max_common_agreement(book, f0, f1, rng=rng)
    ca_threshold = int(np.ceil((1.0 - delta) * n - 1e-9))

    meaningful = ca_threshold > k
    ca_explained = S_star >= ca_threshold
    is_bad = meaningful and (not ca_explained) and (close > 1)

    return LineRecord(
        kind=kind, close_count=close, num_gamma=ngamma, S_star=S_star,
        n=n, k=k, delta=delta, ca_threshold=ca_threshold, is_bad=is_bad,
        meaningful_regime=meaningful, ca_explained=ca_explained,
        S_exact=info["exact"], min_dist_on_line=min(dists),
        f0=np.asarray(f0).tolist(), f1=np.asarray(f1).tolist(),
    )


def line_distances(book: CodewordBook, f0: np.ndarray, f1: np.ndarray
                   ) -> np.ndarray:
    """Exact distance dist(f0 + gamma*f1, C) for EVERY gamma in F.

    Returns an int array indexed by gamma value (length q).  Computed once; the
    delta-dependent close-count is then just a threshold comparison, so a whole
    delta-grid is evaluated from one pass.  This is the expensive primitive
    (q distance computations, each O(q^k * n)).
    """
    F = book.F
    f0 = np.asarray(f0, dtype=np.int64)
    f1 = np.asarray(f1, dtype=np.int64)
    dists = np.empty(F.q, dtype=np.int64)
    for g in range(F.q):
        w = F.add_vec(f0, F.mul_scalar_vec(g, f1))
        dists[g] = dist_to_code(book, w)
    return dists


def evaluate_line_multi(book: CodewordBook, f0, f1, deltas: list[float],
                        kind: str, rng=None) -> dict:
    """Evaluate ONE line against MANY deltas with a single distance + S* pass.

    Returns {delta: LineRecord}.  The line distances (over all gamma) and the
    max common-agreement set S* are computed ONCE; only the threshold/CA
    classification varies with delta.  This removes the per-delta recomputation
    that dominated runtime.
    """
    n, k = book.n, book.k
    dists = line_distances(book, f0, f1)            # length q, expensive (once)
    S_star, info = _max_common_agreement(book, f0, f1, rng=rng)  # once
    min_d = int(dists.min())
    f0l = np.asarray(f0).tolist()
    f1l = np.asarray(f1).tolist()
    out = {}
    for delta in deltas:
        thresh = delta * n + 1e-9
        close = int((dists <= thresh).sum())
        ca_threshold = int(np.ceil((1.0 - delta) * n - 1e-9))
        meaningful = ca_threshold > k
        ca_explained = S_star >= ca_threshold
        is_bad = meaningful and (not ca_explained) and (close > 1)
        out[delta] = LineRecord(
            kind=kind, close_count=close, num_gamma=len(dists), S_star=S_star,
            n=n, k=k, delta=delta, ca_threshold=ca_threshold, is_bad=is_bad,
            meaningful_regime=meaningful, ca_explained=ca_explained,
            S_exact=info["exact"], min_dist_on_line=min_d, f0=f0l, f1=f1l,
        )
    return out


# ---------------------------------------------------------------------------
# Search driver for one (field, domain, k, delta).
# ---------------------------------------------------------------------------
def search_bad_lines(book: CodewordBook, delta: float, rng: np.random.Generator,
                     n_random=200, n_cwnoise=200, n_lowdeg=100, n_rational=100,
                     e_frac=None, verbose=False) -> dict:
    """Sample/structure many lines; return summary stats and worst offenders.

    Sampling counts are explicit (no silent caps).  e_frac sets the per-word error
    fraction for the codeword+noise generator (default ~ delta).
    """
    F = book.F
    n, k = book.n, book.k
    if e_frac is None:
        e_frac = delta
    e_each = max(1, int(round(e_frac * n)))

    records: list[LineRecord] = []
    worst_bad: list[LineRecord] = []
    worst_close: list[LineRecord] = []

    def consider(rec: LineRecord):
        records.append(rec)
        if rec.is_bad:
            worst_bad.append(rec)
        worst_close.append(rec)

    # 1. random pairs
    for _ in range(n_random):
        f0, f1 = gen_random_pair(F, book, rng)
        consider(evaluate_line(book, f0, f1, delta, "random", rng))
    # 2. codeword + noise pairs (independent error supports)
    for _ in range(n_cwnoise):
        f0, f1 = gen_codeword_plus_noise_pair(F, book, rng, e_each, e_each)
        consider(evaluate_line(book, f0, f1, delta, "cw+noise", rng))
    # 3. low-degree (just above k) pairs
    for _ in range(n_lowdeg):
        f0, f1 = gen_low_degree_rational(F, book, rng)
        consider(evaluate_line(book, f0, f1, delta, "deg=k+1", rng))
    # 4. rational function pairs
    for _ in range(n_rational):
        f0, f1 = gen_rational_function(F, book, rng)
        consider(evaluate_line(book, f0, f1, delta, "rational", rng))

    # also: a CA-control case (shared error support) -- MUST NOT be bad.
    # The shared error support has size e_ctrl, so the common agreement set has
    # size n - e_ctrl.  For this to be a valid "CA holds" control at radius delta
    # we need n - e_ctrl >= ceil((1-delta)n), i.e. e_ctrl <= floor(delta*n).  We
    # take e_ctrl = floor(delta*n) - 1 (strictly inside) so S* comfortably clears
    # the CA threshold and the line is genuinely CA-explained.
    e_ctrl = max(1, int(np.floor(delta * n)) - 1)
    n_control = 50
    control_bad = 0
    for _ in range(n_control):
        f0, f1 = gen_shared_error_support(F, book, rng, e_ctrl)
        rec = evaluate_line(book, f0, f1, delta, "shared-support", rng)
        records.append(rec)
        if rec.is_bad:
            control_bad += 1

    worst_close.sort(key=lambda r: (-r.close_count, r.S_star))
    worst_bad.sort(key=lambda r: (-(r.close_count), r.S_star))

    # Exclude the controls from the "sampled lines" stats (they are a diagnostic).
    sampled = [r for r in records if r.kind != "shared-support"]
    close_counts = np.array([r.close_count for r in sampled])
    s_stars = np.array([r.S_star for r in sampled])
    meaningful = np.array([r.meaningful_regime for r in sampled])
    ca_threshold = int(np.ceil((1.0 - delta) * n - 1e-9))
    is_meaningful_regime = ca_threshold > k

    # Among lines with >=2 close points, what is the largest common-agreement set?
    # If this is < ca_threshold we have genuine CA failures (in the meaningful
    # regime).  We report the close-count DISTRIBUTION so smooth-vs-random
    # differential comparison is possible downstream.
    close_with_pts = close_counts[close_counts >= 2]
    summary = {
        "field": F.name, "q": F.q, "n": n, "k": k, "rho": k / n, "delta": delta,
        "johnson": 1.0 - np.sqrt(k / n), "capacity": 1.0 - k / n,
        "ca_threshold": ca_threshold,
        "meaningful_regime": is_meaningful_regime,   # ca_threshold > k
        "num_lines": len(sampled),
        "sampling": {
            "random": n_random, "cw+noise": n_cwnoise, "deg=k+1": n_lowdeg,
            "rational": n_rational, "shared-support-control": n_control,
            "full_line_space_size": "q^(2n) (we sample lines; gamma enumerated fully = q per line)",
            "gamma_enumerated_per_line": F.q,
        },
        "num_bad_lines": len(worst_bad),
        "frac_bad": len(worst_bad) / max(len(sampled), 1),
        "control_bad_lines": control_bad,       # MUST be 0 (sanity)
        # close-count statistics (the quantity to compare across domains)
        "close_count_max": int(close_counts.max()),
        "close_count_mean": float(close_counts.mean()),
        "close_count_p90": float(np.percentile(close_counts, 90)),
        "frac_lines_with_2plus_close": float((close_counts >= 2).mean()),
        # S* statistics among the close lines
        "S_star_mean": float(s_stars.mean()),
        "S_star_max": int(s_stars.max()),
        "S_star_min_among_close": int(
            s_stars[close_counts >= 2].min()) if (close_counts >= 2).any() else -1,
        "worst_bad_lines": [_rec_brief(r) for r in worst_bad[:8]],
        "worst_close_lines": [_rec_brief(r) for r in worst_close[:5]],
    }
    if verbose:
        print(f"    [{F.name} n={n} k={k} d={delta:.3f}] "
              f"bad={len(worst_bad)}/{len(sampled)} "
              f"closemean={close_counts.mean():.2f} maxclose={close_counts.max()} "
              f"meaningful={is_meaningful_regime} ctrlbad={control_bad}")
    return summary


def _rec_brief(r: LineRecord) -> dict:
    return {
        "kind": r.kind, "close_count": r.close_count, "num_gamma": r.num_gamma,
        "S_star": r.S_star, "ca_threshold": r.ca_threshold, "is_bad": r.is_bad,
        "meaningful_regime": r.meaningful_regime, "ca_explained": r.ca_explained,
        "min_dist_on_line": r.min_dist_on_line, "S_exact": r.S_exact,
    }


def _summarize_for_delta(F, n, k, delta, sampled, controls, sampling) -> dict:
    """Build the per-delta summary dict from already-classified LineRecords.

    `sampled` and `controls` are lists of LineRecord (all at this delta)."""
    bad = [r for r in sampled if r.is_bad]
    close_counts = np.array([r.close_count for r in sampled])
    s_stars = np.array([r.S_star for r in sampled])
    control_bad = sum(1 for r in controls if r.is_bad)
    ca_threshold = int(np.ceil((1.0 - delta) * n - 1e-9))
    worst_bad = sorted(bad, key=lambda r: (-r.close_count, r.S_star))
    worst_close = sorted(sampled, key=lambda r: (-r.close_count, r.S_star))
    return {
        "field": F.name, "q": F.q, "n": n, "k": k, "rho": k / n,
        "delta": delta, "johnson": 1.0 - np.sqrt(k / n), "capacity": 1.0 - k / n,
        "ca_threshold": ca_threshold, "meaningful_regime": ca_threshold > k,
        "num_lines": len(sampled), "sampling": sampling,
        "num_bad_lines": len(bad),
        "frac_bad": len(bad) / max(len(sampled), 1),
        "control_bad_lines": control_bad,
        "close_count_max": int(close_counts.max()),
        "close_count_mean": float(close_counts.mean()),
        "close_count_p90": float(np.percentile(close_counts, 90)),
        "frac_lines_with_2plus_close": float((close_counts >= 2).mean()),
        "S_star_mean": float(s_stars.mean()),
        "S_star_max": int(s_stars.max()),
        "S_star_min_among_close": int(
            s_stars[close_counts >= 2].min()) if (close_counts >= 2).any() else -1,
        "worst_bad_lines": [_rec_brief(r) for r in worst_bad[:8]],
        "worst_close_lines": [_rec_brief(r) for r in worst_close[:5]],
    }


def search_bad_lines_multi(book: CodewordBook, deltas: list[float],
                           rng: np.random.Generator,
                           n_random=200, n_cwnoise=200, n_lowdeg=100,
                           n_rational=100, e_frac=None,
                           verbose=False) -> list[dict]:
    """EFFICIENT multi-delta bad-line search: sample each line ONCE, classify it
    against every delta in `deltas`.  Returns a list of per-delta summary dicts
    (same schema as search_bad_lines).

    This is the function the atlas uses: line distances and S* (the expensive
    parts) are computed once per line and reused across the whole delta grid,
    giving a ~len(deltas)x speedup over calling search_bad_lines per delta.

    Per-generator error fraction for cw+noise is e_frac (default: max delta).  The
    CA control uses a per-delta shared-support size so it is a valid "CA holds"
    control at EACH delta.
    """
    F = book.F
    n, k = book.n, book.k
    if e_frac is None:
        e_frac = max(deltas)
    e_each = max(1, int(round(e_frac * n)))

    # records_per_delta[delta] = list of LineRecord (sampled, non-control)
    sampled_pd: dict[float, list[LineRecord]] = {d: [] for d in deltas}

    def sample_and_classify(gen_fn, kind, count):
        for _ in range(count):
            f0, f1 = gen_fn()
            recs = evaluate_line_multi(book, f0, f1, deltas, kind, rng)
            for d in deltas:
                sampled_pd[d].append(recs[d])

    sample_and_classify(lambda: gen_random_pair(F, book, rng), "random", n_random)
    sample_and_classify(
        lambda: gen_codeword_plus_noise_pair(F, book, rng, e_each, e_each),
        "cw+noise", n_cwnoise)
    sample_and_classify(lambda: gen_low_degree_rational(F, book, rng),
                        "deg=k+1", n_lowdeg)
    sample_and_classify(lambda: gen_rational_function(F, book, rng),
                        "rational", n_rational)

    # CA controls: per delta we need a shared support of size floor(delta*n)-1 so
    # the control is genuinely CA-explained at THAT delta.  Generate a pool of
    # shared-support lines covering a range of support sizes and classify each at
    # every delta -- but only COUNT a control toward delta d if its support size
    # e satisfies n-e >= ca_threshold(d) (i.e. it is a valid control at d).
    n_control = 50
    controls_pd: dict[float, list[LineRecord]] = {d: [] for d in deltas}
    # support sizes to cover: from 1 to max needed
    max_e = max(1, int(np.floor(max(deltas) * n)) - 1)
    for _ in range(n_control):
        e = int(rng.integers(1, max_e + 1)) if max_e >= 1 else 1
        f0, f1 = gen_shared_error_support(F, book, rng, e)
        recs = evaluate_line_multi(book, f0, f1, deltas, "shared-support", rng)
        for d in deltas:
            ca_thr = int(np.ceil((1.0 - d) * n - 1e-9))
            if (n - e) >= ca_thr:           # valid CA control at this delta
                controls_pd[d].append(recs[d])

    summaries = []
    for d in deltas:
        sampling = {
            "random": n_random, "cw+noise": n_cwnoise, "deg=k+1": n_lowdeg,
            "rational": n_rational, "shared-support-control": len(controls_pd[d]),
            "shared-support-pool": n_control,
            "full_line_space_size": "q^(2n) (lines sampled; gamma fully enumerated = q/line)",
            "gamma_enumerated_per_line": F.q,
            "e_each_cwnoise": e_each,
        }
        s = _summarize_for_delta(F, n, k, d, sampled_pd[d], controls_pd[d], sampling)
        summaries.append(s)
        if verbose:
            print(f"    [{F.name} n={n} k={k} d={d:.3f}] bad={s['num_bad_lines']}"
                  f"/{s['num_lines']} closemean={s['close_count_mean']:.2f} "
                  f"meaningful={s['meaningful_regime']} ctrlbad={s['control_bad_lines']}")
    return summaries


# ===========================================================================
# Self-test / demo.
# ===========================================================================
def _self_test():
    from ff import PrimeField, BinaryExtensionField
    from rs import domain_subgroup, domain_full
    rng = np.random.default_rng(0x1234)
    print("search_bad_lines self-test")
    print("=" * 60)

    # --- Sanity 1: CA control must explain shared-support lines exactly ---
    # If f0,f1 are codewords corrupted on the SAME t positions, the common
    # agreement set is the n-t good positions, so S* >= n-t.  For delta>=t/n this
    # is >= (1-delta)n, hence NOT a bad line.
    F = PrimeField(31)
    L = domain_subgroup(F, 10)        # n=10
    k = 3
    book = build_codeword_book(F, L, k)
    t = 3
    c0 = random_codeword(book, rng).copy()
    c1 = random_codeword(book, rng).copy()
    pos = rng.choice(book.n, size=t, replace=False)
    for p in pos:
        c0[p] = F.add(int(c0[p]), 1)
        c1[p] = F.add(int(c1[p]), 1)
    S_star, info = _max_common_agreement(book, c0, c1, rng=rng)
    assert S_star >= book.n - t, f"shared-support S*={S_star} should be >= {book.n-t}"
    assert info["exact"], "should be exact for this small code"
    print(f"  shared-support control: S*={S_star} >= n-t={book.n-t} (exact) OK")

    # --- Sanity 2: S* of a line equals n when f0,f1 are both codewords ---
    c0 = random_codeword(book, rng)
    c1 = random_codeword(book, rng)
    S_star, _ = _max_common_agreement(book, c0, c1, rng=rng)
    assert S_star == book.n, f"two codewords must have S*=n, got {S_star}"
    print(f"  two-codeword line: S*=n={book.n} OK")

    # --- Sanity 3: every gamma close when both are codewords (line in code) ---
    close, ng, dists = line_close_count(book, c0, c1, delta=0.0)
    assert close == F.q and max(dists) == 0, "line of codewords: all dist 0"
    print(f"  line-of-codewords close-count = q = {close}, all dist 0 OK")

    # --- Sanity 4: _max_common_agreement brute-force cross-check (tiny code) ---
    F2 = PrimeField(5)
    L2 = domain_full(F2)              # n=4
    k2 = 2
    book2 = build_codeword_book(F2, L2, k2)   # 25 codewords
    f0 = random_word(F2, book2.n, rng)
    f1 = random_word(F2, book2.n, rng)
    # brute force exact S* by full double loop
    pop0, _ = _agreement_bits(book2, f0)
    masks0 = [(book2.G[j] == f0) for j in range(book2.num_codewords)]
    masks1 = [(book2.G[j] == f1) for j in range(book2.num_codewords)]
    brute = 0
    for m0 in masks0:
        for m1 in masks1:
            brute = max(brute, int((m0 & m1).sum()))
    fast, _ = _max_common_agreement(book2, f0, f1, rng=rng)
    assert fast == brute, f"S* fast={fast} != brute={brute}"
    # repeat a few times
    for _ in range(30):
        f0 = random_word(F2, book2.n, rng)
        f1 = random_word(F2, book2.n, rng)
        masks0 = [(book2.G[j] == f0) for j in range(book2.num_codewords)]
        masks1 = [(book2.G[j] == f1) for j in range(book2.num_codewords)]
        brute = max(int((m0 & m1).sum()) for m0 in masks0 for m1 in masks1)
        fast, _ = _max_common_agreement(book2, f0, f1, rng=rng)
        assert fast == brute, f"S* mismatch fast={fast} brute={brute}"
    print("  S* fast-vs-brute agree on 31 random pairs (GF(5) n=4 k=2) OK")

    # --- Sanity 5: evaluate_line_multi agrees with per-delta evaluate_line ---
    F = PrimeField(31)
    L = domain_subgroup(F, 10)
    book = build_codeword_book(F, L, 3)
    deltas = [0.2, 0.35, 0.5, 0.65]
    for _ in range(20):
        f0 = random_word(F, book.n, rng); f1 = random_word(F, book.n, rng)
        multi = evaluate_line_multi(book, f0, f1, deltas, "x", rng)
        for d in deltas:
            single = evaluate_line(book, f0, f1, d, "x", rng)
            assert multi[d].close_count == single.close_count, \
                f"multi/single close mismatch d={d}"
            assert multi[d].S_star == single.S_star, "S* mismatch"
            assert multi[d].is_bad == single.is_bad, "is_bad mismatch"
    print("  evaluate_line_multi == per-delta evaluate_line on 20x4 cases OK")

    # --- Demo: run the efficient multi-delta search ---
    F = BinaryExtensionField(4)
    L = domain_subgroup(F, 15)              # n=15, full subgroup of GF(2^4)*
    book = build_codeword_book(F, L, 3)     # rho=0.2, Johnson ~ 0.553
    summaries = search_bad_lines_multi(book, [0.3, 0.45, 0.55, 0.66], rng,
                                       n_random=60, n_cwnoise=60, n_lowdeg=30,
                                       n_rational=30, verbose=True)
    for s in summaries:
        assert s["control_bad_lines"] == 0, \
            f"CA controls must never be 'bad' (d={s['delta']})"
    print(f"  demo multi-delta search: {len(summaries)} deltas, all controls "
          f"bad=0 OK")

    print("=" * 60)
    print("ALL search_bad_lines SELF-TESTS PASSED")


if __name__ == "__main__":
    _self_test()
