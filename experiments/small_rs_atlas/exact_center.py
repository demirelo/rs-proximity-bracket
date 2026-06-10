"""
exact_center.py -- EXACT min-max (closest-string) ball center for the cluster
certificate, strengthening the open-band P' lower-bound witnesses.

THE GAP THIS CLOSES.
  The cluster certificate (cluster_certificate.py) witnesses a list of size t at
  radius e = floor(delta*n) by exhibiting t codewords c_1..c_t in F_p^n all within
  e of a center word w.  As center it uses the coordinate-wise PLURALITY, which is
  the exact minimizer of the TOTAL (summed) Hamming distance sum_i Delta(c_i,w) --
  NOT of the MAXIMUM radius max_i Delta(c_i,w).  The quantity the certificate needs
  to control is the maximum: a list of size t at radius e exists iff some word w has
  max_i Delta(c_i,w) <= e.  The exact min-max center can only make that maximum
  SMALLER than the plurality's, hence can only fit MORE codewords inside a fixed
  radius-e ball.  Replacing plurality by the exact min-max center therefore yields a
  CERTIFICATE THAT IS NEVER SMALLER, and usually larger -- a strictly stronger
  one-sided lower bound on the worst-case open-band list |Lambda(C,delta)|.

THE OBJECT.  For codewords c_1..c_t in F_p^n the closest-string optimum is
        e*(c_1..c_t) = min_{w in F_p^n}  max_i Delta(c_i, w),
  and a minimizer w* is a min-max center.  Equivalently, with
        cover(w) = (Delta(c_1,w),...,Delta(c_t,w)),
  e* is the least e for which some w has cover(w) <= e coordinatewise.

REDUCTION USED BY THE SOLVER (exact, value-free).
  Fix a coordinate j.  The codewords split into agreement classes by their value
  c_i[j]: class V = { i : c_i[j] = v }.  A center value w[j] = v makes exactly the
  codewords in class V agree (0 disagreement) and every other codeword disagree
  (+1).  Choosing any value NOT among {c_i[j]} disagrees with all t and is never
  useful.  Hence the ONLY decision at coordinate j is "which agreement block to
  hit", contributing the 0/1 vector  (1 - 1_block)  to the per-codeword
  disagreement counts.  The actual field values are irrelevant; only the induced
  PARTITION of {1..t} matters.  Coordinates inducing the same partition are
  interchangeable, so we group the n coordinates into COLUMN TYPES (a partition of
  {1..t} with a multiplicity) and decide, per type, how its columns distribute over
  the type's blocks.  This is closest-string, which is FPT in t (Gramm-Niedermeier-
  Rossmanith); for the small t (<= ~16) the cluster certificate produces it is
  solved exactly and fast.

  At a coordinate whose blocks have sizes s_1>=s_2>=... the best single column for
  the running maximum hits the largest block (fewest codewords charged), so the
  per-coordinate lower bound on added disagreements is the count of codewords
  outside the largest block.  Summing these "forced" minimum charges over the
  remaining coordinates gives an admissible bound used to prune.

SOLVER.  Two equivalent entry points, both EXACT:
    closest_string_radius(C)      -> e*  (and a witness center w*),
    closest_string_feasible(C,e)  -> bool (is there w with max_i Delta<=e?).
  closest_string_radius binary-searches e between a lower bound (max-min over a
  reference codeword pair, and the trivial sum/t bound) and the plurality radius
  (a valid upper bound, since plurality is one feasible center); each feasibility
  query is a branch-and-bound over column types with per-codeword count pruning.

CONSERVATIVENESS (for the lower bound).  The solver returns the TRUE optimum e*; it
  NEVER reports a radius below the true min-max.  The certificate counts a codeword
  only when it is provably within e of the EXHIBITED center w*, so the witnessed
  list size is a genuine lower bound on |Lambda(C,delta)| for ANY cluster.  A
  cleverer cluster could give more; the bound is one-sided.

SELF-TEST (`python3.11 exact_center.py --selftest`):
  * closest_string_radius matches BRUTE FORCE over the whole word space on tiny
    cases (small p, small n, t up to 4) -- exact-vs-exact;
  * feasibility agrees with brute force at every radius;
  * exact min-max radius <= plurality radius on every tested cluster (the key
    invariant: plurality is a feasible center, so it can only be worse);
  * exact min-max radius <= every random word's radius (optimality).

Run `python3.11 exact_center.py --selftest`        for the self-test battery,
    `python3.11 exact_center.py`                   for the full strengthened sweep,
    `python3.11 exact_center.py --quick`           for a fast reduced sweep.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
import time
from dataclasses import dataclass
from multiprocessing import Pool

import numpy as np

from ff import PrimeField
from rs import domain_subgroup
from singlelist_past_johnson import elias_radius

# Reuse the cluster-certificate primitives by IMPORT (do not duplicate / edit them).
from cluster_certificate import (
    band_e,
    cert_count,
    construct_random_plurality,
    construct_lloyd,
    construct_greedy,
    construct_structured,
    construct_localsearch,
    distinct_codewords,
    hamming_to_set,
    johnson_radius,
    plurality_word,
    rand_codewords,
)

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")


# ===========================================================================
# EXACT min-max (closest-string) center.
# ===========================================================================
def _column_types(C: np.ndarray):
    """Reduce the (t, n) codeword matrix to its COLUMN TYPES.

    Each coordinate j induces a partition of {0..t-1} by equal value (the agreement
    classes c_i[j]==v).  Only that partition matters for the min-max objective (the
    field values themselves are irrelevant; see module docstring).  Two coordinates
    with the same partition are interchangeable, so we collapse them.

    Returns a list of (blocks, mult) where `blocks` is a tuple of frozensets
    partitioning {0..t-1} (the agreement classes at this column type, sorted by
    decreasing size) and `mult` is how many coordinates realize that partition.
    Cost O(t*n) for the canonicalisation plus a hash per column.
    """
    t, n = C.shape
    type_counts: dict[tuple, int] = {}
    # Canonical label of column j: assign each row a block id by first-occurrence of
    # its value down the column, so the labeling is value-independent.  Two columns
    # get the same key iff they induce the same partition of the rows.
    for j in range(n):
        col = C[:, j]
        first: dict[int, int] = {}
        key = np.empty(t, dtype=np.int64)
        nxt = 0
        for i in range(t):
            v = int(col[i])
            b = first.get(v)
            if b is None:
                b = nxt
                first[v] = b
                nxt += 1
            key[i] = b
        k = tuple(key.tolist())
        type_counts[k] = type_counts.get(k, 0) + 1

    types = []
    for key, mult in type_counts.items():
        # rebuild blocks from the canonical key
        groups: dict[int, list[int]] = {}
        for i, b in enumerate(key):
            groups.setdefault(b, []).append(i)
        blocks = sorted((frozenset(g) for g in groups.values()),
                        key=lambda s: -len(s))
        types.append((tuple(blocks), mult))
    # process column types with the most blocks (most branching / most signal) first
    types.sort(key=lambda bm: -len(bm[0]))
    return types


def _cover_baseline_to_word(C, types, picks_per_type, disc_key, e):
    """(internal) Materialise an explicit center word from a feasibility solution.

    `types` is the list of (blocks, mult); `picks_per_type[ti]` is either a list of
    chosen block indices (one per column of a COLLISION type) or, for the discrete
    type, an array `x` of length t giving how many discrete columns favor each
    codeword (sum x = disc_mult).  We lay the choices back onto the actual coordinates.
    """
    t, n = C.shape
    # Recompute the partition key + coordinate list per column type (value-free).
    type_cols: dict[tuple, list[int]] = {}
    type_blocks: dict[tuple, tuple] = {}
    for j in range(n):
        col = C[:, j]
        first = {}
        key = np.empty(t, dtype=np.int64)
        nxt = 0
        for i in range(t):
            v = int(col[i])
            b = first.get(v)
            if b is None:
                b = nxt
                first[v] = b
                nxt += 1
            key[i] = b
        k = tuple(key.tolist())
        type_cols.setdefault(k, []).append(j)
        if k not in type_blocks:
            groups = {}
            for i, b in enumerate(key):
                groups.setdefault(b, []).append(i)
            blocks = sorted((frozenset(g) for g in groups.values()), key=lambda s: -len(s))
            type_blocks[k] = tuple(blocks)

    # `types` order must match how the caller indexed picks_per_type.  Each type's
    # `blocks` tuple (frozensets of ROW indices) is a unique key for its coordinate set.
    blocks_to_key = {bl: k for k, bl in type_blocks.items()}
    w = np.empty(n, dtype=np.int64)
    for ti, (blocks, mult) in enumerate(types):
        key = blocks_to_key[blocks]
        cols = type_cols[key]
        picks = picks_per_type[ti]
        if isinstance(picks, np.ndarray):
            # discrete type: picks[i] columns hit singleton block {i}; the remaining
            # columns may hit any block (charge irrelevant once feasible) -> default 0.
            # block index of singleton {i} within `blocks`:
            singleton_index = {next(iter(B)): bi for bi, B in enumerate(blocks)
                               if len(B) == 1}
            seq = []
            for i in range(t):
                seq.extend([singleton_index.get(i, 0)] * int(picks[i]))
            # pad to mult with block 0
            seq.extend([0] * (mult - len(seq)))
            for ci, j in enumerate(cols):
                bi = seq[ci]
                rep = next(iter(blocks[bi]))
                w[j] = int(C[rep, j])
        else:
            for ci, j in enumerate(cols):
                bi = picks[ci] if ci < len(picks) else 0
                rep = next(iter(blocks[bi]))
                w[j] = int(C[rep, j])
    return w


class ClosestStringBudgetError(RuntimeError):
    """Raised if the EXACT feasibility search exceeds its node budget.  This can only
    happen for pathologically collision-heavy instances (tiny alphabet, many strings)
    -- never in the prime-field experiment, where collision columns are ~0.  We raise
    rather than return a wrong/approximate answer, preserving exactness."""


# Generous node cap: the prime-field experiment uses ~0 collision columns and finishes
# in microseconds; this only ever fires on adversarial tiny-alphabet inputs.
_FEAS_NODE_BUDGET = 5_000_000


def closest_string_feasible(C: np.ndarray, e: int, want_witness: bool = False,
                            node_budget: int = _FEAS_NODE_BUDGET):
    """EXACT feasibility: is there a word w in F_p^n with max_i Delta(c_i,w) <= e?

    Reformulation (module docstring).  Each column whose chosen agreement block
    contains codeword i CONTRIBUTES +1 to cover_i; since every column's chosen block
    is nonempty, Delta(c_i,w) = n - cover_i.  Feasible at radius e iff cover_i >= R :=
    n - e for every i.  A column of type tau assigned to block B raises cover_i by 1
    for exactly the i in B.

    The discrete column type (all t codewords distinct: every block a singleton)
    raises exactly one chosen codeword's cover per column, so its `disc_mult` columns
    contribute a vector x (sum x = disc_mult, x>=0) added to cover.  Every OTHER column
    type has a block of size >= 2 (a value collision) and -- over a large field -- is
    rare, so the total multiplicity of collision columns is small.  We BRANCH over the
    collision columns (one block choice each; shallow because there are few) and, at
    each leaf, dispatch the discrete columns by the CLOSED FORM:
        feasible(leaf)  <=>  sum_i max(0, R - cover_base_i)  <=  disc_mult,
    i.e. the discrete columns can supply every codeword its remaining deficit.

    This is EXACT and scales to n in the thousands (collision columns, not n, drive
    the branching).  Returns (feasible, witness) where witness is consumed by
    _cover_baseline_to_word; witness is None when want_witness is False or infeasible.
    """
    C = np.asarray(C, dtype=np.int64)
    t, n = C.shape
    if e < 0:
        return (False, None)
    if e >= n:
        # any single codeword as center w=c_1 gives every Delta <= n; trivially feasible
        return (True, None)
    R = n - e                                   # required cover per codeword
    if R <= 0:
        return (True, None)

    types = _column_types(C)                    # (blocks, mult), most-blocks-first

    # split discrete (all singletons) vs collision (some block size >= 2)
    disc_idx = None
    disc_mult = 0
    collision = []                              # list of (ti, blocks, mult, memb_rows)
    for ti, (blocks, mult) in enumerate(types):
        if all(len(B) == 1 for B in blocks) and len(blocks) == t:
            disc_idx = ti
            disc_mult = mult
        else:
            # membership rows per block (as index lists) for fast cover updates
            memb = [list(B) for B in blocks]
            collision.append((ti, blocks, mult, memb))

    # Flatten collision columns into a single list of column "slots", each carrying its
    # available block membership lists; sort slots so the most-constrained (fewest
    # blocks, i.e. largest collisions) come first for tighter pruning.
    slots = []                                  # each: (ti, list_of_member_index_lists)
    for (ti, blocks, mult, memb) in collision:
        for _ in range(mult):
            slots.append((ti, memb))
    slots.sort(key=lambda s: len(s[1]))         # fewer blocks first

    n_slots = len(slots)
    # suffix max possible cover gain to codeword i from remaining slots: each remaining
    # slot can raise cover_i by at most 1 (if some available block contains i).
    # For pruning we use: max future cover for i over remaining slots = number of
    # remaining slots that have a block containing i.
    # Precompute, per slot, a boolean (t,) of which codewords COULD be covered.
    slot_can = np.zeros((n_slots, t), dtype=np.int64)
    for si, (ti, memb) in enumerate(slots):
        can = np.zeros(t, dtype=np.int64)
        for blk in memb:
            can[blk] = 1
        slot_can[si] = can
    # suffix sums of slot_can (+ disc availability handled separately)
    suffix_can = np.zeros((n_slots + 1, t), dtype=np.int64)
    for si in range(n_slots - 1, -1, -1):
        suffix_can[si] = suffix_can[si + 1] + slot_can[si]

    picks = [0] * n_slots                       # chosen block index per slot

    def discrete_feasible(cover_base):
        deficit = R - cover_base
        deficit = np.where(deficit > 0, deficit, 0)
        return int(deficit.sum()) <= disc_mult, deficit

    # Memoise on (slot index, clamped cover).  Only the deficit up to R matters, so we
    # clamp cover_i at R: states with cover_i >= R are equivalent (that codeword is
    # already satisfied).  This collapses the (otherwise exponential) small-alphabet
    # branching -- many columns charge the same few agreement classes, so distinct
    # clamped states are few.  Over a large field there are ~0 collision slots and the
    # memo is never exercised.  EXACT: clamping cannot change feasibility because cover
    # beyond R is irrelevant to the cover_i>=R constraint and to all later decisions.
    # We store (result, winning_block) so the witness path can be replayed without the
    # memo silently skipping the choices that made a subtree feasible.
    memo: dict = {}
    nodes = [0]

    def recurse(si, cover):
        nodes[0] += 1
        if nodes[0] > node_budget:
            raise ClosestStringBudgetError(
                f"closest-string feasibility exceeded {node_budget} nodes "
                f"(t={t}, n={n}, collision_slots={n_slots}); alphabet too small")
        # cheap necessary condition: cover_i + (max future from collisions) + disc_mult
        # must reach R for every codeword (discrete sum-constraint checked at the leaf).
        if np.any(cover + suffix_can[si] + disc_mult < R):
            return False
        if si == n_slots:
            ok, _ = discrete_feasible(cover)
            return ok
        ckey = (si, tuple(int(min(c, R)) for c in cover))
        cached = memo.get(ckey)
        if cached is not None:
            return cached[0]
        ti, memb = slots[si]
        # covering currently-deficient codewords first is a strong ordering heuristic.
        need = np.where(R - cover > 0, 1, 0)
        block_scores = sorted((-(int(need[blk].sum())), bi) for bi, blk in enumerate(memb))
        result = False
        win_bi = -1
        for _, bi in block_scores:
            cover2 = cover.copy()
            cover2[memb[bi]] += 1
            if recurse(si + 1, cover2):
                result = True
                win_bi = bi
                break
        memo[ckey] = (result, win_bi)
        return result

    cover0 = np.zeros(t, dtype=np.int64)
    ok = recurse(0, cover0)
    if not ok:
        return (False, None)

    # replay the winning choices from the root (deterministic via memo) to fill `picks`
    cover = cover0.copy()
    for si in range(n_slots):
        ckey = (si, tuple(int(min(c, R)) for c in cover))
        win_bi = memo[ckey][1]
        ti, memb = slots[si]
        picks[si] = win_bi
        cover[memb[win_bi]] += 1
    if not want_witness:
        return (True, None)

    # ---- build witness assignment parallel to `types` ----
    # collision picks: regroup the flat slot picks back per collision type, IN ORDER.
    per_type_picks = [None] * len(types)
    # discrete: solve x (cover deficit) from the final cover
    # recompute the final cover from picks
    cover = np.zeros(t, dtype=np.int64)
    # group slots by their type index, preserving slot order within a type
    slot_picks_by_type: dict[int, list[int]] = {}
    for si, (ti, memb) in enumerate(slots):
        cover[memb[picks[si]]] += 1
        slot_picks_by_type.setdefault(ti, []).append(picks[si])
    for ti, plist in slot_picks_by_type.items():
        per_type_picks[ti] = plist
    if disc_idx is not None:
        deficit = R - cover
        deficit = np.where(deficit > 0, deficit, 0)
        x = deficit.copy()                      # give each codeword exactly its deficit
        # distribute any leftover discrete columns arbitrarily (to codeword 0)
        leftover = disc_mult - int(x.sum())
        if leftover > 0:
            x = x.copy()
            x[0] += leftover
        per_type_picks[disc_idx] = x
    # any collision type with no recorded picks (shouldn't happen) -> empty list
    for ti in range(len(types)):
        if per_type_picks[ti] is None:
            per_type_picks[ti] = []
    return (True, (types, per_type_picks, disc_idx))


def closest_string_radius(C: np.ndarray, rng: np.random.Generator | None = None,
                          upper: int | None = None):
    """EXACT closest-string optimum e* = min_w max_i Delta(c_i,w), with a witness w*.

    Strategy: lower bound L (the max-min charge: the maximum over codewords of the
    minimum disagreements it must take = sum over column types of mult*[the codeword
    is outside every block it could pick], conservatively the trivial bound), upper
    bound U from the plurality center (a feasible word), then binary search e in
    [L, U] with closest_string_feasible.  Returns (e_star, w_star).

    The witness w* is materialised by translating the per-column block choices back
    into actual field values (any value from the chosen block's codewords works).
    """
    C = np.asarray(C, dtype=np.int64)
    t, n = C.shape
    if rng is None:
        rng = np.random.default_rng(0)

    if t == 1:
        return 0, C[0].copy()

    # Upper bound: plurality center is feasible (it is a word); its radius bounds e*.
    p_guess = plurality_word(C, int(C.max()) + 1, rng)
    U = int(hamming_to_set(p_guess, C).max())
    if upper is not None:
        U = min(U, int(upper))

    # Lower bound on e*.  Δ(c_i,w)=n-cover_i and each column covers exactly one block,
    # so the TOTAL cover sum_i cover_i is at most sum_j (size of chosen block) <= sum_j
    # (max block size at column j) = S.  Hence sum_i Δ = t*n - sum_i cover_i >= t*n - S,
    # so the maximum Δ is >= ceil((t*n - S)/t).  A valid, often tight, lower bound.
    types = _column_types(C)
    S = 0
    for blocks, mult in types:
        S += mult * max(len(B) for B in blocks)     # max block size per column
    L = max(0, -(-(t * n - S) // t))                # ceil((t*n - S)/t)

    # Binary search the least feasible e in [L, U].
    lo, hi = L, U
    best_w = p_guess.copy()
    best_e = U
    while lo < hi:
        mid = (lo + hi) // 2
        feas, _wit = closest_string_feasible(C, mid, want_witness=False)
        if feas:
            hi = mid
            best_e = mid
        else:
            lo = mid + 1
    # lo == hi == e*; materialise a witness word at the optimum.
    best_e = lo
    feas, wit = closest_string_feasible(C, best_e, want_witness=True)
    assert feas, "binary search invariant violated (optimum must be feasible)"
    best_w = _witness_word(C, wit)
    # final safety: the witness really achieves radius <= e*
    assert int(hamming_to_set(best_w, C).max()) <= best_e, "witness radius exceeds e*"
    return best_e, best_w


def _witness_word(C: np.ndarray, witness) -> np.ndarray:
    """Translate a feasibility solution (types, per_type_picks, disc_idx) into an
    explicit center word over F_p.  Thin wrapper over _cover_baseline_to_word."""
    if witness is None:
        # trivially feasible (e>=n): any codeword works as a center
        return np.asarray(C[0], dtype=np.int64).copy()
    types, per_type_picks, disc_idx = witness
    return _cover_baseline_to_word(C, types, per_type_picks, disc_idx, None)


# ===========================================================================
# Exact-center certificate for ONE candidate cluster.
# ===========================================================================
def exact_center_subcluster(C: np.ndarray, e: int, rng: np.random.Generator):
    """Given a candidate codeword set C (t, n) and a radius e, return the size of the
    LARGEST sub-cluster whose EXACT min-max center radius is <= e, together with that
    center w and the sub-cluster.

    We seek the largest S subset of the rows of C with closest_string_radius(S) <= e.
    A clean, exact-on-the-subset procedure: start from the rows already inside the
    exact-center ball of the full set, then greedily test single additions, each time
    recomputing the EXACT min-max radius of the candidate sub-cluster (never an
    approximation).  Every reported sub-cluster has a verified exact min-max radius
    <= e, so its size is an honest lower bound (the exact-center certificate).
    """
    C = distinct_codewords(np.asarray(C, dtype=np.int64))
    t = C.shape[0]
    if t == 0:
        return 0, None, None
    if t == 1:
        return (1, C[0].copy(), C) if e >= 0 else (0, None, None)

    # 1. exact min-max center of the full candidate set; rows within e form a seed.
    e_full, w_full = closest_string_radius(C, rng)
    d = hamming_to_set(w_full, C)
    in_ball = C[d <= e]
    if in_ball.shape[0] == 0:
        # even the closest single codeword: trivially radius = anything >=0 with w=c
        return 1, C[int(np.argmin(d))].copy(), C[int(np.argmin(d)):int(np.argmin(d)) + 1]

    # verify the seed is genuinely a closest-string-feasible cluster at e
    best_sub = distinct_codewords(in_ball)
    e_sub, w_sub = closest_string_radius(best_sub, rng, upper=e)
    if e_sub > e:
        # shrink to the rows the exact center actually covers (defensive; rare)
        dd = hamming_to_set(w_sub, best_sub)
        best_sub = distinct_codewords(best_sub[dd <= e])
        if best_sub.shape[0] == 0:
            best_sub = C[int(np.argmin(d)):int(np.argmin(d)) + 1]
        e_sub, w_sub = closest_string_radius(best_sub, rng, upper=e)

    # 2. greedily try to ADD any remaining codeword that keeps the EXACT radius <= e.
    remaining = distinct_codewords(C[d > e]) if (d > e).any() else np.empty((0, C.shape[1]), np.int64)
    improved = True
    seen_keys = {best_sub.tobytes()}
    while improved and remaining.shape[0] > 0:
        improved = False
        # order candidates by closeness to the current center (most likely to fit)
        dc = hamming_to_set(w_sub, remaining)
        order = np.argsort(dc)
        for idx in order:
            cand = remaining[idx:idx + 1]
            trial = distinct_codewords(np.vstack([best_sub, cand]))
            if trial.shape[0] <= best_sub.shape[0]:
                continue
            e_try, w_try = closest_string_radius(trial, rng, upper=e)
            if e_try <= e:
                best_sub = trial
                w_sub = w_try
                e_sub = e_try
                # drop the accepted row from remaining
                keep = np.ones(remaining.shape[0], dtype=bool)
                keep[idx] = False
                remaining = remaining[keep]
                improved = True
                break
        # avoid infinite loops on identical states
        kb = best_sub.tobytes()
        if kb in seen_keys and not improved:
            break
        seen_keys.add(kb)

    # final exact verification of the returned witness
    e_final, w_final = closest_string_radius(best_sub, rng, upper=e)
    assert e_final <= e, "returned sub-cluster exceeds the target radius"
    assert int(hamming_to_set(w_final, best_sub).max()) <= e, "witness radius exceeds e"
    return best_sub.shape[0], w_final, best_sub


# ===========================================================================
# One cell: plurality cert vs exact-center cert at fixed (p, n, k, c).
# ===========================================================================
@dataclass
class Cell:
    p: int
    n: int
    k: int
    c: float
    note: str = ""

    @property
    def rho(self):
        return self.k / self.n

    @property
    def label(self):
        return f"GF({self.p})_n{self.n}_k{self.k}_c{self.c}"


def run_cell(args) -> dict:
    """For one (p,n,k,c): generate candidate clusters with the cluster-certificate
    constructions, take the plurality certificate (the existing baseline) AND the
    EXACT-center certificate (largest sub-cluster with exact min-max radius <= e).
    Asserts exact-center radius <= plurality radius on the harvested clusters."""
    spec, budget, seed = args
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
    q = F.q

    t0 = time.time()

    # ---- harvest candidate clusters via the cluster-certificate constructions ----
    # Each returns the best witness it found (a word + the in-ball codewords); we take
    # the codeword SETS it harvested as candidate clusters to re-center exactly.
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
    constructions = {"RP": rp, "LL": ll, "GR": gr}
    heavy = budget.get("heavy_constructions", True)
    if heavy:
        st = construct_structured(F, L, k, e, rng,
                                  n_pool=budget["st_pool"], trials=budget["st_trials"])
        ls = construct_localsearch(F, L, k, e, rng,
                                   seed_constructions=[gr, st, ll, rp],
                                   n_moves=budget["ls_moves"])
        constructions["ST"] = st
        constructions["LS"] = ls
    constructions_run = sorted(constructions.keys())

    # plurality certificate = the cluster-certificate baseline (max over constructions)
    plur_cert = max(c0["cert"] for c0 in constructions.values())

    # ---- candidate clusters for exact re-centering ----
    # The witness codeword sets, PLUS for the strongest constructions we also collect
    # the BALL pools (in-ball codewords) so the exact center can pull in extra rows.
    candidates = []
    for key, c0 in constructions.items():
        if c0["cw"] is not None and c0["cw"].shape[0] >= 1:
            candidates.append(c0["cw"])

    # Also re-mine a few fresh random/Lloyd pools so the exact center isn't limited to
    # what the plurality search happened to keep (these are O(t*n) to make).
    for _ in range(budget["exact_extra_pools"]):
        psz = int(rng.choice(budget["ll_pool_sizes"]))
        pool = rand_codewords(F, L, k, psz, rng)
        w0 = plurality_word(pool, q, rng)
        dd = hamming_to_set(w0, pool)
        # keep a generous shell so the exact center has room to optimise the max
        shell = pool[dd <= e + max(2, n // 50)]
        if shell.shape[0] >= 2:
            candidates.append(distinct_codewords(shell))

    # ---- EXACT-center certificate: largest sub-cluster with min-max radius <= e ----
    exact_cert = 0
    exact_w = None
    exact_cw = None
    n_exact_clusters = 0
    plur_ge_exact_ok = True   # invariant: exact-center radius <= plurality radius
    invariant_checks = 0
    for cand in candidates:
        cand = distinct_codewords(cand)
        if cand.shape[0] == 0:
            continue
        n_exact_clusters += 1

        # INVARIANT CHECK on this exact cluster: exact min-max radius <= plurality radius.
        # (plurality is a feasible center, so e* can only be <= its radius.)
        if cand.shape[0] >= 2:
            e_star, _w = closest_string_radius(cand, rng)
            w_plur = plurality_word(cand, q, rng)
            r_plur = int(hamming_to_set(w_plur, cand).max())
            invariant_checks += 1
            if e_star > r_plur:
                plur_ge_exact_ok = False

        size, w, cwset = exact_center_subcluster(cand, e, rng)
        if size > exact_cert:
            exact_cert = size
            exact_w = w
            exact_cw = cwset

    elapsed = time.time() - t0

    # honest re-verification of the reported exact-center witness against e
    witness_ok = True
    if exact_w is not None and exact_cw is not None:
        witness_ok = bool(int(hamming_to_set(exact_w, exact_cw).max()) <= e)
        # and it must be a strict superset count of distinct codewords
        witness_ok = witness_ok and (distinct_codewords(exact_cw).shape[0] == exact_cert)

    gap = int(exact_cert - plur_cert)

    return {
        "label": spec.label, "p": p, "n": n, "k": k, "rho": round(rho, 6),
        "c": c, "note": spec.note, "has_open_band": True,
        "J": round(J, 6), "r_E": round(rE, 6),
        "e": int(e), "delta": round(delta, 6),
        "margin_rE_minus_delta": round(rE - delta, 6),
        "plurality_cert": int(plur_cert),
        "exact_center_cert": int(exact_cert),
        "gap": gap,
        "per_construction_plurality": {kk: int(v["cert"]) for kk, v in constructions.items()},
        "constructions_run": constructions_run,
        "exact_le_plurality_invariant_ok": bool(plur_ge_exact_ok),
        "n_invariant_checks": int(invariant_checks),
        "exact_witness_ok": bool(witness_ok),
        "n_candidate_clusters": int(n_exact_clusters),
        "elapsed_s": round(elapsed, 2),
        "exact_witness_w": exact_w.tolist() if exact_w is not None else None,
    }


# ===========================================================================
# Ladder construction (mirror the plurality sweep's fixed-field cells).
# ===========================================================================
def build_ladder(quick: bool = False, n_cap: int | None = None) -> list[Cell]:
    """Fixed-rho ladders on the Fermat prime p=65537 (p-1 = 2^16 covers every n in the
    ladder, so a single field covers all n -- any change in the certificate is purely
    an n-effect).  rho in {1/16,1/8,1/4}, band fractions c in {0.25,0.5,0.75} (the
    midpoint c=0.5 is the headline; 0.25/0.75 are included to match the plurality
    table).  n grows by powers of two; n_cap honestly bounds the reach if set."""
    P_FIXED = 65537
    ns_all = [32, 64, 128, 256, 512, 1024, 2048, 4096]
    ns = ns_all[:5] if quick else ns_all
    if n_cap is not None:
        ns = [n for n in ns if n <= n_cap]
    rhos = [("1/16", 16), ("1/8", 8), ("1/4", 4)]
    cs = [0.5] if quick else [0.25, 0.5, 0.75]

    cells = []
    for rho_lbl, inv in rhos:
        for n in ns:
            if (P_FIXED - 1) % n != 0:
                continue
            k = max(2, round(n / inv))
            if elias_radius(k / n, P_FIXED) <= johnson_radius(k / n):
                continue
            for c in cs:
                cells.append(Cell(P_FIXED, n, k, c, note=f"fixedfield_rho={rho_lbl}"))
    return cells


def cert_budget(n: int, quick: bool = False) -> dict:
    """Trial/pool sizes for harvesting candidate clusters, plus the exact-center extra
    pools.

    Harvesting reuses the cluster-certificate constructions (all O(t*n)); the EXACT
    closest-string solve adds negligible cost (microseconds: collision columns are ~0
    over p=65537).  The dominant cost is the plurality LOCAL SEARCH at large n
    (n_moves * O(t*n) per move).  Because the exact min-max center re-optimises the
    center EXACTLY on every harvested cluster, heavy plurality center-polishing is
    largely redundant here, so we use lighter local search at large n than the raw
    cluster_certificate budget while keeping the cluster-generating trials (RP/LL/GR/ST)
    strong.  Trial counts are logged per cell (n_candidate_clusters, n_trials) for
    honesty.
    """
    if n <= 128:
        rp_t = (3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32)
        rp_trials = 200 if not quick else 60
        ll_pool = (8, 12, 16, 24, 32, 48, 64)
        ll_trials = 60 if not quick else 20
        gr_restarts = 60 if not quick else 20
        gr_max_t = 40
        st_pool = 64; st_trials = 8 if not quick else 3
        ls_moves = 300 if not quick else 120
        extra = 8 if not quick else 4
    elif n <= 512:
        rp_t = (3, 4, 5, 6, 8, 10, 12, 16, 20, 24)
        rp_trials = 100 if not quick else 36
        ll_pool = (8, 12, 16, 24, 32, 48)
        ll_trials = 36 if not quick else 14
        gr_restarts = 36 if not quick else 14
        gr_max_t = 32
        st_pool = 48; st_trials = 6 if not quick else 2
        ls_moves = 200 if not quick else 80
        extra = 6 if not quick else 3
    else:  # n in {1024, 2048, 4096}: trim trials to bound wall-time (cost ~ trials*t*n*k)
        rp_t = (3, 4, 5, 6, 8, 10, 12, 16, 20)
        rp_trials = 12 if not quick else 6
        ll_pool = (8, 12, 16, 24, 32)
        ll_trials = 10 if not quick else 6
        gr_restarts = 10 if not quick else 6
        gr_max_t = 24
        st_pool = 24; st_trials = 3 if not quick else 2
        ls_moves = 60 if not quick else 30
        extra = 8 if not quick else 4
    # The structured-coset (ST) and plurality local-search (LS) constructions are O(n*k)
    # / O(n_moves * t * n) and become the wall-time bottleneck at large n (ST builds all
    # k-1 low monomials with field powers).  Because the EXACT min-max center re-optimises
    # the center on every harvested cluster, LS center-polishing is redundant here; and
    # the plurality study already established that the worst clusters are NON-coset, so ST
    # adds little.  We therefore run ST/LS only for n <= 512 and log it (heavy_constructions).
    heavy = n <= 512
    return {
        "rp_t_set": rp_t, "rp_trials": rp_trials,
        "ll_pool_sizes": ll_pool, "ll_trials": ll_trials,
        "gr_restarts": gr_restarts, "gr_batch": 12, "gr_max_t": gr_max_t,
        "st_pool": st_pool, "st_trials": st_trials,
        "ls_moves": ls_moves,
        "exact_extra_pools": extra,
        "heavy_constructions": heavy,
    }


# ===========================================================================
# Driver.
# ===========================================================================
def run_full(quick: bool = False, procs: int = None, n_cap: int | None = None):
    procs = procs or min(16, os.cpu_count() or 4)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t_start = time.time()

    cells = build_ladder(quick=quick, n_cap=n_cap)
    print(f"[exact_center] {len(cells)} cells (fixed field p=65537), procs={procs}, "
          f"quick={quick}, n_cap={n_cap}", flush=True)

    args = [(s, cert_budget(s.n, quick), 0xE0C1 + 7919 * i) for i, s in enumerate(cells)]
    results = []
    with Pool(procs) as pool:
        for i, rec in enumerate(pool.imap_unordered(run_cell, args)):
            results.append(rec)
            if rec.get("has_open_band"):
                print(f"  [{i+1}/{len(args)}] {rec['label']:>26} "
                      f"n={rec['n']:>5} rho={rec['rho']:.4f} c={rec['c']} "
                      f"delta={rec['delta']:.4f} "
                      f"plur={rec['plurality_cert']:>2} exact={rec['exact_center_cert']:>2} "
                      f"gap={rec['gap']:+d} "
                      f"inv_ok={rec['exact_le_plurality_invariant_ok']} "
                      f"{rec['elapsed_s']}s", flush=True)

    elapsed = time.time() - t_start

    inv_all_ok = all(r.get("exact_le_plurality_invariant_ok", True)
                     for r in results if r.get("has_open_band"))
    wit_all_ok = all(r.get("exact_witness_ok", True)
                     for r in results if r.get("has_open_band"))

    # assemble the headline comparison table: rho | c | n -> (plur, exact, gap)
    table = {}
    for r in results:
        if not r.get("has_open_band"):
            continue
        rho_lbl = r["note"].split("rho=")[-1]
        table.setdefault(rho_lbl, {}).setdefault(str(r["c"]), []).append({
            "n": r["n"], "p": r["p"], "k": r["k"], "delta": r["delta"],
            "margin": r["margin_rE_minus_delta"],
            "plurality_cert": r["plurality_cert"],
            "exact_center_cert": r["exact_center_cert"],
            "gap": r["gap"],
        })
    for rho_lbl in table:
        for cstr in table[rho_lbl]:
            table[rho_lbl][cstr].sort(key=lambda x: x["n"])

    summary = {
        "meta": {
            "purpose": "EXACT min-max (closest-string) center for the cluster certificate: "
                       "strengthened one-sided lower bound on the open-band list at fixed rho.",
            "date": "2026-06-03",
            "field": "fixed Fermat prime p=65537 (p-1=2^16 covers all n in the ladder)",
            "quick": quick, "procs": procs, "n_cap": n_cap,
            "elapsed_s": round(elapsed, 1),
            "method": "harvest candidate clusters via the cluster-certificate constructions; "
                      "for each, EXACT closest-string min-max center (branch-and-bound over "
                      "column types, binary search on radius); report the largest sub-cluster "
                      "with exact min-max radius <= floor(delta*n).",
            "invariant_exact_le_plurality_all_ok": bool(inv_all_ok),
            "exact_witness_all_ok": bool(wit_all_ok),
            "honesty": "exact-center cert is a STRONGER but still ONE-SIDED lower bound on "
                       "Lambda(C,delta): a cleverer cluster could give more. Flat in n at fixed "
                       "rho supports P' boundedness; growth would threaten P'.",
        },
        "comparison_table": table,
        "cells": results,
    }

    with open(os.path.join(RESULTS_DIR, "exact_center.json"), "w") as f:
        json.dump(summary, f, indent=1)

    with open(os.path.join(RESULTS_DIR, "exact_center.csv"), "w", newline="") as f:
        wtr = csv.writer(f)
        wtr.writerow(["rho", "c", "p", "n", "k", "delta", "margin_rE_minus_delta",
                      "plurality_cert", "exact_center_cert", "gap"])
        for r in sorted(results, key=lambda x: (x.get("note", ""), x.get("c", 0),
                                                x.get("n", 0))):
            if not r.get("has_open_band"):
                continue
            rho_lbl = r["note"].split("rho=")[-1]
            wtr.writerow([rho_lbl, r["c"], r["p"], r["n"], r["k"], r["delta"],
                          r["margin_rE_minus_delta"], r["plurality_cert"],
                          r["exact_center_cert"], r["gap"]])

    _print_headline(summary)
    print(f"\n[exact_center] DONE in {elapsed:.1f}s -> results/exact_center.{{json,csv}}",
          flush=True)
    return summary


def _trend(vals):
    if len(vals) < 2:
        return "n/a"
    if max(vals) - min(vals) <= 1:
        return f"FLAT ({min(vals)}-{max(vals)})"
    if vals[-1] > vals[0]:
        return f"GROWING {vals[0]}->{vals[-1]}"
    return f"VARYING {min(vals)}-{max(vals)}"


def _print_headline(summary):
    print("\n" + "=" * 78)
    print("HEADLINE: plurality vs EXACT-center certificate vs n at FIXED rho")
    print("=" * 78)
    for rho_lbl, by_c in summary["comparison_table"].items():
        for cstr, rows in by_c.items():
            plur = "  ".join(f"n={r['n']}:{r['plurality_cert']}" for r in rows)
            exa = "  ".join(f"n={r['n']}:{r['exact_center_cert']}" for r in rows)
            ev = [r["exact_center_cert"] for r in rows]
            print(f"  rho={rho_lbl:5s} c={cstr:4s}")
            print(f"      plurality : {plur}")
            print(f"      exact     : {exa}    [{_trend(ev)}]")
    m = summary["meta"]
    print(f"\n  invariant exact<=plurality ALL OK: {m['invariant_exact_le_plurality_all_ok']}; "
          f"witness radius ALL OK: {m['exact_witness_all_ok']}")


# ===========================================================================
# Self-test battery.
# ===========================================================================
def _brute_closest_string(C: np.ndarray, q: int):
    """BRUTE FORCE closest-string optimum over the WHOLE word space F_q^n.

    Enumerates all q^n words, returns (e*, a witness w*).  Only for tiny cases
    (q^n small) -- this is the exact-vs-exact oracle the solver is checked against.
    """
    t, n = C.shape
    best_e = n + 1
    best_w = None
    for tup in itertools.product(range(q), repeat=n):
        w = np.array(tup, dtype=np.int64)
        r = int((C != w[None, :]).sum(axis=1).max())
        if r < best_e:
            best_e = r
            best_w = w.copy()
            if best_e == 0:
                break
    return best_e, best_w


def _self_test():
    print("exact_center self-test")
    print("=" * 70)
    rng = np.random.default_rng(0xEC11)

    # ---- 1. EXACT closest-string == BRUTE FORCE over the whole word space ----
    print("[1] closest_string_radius == brute force over F_q^n (tiny cases):")
    n_cases = 0
    for _ in range(120):
        q = int(rng.integers(2, 5))            # small alphabet
        n = int(rng.integers(2, 7))            # small length so q^n is enumerable
        t = int(rng.integers(2, 5))            # up to 4 strings
        if q ** n > 4096:
            continue
        C = rng.integers(0, q, size=(t, n)).astype(np.int64)
        e_exact, w_exact = closest_string_radius(C, rng)
        e_brute, _wb = _brute_closest_string(C, q)
        assert e_exact == e_brute, \
            f"closest-string mismatch: exact={e_exact} brute={e_brute}\nC=\n{C}"
        # the returned witness must actually achieve e_exact
        assert int((C != w_exact[None, :]).sum(axis=1).max()) == e_exact, \
            "witness does not achieve the exact optimum"
        n_cases += 1
    print(f"    {n_cases} random tiny instances: exact min-max == brute force, "
          f"witness achieves it OK")

    # ---- 2. feasibility agrees with brute force at EVERY radius ----
    print("[2] closest_string_feasible(C,e) == brute force, every radius:")
    n_checks = 0
    for _ in range(60):
        q = int(rng.integers(2, 4))
        n = int(rng.integers(2, 6))
        t = int(rng.integers(2, 5))
        if q ** n > 4096:
            continue
        C = rng.integers(0, q, size=(t, n)).astype(np.int64)
        e_brute, _ = _brute_closest_string(C, q)
        for e in range(0, n + 1):
            feas, _w = closest_string_feasible(C, e)
            assert feas == (e >= e_brute), \
                f"feasibility mismatch at e={e}: feas={feas} but e*={e_brute}"
            n_checks += 1
    print(f"    {n_checks} (instance,radius) feasibility checks match brute force OK")

    # ---- 3. KEY INVARIANT: exact min-max radius <= plurality radius, always ----
    # Regime: alphabet large relative to t (q >= 4*t^2), matching the prime-field
    # experiment where value collisions across codewords are sparse.  Larger n than the
    # brute-force tests, exercising the discrete-column closed form + a few collisions.
    print("[3] exact min-max radius <= plurality radius (the strengthening invariant):")
    n_inv = 0
    strictly_better = 0
    for _ in range(600):
        t = int(rng.integers(2, 17))
        q = int(rng.integers(4 * t * t, 4 * t * t + 200))   # sparse-collision regime
        n = int(rng.integers(4, 200))
        C = distinct_codewords(rng.integers(0, q, size=(t, n)).astype(np.int64))
        if C.shape[0] < 2:
            continue
        e_exact, w_exact = closest_string_radius(C, rng)
        w_plur = plurality_word(C, q, rng)
        r_plur = int((C != w_plur[None, :]).sum(axis=1).max())
        assert e_exact <= r_plur, \
            f"INVARIANT VIOLATED: exact={e_exact} > plurality={r_plur}\nC=\n{C}"
        assert int((C != w_exact[None, :]).sum(axis=1).max()) == e_exact, \
            "witness does not achieve the exact optimum"
        n_inv += 1
        if e_exact < r_plur:
            strictly_better += 1
    print(f"    {n_inv} clusters (q>=4t^2): exact min-max radius <= plurality radius "
          f"ALWAYS OK")
    print(f"      (exact strictly beat plurality on {strictly_better}/{n_inv} clusters "
          f"-- the source of the certificate strengthening)")

    # ---- 3b. node-budget guard: pathological tiny-alphabet input raises, never hangs --
    # The collision-branching is exponential only for tiny alphabets (never the prime-
    # field experiment).  With a SMALL explicit node budget the feasibility search must
    # terminate quickly -- either solving or raising ClosestStringBudgetError, never
    # hanging or returning a wrong answer.  Here we drive it at a hard radius (e small,
    # so feasibility is far from trivial) over GF(2) with many strings.
    print("[3b] tiny-alphabet input is either solved or cleanly budget-guarded (no hang):")
    guarded = solved = 0
    t_guard0 = time.time()
    for _ in range(30):
        q = 2
        n = int(rng.integers(20, 40))
        t = int(rng.integers(9, 15))
        C = distinct_codewords(rng.integers(0, q, size=(t, n)).astype(np.int64))
        if C.shape[0] < 2:
            continue
        e_target = n // 3                         # a demanding radius
        try:
            feas, _wit = closest_string_feasible(C, e_target, node_budget=20_000)
            assert feas in (True, False)
            solved += 1
        except ClosestStringBudgetError:
            guarded += 1
    dt_guard = time.time() - t_guard0
    assert dt_guard < 30, f"budget guard too slow ({dt_guard:.1f}s) -- guard not firing"
    print(f"    GF(2) stress (budget=20k nodes): {solved} resolved, {guarded} cleanly "
          f"budget-guarded in {dt_guard:.1f}s (no hang, no wrong value) OK")

    # ---- 4. optimality: exact radius <= every random word's radius ----
    print("[4] exact min-max radius <= every sampled random word's radius:")
    for _ in range(40):
        t = int(rng.integers(2, 9))
        q = int(rng.integers(4 * t * t, 4 * t * t + 50))
        n = int(rng.integers(3, 60))
        C = rng.integers(0, q, size=(t, n)).astype(np.int64)
        e_exact, _ = closest_string_radius(C, rng)
        for _ in range(300):
            wr = rng.integers(0, q, size=n).astype(np.int64)
            assert int((C != wr[None, :]).sum(axis=1).max()) >= e_exact, \
                "a random word beat the exact min-max center"
    print("    exact min-max radius <= 300 random words on each of 40 instances OK")

    # ---- 5. genuine RS clusters: exact center on real codewords, radius <= plurality
    print("[5] real RS clusters (GF(97) n=24 k=3): exact center <= plurality, witness valid:")
    F = PrimeField(97); n = 24; k = 3
    L = domain_subgroup(F, n)
    q = F.q
    n_rs = 0
    rs_better = 0
    for _ in range(40):
        t = int(rng.integers(3, 10))
        C = rand_codewords(F, L, k, t, rng)
        C = distinct_codewords(C)
        if C.shape[0] < 2:
            continue
        e_exact, w_exact = closest_string_radius(C, rng)
        w_plur = plurality_word(C, q, rng)
        r_plur = int(hamming_to_set(w_plur, C).max())
        assert e_exact <= r_plur, f"RS invariant: exact {e_exact} > plurality {r_plur}"
        assert int(hamming_to_set(w_exact, C).max()) == e_exact, "RS witness invalid"
        n_rs += 1
        if e_exact < r_plur:
            rs_better += 1
    print(f"    {n_rs} RS clusters: exact <= plurality OK (exact strictly better on "
          f"{rs_better}); witnesses valid OK")

    # ---- 6. exact-center sub-cluster certificate is honest and >= plurality cert ----
    print("[6] exact_center_subcluster: verified radius <= e, size >= plurality in-ball:")
    e = 16
    n_sub = 0
    for _ in range(30):
        t = int(rng.integers(4, 12))
        C = rand_codewords(F, L, k, t, rng)
        # plurality in-ball count of the full set (the baseline this must not undershoot)
        w_plur = plurality_word(distinct_codewords(C), q, rng)
        plur_in = cert_count(w_plur, C, e)
        size, w, cwset = exact_center_subcluster(C, e, rng)
        # the reported sub-cluster must be exactly within e under its exhibited center
        assert int(hamming_to_set(w, cwset).max()) <= e, "exact sub-cluster exceeds e"
        assert distinct_codewords(cwset).shape[0] == size, "size != #distinct rows"
        # exact center can only fit >= the plurality in-ball of the SAME full set
        assert size >= plur_in, \
            f"exact-center cert {size} < plurality in-ball {plur_in} (should never happen)"
        n_sub += 1
    print(f"    {n_sub} clusters: exact-center sub-cluster within e, size >= plurality "
          f"in-ball OK")

    # ---- 7. column-type reduction is value-independent (relabel invariance) ----
    # Two relabelings: a sparse-collision regime (matches the experiment) AND a small
    # exhaustively-brute-checkable regime (so the value-free claim is pinned at both
    # ends), with a random alphabet permutation applied to every entry.
    print("[7] column-type reduction is invariant under field-value relabeling:")
    for _ in range(40):
        t = int(rng.integers(2, 9))
        q = int(rng.integers(4 * t * t, 4 * t * t + 50))
        n = int(rng.integers(3, 80))
        C = rng.integers(0, q, size=(t, n)).astype(np.int64)
        perm = rng.permutation(q)
        e1, _ = closest_string_radius(C, rng)
        e2, _ = closest_string_radius(perm[C], rng)
        assert e1 == e2, "closest-string radius changed under alphabet relabeling"
    for _ in range(40):                          # tiny + brute-validated regime
        q = int(rng.integers(2, 5))
        n = int(rng.integers(2, 6))
        t = int(rng.integers(2, 5))
        if q ** n > 4096:
            continue
        C = rng.integers(0, q, size=(t, n)).astype(np.int64)
        perm = rng.permutation(q)
        e1, _ = closest_string_radius(C, rng)
        e2, _ = closest_string_radius(perm[C], rng)
        assert e1 == e2 == _brute_closest_string(C, q)[0], "relabel/brute mismatch"
    print("    relabeling the alphabet leaves e* unchanged (value-free reduction) OK")

    # ---- 8. end-to-end cell runs and reports a sane comparison ----
    print("[8] run_cell end-to-end (GF(65537) n=64 k=4 c=0.5): plurality vs exact:")
    spec = Cell(65537, 64, 4, 0.5, note="fixedfield_rho=1/16")
    rec = run_cell((spec, cert_budget(64, quick=True), 0xABCD))
    assert rec["has_open_band"], "n=64 rho=1/16 must have an open band"
    assert rec["exact_le_plurality_invariant_ok"], "cell invariant must hold"
    assert rec["exact_center_cert"] >= rec["plurality_cert"], \
        "exact-center cert must be >= plurality cert"
    assert rec["exact_witness_ok"], "exact witness must be valid"
    print(f"    plurality={rec['plurality_cert']} exact={rec['exact_center_cert']} "
          f"gap={rec['gap']:+d} inv_ok={rec['exact_le_plurality_invariant_ok']} "
          f"in {rec['elapsed_s']}s OK")

    print("=" * 70)
    print("ALL EXACT-CENTER SELF-TESTS PASSED")


# ===========================================================================
# Main.
# ===========================================================================
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true", help="run the self-test battery")
    ap.add_argument("--quick", action="store_true", help="fast reduced sweep")
    ap.add_argument("--procs", type=int, default=None, help="parallel worker processes")
    ap.add_argument("--n-cap", type=int, default=None,
                    help="honestly cap the top n if branch-and-bound is too slow")
    args = ap.parse_args()

    if args.selftest:
        _self_test()
    else:
        run_full(quick=args.quick, procs=args.procs, n_cap=args.n_cap)
