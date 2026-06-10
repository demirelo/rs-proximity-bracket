"""
search_interleaved.py -- list sizes for interleaved Reed-Solomon codes.

The m-interleaved code C^{equiv m} (a.k.a. the "equivalent code") has as codewords
the m-tuples (c_1,...,c_m) of codewords c_j in C.  The metric is COLUMN-wise: the
distance between two m-tuples W=(w_1,...,w_m) and (c_1,...,c_m) is

    d(W, c) = #{ i in L : exists j with w_j[i] != c_j[i] }
            = n - #{ i in L : w_j[i] = c_j[i] for ALL j }.

So a column i is an "agreement column" iff every coordinate agrees there.  The
list at radius delta for a target W is

    Lambda(C^{equiv m}, W, delta) = { (c_1,...,c_m) : d(W,c) <= delta*n }.

Equivalently (c_1,...,c_m) is in the list iff there is a set S, |S| >= (1-delta)n,
on which every c_j agrees with w_j simultaneously.  This is the m-fold analogue of
the common-agreement-set computation in search_bad_lines.py.

Large interleaved lists hurt SNARK soundness (the soundness error carries a term
|Lambda(C^{equiv m},delta)| / |F|).  We search for target words W maximizing the
list size, for small m (2,3), and report max list size vs delta.

EXACT list-size computation
---------------------------
For each coordinate j, precompute the agreement mask A_j(c_j) = {i : w_j[i]=c_j[i]}
of every codeword c_j with the target w_j (bit-packed).  Then

    list size = #{ (c_1,...,c_m) : popcount( A_1(c_1) & ... & A_m(c_m) ) >= (1-delta)n }.

We enumerate codewords for coordinate 1, prune by popcount, and for surviving
prefixes intersect with coordinate-2 masks, etc.  For m=2,3 and the atlas field
sizes this is a handful of vectorized passes.  Branch-and-bound: a prefix whose
partial-intersection popcount already < (1-delta)n can be pruned (intersection
only shrinks), making the search fast.

We count the list size for SEVERAL delta values at once (a list is monotone in
delta), reusing the agreement masks.
"""

from __future__ import annotations

import numpy as np

from ff import FiniteField
from rs import (CodewordBook, build_codeword_book, encode, random_word,
                random_codeword, codeword_plus_noise)
from search_bad_lines import _agreement_bits, _POPCOUNT8


def _popcount_packed(packed_rows: np.ndarray) -> np.ndarray:
    """Popcount each row of a (R, B) uint8 packed-bit matrix -> (R,) int."""
    return _POPCOUNT8[packed_rows].sum(axis=1)


def interleaved_list_sizes(book: CodewordBook, target: list[np.ndarray],
                           deltas: list[float],
                           hard_cap_prefix: int = 4_000_000
                           ) -> tuple[dict, dict]:
    """Exact list sizes |Lambda(C^{equiv m}, target, delta)| for each delta.

    target is a list of m words (w_1,...,w_m), each length n.  Returns
    (sizes, info) where sizes[delta] is the integer list size and info documents
    exactness.  Also returns, for the smallest delta, a witness count.

    Algorithm (m up to 3): branch-and-bound over codeword tuples using bit-packed
    agreement masks and the requirement popcount(intersection) >= tau where
    tau = ceil((1-delta_max)*n) for pruning (the loosest, i.e. smallest tau, so we
    keep every tuple that could be in ANY requested list, then bucket by delta).
    """
    m = len(target)
    n = book.n
    ncw = book.num_codewords
    assert 1 <= m <= 3, "this searcher supports m in {1,2,3}"

    # Agreement masks (bit-packed) and popcounts for each coordinate.
    packs = []
    pops = []
    for j in range(m):
        pj, packj = _agreement_bits(book, target[j])
        packs.append(packj)
        pops.append(pj)

    # tau for each delta: a tuple is in Lambda(delta) iff intersection popcount
    # >= ceil((1-delta)n).  Smallest tau corresponds to largest delta.
    taus = {d: int(np.ceil((1.0 - d) * n - 1e-9)) for d in deltas}
    tau_min = min(taus.values())          # loosest requirement, for pruning

    # Histogram of intersection popcounts over ALL in-list tuples; from this we
    # read off every delta's list size as the count of tuples with popcount>=tau.
    inter_hist = np.zeros(n + 1, dtype=np.int64)
    exact = True
    examined_prefixes = 0

    if m == 1:
        # list size = #{c : popcount(A(c)) >= tau}
        for d, tau in taus.items():
            pass
        cnt = np.bincount(pops[0], minlength=n + 1)
        # cumulative from top
        for v in range(n + 1):
            inter_hist[v] = cnt[v]
    else:
        # Sort coordinate-1 codewords by popcount desc; prune when pop1 < tau_min.
        order1 = np.argsort(-pops[0], kind="stable")
        pop1s = pops[0][order1]
        pack1s = packs[0][order1]

        # Pre-sort coordinate-2 (and 3) by popcount desc for prefix pruning.
        order2 = np.argsort(-pops[1], kind="stable")
        pop2s = pops[1][order2]
        pack2s = packs[1][order2]
        if m == 3:
            order3 = np.argsort(-pops[2], kind="stable")
            pop3s = pops[2][order3]
            pack3s = packs[2][order3]

        for i1 in range(ncw):
            if pop1s[i1] < tau_min:
                break                      # remaining c1 cannot reach tau_min
            if examined_prefixes >= hard_cap_prefix:
                exact = False
                break
            examined_prefixes += 1
            row1 = pack1s[i1]
            # Intersect with all c2 whose pop2 >= tau_min (prefix of sorted).
            n2 = int(np.searchsorted(-pop2s, -tau_min, side="right"))
            if n2 == 0:
                continue
            and12 = np.bitwise_and(pack2s[:n2], row1[None, :])    # (n2, B)
            pc12 = _popcount_packed(and12)                       # (n2,)
            keep2 = pc12 >= tau_min
            if not keep2.any():
                continue
            if m == 2:
                vals = pc12[keep2]
                # accumulate into histogram
                hh = np.bincount(vals, minlength=n + 1)
                inter_hist += hh
            else:  # m == 3
                and12_keep = and12[keep2]                          # (k2, B)
                pc12_keep = pc12[keep2]
                # For each surviving (c1,c2) prefix, intersect with c3.
                n3 = int(np.searchsorted(-pop3s, -tau_min, side="right"))
                if n3 == 0:
                    continue
                pack3_pref = pack3s[:n3]                           # (n3, B)
                for r in range(and12_keep.shape[0]):
                    rowp = and12_keep[r]
                    if pc12_keep[r] < tau_min:
                        continue
                    and123 = np.bitwise_and(pack3_pref, rowp[None, :])
                    pc123 = _popcount_packed(and123)
                    keep3 = pc123 >= tau_min
                    if keep3.any():
                        hh = np.bincount(pc123[keep3], minlength=n + 1)
                        inter_hist += hh

    # Read off list sizes: list(delta) = sum of inter_hist[v] for v >= tau(delta).
    # inter_hist currently holds, for each intersection-popcount value v, the
    # number of in-(tau_min)-list tuples with exactly that popcount.
    suffix = np.cumsum(inter_hist[::-1])[::-1]    # suffix[v] = #tuples with pop>=v
    sizes = {}
    for d, tau in taus.items():
        tau_c = min(max(tau, 0), n)
        sizes[d] = int(suffix[tau_c])

    info = {
        "exact": exact, "m": m, "num_codewords": ncw,
        "examined_prefixes": examined_prefixes,
        "tau_per_delta": taus,
    }
    return sizes, info


# ---------------------------------------------------------------------------
# Target-word generators for the interleaved search.
# ---------------------------------------------------------------------------
def gen_target_random(F, book, rng, m):
    return [random_word(F, book.n, rng) for _ in range(m)]


def gen_target_cwnoise(F, book, rng, m, e):
    return [codeword_plus_noise(book, rng, e) for _ in range(m)]


def gen_target_shared_support(F, book, rng, m, e):
    """All m coordinates are codewords corrupted on the SAME positions.  This
    maximizes the common agreement set and should produce the LARGEST lists
    (every coordinate can independently pick among codewords agreeing on the n-e
    good columns) -- the structured worst case for list size."""
    cws = [random_codeword(book, rng).copy() for _ in range(m)]
    pos = rng.choice(book.n, size=e, replace=False)
    for c in cws:
        for p in pos:
            c[p] = F.add(int(c[p]), int(rng.integers(1, F.q)))
    return cws


# ---------------------------------------------------------------------------
# Search driver.
# ---------------------------------------------------------------------------
def search_interleaved(book: CodewordBook, m: int, deltas: list[float],
                       rng: np.random.Generator,
                       n_random=60, n_cwnoise=60, n_shared=60,
                       e_frac=None, verbose=False) -> dict:
    """Search target words maximizing interleaved list size for given m.

    Returns a summary with, per delta, the max list size found and which
    generator produced it.  Sampling counts are explicit.
    """
    F = book.F
    n, k = book.n, book.k
    if e_frac is None:
        e_frac = max(deltas) * 0.9
    e_each = max(1, int(round(e_frac * n)))

    # Baseline: the all-(uncorrupted-codeword) target has list size... we also
    # track the GENERIC list size for a random target as a reference.
    best_by_delta = {d: {"size": 0, "kind": None, "info": None} for d in deltas}
    sum_sizes = {d: 0 for d in deltas}
    count = 0

    def consider(target, kind):
        nonlocal count
        sizes, info = interleaved_list_sizes(book, target, deltas)
        count += 1
        for d in deltas:
            sum_sizes[d] += sizes[d]
            if sizes[d] > best_by_delta[d]["size"]:
                best_by_delta[d] = {"size": sizes[d], "kind": kind,
                                    "info": {"exact": info["exact"]}}
        return sizes

    for _ in range(n_random):
        consider(gen_target_random(F, book, rng, m), "random")
    for _ in range(n_cwnoise):
        consider(gen_target_cwnoise(F, book, rng, m, e_each), "cw+noise")
    for _ in range(n_shared):
        consider(gen_target_shared_support(F, book, rng, m, e_each), "shared-support")

    # The single most extremal target we can write down analytically: a codeword
    # tuple corrupted on a SHARED support of size e where n-e == ceil((1-d)n).
    # For each delta, that gives a guaranteed-large list; record it explicitly.
    structured = {}
    for d in deltas:
        e = max(0, n - int(np.ceil((1.0 - d) * n - 1e-9)))   # exactly on boundary
        if e < 1:
            structured[d] = None
            continue
        tgt = gen_target_shared_support(F, book, rng, m, e)
        sizes, info = interleaved_list_sizes(book, tgt, [d])
        structured[d] = {"e": e, "size": sizes[d], "exact": info["exact"]}
        if sizes[d] > best_by_delta[d]["size"]:
            best_by_delta[d] = {"size": sizes[d], "kind": "structured-boundary",
                                "info": {"exact": info["exact"]}}

    summary = {
        "field": F.name, "q": F.q, "n": n, "k": k, "m": m, "rho": k / n,
        "johnson": 1.0 - np.sqrt(k / n), "capacity": 1.0 - k / n,
        "deltas": list(deltas),
        "sampling": {"random": n_random, "cw+noise": n_cwnoise,
                     "shared-support": n_shared,
                     "e_each": e_each,
                     "full_space_size": f"q^(k*m) targets implicitly = {F.q}^({k}*{m}); "
                                        f"we sample targets, list computed EXACTLY over all "
                                        f"q^(k*m)={F.q**(k*m)} interleaved codewords"},
        "interleaved_codewords_total": F.q ** (k * m),
        "max_list_by_delta": {str(d): best_by_delta[d] for d in deltas},
        "mean_list_by_delta": {str(d): sum_sizes[d] / max(count, 1) for d in deltas},
        "structured_boundary_target": {str(d): structured[d] for d in deltas},
    }
    if verbose:
        ds = " ".join(f"d={d}:{best_by_delta[d]['size']}" for d in deltas)
        print(f"    [{F.name} n={n} k={k} m={m}] max list sizes: {ds}")
    return summary


# ===========================================================================
# Self-test.
# ===========================================================================
def _brute_list_size(book, target, delta):
    """Independent brute-force list size: enumerate ALL interleaved codeword
    tuples and count those within delta in the columnwise metric.  Only used in
    tests on tiny codes (q^(k*m) small)."""
    m = len(target)
    n = book.n
    tau = int(np.ceil((1.0 - delta) * n - 1e-9))
    G = book.G
    ncw = book.num_codewords
    # masks per coordinate
    masks = [[(G[c] == target[j]) for c in range(ncw)] for j in range(m)]
    import itertools
    cnt = 0
    for tup in itertools.product(range(ncw), repeat=m):
        inter = masks[0][tup[0]].copy()
        for j in range(1, m):
            inter = inter & masks[j][tup[j]]
        if int(inter.sum()) >= tau:
            cnt += 1
    return cnt


def _self_test():
    from ff import PrimeField, BinaryExtensionField
    from rs import domain_full, domain_subgroup
    rng = np.random.default_rng(0x5151)
    print("search_interleaved self-test")
    print("=" * 60)

    # --- Brute-force cross-check on a TINY code, m=2 and m=3 ---
    F = PrimeField(5)
    L = domain_full(F)              # n=4
    book = build_codeword_book(F, L, 2)     # 25 codewords
    for m in [2, 3]:
        for _ in range(15):
            target = [random_word(F, book.n, rng) for _ in range(m)]
            for delta in [0.25, 0.5, 0.75, 1.0]:
                fast, info = interleaved_list_sizes(book, target, [delta])
                brute = _brute_list_size(book, target, delta)
                assert fast[delta] == brute, (
                    f"m={m} delta={delta}: fast={fast[delta]} brute={brute}")
                assert info["exact"]
        print(f"  GF(5) n=4 k=2 m={m}: fast list size == brute on 15x4 cases OK")

    # --- Another tiny code GF(2^2) (q=4), n=3, k=2, m=2 ---
    F = BinaryExtensionField(2)
    L = domain_full(F)              # n=3
    book = build_codeword_book(F, L, 2)     # 16 codewords
    for _ in range(15):
        target = [random_word(F, book.n, rng) for _ in range(2)]
        for delta in [1/3, 2/3, 1.0]:
            fast, _ = interleaved_list_sizes(book, target, [delta])
            brute = _brute_list_size(book, target, delta)
            assert fast[delta] == brute, f"GF4 m=2 d={delta}: {fast[delta]} vs {brute}"
    print("  GF(2^2) n=3 k=2 m=2: fast == brute on 15x3 cases OK")

    # --- Structural sanity: list at delta=1 is ALL interleaved codewords ---
    F = PrimeField(7)
    L = domain_full(F)              # n=6
    book = build_codeword_book(F, L, 2)
    target = [random_word(F, book.n, rng) for _ in range(2)]
    sizes, _ = interleaved_list_sizes(book, target, [1.0])
    assert sizes[1.0] == book.num_codewords ** 2, \
        f"delta=1 list must be all {book.num_codewords**2} tuples, got {sizes[1.0]}"
    print(f"  delta=1 list size = (q^k)^2 = {sizes[1.0]} OK")

    # --- Structural sanity: a codeword-tuple target has list size >= 1 (itself)
    #     at delta=0, and exactly the tuples agreeing everywhere. ---
    c0 = random_codeword(book, rng)
    c1 = random_codeword(book, rng)
    sizes, _ = interleaved_list_sizes(book, [c0, c1], [0.0])
    assert sizes[0.0] == 1, f"a codeword tuple has exactly itself at delta=0, got {sizes[0.0]}"
    print(f"  codeword-tuple target: list size at delta=0 is 1 OK")

    # --- Demo search ---
    F = BinaryExtensionField(4)
    L = domain_subgroup(F, 15)
    book = build_codeword_book(F, L, 3)     # rho=0.2
    deltas = [0.4, 0.5, 0.6]
    for m in [2, 3]:
        summary = search_interleaved(book, m, deltas, rng,
                                     n_random=20, n_cwnoise=20, n_shared=20,
                                     verbose=True)
        # list size is monotone nondecreasing in delta
        szs = [summary["max_list_by_delta"][str(d)]["size"] for d in deltas]
        assert szs == sorted(szs), f"list size not monotone in delta: {szs}"
    print("  demo m=2,3 search: list sizes monotone in delta OK")

    print("=" * 60)
    print("ALL search_interleaved SELF-TESTS PASSED")


if __name__ == "__main__":
    _self_test()
