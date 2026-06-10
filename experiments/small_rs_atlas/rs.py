"""
rs.py -- Reed-Solomon codes over small finite fields, with EXACT distance.

C = RS[F, L, k] = { (p(x))_{x in L} : p in F[X], deg p < k }, |L| = n, rate rho=k/n.
There are exactly q^k codewords (one per coefficient vector in F^k).

Correctness-first philosophy
----------------------------
We never decode.  To compute dist(w, C) we ENUMERATE all q^k codewords and take
the minimum Hamming distance.  This is exponential in k but we deliberately keep
q^k <= ~2e6, which is fine for the atlas (q<=31,k<=4; q<=16,k<=5; GF(2^m) m<=6).
This sidesteps every decoder/list-decoder bug -- the distances are ground truth.

Vectorization
-------------
The codeword matrix G has shape (q^k, n): row j is the evaluation of the j-th
polynomial on L.  We build it once (the Vandermonde-style generator applied to all
coefficient vectors) and reuse it.  Then:

    dist_to_code(w)            : (G != w).sum(axis=1).min()     -- one numpy pass
    counts on a line f0+gamma f1 : reuse the same G for every gamma

For larger q^k we chunk the codeword matrix to bound memory (see CodewordBook).

Encoding uses Horner's rule in the field's vectorized ops so it works uniformly
for prime and GF(2^m) fields.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np

from ff import FiniteField, _divisors


def _codeword_dtype(q: int):
    """Smallest unsigned int dtype that faithfully stores field elements in [0,q).

    q<=256 -> uint8 (8x less RAM than int64, the Wave-1 common case); q<=65536 ->
    uint16 (Wave-2 primes 257..1031, GF(2^7), GF(2^8)).  Comparisons (G==w/G!=w)
    against an int64 query word are exact for any unsigned width, so distance and
    agreement kernels are unaffected.
    """
    if q <= 256:
        return np.uint8
    if q <= 65536:
        return np.uint16
    return np.int64


# ---------------------------------------------------------------------------
# Domain constructors.
# ---------------------------------------------------------------------------
def domain_subgroup(F: FiniteField, n: int) -> np.ndarray:
    """Smooth domain: the order-n multiplicative subgroup H_n <= F* (n | q-1)."""
    return F.subgroup(n)


def domain_coset(F: FiniteField, n: int, shift: int) -> np.ndarray:
    """Smooth coset domain: shift * H_n."""
    return F.coset(n, shift)


def domain_full(F: FiniteField) -> np.ndarray:
    """Full multiplicative domain L = F* (n = q-1)."""
    return F.nonzero_elements()


def domain_random(F: FiniteField, n: int, rng: np.random.Generator,
                  include_zero: bool = False) -> np.ndarray:
    """Random size-n subset of F* (or F if include_zero)."""
    return F.random_subset(n, rng, include_zero=include_zero)


# ---------------------------------------------------------------------------
# Encoding: poly coefficients -> codeword on L.
# ---------------------------------------------------------------------------
def encode(F: FiniteField, L: np.ndarray, coeffs: np.ndarray) -> np.ndarray:
    """Evaluate polynomial with given coeffs (low->high degree) on all of L.

    coeffs has length k (deg < k).  Uses Horner in the field's vectorized ops.
    Returns a length-n codeword (field elements).
    """
    coeffs = np.asarray(coeffs, dtype=np.int64)
    L = np.asarray(L, dtype=np.int64)
    acc = np.full(len(L), int(coeffs[-1]), dtype=np.int64)  # highest coeff
    for c in coeffs[-2::-1]:                                  # down to constant
        acc = F.add_vec(F.mul_vec(acc, L), int(c))
    return acc


# ---------------------------------------------------------------------------
# Codeword book: all q^k codewords as a (possibly chunked) matrix.
# ---------------------------------------------------------------------------
def _all_coeff_vectors(q: int, k: int) -> np.ndarray:
    """All q^k coefficient vectors in F^k, shape (q^k, k), low->high degree.

    Generated with a mixed-radix odometer via numpy meshgrid-style construction.
    """
    # Build via repeated tiling: column j cycles with period q^j.
    total = q ** k
    out = np.empty((total, k), dtype=np.int64)
    for j in range(k):
        period = q ** j
        block = q ** (j + 1)
        pattern = np.repeat(np.arange(q, dtype=np.int64), period)
        out[:, j] = np.tile(pattern, total // block)
    return out


@dataclass
class CodewordBook:
    """All q^k codewords of RS[F,L,k] for fast exact distance computations.

    If q^k * n is small enough we store the full (q^k, n) matrix `G`.  Otherwise
    we keep coefficient vectors and regenerate codeword chunks on demand.
    """
    F: FiniteField
    L: np.ndarray
    k: int
    G: np.ndarray | None          # (q^k, n) full matrix, or None if chunked
    coeffs: np.ndarray            # (q^k, k) all coefficient vectors
    chunk: int                    # chunk size (rows) when G is None

    @property
    def n(self) -> int:
        return len(self.L)

    @property
    def num_codewords(self) -> int:
        return self.F.q ** self.k

    def iter_chunks(self):
        """Yield (start, codeword_matrix_chunk) covering all codewords."""
        if self.G is not None:
            yield 0, self.G
            return
        F, L, k = self.F, self.L, self.k
        Larr = np.asarray(L, dtype=np.int64)
        total = self.num_codewords
        for start in range(0, total, self.chunk):
            block = self.coeffs[start:start + self.chunk]   # (b, k)
            # Horner over the whole block at once: (b,) acc broadcast against L.
            acc = np.broadcast_to(block[:, -1:], (block.shape[0], len(Larr))).astype(np.int64)
            acc = acc.copy()
            for j in range(k - 2, -1, -1):
                acc = F.add_vec(F.mul_vec(acc, Larr[None, :]), block[:, j:j+1])
            yield start, acc.astype(_codeword_dtype(self.F.q))   # width matches full path


def build_codeword_book(F: FiniteField, L: np.ndarray, k: int,
                        max_full_cells: int = 60_000_000,
                        chunk: int = 50_000) -> CodewordBook:
    """Construct the codeword book.

    If q^k * n <= max_full_cells we materialize the full (q^k, n) matrix once
    (fast, the common case for the atlas).  Otherwise we store coefficient
    vectors and stream codeword chunks of `chunk` rows.
    """
    q = F.q
    n = len(L)
    total = q ** k
    coeffs = _all_coeff_vectors(q, k)
    Larr = np.asarray(L, dtype=np.int64)

    if total * n <= max_full_cells:
        # Materialize all codewords: Horner over the full coeff matrix.
        acc = np.broadcast_to(coeffs[:, -1:], (total, n)).astype(np.int64).copy()
        for j in range(k - 2, -1, -1):
            acc = F.add_vec(F.mul_vec(acc, Larr[None, :]), coeffs[:, j:j+1])
        # Store codewords in the smallest faithful unsigned dtype (uint8 for q<=256,
        # uint16 up to 65536) -> much less RAM than int64 at the same comparison
        # speed.  Distance/agreement kernels compare against an int64 query word,
        # which is exact for any unsigned width, so correctness is preserved.
        assert q <= 65536, "codeword storage supports q<=65536 (uint16)"
        return CodewordBook(F=F, L=Larr, k=k, G=acc.astype(_codeword_dtype(q)),
                            coeffs=coeffs, chunk=chunk)
    else:
        return CodewordBook(F=F, L=Larr, k=k, G=None, coeffs=coeffs, chunk=chunk)


# ---------------------------------------------------------------------------
# Exact distance and agreement computations.
# ---------------------------------------------------------------------------
def dist_to_code(book: CodewordBook, w: np.ndarray) -> int:
    """Exact Hamming distance from word w to the code (min over all codewords)."""
    w = np.asarray(w, dtype=np.int64)
    best = book.n + 1
    for _, G in book.iter_chunks():
        d = (G != w[None, :]).sum(axis=1)
        m = int(d.min())
        if m < best:
            best = m
            if best == 0:
                break
    return best


def nearest_codeword(book: CodewordBook, w: np.ndarray) -> tuple[int, np.ndarray]:
    """Return (distance, a nearest codeword) for word w."""
    w = np.asarray(w, dtype=np.int64)
    best = book.n + 1
    best_cw = None
    for start, G in book.iter_chunks():
        d = (G != w[None, :]).sum(axis=1)
        j = int(d.argmin())
        if int(d[j]) < best:
            best = int(d[j])
            best_cw = G[j].copy()
    return best, best_cw


def agreement_mask_best(book: CodewordBook, w: np.ndarray) -> np.ndarray:
    """Boolean mask (length n) of agreement positions for a *nearest* codeword.

    Ties broken arbitrarily (first argmin).  Used as a quick agreement-set proxy.
    """
    _, cw = nearest_codeword(book, w)
    return (cw == np.asarray(w, dtype=np.int64))


def is_delta_close(book: CodewordBook, w: np.ndarray, delta: float) -> bool:
    """True iff dist(w, C)/n <= delta (within floating tolerance)."""
    return dist_to_code(book, w) <= delta * book.n + 1e-9


# ---------------------------------------------------------------------------
# Min-distance / MDS check.
# ---------------------------------------------------------------------------
def min_distance(book: CodewordBook) -> int:
    """Exact minimum distance of the code = min weight of a nonzero codeword.

    For a LINEAR code, min distance = min Hamming weight over nonzero codewords.
    We compute weights of all codewords (distance from the all-zero word) and take
    the min over nonzero ones.  RS is MDS so this must equal n-k+1.
    """
    zero = np.zeros(book.n, dtype=np.int64)
    best = book.n + 1
    for _, G in book.iter_chunks():
        w = (G != zero[None, :]).sum(axis=1)   # Hamming weight of each codeword
        nz = w[w > 0]
        if nz.size:
            m = int(nz.min())
            if m < best:
                best = m
    return best


def assert_mds(book: CodewordBook):
    """Assert the code is MDS: min distance == n - k + 1."""
    d = min_distance(book)
    expected = book.n - book.k + 1
    assert d == expected, (
        f"MDS check FAILED for {book.F.name} n={book.n} k={book.k}: "
        f"min-dist={d}, expected n-k+1={expected}"
    )
    return d


# ---------------------------------------------------------------------------
# Random words / structured words helpers (used by the searches).
# ---------------------------------------------------------------------------
def random_word(F: FiniteField, n: int, rng: np.random.Generator) -> np.ndarray:
    """Uniform random word in F^n."""
    return rng.integers(0, F.q, size=n).astype(np.int64)


def random_codeword(book: CodewordBook, rng: np.random.Generator) -> np.ndarray:
    """A uniformly random codeword (random coefficient vector)."""
    c = rng.integers(0, book.F.q, size=book.k).astype(np.int64)
    return encode(book.F, book.L, c)


def codeword_plus_noise(book: CodewordBook, rng: np.random.Generator,
                        num_errors: int) -> np.ndarray:
    """A codeword with `num_errors` positions perturbed to random other values."""
    F = book.F
    w = random_codeword(book, rng).copy()
    pos = rng.choice(book.n, size=num_errors, replace=False)
    for p in pos:
        delta = int(rng.integers(1, F.q))   # nonzero offset guarantees a change
        w[p] = F.add(int(w[p]), delta)
    return w


# ===========================================================================
# Self-test battery.
# ===========================================================================
def _brute_min_distance_linear_check(F, L, k):
    """Independent brute check used only in tests for very small codes."""
    book = build_codeword_book(F, L, k)
    return min_distance(book)


def _self_test():
    from ff import PrimeField, BinaryExtensionField
    rng = np.random.default_rng(0xBEEF)
    print("RS code self-test")
    print("=" * 60)

    # ---- 1. Encoding correctness against direct Horner evaluation ----
    print("Encoding vs direct evaluation:")
    for F in [PrimeField(7), PrimeField(31), BinaryExtensionField(4)]:
        L = domain_full(F)
        k = 3
        coeffs = rng.integers(0, F.q, size=k).astype(np.int64)
        cw = encode(F, L, coeffs)
        # direct: p(x) = sum coeffs[i] x^i, computed scalar-wise via field ops
        for idx in range(0, len(L), max(1, len(L) // 7)):
            x = int(L[idx])
            val = 0
            for i, c in enumerate(coeffs):
                val = F.add(val, F.mul(int(c), F.pow(x, i)))
            assert cw[idx] == val, f"{F.name}: encode mismatch at x={x}"
        print(f"  {F.name}: OK")

    # ---- 2. q^k codewords, all distinct as coeff vectors ----
    print("Codeword book size:")
    for F, k in [(PrimeField(5), 3), (BinaryExtensionField(4), 3)]:
        L = domain_full(F)
        book = build_codeword_book(F, L, k)
        assert book.num_codewords == F.q ** k
        # codewords are distinct because RS is injective (n >= k, deg<k):
        rows = set()
        for _, G in book.iter_chunks():
            for r in G:
                rows.add(tuple(r.tolist()))
        assert len(rows) == F.q ** k, f"{F.name}: codewords not distinct"
        print(f"  {F.name} k={k}: {book.num_codewords} distinct codewords OK")

    # ---- 3. MDS min-distance check across domains/fields ----
    print("MDS min-distance (must equal n-k+1):")
    cases = [
        (PrimeField(7),  domain_full(PrimeField(7)),            3),   # n=6,k=3 -> d=4
        (PrimeField(11), domain_subgroup(PrimeField(11), 5),    2),   # n=5,k=2 -> d=4
        (PrimeField(11), domain_subgroup(PrimeField(11), 10),   3),   # n=10,k=3-> d=8
        (PrimeField(31), domain_subgroup(PrimeField(31), 6),    3),   # n=6,k=3 -> d=4
        (PrimeField(31), domain_subgroup(PrimeField(31), 10),   3),   # n=10
        (PrimeField(31), domain_subgroup(PrimeField(31), 15),   4),   # n=15
        (PrimeField(31), domain_coset(PrimeField(31), 10, 3),   3),   # coset n=10
        (BinaryExtensionField(4), domain_subgroup(BinaryExtensionField(4), 5), 2),   # n=5
        (BinaryExtensionField(4), domain_subgroup(BinaryExtensionField(4), 15), 4),  # n=15
        (BinaryExtensionField(4), domain_full(BinaryExtensionField(4)), 4),          # n=15
        (BinaryExtensionField(5), domain_subgroup(BinaryExtensionField(5), 31), 4),  # n=31
        (BinaryExtensionField(6), domain_subgroup(BinaryExtensionField(6), 9), 3),   # n=9
        (BinaryExtensionField(6), domain_subgroup(BinaryExtensionField(6), 21), 4),  # n=21
    ]
    for F, L, k in cases:
        book = build_codeword_book(F, L, k)
        d = assert_mds(book)
        rho = k / len(L)
        print(f"  {F.name:>10}  n={len(L):>2} k={k}  rho={rho:.3f}  "
              f"min-dist={d}  (n-k+1={len(L)-k+1})  OK")

    # ---- 4. Random-subset domain is also MDS (RS is MDS for ANY domain) ----
    print("Random-subset domain MDS:")
    F = PrimeField(31)
    L = domain_random(F, 12, rng)
    book = build_codeword_book(F, L, 3)
    d = assert_mds(book)
    print(f"  {F.name} random n=12 k=3: min-dist={d} (expected {12-3+1}) OK")

    # ---- 5. dist_to_code: a codeword has distance 0; codeword+e errors has dist e
    #         (for e < d/2 the nearest codeword is unique = the original) ----
    print("dist_to_code sanity:")
    F = PrimeField(31)
    L = domain_subgroup(F, 15)
    k = 4
    book = build_codeword_book(F, L, k)      # d = 15-4+1 = 12, unique dec radius 5
    cw = random_codeword(book, rng)
    assert dist_to_code(book, cw) == 0, "codeword must have distance 0"
    for e in [1, 2, 3, 4, 5]:
        w = cw.copy()
        pos = rng.choice(book.n, size=e, replace=False)
        for p in pos:
            w[p] = F.add(int(w[p]), int(rng.integers(1, F.q)))
        dd = dist_to_code(book, w)
        assert dd == e, f"within unique-dec radius, dist must equal #errors: {dd} vs {e}"
    print(f"  {F.name} n=15 k=4 (d=12): codeword dist=0; e in 1..5 errors -> dist=e OK")

    # ---- 6. chunked vs full codeword book agree on distance ----
    print("Chunked vs full agreement:")
    F = PrimeField(31)
    L = domain_subgroup(F, 15)
    k = 4
    full = build_codeword_book(F, L, k)                          # full matrix
    chunked = build_codeword_book(F, L, k, max_full_cells=1)     # force chunking
    assert chunked.G is None
    for _ in range(20):
        w = random_word(F, book.n, rng)
        assert dist_to_code(full, w) == dist_to_code(chunked, w), "chunk mismatch"
    assert min_distance(full) == min_distance(chunked)
    print("  chunked path matches full path on 20 random words + min-dist OK")

    print("=" * 60)
    print("ALL RS SELF-TESTS PASSED")


if __name__ == "__main__":
    _self_test()
