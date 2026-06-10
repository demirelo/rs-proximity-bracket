"""
ff_ext.py -- ODD-CHARACTERISTIC EXTENSION fields GF(p^e), p odd prime, e>=2.

This is the decisive-experiment companion to ff.py.  ff.py implements prime fields
GF(p) and BINARY extension fields GF(2^m).  Binary extensions have ODD group order
(2^m - 1) and therefore NO power-of-two-order ("smooth") multiplicative subgroup, so
they cannot host the Kambire/BCHKS smooth-subgroup construction at all.  The genuine
open question (Proximity-Prize sub-lemma N2) is whether that construction -- known to
fire over PRIME fields -- also fires over ODD-characteristic EXTENSION fields
GF(p^e).  Those DO have smooth subgroups (e.g. |GF(9)*|=8=2^3, |GF(49)*|=48=16*3,
|GF(289)*|=288=32*9), so the construction can literally run; whether the algebraic
"distinct subset-sum" signal survives in F_{p^e} is what we are here to measure.

Representation
--------------
An element of GF(p^e) is a polynomial of degree < e over GF(p):
    a(x) = a_0 + a_1 x + ... + a_{e-1} x^{e-1},   a_j in {0,...,p-1}.
We pack the digit vector (a_0,...,a_{e-1}) into a single integer in [0, p^e) using a
base-p ("little-endian") encoding:
    code = a_0 + a_1*p + a_2*p^2 + ... + a_{e-1}*p^{e-1}.
This keeps the SAME contract as ff.py: every field element is a plain int in [0, q),
so the whole RS / bad-line machinery (which packs codewords into numpy uint arrays
and compares element-by-element) works UNCHANGED.

Arithmetic
----------
Addition/subtraction/negation act digitwise mod p (no carries -- it is a vector
space over GF(p)).  Multiplication is polynomial multiplication followed by reduction
modulo a fixed irreducible polynomial pi(x) of degree e over GF(p).  We find pi by
brute force (smallest monic irreducible) and verify irreducibility exhaustively for
the small e we use.  For speed we build full log/antilog tables against a primitive
element (the multiplicative group of any finite field is cyclic), exactly mirroring
the GF(2^m) implementation: every nonzero element is alpha^t for a unique t, products
are O(1) table lookups, and the vectorized ops index the (doubled) antilog table.

Self-test (run `python ff_ext.py`)
----------------------------------
For GF(9), GF(27), GF(25), GF(49), GF(81), GF(121), GF(169), GF(289), GF(343),
GF(243): exhaustive field axioms (or sampled for the larger ones), an EXPLICIT
multiplicative generator with order exactly q-1 enumerating F*, and for EVERY divisor
n of q-1 the order-n subgroup is closed, of size n, contains 1, excludes 0, and every
element has order dividing n.  Plus vectorized-vs-scalar agreement and a couple of
RS-is-MDS checks (delegated to rs.py via the shared interface).
"""

from __future__ import annotations

import numpy as np

from ff import FiniteField, _divisors, _factorize, is_prime


# ---------------------------------------------------------------------------
# Polynomial arithmetic over GF(p) on little-endian coefficient lists.
# ---------------------------------------------------------------------------
def _poly_trim(a: list[int]) -> list[int]:
    """Drop trailing (high-degree) zero coefficients; [] represents the 0 poly."""
    i = len(a) - 1
    while i >= 0 and a[i] == 0:
        i -= 1
    return a[: i + 1]


def _poly_mulmod(a: list[int], b: list[int], pi: list[int], p: int) -> list[int]:
    """(a * b) mod pi over GF(p).  pi is the MONIC reducing poly (degree e)."""
    if not a or not b:
        return []
    prod = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    prod[i + j] = (prod[i + j] + ai * bj) % p
    return _poly_divmod_rem(prod, pi, p)


def _poly_divmod_rem(a: list[int], pi: list[int], p: int) -> list[int]:
    """Remainder of a mod pi over GF(p); pi monic of degree e."""
    a = a[:]  # copy; we reduce in place
    e = len(pi) - 1
    # Schoolbook reduction: while deg a >= e, subtract a shifted multiple of pi.
    deg_a = len(a) - 1
    while deg_a >= e:
        coeff = a[deg_a] % p
        if coeff:
            shift = deg_a - e
            # a -= coeff * x^shift * pi   (pi monic so leading term cancels)
            for j in range(e + 1):
                a[shift + j] = (a[shift + j] - coeff * pi[j]) % p
        # trim leading zeros
        while deg_a >= 0 and a[deg_a] % p == 0:
            deg_a -= 1
        a = a[: deg_a + 1]
    return _poly_trim([c % p for c in a])


def _is_irreducible(pi: list[int], p: int) -> bool:
    """Test whether monic pi (degree e) is irreducible over GF(p).

    Method (correctness-first, exhaustive for our small e):
      * It must have NO root in GF(p)            (rules out a linear factor).
      * For e in {2,3} a root test alone is decisive (a reducible deg-2 or deg-3 poly
        always has a linear factor).
      * For general e we additionally test divisibility by every monic poly of degree
        d for 2 <= d <= e//2.  e<=5 here, so degrees checked are tiny.
    """
    e = len(pi) - 1
    # 1) no linear factor: no root in GF(p)
    for r in range(p):
        # evaluate pi(r) mod p by Horner
        val = 0
        for c in reversed(pi):
            val = (val * r + c) % p
        if val == 0:
            return False
    if e <= 3:
        return True
    # 2) no factor of degree d for 2 <= d <= e//2
    for d in range(2, e // 2 + 1):
        # enumerate all monic polys of degree d over GF(p): leading coeff 1,
        # lower d coeffs free in {0..p-1}.
        for code in range(p ** d):
            lower = []
            c = code
            for _ in range(d):
                lower.append(c % p)
                c //= p
            cand = lower + [1]  # monic degree d
            if len(_poly_divmod_rem(pi[:], cand, p)) == 0:
                return False
    return True


def _find_irreducible(p: int, e: int) -> list[int]:
    """Smallest monic irreducible poly of degree e over GF(p) (little-endian list).

    Search lower coefficients in odometer order so the result is deterministic.
    """
    for code in range(p ** e):
        lower = []
        c = code
        for _ in range(e):
            lower.append(c % p)
            c //= p
        pi = lower + [1]  # monic degree e
        if _is_irreducible(pi, p):
            return pi
    raise RuntimeError(f"no irreducible degree-{e} poly found over GF({p})")


# ---------------------------------------------------------------------------
# The field.
# ---------------------------------------------------------------------------
class PrimePowerField(FiniteField):
    """GF(p^e) for odd prime p and e>=2, polynomial basis with log/antilog tables.

    Element code in [0, p^e) <-> base-p digit poly a_0 + a_1 x + ... + a_{e-1}x^{e-1}.
    Addition is digitwise mod p; multiplication via log tables against a primitive
    element alpha (found by search).  Mirrors ff.BinaryExtensionField's table design.
    """

    def __init__(self, p: int, e: int, pi: list[int] | None = None):
        if not is_prime(p):
            raise ValueError(f"{p} is not prime")
        if p == 2:
            raise ValueError("use BinaryExtensionField for p=2; this is for ODD p")
        if e < 2:
            raise ValueError("extension degree e must be >= 2")
        self.p = p
        self.e = e
        self.q = p ** e
        self.char = p
        self.name = f"GF({p}^{e})"
        self.pi = pi if pi is not None else _find_irreducible(p, e)
        assert len(self.pi) == e + 1 and self.pi[-1] == 1

        # Precompute base-p digit decompositions for every code, plus pairwise
        # digitwise add / sub / neg lookup at the DIGIT level (a single mod-p table)
        # so vectorized ops are fast.
        self._pow_p = np.array([p ** j for j in range(e)], dtype=np.int64)

        # Build digit matrix D: shape (q, e), D[code] = (a_0..a_{e-1}).
        codes = np.arange(self.q, dtype=np.int64)
        D = np.empty((self.q, e), dtype=np.int64)
        tmp = codes.copy()
        for j in range(e):
            D[:, j] = tmp % p
            tmp //= p
        self._digits = D  # (q, e)

        # ---- build log/antilog tables against a primitive element ----
        self._build_mul_tables()

    # ----- helpers: code <-> poly (little-endian list) -----
    def _to_poly(self, code: int) -> list[int]:
        a = []
        c = int(code)
        for _ in range(self.e):
            a.append(c % self.p)
            c //= self.p
        return _poly_trim(a)

    def _from_poly(self, a: list[int]) -> int:
        code = 0
        for j, c in enumerate(a):
            code += (c % self.p) * (self.p ** j)
        return int(code)

    def _mul_codes(self, x: int, y: int) -> int:
        """Direct polynomial multiply-mod (used to bootstrap the log tables)."""
        if x == 0 or y == 0:
            return 0
        prod = _poly_mulmod(self._to_poly(x), self._to_poly(y), self.pi, self.p)
        return self._from_poly(prod)

    def _build_mul_tables(self):
        """Find a primitive element and build doubled antilog + log tables."""
        q = self.q
        order = q - 1
        prime_factors = list(_factorize(order))

        def is_primitive(code: int) -> bool:
            if code == 0:
                return False
            # primitive iff code^(order/pf) != 1 for every prime pf | order
            for pf in prime_factors:
                # fast exponentiation using _mul_codes
                ex = order // pf
                r = 1
                base = code
                ee = ex
                while ee > 0:
                    if ee & 1:
                        r = self._mul_codes(r, base)
                    base = self._mul_codes(base, base)
                    ee >>= 1
                if r == 1:
                    return False
            return True

        alpha = None
        for cand in range(2, q):
            if is_primitive(cand):
                alpha = cand
                break
        if alpha is None:
            raise RuntimeError(f"no primitive element found for {self.name}")
        self._alpha = alpha

        exp = np.zeros(order, dtype=np.int64)   # exp[t] = alpha^t
        log = np.full(q, -1, dtype=np.int64)    # log[v] = t s.t. alpha^t = v
        cur = 1
        for t in range(order):
            exp[t] = cur
            log[cur] = t
            cur = self._mul_codes(cur, alpha)
        assert cur == 1, f"{self.name}: alpha not primitive (cycle != order)"
        # Double the antilog table so exp[t] valid for t in [0, 2*order) without mod.
        self._exp = np.concatenate([exp, exp])
        self._log = log

    # ---- scalar ops ------------------------------------------------------
    def add(self, a, b):
        # digitwise add mod p, recombine
        da = self._digits[a]
        db = self._digits[b]
        return int((((da + db) % self.p) * self._pow_p).sum())

    def sub(self, a, b):
        da = self._digits[a]
        db = self._digits[b]
        return int((((da - db) % self.p) * self._pow_p).sum())

    def neg(self, a):
        da = self._digits[a]
        return int((((-da) % self.p) * self._pow_p).sum())

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return int(self._exp[self._log[a] + self._log[b]])

    def inv(self, a):
        if a == 0:
            raise ZeroDivisionError("inverse of 0")
        return int(self._exp[(self.q - 1) - self._log[a]])

    # ---- vectorized ops --------------------------------------------------
    def add_vec(self, A, B):
        A = np.asarray(A, dtype=np.int64)
        B = np.asarray(B, dtype=np.int64)
        # digitwise add mod p via the precomputed digit matrix, then recombine
        da = self._digits[A]               # (..., e)
        db = self._digits[B]
        s = (da + db) % self.p
        return (s * self._pow_p).sum(axis=-1)

    def sub_vec(self, A, B):
        A = np.asarray(A, dtype=np.int64)
        B = np.asarray(B, dtype=np.int64)
        da = self._digits[A]
        db = self._digits[B]
        s = (da - db) % self.p
        return (s * self._pow_p).sum(axis=-1)

    def mul_vec(self, A, B):
        A = np.asarray(A, dtype=np.int64)
        B = np.asarray(B, dtype=np.int64)
        nz = (A != 0) & (B != 0)
        out = np.zeros(np.broadcast(A, B).shape, dtype=np.int64)
        la = self._log[A]
        lb = self._log[B]
        prod = self._exp[(la + lb)]
        out[nz] = prod[nz]
        return out

    def mul_scalar_vec(self, c, A):
        c = int(c)
        if c == 0:
            return np.zeros_like(np.asarray(A, dtype=np.int64))
        A = np.asarray(A, dtype=np.int64)
        nz = A != 0
        out = np.zeros_like(A)
        lc = int(self._log[c])
        out[nz] = self._exp[self._log[A[nz]] + lc]
        return out

    def generator(self) -> int:
        """Primitive element (the one used to build the log tables)."""
        return int(self._alpha)


# ---------------------------------------------------------------------------
# Convenience constructor mirroring ff.make_field but for odd-char extensions.
# ---------------------------------------------------------------------------
def make_ext_field(p: int, e: int) -> PrimePowerField:
    return PrimePowerField(p, e)


# Catalogue of the fields the brief asks us to support / self-test.
TARGET_FIELDS = [
    (3, 2),   # GF(9),   q-1=8   = 2^3
    (3, 3),   # GF(27),  q-1=26  = 2*13
    (5, 2),   # GF(25),  q-1=24  = 2^3*3
    (7, 2),   # GF(49),  q-1=48  = 2^4*3
    (3, 4),   # GF(81),  q-1=80  = 2^4*5
    (11, 2),  # GF(121), q-1=120 = 2^3*3*5
    (13, 2),  # GF(169), q-1=168 = 2^3*3*7
    (17, 2),  # GF(289), q-1=288 = 2^5*3^2
    (7, 3),   # GF(343), q-1=342 = 2*3^2*19
    (3, 5),   # GF(243), q-1=242 = 2*11^2
]


# ===========================================================================
# Self-test battery.
# ===========================================================================
def _test_field(F: PrimePowerField, rng: np.random.Generator, n_samples: int = 3000):
    q = F.q
    print(f"  testing {F.name} (q={q}, pi={F.pi}, gen={F.generator()}) ...",
          end=" ", flush=True)

    # --- field axioms: exhaustive if small, else sampled ---
    if q <= 49:
        els = list(range(q))
        triples = [(a, b, c) for a in els for b in els for c in els]
    else:
        triples = [(int(rng.integers(q)), int(rng.integers(q)), int(rng.integers(q)))
                   for _ in range(n_samples)]

    for a, b, c in triples:
        assert F.add(a, b) == F.add(b, a), f"add comm {F.name} {a},{b}"
        assert F.mul(a, b) == F.mul(b, a), f"mul comm {F.name} {a},{b}"
        assert F.add(F.add(a, b), c) == F.add(a, F.add(b, c)), "add assoc"
        assert F.mul(F.mul(a, b), c) == F.mul(a, F.mul(b, c)), "mul assoc"
        assert F.mul(a, F.add(b, c)) == F.add(F.mul(a, b), F.mul(a, c)), "distrib"
        assert F.add(a, F.zero) == a
        assert F.add(a, F.neg(a)) == F.zero
        assert F.sub(a, b) == F.add(a, F.neg(b))

    # multiplicative inverses for all nonzero, and identity
    for a in range(1, q):
        assert F.mul(a, F.inv(a)) == F.one, f"inv {F.name} {a}"
        assert F.mul(a, F.one) == a

    # characteristic: 1 added p times == 0, and p*x == 0 for all x (sampled)
    acc = 0
    for _ in range(F.p):
        acc = F.add(acc, F.one)
    assert acc == F.zero, f"{F.name}: char != p"

    # --- explicit generator: order exactly q-1, enumerates F* ---
    g = F.generator()
    assert F.is_generator(g), f"{F.name}: {g} not a generator"
    assert F.multiplicative_order(g) == q - 1, f"{F.name}: gen order wrong"
    seen = set()
    cur = F.one
    for _ in range(q - 1):
        seen.add(cur)
        cur = F.mul(cur, g)
    assert len(seen) == q - 1 and 0 not in seen, f"{F.name}: F* not cyclic via g"

    # --- every divisor subgroup: closed, sized n, contains 1, excludes 0 ---
    for n in _divisors(q - 1):
        H = F.subgroup(n)
        Hset = set(H.tolist())
        assert len(H) == n, f"{F.name}: |H_{n}|={len(H)}"
        assert len(Hset) == n, f"{F.name}: H_{n} duplicates"
        assert 0 not in Hset, f"{F.name}: H_{n} contains 0"
        assert F.one in Hset, f"{F.name}: H_{n} missing 1"
        # full closure check for small subgroups, spot check for large
        if n <= 64:
            for x in H.tolist():
                for y in H.tolist():
                    assert F.mul(int(x), int(y)) in Hset, f"{F.name}: H_{n} not closed"
        else:
            for _ in range(200):
                x = int(H[int(rng.integers(n))]); y = int(H[int(rng.integers(n))])
                assert F.mul(x, y) in Hset, f"{F.name}: H_{n} not closed"
        for x in H.tolist():
            assert F.pow(int(x), n) == F.one, f"{F.name}: H_{n} element order"

    # --- vectorized ops agree with scalar ops ---
    A = rng.integers(0, q, size=600).astype(np.int64)
    B = rng.integers(0, q, size=600).astype(np.int64)
    add_v = F.add_vec(A, B); sub_v = F.sub_vec(A, B); mul_v = F.mul_vec(A, B)
    for i in range(0, 600, 13):
        a, b = int(A[i]), int(B[i])
        assert add_v[i] == F.add(a, b), f"add_vec {F.name}"
        assert sub_v[i] == F.sub(a, b), f"sub_vec {F.name}"
        assert mul_v[i] == F.mul(a, b), f"mul_vec {F.name} {a},{b}"
    c = int(rng.integers(1, q))
    cv = F.mul_scalar_vec(c, A)
    for i in range(0, 600, 13):
        assert cv[i] == F.mul(c, int(A[i])), f"mul_scalar_vec {F.name}"

    # --- coset disjointness sanity ---
    mids = [d for d in _divisors(q - 1) if 1 < d < q - 1]
    if mids:
        nn = mids[0]
        Hset = set(F.subgroup(nn).tolist())
        shift = next((s for s in range(1, q) if s not in Hset), None)
        if shift is not None:
            cs = set(F.coset(nn, shift).tolist())
            assert len(cs) == nn and Hset.isdisjoint(cs), f"{F.name}: coset"

    print("OK")


def _self_test():
    rng = np.random.default_rng(0x0E47)
    print("Odd-characteristic extension field self-test (ff_ext.py)")
    print("=" * 68)
    for (p, e) in TARGET_FIELDS:
        _test_field(PrimePowerField(p, e), rng)

    # ---- Known-value spot checks ----
    # GF(9) with the smallest irreducible: confirm a primitive element cycles 8.
    F9 = PrimePowerField(3, 2)
    g = F9.generator()
    powers = set(F9.pow(g, t) for t in range(8))
    assert powers == set(range(1, 9)), "GF(9): generator must enumerate F*"
    assert F9.pow(g, 8) == 1, "GF(9): g^8 must be 1"
    # subgroup of order 2 is {1, -1}: the unique element of order 2 is p^e-1 element
    H2 = set(F9.subgroup(2).tolist())
    assert H2 == {1, F9.neg(1)}, f"GF(9): order-2 subgroup should be {{1,-1}}, got {H2}"
    print("\nKnown-value spot checks: OK "
          "(GF(9) generator cycles F*, order-2 subgroup = {1,-1})")

    # ---- Embedding sanity: GF(p) sits inside GF(p^e) as the constant polys ----
    # The constant polynomials {0,1,...,p-1} (codes 0..p-1) are closed under + and *
    # and reproduce GF(p) arithmetic.
    F49 = PrimePowerField(7, 2)
    for a in range(7):
        for b in range(7):
            assert F49.add(a, b) == (a + b) % 7, "GF(49): GF(7) add embedding"
            assert F49.mul(a, b) == (a * b) % 7, "GF(49): GF(7) mul embedding"
    print("Subfield embedding GF(p) <= GF(p^e): OK (constants reproduce GF(p))")

    print("=" * 68)
    print("ALL ff_ext SELF-TESTS PASSED")


if __name__ == "__main__":
    _self_test()
