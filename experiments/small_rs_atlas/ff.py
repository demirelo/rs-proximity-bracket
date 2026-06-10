"""
ff.py -- Finite field arithmetic for the small-field Reed-Solomon proximity atlas.

We hand-roll finite fields because `galois`/sage are unavailable. Two families are
provided, behind a common interface:

  * Prime fields GF(p): elements are integers 0..p-1, arithmetic is mod p.
  * Binary extension fields GF(2^m): elements are integers 0..2^m-1 interpreted as
    polynomials over GF(2) in a fixed irreducible-polynomial (Conway-ish) basis.
    Multiplication uses precomputed log/antilog tables built from a primitive
    element, so element x is represented by its integer value and products are
    O(1) table lookups.

Design goals for the rest of the atlas:
  * Field elements are *plain Python ints in [0, q)*.  This is deliberate: it lets
    us pack whole codewords into numpy uint arrays and vectorize Hamming-distance
    computations.  Each field therefore exposes BOTH scalar ops (for clarity /
    testing) and *vectorized* numpy ops (the workhorse for the searches):
        add_vec(a, b), mul_vec(a, b), mul_scalar_vec(c, A), ...
    operating elementwise on numpy int arrays whose entries are field elements.

  * Every field knows:
        q                : field size
        char             : characteristic
        elements()       : iterator/array over all elements
        zero, one        : identities
        generator()      : a multiplicative generator (primitive element) of F*
        subgroup(n)      : the order-n multiplicative subgroup (n | q-1), as a
                           numpy array of its elements (powers of g^((q-1)/n))
        coset(n, shift)  : a multiplicative coset shift * H_n

The multiplicative group of any finite field is cyclic; the self-test verifies we
actually found a generator (its order is exactly q-1) and checks the field axioms
(associativity, distributivity, inverses) on random samples plus exhaustively on
the small fields.

Run `python ff.py` to execute the self-test battery.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Integer factorization helpers (stdlib-only; q-1 is tiny for our fields).
# ---------------------------------------------------------------------------
def _factorize(n: int) -> dict[int, int]:
    """Trial-division prime factorization of a small positive integer."""
    factors: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def _divisors(n: int) -> list[int]:
    """All positive divisors of n, sorted ascending."""
    divs = []
    d = 1
    while d * d <= n:
        if n % d == 0:
            divs.append(d)
            if d != n // d:
                divs.append(n // d)
        d += 1
    return sorted(divs)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


# ---------------------------------------------------------------------------
# Common base class.
# ---------------------------------------------------------------------------
class FiniteField:
    """Abstract base.  Subclasses must set q, char and implement the ops below."""

    q: int
    char: int
    name: str

    # ---- scalar ops (subclasses implement) -------------------------------
    def add(self, a: int, b: int) -> int:
        raise NotImplementedError

    def sub(self, a: int, b: int) -> int:
        raise NotImplementedError

    def mul(self, a: int, b: int) -> int:
        raise NotImplementedError

    def inv(self, a: int) -> int:
        raise NotImplementedError

    def neg(self, a: int) -> int:
        raise NotImplementedError

    def pow(self, a: int, e: int) -> int:
        """Exponentiation by squaring using the field's own mul.

        Handles negative exponents via inversion.  e==0 -> one.
        """
        if e < 0:
            a = self.inv(a)
            e = -e
        result = self.one
        base = a
        while e > 0:
            if e & 1:
                result = self.mul(result, base)
            base = self.mul(base, base)
            e >>= 1
        return result

    # ---- vectorized ops (subclasses implement) ---------------------------
    def add_vec(self, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def sub_vec(self, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def mul_vec(self, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def mul_scalar_vec(self, c: int, A: np.ndarray) -> np.ndarray:
        """Multiply every entry of array A by the scalar field element c."""
        raise NotImplementedError

    # ---- structure -------------------------------------------------------
    @property
    def zero(self) -> int:
        return 0

    @property
    def one(self) -> int:
        return 1

    def elements(self) -> np.ndarray:
        """All q field elements as a uint array [0, 1, ..., q-1]."""
        return np.arange(self.q, dtype=np.int64)

    def nonzero_elements(self) -> np.ndarray:
        return np.arange(1, self.q, dtype=np.int64)

    # ---- multiplicative-group machinery ---------------------------------
    def multiplicative_order(self, a: int) -> int:
        """Order of a in F* (a != 0).  Found via divisors of q-1."""
        if a == 0:
            raise ValueError("0 has no multiplicative order")
        order = self.q - 1
        for d in _divisors(order):
            if self.pow(a, d) == self.one:
                return d
        return order  # unreachable for a valid field element

    def is_generator(self, g: int) -> bool:
        """True iff g generates all of F* (i.e. has order q-1)."""
        if g == 0:
            return False
        order = self.q - 1
        # g is primitive iff g^(order/p) != 1 for every prime p | order.
        for p in _factorize(order):
            if self.pow(g, order // p) == self.one:
                return False
        return True

    _gen_cache: int | None = None

    def generator(self) -> int:
        """Return a primitive element (generator) of F*, searching upward."""
        if self._gen_cache is not None:
            return self._gen_cache
        for g in range(2, self.q):
            if self.is_generator(g):
                self._gen_cache = g
                return g
        # q == 2: F* = {1}, trivially generated by 1.
        if self.q == 2:
            self._gen_cache = 1
            return 1
        raise RuntimeError(f"no generator found for {self.name}")

    def subgroup(self, n: int) -> np.ndarray:
        """The unique order-n multiplicative subgroup H_n <= F*, as an array.

        Requires n | (q-1).  Returned in the natural cyclic order
        [1, h, h^2, ..., h^{n-1}] where h = g^((q-1)/n) has order n.
        """
        if (self.q - 1) % n != 0:
            raise ValueError(f"order {n} does not divide |F*|={self.q - 1}")
        g = self.generator()
        h = self.pow(g, (self.q - 1) // n)
        elems = np.empty(n, dtype=np.int64)
        cur = self.one
        for i in range(n):
            elems[i] = cur
            cur = self.mul(cur, h)
        return elems

    def coset(self, n: int, shift: int) -> np.ndarray:
        """The multiplicative coset shift * H_n (shift != 0).

        Same cyclic ordering as subgroup(n), scaled by `shift`.
        """
        if shift == 0:
            raise ValueError("coset shift must be nonzero")
        H = self.subgroup(n)
        return self.mul_scalar_vec(shift, H)

    def random_subset(self, n: int, rng: np.random.Generator,
                      include_zero: bool = False) -> np.ndarray:
        """A random size-n subset of field elements (distinct).

        By default drawn from F* (nonzero), matching RS evaluation-domain
        convention; set include_zero=True to draw from all of F.
        """
        pool = self.elements() if include_zero else self.nonzero_elements()
        if n > len(pool):
            raise ValueError(f"cannot pick {n} distinct from pool of {len(pool)}")
        idx = rng.choice(len(pool), size=n, replace=False)
        return np.sort(pool[idx])

    def __repr__(self) -> str:
        return f"<{self.name}>"


# ---------------------------------------------------------------------------
# Prime fields GF(p).
# ---------------------------------------------------------------------------
class PrimeField(FiniteField):
    """GF(p) for prime p.  Elements are ints mod p."""

    def __init__(self, p: int):
        if not is_prime(p):
            raise ValueError(f"{p} is not prime")
        self.q = p
        self.char = p
        self.name = f"GF({p})"
        # Precompute multiplicative inverses for all nonzero elements.
        self._inv = np.zeros(p, dtype=np.int64)
        for a in range(1, p):
            self._inv[a] = pow(a, p - 2, p)  # Fermat

    # scalar ops
    def add(self, a, b):
        return (a + b) % self.q

    def sub(self, a, b):
        return (a - b) % self.q

    def mul(self, a, b):
        return (a * b) % self.q

    def neg(self, a):
        return (-a) % self.q

    def inv(self, a):
        if a == 0:
            raise ZeroDivisionError("inverse of 0")
        return int(self._inv[a])

    # vectorized ops -- arithmetic on int64 arrays then reduce mod p.
    def add_vec(self, A, B):
        return (A + B) % self.q

    def sub_vec(self, A, B):
        return (A - B) % self.q

    def mul_vec(self, A, B):
        return (A * B) % self.q

    def mul_scalar_vec(self, c, A):
        return (int(c) * A) % self.q


# ---------------------------------------------------------------------------
# Binary extension fields GF(2^m) via log/antilog tables.
# ---------------------------------------------------------------------------
# Fixed primitive (irreducible) polynomials, given as the low bits of x^m + ...
# (the x^m term is implicit).  These are standard primitive polynomials so that
# x (i.e. the integer 2) is itself a primitive element.
_GF2M_POLY = {
    2: 0b111,        # x^2 + x + 1
    3: 0b1011,       # x^3 + x + 1
    4: 0b10011,      # x^4 + x + 1
    5: 0b100101,     # x^5 + x^2 + 1
    6: 0b1000011,    # x^6 + x + 1
    7: 0b10001001,   # x^7 + x^3 + 1
    8: 0b100011101,  # x^8 + x^4 + x^3 + x^2 + 1
}


class BinaryExtensionField(FiniteField):
    """GF(2^m), polynomial basis, log/antilog tables.

    Element i in [0, 2^m) encodes the polynomial sum_j bit_j(i) x^j.  Addition is
    XOR.  Multiplication uses logs base the primitive element alpha = x (=2):
    every nonzero element equals alpha^k for a unique k in [0, 2^m-1).
    """

    def __init__(self, m: int):
        if m not in _GF2M_POLY:
            raise ValueError(f"GF(2^{m}) not configured (have {sorted(_GF2M_POLY)})")
        self.m = m
        self.q = 1 << m
        self.char = 2
        self.name = f"GF(2^{m})"
        self._poly = _GF2M_POLY[m]

        # Build antilog (exp) and log tables using alpha = x = 2 as primitive root.
        q = self.q
        exp = np.zeros(q - 1, dtype=np.int64)   # exp[k] = alpha^k, k in [0,q-1)
        log = np.full(q, -1, dtype=np.int64)    # log[v] = k s.t. alpha^k = v
        x = 1
        for k in range(q - 1):
            exp[k] = x
            log[x] = k
            # multiply x by alpha (=2): shift left, reduce by poly if degree m.
            x <<= 1
            if x & q:                # bit m set -> reduce
                x ^= self._poly
        # Sanity: alpha must be primitive, i.e. we cycled through all q-1 nonzero.
        assert x == 1, f"alpha=2 is not primitive for GF(2^{m}); fix the polynomial"

        # Double the exp table so exp[k] is valid for k in [0, 2(q-1)) without mod.
        self._exp = np.concatenate([exp, exp])
        self._log = log

    # ---- scalar ops ------------------------------------------------------
    def add(self, a, b):
        return a ^ b

    def sub(self, a, b):
        return a ^ b           # char 2

    def neg(self, a):
        return a               # char 2: -a == a

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
        return np.bitwise_xor(A, B)

    def sub_vec(self, A, B):
        return np.bitwise_xor(A, B)

    def mul_vec(self, A, B):
        """Elementwise product of two field-element arrays via log/antilog.

        Zeros in either operand must map to zero output; we handle them with a
        mask since log[0] is undefined.
        """
        A = np.asarray(A, dtype=np.int64)
        B = np.asarray(B, dtype=np.int64)
        nz = (A != 0) & (B != 0)
        out = np.zeros(np.broadcast(A, B).shape, dtype=np.int64)
        la = self._log[A]
        lb = self._log[B]
        # exp table is doubled, so la+lb in [0, 2(q-1)-2] indexes safely.
        prod = self._exp[(la + lb)]
        out[nz] = prod[nz]
        return out

    def mul_scalar_vec(self, c, A):
        c = int(c)
        if c == 0:
            return np.zeros_like(A)
        A = np.asarray(A, dtype=np.int64)
        nz = A != 0
        out = np.zeros_like(A)
        lc = int(self._log[c])
        out[nz] = self._exp[self._log[A[nz]] + lc]
        return out


# ---------------------------------------------------------------------------
# Convenience constructor.
# ---------------------------------------------------------------------------
def make_field(spec: str) -> FiniteField:
    """Build a field from a string spec.

    Examples: "GF(31)", "GF(2^4)", "31", "2^5".
    """
    s = spec.strip().upper().replace("GF(", "").replace(")", "")
    if "^" in s:
        base, exp = s.split("^")
        if int(base) != 2:
            raise ValueError("only GF(2^m) extension fields are implemented")
        return BinaryExtensionField(int(exp))
    return PrimeField(int(s))


# ===========================================================================
# Self-test battery.
# ===========================================================================
def _test_field(F: FiniteField, rng: np.random.Generator, n_samples: int = 2000):
    """Check field axioms, generator, subgroup structure for one field."""
    q = F.q
    print(f"  testing {F.name} (q={q}) ...", end=" ", flush=True)

    # --- exhaustive on small fields, else random sampling ---
    if q <= 64:
        els = list(range(q))
        triples = [(a, b, c) for a in els for b in els for c in els]
        # cap exhaustive triples for q near 64 (64^3 = 262144, fine)
    else:
        triples = [(int(rng.integers(q)), int(rng.integers(q)), int(rng.integers(q)))
                   for _ in range(n_samples)]

    for a, b, c in triples:
        # commutativity
        assert F.add(a, b) == F.add(b, a), f"add comm {F.name} {a},{b}"
        assert F.mul(a, b) == F.mul(b, a), f"mul comm {F.name} {a},{b}"
        # associativity
        assert F.add(F.add(a, b), c) == F.add(a, F.add(b, c)), f"add assoc {a},{b},{c}"
        assert F.mul(F.mul(a, b), c) == F.mul(a, F.mul(b, c)), f"mul assoc {a},{b},{c}"
        # distributivity
        assert F.mul(a, F.add(b, c)) == F.add(F.mul(a, b), F.mul(a, c)), \
            f"distrib {F.name} {a},{b},{c}"
        # additive inverse / identity
        assert F.add(a, F.zero) == a
        assert F.add(a, F.neg(a)) == F.zero
        assert F.sub(a, b) == F.add(a, F.neg(b))

    # multiplicative inverses for all nonzero
    for a in range(1, q):
        assert F.mul(a, F.inv(a)) == F.one, f"inv {F.name} {a}"
        assert F.mul(a, F.one) == a

    # --- generator: order exactly q-1, and it enumerates F* ---
    g = F.generator()
    assert F.is_generator(g), f"{F.name}: {g} not a generator"
    assert F.multiplicative_order(g) == q - 1, f"{F.name}: gen order wrong"
    seen = set()
    cur = F.one
    for _ in range(q - 1):
        seen.add(cur)
        cur = F.mul(cur, g)
    assert len(seen) == q - 1 and 0 not in seen, f"{F.name}: F* not cyclic via g"

    # --- subgroups for every divisor n of q-1: closure, order, MDS-relevant ---
    for n in _divisors(q - 1):
        H = F.subgroup(n)
        assert len(H) == n, f"{F.name}: |H_{n}|={len(H)}"
        assert len(set(H.tolist())) == n, f"{F.name}: H_{n} has duplicates"
        assert 0 not in set(H.tolist()), f"{F.name}: H_{n} contains 0"
        assert F.one in set(H.tolist()), f"{F.name}: H_{n} missing identity"
        # closure under multiplication (spot check: product of first two in H)
        if n >= 2:
            prod = F.mul(int(H[1]), int(H[1]))
            assert prod in set(H.tolist()), f"{F.name}: H_{n} not closed"
        # every element has order dividing n
        for x in H.tolist():
            if x != F.one:
                assert F.pow(int(x), n) == F.one, f"{F.name}: H_{n} elt order"

    # --- vectorized ops agree with scalar ops on random arrays ---
    A = rng.integers(0, q, size=500).astype(np.int64)
    B = rng.integers(0, q, size=500).astype(np.int64)
    add_v = F.add_vec(A, B)
    sub_v = F.sub_vec(A, B)
    mul_v = F.mul_vec(A, B)
    for i in range(0, 500, 17):
        a, b = int(A[i]), int(B[i])
        assert add_v[i] == F.add(a, b), f"add_vec {F.name} {a},{b}"
        assert sub_v[i] == F.sub(a, b), f"sub_vec {F.name} {a},{b}"
        assert mul_v[i] == F.mul(a, b), f"mul_vec {F.name} {a},{b}: {mul_v[i]} vs {F.mul(a,b)}"
    # scalar-times-vector
    c = int(rng.integers(0, q))
    cv = F.mul_scalar_vec(c, A)
    for i in range(0, 500, 17):
        assert cv[i] == F.mul(c, int(A[i])), f"mul_scalar_vec {F.name}"

    # --- coset sanity: shift * H_n is disjoint from H_n unless shift in H_n ---
    if q - 1 >= 4:
        n = [d for d in _divisors(q - 1) if 1 < d < q - 1]
        if n:
            nn = n[0]
            H = set(F.subgroup(nn).tolist())
            # pick a shift outside H
            shift = next((s for s in range(1, q) if s not in H), None)
            if shift is not None:
                cs = set(F.coset(nn, shift).tolist())
                assert len(cs) == nn, f"{F.name}: coset size"
                assert H.isdisjoint(cs), f"{F.name}: nontrivial coset not disjoint"

    print("OK")


def _self_test():
    rng = np.random.default_rng(0xC0DE)
    print("Finite field self-test")
    print("=" * 60)

    print("Prime fields:")
    for p in [2, 3, 5, 7, 11, 13, 17, 31, 127, 257]:
        _test_field(PrimeField(p), rng)

    print("Binary extension fields:")
    for m in [2, 3, 4, 5, 6, 7, 8]:
        _test_field(BinaryExtensionField(m), rng)

    # Cross-check: GF(2^m) distributivity already covered; verify a known product.
    F16 = BinaryExtensionField(4)
    # In GF(2^4) with x^4+x+1: alpha=2, alpha^4 = alpha+1 = 3.
    assert F16.pow(2, 4) == 3, f"GF(2^4) alpha^4 should be 3, got {F16.pow(2,4)}"
    # alpha^15 = 1
    assert F16.pow(2, 15) == 1
    print("\nKnown-value spot checks: OK (GF(2^4) alpha^4=3, alpha^15=1)")

    # make_field parser
    assert make_field("GF(31)").q == 31
    assert make_field("GF(2^5)").q == 32
    assert make_field("2^6").q == 64
    print("make_field parser: OK")

    print("=" * 60)
    print("ALL FIELD SELF-TESTS PASSED")


if __name__ == "__main__":
    _self_test()
