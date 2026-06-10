"""
n2_char0_count.py -- the FIELD-AGNOSTIC generic distinct-bad-scalar count, computed in
CHARACTERISTIC ZERO (the cyclotomic ring Z[xi_s]), plus the smart DP there.

THE POINT (why this settles 'extension vs prime' for general r)
---------------------------------------------------------------
The firing condition e_2=0 and the bad scalar lambda=e_1 are computed PURELY from the
order-s subgroup H = <xi>, xi a primitive s-th root of unity.  Over ANY field of
characteristic p in which xi exists, e_1 and e_2 are the images of FIXED elements of the
cyclotomic ring Z[xi_s] under the reduction Z[xi_s] -> F.  Two firing subsets collide
(same lambda) in F iff p divides the (integer) resultant-type quantity measuring their
difference in Z[xi_s].  Hence:

  * The GENERIC count (large characteristic, avoiding the finitely many 'bad' primes that
    cause extra collisions) equals the CHARACTERISTIC-ZERO count N0(s,r) := #distinct
    e_1 over firing r-subsets of the abstract s-th roots of unity in C.
  * N0(s,r) does NOT depend on the field, the extension degree e, or even whether xi lands
    in F_p or F_{p^e}.  A primitive s-th root exists in F_p iff p==1 (mod s), and GENUINELY
    in F_{p^2}\\F_p iff p==-1 (mod s); in BOTH cases the reduction Z[xi_s]->F is injective on
    the relevant differences once p is large, so BOTH realize exactly N0(s,r).

So 'does the general-r distinct-bad-scalar count match between prime and genuine
extension?' has the answer: YES, both equal N0(s,r), EXCEPT for finite-field saturation
(when q < N0, the count is capped at ~q on EITHER side -- a field-SIZE effect, not an
extension effect).  This file computes N0(s,r) exactly and the per-field realized count,
and the finite-field experiment confirms genuine GF(p^2) hits N0 when q is large enough.

s = 2^t (smooth) makes the cyclotomic arithmetic trivial: Phi_{2^t}(x) = x^{s/2}+1, so
xi^{s/2} = -1 and every power xi^j is a SIGNED unit vector in the basis {1,xi,..,xi^{s/2-1}}:
   xi^j =  +e_j         for 0 <= j < s/2
   xi^j =  -e_{j-s/2}   for s/2 <= j < s.
Thus e_1 is an integer vector of length s/2 with entries in [-r, r]; e_2 likewise.  Exact.
"""

from __future__ import annotations

import itertools
from math import comb


# ---------------------------------------------------------------------------
# Cyclotomic arithmetic for s = 2^t : Z[xi]/(xi^{s/2}+1).  Elements = int vectors
# of length d = s/2.  Multiplication is negacyclic convolution (x^d = -1).
# ---------------------------------------------------------------------------
def _is_pow2(n: int) -> bool:
    return n >= 1 and (n & (n - 1)) == 0


def power_vectors_pow2(s: int):
    """xi^j (j=0..s-1) as signed unit vectors of length d=s/2 in Z[xi]/(xi^{d}+1)."""
    assert _is_pow2(s) and s >= 4
    d = s // 2
    P = []
    for j in range(s):
        v = [0] * d
        if j < d:
            v[j] = 1
        else:
            v[j - d] = -1
        P.append(tuple(v))
    return P, d


def negacyclic_mul(a, b, d):
    """(a*b) mod (x^d + 1), integer vectors length d."""
    out = [0] * d
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for k, bk in enumerate(b):
            if bk == 0:
                continue
            idx = i + k
            if idx < d:
                out[idx] += ai * bk
            else:
                out[idx - d] -= ai * bk   # x^d = -1
    return out


# ---------------------------------------------------------------------------
# Char-0 brute distinct-lambda (exact), for validation / small (s,r).
# ---------------------------------------------------------------------------
def char0_distinct_brute(s: int, r: int, cap: int = 3_000_000):
    """EXACT char-0 distinct-lambda over firing r-subsets, brute force.
    Returns (distinct, num_firing) or (None,None) if binom(s,r)>cap."""
    if comb(s, r) > cap:
        return None, None
    P, d = power_vectors_pow2(s)
    lam = set()
    nf = 0
    for c in itertools.combinations(range(s), r):
        e1 = [0] * d
        e2 = [0] * d
        for j in c:
            pj = P[j]
            # e2 += xi^j * e1  (negacyclic), then e1 += xi^j
            prod = negacyclic_mul(pj, e1, d)
            for i in range(d):
                e2[i] += prod[i]
            for i in range(d):
                e1[i] += pj[i]
        if not any(e2):
            nf += 1
            lam.add(tuple(e1))
    return len(lam), nf


# ---------------------------------------------------------------------------
# Char-0 SMART DP (no binom(s,r)).  State = (e1_vec, e2_vec) as tuples; width is the
# number of distinct reachable integer (e1,e2) pairs at each size.  In char 0 there is
# NO field saturation, so this is the pure combinatorial count -- the width can grow,
# but for our (s,r) it stays modest because the entries are bounded in [-r,r].
# ---------------------------------------------------------------------------
def char0_distinct_dp(s: int, r: int, state_cap: int = 8_000_000):
    """EXACT char-0 distinct-lambda via Newton DP on integer (e1,e2) vectors.
    Returns (distinct, max_width) or (None, max_width) if a layer exceeds state_cap."""
    P, d = power_vectors_pow2(s)
    z = (0,) * d
    dp = [set() for _ in range(r + 1)]
    dp[0].add((z, z))
    max_w = 1
    for j in range(s):
        pj = P[j]
        top = max((k for k in range(r + 1) if dp[k]), default=0)
        for k in range(min(r - 1, top), -1, -1):
            cur = dp[k]
            if not cur:
                continue
            dst = dp[k + 1]
            for (e1, e2) in cur:
                prod = negacyclic_mul(pj, list(e1), d)
                ne1 = tuple(e1[i] + pj[i] for i in range(d))
                ne2 = tuple(e2[i] + prod[i] for i in range(d))
                dst.add((ne1, ne2))
            if len(dst) > state_cap:
                return None, len(dst)
        m = max(len(x) for x in dp)
        if m > max_w:
            max_w = m
    firing = {e1 for (e1, e2) in dp[r] if not any(e2)}
    return len(firing), max_w


if __name__ == "__main__":
    import time
    print("char-0 exact firing count N0(s,r): DP vs brute, and vs the r=4 formula")
    print("=" * 84)
    print("  (s/2-1)^2 is the KNOWN generic r=4 count; binom(s/2,r) is Kambire's reference")
    print("-" * 84)
    # validate DP==brute on small cases
    for s, r in [(8, 4), (16, 4), (16, 5), (16, 6), (8, 3), (16, 3)]:
        bf, nfb = char0_distinct_brute(s, r)
        dp, w = char0_distinct_dp(s, r)
        tag = "OK" if bf == dp else "*** MISMATCH ***"
        print(f"  s={s:>3} r={r}: brute={bf} dp={dp} {tag}  fire={nfb}  "
              f"(s/2-1)^2={(s//2-1)**2 if r==4 else '-'}  binom(s/2,r)={comb(s//2,r)}")
    print("-" * 84)
    # the law, DP only (push r), for s up to 64
    print("  N0(s,r) law (char-0 DP):")
    for s in [8, 16, 32, 64]:
        for r in range(3, s // 2 + 1):
            t = time.time()
            dp, w = char0_distinct_dp(s, r)
            dt = time.time() - t
            if dp is None:
                print(f"    s={s:>3} r={r:>2}: DP width {w} > cap, skipped")
                break
            print(f"    s={s:>3} r={r:>2} (rho={(r-2)/s:.4f}): N0={dp:>10}  "
                  f"binom(s/2,r)={comb(s//2,r):>14}  width={w:>9}  [{dt:.1f}s]")
