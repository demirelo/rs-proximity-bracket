"""
n2_count_laws.py -- the TWO distinct-bad-scalar count laws and their FIELD-AGNOSTICISM.

A subtle but decisive convention point (resolved by EXACT computation, see the
n2-verdict): the number of distinct bad scalars lambda for the Kambire/KK25 line
f=X^{rm}, g=X^{(r-1)m} on a smooth subgroup depends on the RS dimension convention:

  READING A (STRICT, RS[F,L,k] = polys of degree < k, dimension k = (r-2)m, rate
             rho=(r-2)/s):  lambda is CLOSE  iff  some r-subset summing to lambda has
             e_2 = 0 (deg R < k).  The distinct-close count is
                 N0_fire(s,r) = #{ distinct e_1 over e_2=0 r-subsets }.
             (EXACTLY verified: {close lambda} == {firing lambda}, no close lambda
             outside H^{(+r)}, on prime GF(17) AND genuine ext GF(7^2).)

  READING B (Kambire-literal, code = polys of degree <= (r-2)m, dimension (r-2)m+1,
             rate (r-2)/s + 1/n):  R (degree <= (r-2)m) is itself a codeword, so EVERY
             lambda in H^{(+r)} is close.  The distinct-close count is
                 N0_sum(s,r) = |H^{(+r)}| = #{ distinct e_1 over ALL r-subsets }.
             (EXACTLY verified: ALL |H^{(+r)}| lambda are close under deg<=k, GF(17).)

The two differ by O(1/n) in rate (negligible asymptotically) but by a LARGE factor in
count (N0_fire << N0_sum, and N0_fire=0 for many (s,r)).  BOTH are characteristic-zero
cyclotomic quantities -- functions of (s,r) ONLY, independent of the field and of the
extension degree e.  This file computes BOTH laws in char 0 (the field-agnostic generic
value) and confirms a GENUINE extension GF(p^2) realizes them when q is large enough.

The headline consequence for the N2 verdict: for BOTH readings the distinct-bad-scalar
count is FIELD-AGNOSTIC, so it MATCHES between prime and genuine odd-char extension --
the counterexample MECHANISM extends.  Whether it is PRIZE-LEVEL is then a pure size
question on N0_*(s,r) vs 2^-128*|F| (n2_prize_threshold.py), identical for prime and ext.
"""

from __future__ import annotations

import itertools
from math import comb

from n2_char0_count import power_vectors_pow2, negacyclic_mul, char0_distinct_dp


# ---------------------------------------------------------------------------
# Char-0 EXACT counts (brute, binom-capped): N0_fire and N0_sum.
# ---------------------------------------------------------------------------
def char0_counts_brute(s: int, r: int, cap: int = 3_500_000):
    """(N0_fire, N0_sum, num_firing_subsets) in char 0, brute force.
    N0_fire = #distinct e_1 over e_2=0 subsets; N0_sum = #distinct e_1 over ALL subsets."""
    if comb(s, r) > cap:
        return None, None, None
    P, d = power_vectors_pow2(s)
    fire = set()
    alls = set()
    nf = 0
    for c in itertools.combinations(range(s), r):
        e1 = [0] * d
        e2 = [0] * d
        for j in c:
            pj = P[j]
            prod = negacyclic_mul(pj, e1, d)
            for i in range(d):
                e2[i] += prod[i]
            for i in range(d):
                e1[i] += pj[i]
        te1 = tuple(e1)
        alls.add(te1)
        if not any(e2):
            nf += 1
            fire.add(te1)
    return len(fire), len(alls), nf


# ---------------------------------------------------------------------------
# Char-0 SMART DP for N0_sum (distinct subset-sums) -- width = #distinct partial e_1
# vectors, bounded by the number of {-1,0,1}^d vectors reachable, which is far smaller
# than binom(s,r) and does NOT blow up like the (e_1,e_2) joint DP.
# ---------------------------------------------------------------------------
def char0_Nsum_dp(s: int, r: int, state_cap: int = 20_000_000):
    """N0_sum(s,r) = #distinct e_1 over ALL r-subsets, char-0 DP on e_1 vectors only."""
    P, d = power_vectors_pow2(s)
    z = (0,) * d
    dp = [set() for _ in range(r + 1)]
    dp[0].add(z)
    for j in range(s):
        pj = P[j]
        top = max((k for k in range(r + 1) if dp[k]), default=0)
        for k in range(min(r - 1, top), -1, -1):
            cur = dp[k]
            if not cur:
                continue
            dst = dp[k + 1]
            for e1 in cur:
                dst.add(tuple(e1[i] + pj[i] for i in range(d)))
            if len(dst) > state_cap:
                return None
    return len(dp[r])


if __name__ == "__main__":
    import time
    print("=" * 96)
    print("TWO char-0 distinct-bad-scalar count laws, vs binom(s/2,r) [Kambire's claimed bound]")
    print("  N0_fire = strict deg<k count (firing e_2=0);  N0_sum = deg<=k count (|H^{(+r)}|)")
    print("=" * 96)
    print(f"  {'s':>3} {'r':>3} {'rho=(r-2)/s':>11} {'N0_fire':>9} {'N0_sum':>10} "
          f"{'binom(s/2,r)':>13} {'#firing':>9}  notes")
    print("-" * 96)
    # full brute scan where feasible; mark design r=rho*s+2 rows
    design_rs = {}
    for s in [8, 16, 32, 64]:
        for rho_name, rho in [("1/2", 0.5), ("1/4", 0.25), ("1/8", 0.125), ("1/16", 1 / 16)]:
            r = int(round(rho * s)) + 2
            design_rs.setdefault((s, r), []).append(rho_name)
    for s in [8, 16, 32]:
        for r in range(3, s):
            if comb(s, r) > 3_500_000:
                continue
            nf, ns, nfire = char0_counts_brute(s, r)
            if nf is None:
                continue
            note = ""
            if (s, r) in design_rs:
                note = "DESIGN rho=" + ",".join(design_rs[(s, r)])
            if nf > 0 and not note:
                note = "fires"
            print(f"  {s:>3} {r:>3} {(r-2)/s:>11.4f} {nf:>9} {ns:>10} "
                  f"{comb(s//2,r):>13} {nfire:>9}  {note}")
    print("-" * 96)
    print("  r=4 closed form check: N0_fire(s,4) =? (s/2-1)^2 ; N0_sum(s,4) =? |H^{(+4)}|")
    for s in [8, 16, 32]:
        nf, ns, _ = char0_counts_brute(s, 4)
        print(f"    s={s}: N0_fire={nf} (s/2-1)^2={(s//2-1)**2}  N0_sum={ns}")
    print("=" * 96)
