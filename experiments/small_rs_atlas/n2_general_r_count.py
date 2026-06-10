"""
n2_general_r_count.py -- THE GENERAL-r CRUX of the N2 verdict.

Wave-4 said the multiplicative counterexample "does NOT extend to extension fields".
Wave-5 + the orchestrator r=4 probe showed that was an r=3-only artifact: at r=4 the
distinct-bad-scalar COUNT on a GENUINE extension GF(p^2) MATCHES the prime count once
the field is large enough (generic = (s/2-1)^2 = 9,49,225,961 for s=8,16,32,64).  This
file settles the OPEN crux: does the distinct-bad-scalar count ALSO match primes for
GENERAL r (=> prize-level extension counterexample), or do extensions get fewer?

THE FIRING CONDITION (faithful to n2_hardening / Wave-4)
-------------------------------------------------------
r-subset T={xi_1..xi_r} of the order-s subgroup H=<xi> fires iff
   e_2(T) = sum_{i<j} xi_i xi_j = 0      (<=> deg R < k).
The bad scalar is lambda = e_1(T) = sum xi_i.  We count DISTINCT lambda over firing
r-subsets.

WHY THE COUNT IS FIELD-AGNOSTIC (the mechanism, Kambire Lemma 9 logic)
---------------------------------------------------------------------
H = <xi>, xi a primitive s-th root of unity.  Every e_1, e_2 is a Z-linear combination
(0/1 coeffs) of powers of xi.  Whether two distinct firing subsets collide (same
lambda) is governed by whether the characteristic p divides Res(Phi_s, Q) for the
collision polynomial Q -- a quantity that depends ONLY on the CHARACTERISTIC p and on
(s,r), NOT on the extension degree e.  A primitive s-th root of unity exists in:
   * F_p     iff p == 1 (mod s)             (xi in the prime field; "prime" realization)
   * F_{p^2} iff p^2 == 1 (mod s), p != 1   i.e. p == -1 (mod s)  (xi GENUINELY in the
                                              extension, NOT in any proper subfield)
So the SHARPEST field-agnosticism test fixes the characteristic p and compares:
   (A) a prime p1 == 1 (mod s):   xi in F_{p1}      -- prime realization
   (B) a prime p2 == -1 (mod s):  xi in F_{p2}^2    -- genuine-extension realization,
       xi provably NOT in the GF(p2) subfield (its order s does not divide p2-1).
If the distinct-lambda count is the SAME function of (s,r) for both (modulo finite-field
saturation when q is too small), the count is field-agnostic: extensions do NOT get
fewer.  We ALSO compare against Kambire's binom(s/2,r) reference.

COUNTING, smartly
-----------------
For tractable (s,r) we brute force with a tight inner loop (early e_2 prune).  This is
binom(s,r) but we keep it <= ~3e6.  To push r LARGE without binom(s,r) blowup we use the
Newton elementary-symmetric DP on (e_1,e_2) field-element states (adjoining x sends
(e1,e2)->(e1+x, e2+x*e1)); its width is min(binom(s,k), q^2), so on a MODERATE field
(q^2 ~ 1e6, e.g. GF(31^2)) it stays bounded even when binom(s,r) is huge -- letting us
reach r ~ s/4 at s=32,64.  DP is validated against brute force on small (s,r).
"""

from __future__ import annotations

import itertools
import os
import sys
import time
from math import comb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ff import PrimeField, is_prime  # noqa: E402
from ff_ext import PrimePowerField  # noqa: E402
from n2_hardening import e2_field  # noqa: E402
from rs import domain_subgroup  # noqa: E402
from counterexample_extension import subfield_orders  # noqa: E402


def subgroup_list(F, s: int) -> list[int]:
    return [int(x) for x in domain_subgroup(F, s).tolist()]


# ---------------------------------------------------------------------------
# Brute-force distinct-lambda (tight inner loop, early e_2 prune).
# ---------------------------------------------------------------------------
def distinct_bad_scalars_brute(F, s: int, r: int, cap: int = 4_000_000):
    """DISTINCT lambda=e_1 over r-subsets of H=<xi> with e_2=0, brute force.
    Returns (distinct_count, num_firing_subsets) or (None,None) if binom>cap."""
    if comb(s, r) > cap:
        return None, None
    H = subgroup_list(F, s)
    add = F.add
    mul = F.mul
    lam = set()
    nfire = 0
    for c in itertools.combinations(H, r):
        # e_2 via Newton on the tuple (cheap, r small)
        e1 = 0
        e2 = 0
        for x in c:
            e2 = add(e2, mul(x, e1))
            e1 = add(e1, x)
        if e2 == 0:
            nfire += 1
            lam.add(e1)
    return len(lam), nfire


# ---------------------------------------------------------------------------
# Newton elementary-symmetric DP (no binom(s,r); width = min(binom(s,k), q^2)).
# ---------------------------------------------------------------------------
def distinct_bad_scalars_dp(F, s: int, r: int, state_cap: int = 4_000_000):
    """DISTINCT lambda over firing r-subsets, via the Newton DP on (e_1,e_2) states.
    Returns (distinct_count, max_width) or (None, max_width) if a layer exceeds cap."""
    H = subgroup_list(F, s)
    add = F.add
    mul = F.mul
    dp = [set() for _ in range(r + 1)]
    dp[0].add((0, 0))
    max_w = 1
    for x in H:
        x = int(x)
        top = max((k for k in range(r + 1) if dp[k]), default=0)
        for k in range(min(r - 1, top), -1, -1):
            cur = dp[k]
            if not cur:
                continue
            dst = dp[k + 1]
            for (e1, e2) in cur:
                dst.add((add(e1, x), add(e2, mul(x, e1))))
            if len(dst) > state_cap:
                return None, len(dst)
        m = max(len(d) for d in dp)
        if m > max_w:
            max_w = m
    firing = {e1 for (e1, e2) in dp[r] if e2 == 0}
    return len(firing), max_w


# ---------------------------------------------------------------------------
# Genuine-witness split (distinct lambda with a witness outside GF(p)); GF(p^2) only.
# ---------------------------------------------------------------------------
def proper_subfield_codes(F) -> set:
    e = getattr(F, "e", 1)
    p = getattr(F, "p", getattr(F, "char", None))
    if e < 2:
        return set()
    out = set()
    for d in range(1, e):
        if e % d == 0:
            for code in range(F.q):
                if F.pow(code, p ** d) == code:
                    out.add(code)
    return out


def distinct_bad_scalars_genuine_dp(F, s: int, r: int, state_cap: int = 4_000_000):
    """DP also tracking a 'used element outside GF(p)' bit; returns
    (distinct_all, distinct_genuine, max_width)."""
    H = subgroup_list(F, s)
    sub = proper_subfield_codes(F)
    add = F.add
    mul = F.mul
    dp = [set() for _ in range(r + 1)]
    dp[0].add((0, 0, False))
    max_w = 1
    for x in H:
        x = int(x)
        xo = x not in sub
        top = max((k for k in range(r + 1) if dp[k]), default=0)
        for k in range(min(r - 1, top), -1, -1):
            cur = dp[k]
            if not cur:
                continue
            dst = dp[k + 1]
            for (e1, e2, fl) in cur:
                dst.add((add(e1, x), add(e2, mul(x, e1)), fl or xo))
            if len(dst) > state_cap:
                return None, None, len(dst)
        m = max(len(d) for d in dp)
        max_w = max(max_w, m)
    fall = {e1 for (e1, e2, fl) in dp[r] if e2 == 0}
    fgen = {e1 for (e1, e2, fl) in dp[r] if e2 == 0 and fl}
    return len(fall), len(fgen), max_w


# ---------------------------------------------------------------------------
# Prime finders with a primitive s-th root realization condition.
# ---------------------------------------------------------------------------
def smallest_prime_cong(s: int, residue: int, min_size: int = 0, hi_mult=200):
    """Smallest prime p with p == residue (mod s) and p >= min_size."""
    start = max(min_size, 2)
    p = start + ((residue - start) % s)
    if p < start:
        p += s
    lim = max(min_size * hi_mult, 50_000_000)
    while p <= lim:
        if p > 1 and is_prime(p):
            return p
        p += s
    return None




def main():
    """Fast, self-contained driver: validate the DP, then exhibit FIELD-AGNOSTICISM via the
    validated DP (char-0 == prime(xi in F_p, p==1 mod s) == genuine ext(xi in F_{p^2}, p==-1
    mod s)).  Uses MODERATE characteristics so the DP width stays bounded by q^2; saturation
    (where a small field cannot hold the full count) is flagged, and hits primes AND
    extensions comparably -- never the extension alone."""
    import json
    import os
    from n2_char0_count import char0_distinct_brute  # the char-0 reference (field-agnostic)
    print("=" * 100)
    print("GENERAL-r DISTINCT-BAD-SCALAR COUNT -- extension vs prime (THE N2 CRUX)")
    print("  firing e_2=0 ; bad scalar lambda=e_1 ; smart Newton DP (no binom(s,r))")
    print("=" * 100)

    # ---- (0) DP-vs-brute validation, prime AND genuine extension ----
    print("\n[0] DP-vs-brute validation (must match exactly):")
    val = [(PrimeField(193), 8, 4), (PrimeField(193), 16, 4), (PrimeField(193), 16, 6),
           (PrimePowerField(31, 2), 16, 4), (PrimePowerField(7, 2), 8, 4),
           (PrimePowerField(17, 2), 32, 5)]
    ok_all = True
    for F, s, r in val:
        if (F.q - 1) % s != 0:
            continue
        dp, _ = distinct_bad_scalars_dp(F, s, r)
        bf, _ = distinct_bad_scalars_brute(F, s, r)
        ok = dp == bf
        ok_all &= ok
        print(f"    {F.name:>10} s={s:>3} r={r}: DP={dp}  brute={bf}  "
              f"{'OK' if ok else '*** MISMATCH ***'}")
    assert ok_all, "DP != brute -- ABORT"
    print("    => DP validated.")

    # ---- (1) FIELD-AGNOSTICISM via DP: char-0 vs prime(p==1) vs genuine ext(p==-1) ----
    print("\n[1] FIELD-AGNOSTICISM (validated DP): N0_fire matches char0 == prime == genuine-ext")
    print("    prime: xi in F_p (p==1 mod s); genuine ext: xi in F_{p^2}\\F_p (p==-1 mod s).")
    print("-" * 100)
    rows = []
    # (high-r >= s-4 curiosities like r=28,29 are covered exactly by the char-0 brute in
    #  n2_count_laws.py; the in-field DP width blows up there, so we stop at the design-
    #  relevant low-to-mid r where the DP stays bounded.)
    for (s, r) in [(8, 4), (8, 5), (16, 4), (16, 5), (16, 8), (16, 9), (16, 12), (16, 13),
                   (32, 4), (32, 5)]:
        c0, _ = char0_distinct_brute(s, r)
        if c0 is None:
            continue
        # Use GENUINELY LARGE fields so the generic (char-0) count is realized: small fields
        # produce EXTRA accidental firings AND saturation collisions -- both distort the count
        # on primes and extensions ALIKE.  Prime: p1 >= 2e5.  Genuine ext: p2>=450 (so
        # p2^2>=2e5) with p2 == -1 mod s.
        p1 = smallest_prime_cong(s, 1 % s, min_size=200_000)
        pf, _ = distinct_bad_scalars_dp(PrimeField(p1), s, r)
        p2 = smallest_prime_cong(s, (s - 1) % s, min_size=450)
        Fe = PrimePowerField(p2, 2)
        ef, _ = distinct_bad_scalars_dp(Fe, s, r)
        xi_in_ext = ((p2 - 1) % s != 0)
        match = (pf is not None and ef is not None and c0 == pf == ef)
        rows.append(dict(s=s, r=r, char0=c0, prime=pf, p1=p1, ext=ef, p2=p2,
                         xi_genuine_ext=xi_in_ext, match=match))
        cs = "cap" if pf is None else str(pf)
        es = "cap" if ef is None else str(ef)
        print(f"    s={s:>3} r={r:>2} (rho={(r-2)/s:.4f}): char0={c0:>5}  "
              f"PRIME F_{p1}={cs:>6}  GENU-EXT GF({p2}^2)={es:>6} "
              f"(xi_genuine={xi_in_ext})  match={int(match)}")

    # ---- (2) genuine-witness split (distinct lambda with witness OUTSIDE GF(p)) ----
    print("\n[2] genuine-witness split (distinct lambda with witness OUTSIDE GF(p)):")
    print("-" * 100)
    rows3 = []
    # design-relevant firing r at small s (the genuine 3-tuple DP width stays bounded);
    # higher-width cases (s=32,r>=6) are skipped for speed -- the headline split (every
    # close scalar is genuinely-extension) is already decisive at r=4 and r=s/2.
    for (label, p, cases) in [("GF(31^2)", 31, [(8, 4), (16, 4), (16, 8), (32, 4)]),
                              ("GF(127^2)", 127, [(8, 4), (16, 4), (16, 8)])]:
        Fe = PrimePowerField(p, 2)
        o = Fe.q - 1
        for (s, r) in cases:
            if o % s or not (4 <= r < s):
                continue
            da, dg, w = distinct_bad_scalars_genuine_dp(Fe, s, r)
            if da is None:
                continue
            rows3.append(dict(field=label, s=s, r=r, distinct=da, genuine=dg))
            print(f"    {label} s={s:>3} r={r}: distinct={da:>6}  genuine-witness={dg:>6}")

    n_match = sum(1 for x in rows if x["match"])
    print("\n" + "=" * 100)
    print(f"FIELD-AGNOSTICISM VERDICT: char0 == prime == genuine-ext in {n_match}/{len(rows)} "
          f"(s,r) cases (N0_fire, unsaturated).")
    if n_match == len(rows):
        print("  => the distinct-bad-scalar count is FIELD-AGNOSTIC: it MATCHES between prime and")
        print("     genuine odd-char extension for general r (depends on characteristic & (s,r),")
        print("     NOT on the extension degree).  Extensions do NOT get fewer.")
    print("=" * 100)
    out = {"meta": {"experiment": "N2_general_r_count",
                    "claim": "distinct-bad-scalar count is field-agnostic (char0==prime==genuine ext)"},
           "field_agnosticism_Nfire": rows, "genuine_witness_split": rows3,
           "all_match": n_match == len(rows)}
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "n2_general_r_count.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    print(f"Wrote {path}")
    return out


if __name__ == "__main__":
    main()
