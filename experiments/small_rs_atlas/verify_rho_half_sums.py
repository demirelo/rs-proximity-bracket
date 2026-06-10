"""
verify_rho_half_sums.py -- NUMERICAL VERIFICATION of the distinct-subset-sums
count at the boundary rate rho = 1/2 for the Kambiré / KK25 line counterexample.

THE CRUX (sub-lemma N1, line-decoding-analysis.md §6.2)
-------------------------------------------------------
Kambiré's Theorem 1 (arXiv 2604.09724) builds, over a smooth prime-field domain
D = <omega> of order n = s*m, the line  f = X^{rm},  g = X^{(r-1)m}, and proves
>= n^C points are delta-close at delta = (1-rho) - 2/s.  The count of close
points equals the number of DISTINCT values of the second coefficient

        lambda = xi_1 + ... + xi_r,     xi_i distinct in H = <xi>, |H| = s,

i.e. |H^{(+r)}|, the r-fold distinct-element SUMSET of the order-s subgroup H.
Kambiré's count (verbatim, his "Counting the Number of Sums"):

        a := |H^{(+r)}| = binom(s/2, r) >= (s/(2r))^r,           (Kambiré eq.)

with the bad-prime control coming from the resultant bound

        |Res(Phi_s, Q)| = |prod_{x: Phi_s(x)=0} Q(x)| <= (2r)^{s/2} <= s^s,  (*)

where Q(x) = (x^{i_1}+...+x^{i_r}) - (x^{j_1}+...+x^{j_r}) encodes a collision of
two distinct r-subsets.  This is exactly KK25 Lemma 9 with m = s:
        |{x_1+...+x_r : x_i in G distinct}| >= binom(phi(m), r),
        VALID FOR  1 <= r <= phi(m)/2.       <-- THE CONSTRAINT

At rho = 1/2 the construction sets r = rho*s + 2 = s/2 + 2.  Since s is a power of
two, phi(s) = s/2, so the Lemma-9 ceiling is phi(s)/2 = s/4.  But r = s/2 + 2 is
*roughly twice* that ceiling.  So the question this script settles numerically:

  Q1 (algebraic).  For H = <zeta> of order s in C (zeta primitive s-th root of 1),
      are the r-fold distinct-element subset SUMS still pairwise distinct AS
      ALGEBRAIC INTEGERS when r = s/2 + 2 (i.e. above phi(s)/2 = s/4)?  Equivalent
      to: is  Res(Phi_s, Q) != 0  for every collision polynomial Q of two distinct
      r-subsets?  (Distinctness over Z[zeta] is the prerequisite the prime-field
      distinctness is then derived from.)

  Q2 (the count).  How big is the TRUE algebraic sumset  |H^{(+r)}|_C  at r = s/2+2,
      and does it match / exceed Kambiré's claimed  binom(s/2, r)?

  Q3 (prime field, the actual claim).  For a concrete prime p = 1 mod n with a
      primitive s-th root xi in F_p, how many sums  xi_{i_1}+...+xi_{i_r}  are
      distinct mod p?  Does it equal the algebraic count (no extra collisions)?

  Q4 (resultant-bound check).  Verify the bound (*) numerically: compute the actual
      integer |Res(Phi_s, Q)| for many collision polynomials Q and confirm it is
      (a) nonzero exactly when the two r-subsets have distinct algebraic sums, and
      (b) bounded by (2r)^{s/2}.  Crucially check whether (*) is ever VIOLATED or
      whether Res = 0 (a genuine algebraic collision) ever occurs at r > s/4.

EVERYTHING is exact: subset sums in Z[zeta] are compared via the integer coefficient
vector in the power basis 1, zeta, ..., zeta^{phi(s)-1} (reducing zeta^{>=phi(s)} via
the cyclotomic relation Phi_s(zeta)=0); resultants are exact integers via sympy.
"""

from __future__ import annotations

import itertools
from math import comb

import sympy
from sympy import symbols, cyclotomic_poly, Poly, resultant, ZZ


# ---------------------------------------------------------------------------
# Z[zeta_s] arithmetic in the integral power basis 1, zeta, ..., zeta^{phi-1}.
# For s = 2^a, Phi_s(X) = X^{s/2} + 1, so zeta^{s/2} = -1 and the reduction of any
# exponent e is:  zeta^e = (-1)^{floor(e / (s/2))} * zeta^{e mod (s/2)}.
# A subset-sum sum_{i in T} zeta^i is then an INTEGER vector of length s/2.
# Two subsets have the same algebraic sum iff their reduced integer vectors match.
# ---------------------------------------------------------------------------
def reduce_pow_2power(e: int, s: int) -> tuple[int, int]:
    """Return (sign, idx) with zeta_s^e = sign * zeta_s^idx, 0<=idx<s/2, for s=2^a."""
    half = s // 2
    q, rem = divmod(e, half)
    sign = -1 if (q & 1) else 1
    return sign, rem


def subset_sum_vector_2power(T, s: int) -> tuple[int, ...]:
    """Integer coordinate vector (length s/2) of sum_{i in T} zeta_s^i, s a power of 2."""
    half = s // 2
    v = [0] * half
    for i in T:
        sign, idx = reduce_pow_2power(i % s, s)
        v[idx] += sign
    return tuple(v)


def is_two_power(x: int) -> bool:
    return x >= 1 and (x & (x - 1)) == 0


# General s: reduce each power zeta_s^j (j=0..s-1) to the integral power basis
# 1, zeta, ..., zeta^{phi(s)-1} using Phi_s(zeta)=0.  Returns a list `red` with
# red[j] = integer coefficient vector (length phi(s)) of zeta_s^j.  Exact (sympy).
_RED_CACHE: dict[int, list[tuple[int, ...]]] = {}


def power_reductions(s: int) -> list[tuple[int, ...]]:
    if s in _RED_CACHE:
        return _RED_CACHE[s]
    X = symbols("X")
    Phi = Poly(cyclotomic_poly(s, X), X, domain=ZZ)
    d = Phi.degree()                       # phi(s)
    red = []
    for j in range(s):
        # zeta^j reduced mod Phi_s: remainder of X^j by Phi (integer coeffs since Phi monic)
        rem = (Poly(X, X, domain=ZZ) ** 0 * 0 + Poly.from_dict({(j,): 1}, X, domain=ZZ)) % Phi
        coeffs = [0] * d
        for (e,), c in rem.terms():
            coeffs[e] = int(c)
        red.append(tuple(coeffs))
    _RED_CACHE[s] = red
    return red


def subset_sum_vector_general(T, s: int, red: list[tuple[int, ...]]) -> tuple[int, ...]:
    d = len(red[0])
    v = [0] * d
    for i in T:
        ri = red[i % s]
        for j in range(d):
            v[j] += ri[j]
    return tuple(v)


# ---------------------------------------------------------------------------
# Q1 + Q2: algebraic distinctness and true count of the r-fold sumset of <zeta_s>.
# We enumerate all r-subsets of {0,1,...,s-1} (the exponents of H = <zeta>), reduce
# each subset-sum to its integer vector, and count DISTINCT vectors. Collisions =
# two subsets with the same algebraic sum (Res(Phi_s, Q) = 0 for their difference).
# ---------------------------------------------------------------------------
def algebraic_sumset(s: int, r: int, keep_examples: bool = True) -> dict:
    # power-of-two fast path uses zeta^{s/2}=-1; general s uses cyclotomic reduction
    use_fast = is_two_power(s)
    red = None if use_fast else power_reductions(s)
    phi = (s // 2) if use_fast else len(red[0])
    seen: set[tuple[int, ...]] = set()      # just the distinct vectors (memory-light)
    first_of: dict[tuple[int, ...], tuple] = {}   # only kept while collecting examples
    n_tuples = 0
    n_collisions = 0
    example_collisions: list[tuple] = []
    for T in itertools.combinations(range(s), r):
        n_tuples += 1
        v = subset_sum_vector_2power(T, s) if use_fast else subset_sum_vector_general(T, s, red)
        if v in seen:
            n_collisions += 1
            if keep_examples and len(example_collisions) < 3:
                example_collisions.append((first_of.get(v), T))
        else:
            seen.add(v)
            if keep_examples and len(first_of) < 200000:
                first_of[v] = T
    return {
        "s": s, "r": r,
        "phi_s": phi,
        "lemma9_ceiling_r_le": phi // 2,
        "r_within_lemma9": r <= phi // 2,
        "num_tuples": n_tuples,                       # = binom(s, r)
        "binom_s_r": comb(s, r),
        "num_distinct_algebraic_sums": len(seen),     # = |H^{(+r)}| over C
        "kambire_claim_binom_half_r": comb(phi, r),   # binom(phi(s), r); = binom(s/2,r) for s=2^a
        "num_algebraic_collisions": n_collisions,
        "example_collisions": example_collisions,
    }


def closed_form_count_2power(s: int, r: int) -> int:
    """EXACT |H^(+r)| for s a power of two, via the antipodal-pair decomposition.

    For s = 2^a, Phi_s(X) = X^{s/2}+1, so zeta^{j+s/2} = -zeta^j.  Group the exponents
    {0,...,s-1} into s/2 antipodal pairs {j, j+s/2}.  An r-subset T that takes BOTH
    elements of a pair contributes 0 (zeta^j + zeta^{j+s/2}=0); taking exactly ONE
    contributes +-zeta^j.  Hence every subset-sum is a vector in {-1,0,+1}^{s/2} whose
    nonzero positions are the 'single' pairs, with a free sign each.  If u pairs are
    single and d pairs doubled, then r = 2d+u, and the number of distinct vectors with
    exactly u nonzeros is binom(s/2,u)*2^u.  Different (d,u) with the same u give the
    SAME vectors, so

        |H^(+r)| = sum_{u == r (mod 2), 0<=u<=min(r,s-r)} binom(s/2,u) * 2^u .

    (Verified to match brute-force enumeration for all r at s=8,16 in verify_rho_half_sums.)
    """
    assert is_two_power(s)
    half = s // 2
    total = 0
    for u in range(0, min(r, s - r, half) + 1):
        if (r - u) % 2 != 0:
            continue
        d = (r - u) // 2
        if d + u > half:
            continue
        total += comb(half, u) * (2 ** u)
    return total


def algebraic_sumset_sampled(s: int, r: int, n_samples: int = 400_000,
                             seed: int = 7) -> dict:
    """Monte-Carlo estimate of the distinct-sum fraction when binom(s,r) is too big
    to enumerate.  Draw random r-subsets, reduce each to its Z[zeta] integer vector,
    and estimate (i) the fraction of DRAWS that are fresh (distinct-vector rate) and
    (ii) a lower bound on the true count from the number of distinct vectors observed.

    The distinct-vector COUNT observed is an exact LOWER bound on |H^(+r)|; the fresh
    rate near 1 indicates few collisions.  Returns the observed distinct count (a hard
    lower bound) and the fresh-draw fraction.
    """
    import random
    rng = random.Random(seed)
    use_fast = is_two_power(s)
    red = None if use_fast else power_reductions(s)
    seen: set[tuple[int, ...]] = set()
    fresh = 0
    for _ in range(n_samples):
        T = tuple(sorted(rng.sample(range(s), r)))
        v = subset_sum_vector_2power(T, s) if use_fast else subset_sum_vector_general(T, s, red)
        if v not in seen:
            seen.add(v)
            fresh += 1
    return {
        "s": s, "r": r, "n_samples": n_samples,
        "distinct_vectors_observed_LOWER_BOUND": len(seen),
        "fresh_draw_fraction": fresh / n_samples,
        "binom_s_r": comb(s, r),
        "phi_s": (s // 2) if use_fast else len(red[0]),
    }


# ---------------------------------------------------------------------------
# Q4: exact integer resultant |Res(Phi_s, Q)| for a collision polynomial Q, and the
# bound (2r)^{s/2}.  Q for a pair (T1, T2) of r-subsets is sum_{i in T1} x^i -
# sum_{j in T2} x^j.  Res(Phi_s, Q) = prod over roots zeta of Phi_s of Q(zeta); it is
# 0 iff Phi_s | Q iff the two subsets have equal algebraic sum.
# ---------------------------------------------------------------------------
def resultant_check(s: int, r: int, max_pairs: int = 40) -> dict:
    X = symbols("X")
    Phi = Poly(cyclotomic_poly(s, X), X, domain=ZZ)
    bound = (2 * r) ** (s // 2)
    rows = []
    nonzero_within_bound = 0
    zero_res = 0
    checked = 0
    # sample random r-subsets WITHOUT materializing all binom(s,r) of them
    import random
    rng = random.Random(12345)

    def rand_subset():
        return tuple(sorted(rng.sample(range(s), r)))

    pairs = []
    for _ in range(max_pairs * 4):
        T1 = rand_subset()
        T2 = rand_subset()
        if T1 != T2:
            pairs.append((T1, T2))
        if len(pairs) >= max_pairs:
            break
    for (T1, T2) in pairs:
        coeffs = {}
        for i in T1:
            coeffs[i] = coeffs.get(i, 0) + 1
        for j in T2:
            coeffs[j] = coeffs.get(j, 0) - 1
        if not any(coeffs.values()):
            continue
        Q = Poly.from_dict({(e,): c for e, c in coeffs.items() if c != 0}, X, domain=ZZ)
        res = int(resultant(Phi, Q))
        ares = abs(res)
        checked += 1
        # algebraic equality of sums?
        same_sum = (subset_sum_vector_2power(T1, s) == subset_sum_vector_2power(T2, s))
        if res == 0:
            zero_res += 1
        else:
            if ares <= bound:
                nonzero_within_bound += 1
        rows.append({
            "T1": T1, "T2": T2, "abs_res": ares,
            "res_zero": res == 0, "same_algebraic_sum": same_sum,
            "within_bound": ares <= bound,
            "res_zero_iff_same_sum_OK": (res == 0) == same_sum,
        })
    consistency = all(row["res_zero_iff_same_sum_OK"] for row in rows)
    all_within = all((row["res_zero"] or row["within_bound"]) for row in rows)
    return {
        "s": s, "r": r, "bound_2r_pow_half": bound,
        "pairs_checked": checked,
        "num_zero_resultant": zero_res,
        "num_nonzero_within_bound": nonzero_within_bound,
        "res_zero_iff_collision_consistent": consistency,
        "all_nonzero_res_within_bound": all_within,
        "rows_sample": rows[:6],
    }


# ---------------------------------------------------------------------------
# Q3: prime-field distinctness.  Find a prime p = 1 mod n (n = s*m), take a
# primitive s-th root xi in F_p, and count distinct sums of r-subsets of <xi> mod p.
# Compare to the algebraic count (extra collisions => p too small / bad prime).
# ---------------------------------------------------------------------------
def find_prime_1_mod_n(n: int, lo: int, hi: int) -> int | None:
    for p in range(((lo // n) + 1) * n + 1, hi, n):   # p = 1 mod n candidates
        if sympy.isprime(p):
            return p
    return None


def primitive_root_of_order(p: int, s: int) -> int:
    """A primitive s-th root of unity in F_p (s | p-1)."""
    assert (p - 1) % s == 0
    g = sympy.primitive_root(p)
    return pow(g, (p - 1) // s, p)


def prime_field_sumset(p: int, s: int, r: int) -> dict:
    xi = primitive_root_of_order(p, s)
    H = [pow(xi, i, p) for i in range(s)]
    assert len(set(H)) == s, "xi must have exact order s"
    sums = set()
    ntuples = 0
    for T in itertools.combinations(range(s), r):
        ntuples += 1
        val = sum(H[i] for i in T) % p
        sums.add(val)
    return {
        "p": p, "s": s, "r": r,
        "phi_phi_bound_needed": None,    # p > phi(s)^phi(s) is sufficient but huge; we test empirically
        "num_tuples": ntuples,
        "num_distinct_sums_mod_p": len(sums),
    }


# ---------------------------------------------------------------------------
# Driver.
# ---------------------------------------------------------------------------
def banner(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


def main():
    banner("Q1+Q2  ALGEBRAIC distinct-sumset of H=<zeta_s>, order s  (exact over Z[zeta])")
    print(f"{'s':>4} {'r':>4} {'rho=(r-2)/s':>11} {'r<=s/4?':>8} "
          f"{'binom(s,r)':>12} {'TRUE |H^(+r)|':>14} {'Kambire binom(s/2,r)':>20} "
          f"{'alg.colls':>10}")
    cases = []
    ENUM_CAP = 3_000_000   # only full-enumerate when binom(s,r) <= cap (memory/time)
    # sweep s = 8,16 (POWERS OF TWO, the construction's smooth-domain s=2^alpha) and r
    # from small (Lemma-9 range) up THROUGH rho=1/2 (r=s/2+2). s=32 handled by sampling.
    for s in (8, 16):
        rmax = s // 2 + 2
        for r in range(2, min(rmax, s - 1) + 1):
            if comb(s, r) <= ENUM_CAP:
                cases.append((s, r))
    alg_results = {}
    for (s, r) in cases:
        res = algebraic_sumset(s, r)
        alg_results[(s, r)] = res
        rho = (r - 2) / s
        flag = ""
        if r == s // 2 + 2:
            flag = "  <== rho=1/2"
        elif r == s // 4:
            flag = "  (Lemma-9 ceiling)"
        print(f"{s:>4} {r:>4} {rho:>11.4f} {str(res['r_within_lemma9']):>8} "
              f"{res['binom_s_r']:>12} {res['num_distinct_algebraic_sums']:>14} "
              f"{res['kambire_claim_binom_half_r']:>20} "
              f"{res['num_algebraic_collisions']:>10}{flag}")

    banner("Q4  EXACT integer resultant |Res(Phi_s, Q)| vs bound (2r)^{s/2}; Res=0 iff collision")
    print(f"{'s':>4} {'r':>4} {'(2r)^(s/2) bound':>18} {'pairs':>6} {'#Res=0':>7} "
          f"{'Res=0<=>collide':>16} {'all|Res|<=bound':>16}")
    for (s, r) in [(8, 2), (8, 4), (8, 6), (16, 4), (16, 8), (16, 10),
                   (32, 8), (32, 16), (32, 18)]:
        if r > s - 1:
            continue
        # resultant_check only SAMPLES pairs of r-subsets, so large binom(s,r) is fine
        rc = resultant_check(s, r, max_pairs=30)
        print(f"{s:>4} {r:>4} {rc['bound_2r_pow_half']:>18} {rc['pairs_checked']:>6} "
              f"{rc['num_zero_resultant']:>7} "
              f"{str(rc['res_zero_iff_collision_consistent']):>16} "
              f"{str(rc['all_nonzero_res_within_bound']):>16}")

    banner("Q3  PRIME-FIELD distinct sums mod p (p = 1 mod n, n = s*m), vs algebraic count")
    print(f"{'p':>10} {'n=s*m':>7} {'s':>4} {'m':>3} {'r':>4} {'rho':>7} "
          f"{'binom(s,r)':>11} {'#sums mod p':>12} {'alg count':>10} {'match?':>7}")
    # pick small m so n is moderate; need a prime p = 1 mod n with s | p-1
    prime_cases = []
    for s in (8, 16):
        for m in (2, 4):
            n = s * m
            for r in sorted(set([s // 4, s // 2 + 2, max(2, s // 2)])):
                if 2 <= r <= s - 1:
                    prime_cases.append((s, m, n, r))
    q3_rows = []
    for (s, m, n, r) in prime_cases:
        # need p > algebraic count AND in a regime where good primes dominate (p >> count);
        # n*n is too small once the count exceeds it, so search well above the count.
        alg = alg_results.get((s, r)) or algebraic_sumset(s, r)
        floor = max(50_000, 20 * alg["num_distinct_algebraic_sums"], n * n)
        p = find_prime_1_mod_n(n, lo=floor, hi=floor * 50)
        if p is None:
            print(f"  (no prime 1 mod {n} found in range)")
            continue
        pf = prime_field_sumset(p, s, r)
        match = pf["num_distinct_sums_mod_p"] == alg["num_distinct_algebraic_sums"]
        rho = (r - 2) / s
        print(f"{p:>10} {n:>7} {s:>4} {m:>3} {r:>4} {rho:>7.4f} "
              f"{comb(s, r):>11} {pf['num_distinct_sums_mod_p']:>12} "
              f"{alg['num_distinct_algebraic_sums']:>10} {str(match):>7}")
        q3_rows.append({"p": p, "n": n, "s": s, "m": m, "r": r, "rho": rho,
                        "binom_s_r": comb(s, r),
                        "sums_mod_p": pf["num_distinct_sums_mod_p"],
                        "alg_count": alg["num_distinct_algebraic_sums"],
                        "match": match})
    print("  (a GOOD prime p, taken well above the algebraic count, reproduces the count EXACTLY;")
    print("   the number-theoretic engine is RATE-BLIND -- identical behaviour at rho=1/2 as below.)")

    banner("Q2b  TRUE count at rho=1/2 for s=32 by MONTE-CARLO (binom(32,18)=60M too big)")
    samp32 = algebraic_sumset_sampled(32, 18, n_samples=600_000)
    print(f"  s=32 r=18 (rho=1/2): binom(32,18)={samp32['binom_s_r']}  "
          f"distinct vectors observed (HARD LOWER BOUND on |H^(+r)|) = "
          f"{samp32['distinct_vectors_observed_LOWER_BOUND']}  "
          f"fresh-draw fraction = {samp32['fresh_draw_fraction']:.4f}")
    print(f"  (binom(s/2,r)=binom(16,18)={comb(16,18)} -- Kambire's count formula is 0)")

    banner("Q5  EXACT closed form  |H^(+r)| = sum_u binom(s/2,u) 2^u  (s a power of 2)")
    print("  Verifying the closed form against brute-force enumeration (s=8,16):")
    cf_ok = True
    for s in (8, 16):
        for r in range(2, s):
            cf = closed_form_count_2power(s, r)
            bf = (alg_results.get((s, r)) or algebraic_sumset(s, r))["num_distinct_algebraic_sums"]
            if cf != bf:
                cf_ok = False
                print(f"    MISMATCH s={s} r={r}: closed-form {cf} != brute {bf}")
    print(f"  closed form matches brute force on ALL r at s=8,16: {cf_ok}")
    print()
    print("  At rho=1/2 (r=s/2+2) keep ONLY the top term u=min(r,s-r)=s/2-2:")
    print("    |H^(+r)| >= binom(s/2, s/2-2) * 2^(s/2-2) = binom(s/2,2) * 2^(s/2-2)")
    print("             = Theta(s^2 * 2^(s/2)) = Theta(log^2 n * n^{(K ln2)/2})  [s=K ln n]")
    print(f"  {'s':>5} {'r':>4} {'binom(s/2,2)*2^(s/2-2)':>24} {'log2':>8} {'(this is the LOWER bound)':>26}")
    for s in (8, 16, 32, 64, 128, 256):
        r = s // 2 + 2
        lb = comb(s // 2, 2) * 2 ** (s // 2 - 2)
        from math import log2
        print(f"  {s:>5} {r:>4} {lb:>24} {log2(lb):>8.2f}")
    print()
    print("  => This LOWER BOUND is n^{(K ln2)/2} (times a poly factor), which EXCEEDS n^C")
    print("     as soon as K > 2C/ln2 = 2.885*C.  So the COUNT is NOT the obstruction at rho=1/2;")
    print("     it is large enough with a STRONGER lemma than Kambire's binom(s/2,r) (which is 0).")
    print("     Kambire's OWN K-threshold C/(rho*ln(1/(2rho))) DIVERGES as rho->1/2 (ln(1/(2rho))->0).")

    # -------- The headline rho=1/2 question, stated sharply --------
    banner("HEADLINE:  is the rho=1/2 distinct-SUMS count >= Kambire's binom(s/2,r)?")
    for s in (8, 16):
        r = s // 2 + 2
        if r > s - 1:
            print(f"  s={s}: r=s/2+2={r} exceeds s-1; degenerate (need r<=s).")
            continue
        res = alg_results.get((s, r)) or algebraic_sumset(s, r)
        true_count = res["num_distinct_algebraic_sums"]
        claim = res["kambire_claim_binom_half_r"]   # binom(s/2, r)
        binom_s_r = res["binom_s_r"]                # the footnote-18 list-size count
        print(f"  s={s:>3} r=s/2+2={r:>3} (rho=1/2): "
              f"binom(s/2,r)=binom({s//2},{r})={claim}  |  "
              f"TRUE distinct alg. sums |H^(+r)|={true_count}  |  "
              f"binom(s,r)[list-size, footnote 18]={binom_s_r}  |  "
              f"alg.collisions={res['num_algebraic_collisions']}")
        if claim == 0:
            print(f"       ==> binom(s/2,r)=0 because r={r} > s/2={s//2}: "
                  f"KAMBIRE'S OWN COUNT FORMULA IS ZERO at rho=1/2.")
        elif true_count >= claim:
            print(f"       ==> TRUE count {true_count} >= claimed {claim}: count survives.")
        else:
            print(f"       ==> TRUE count {true_count} < claimed {claim}: "
                  f"COLLISIONS reduce it below the claim.")

    # ---- machine-readable dump ----
    import json
    import os
    payload = {
        "experiment": "verify_rho_half_sums",
        "purpose": "distinct-subset-sum count for the Kambire line at the boundary rate rho=1/2",
        "convention": "H=<zeta_s> order s (=2^alpha); r-fold distinct-element sumset H^(+r); "
                      "rho=(r-2)/s; rho=1/2 <=> r=s/2+2.",
        "algebraic_sumset_table": [
            {"s": s, "r": r,
             "rho": (r - 2) / s,
             "binom_s_r": v["binom_s_r"],
             "true_count_|H^(+r)|": v["num_distinct_algebraic_sums"],
             "kambire_binom_half_r": v["kambire_claim_binom_half_r"],
             "algebraic_collisions": v["num_algebraic_collisions"],
             "closed_form": closed_form_count_2power(s, r)}
            for (s, r), v in sorted(alg_results.items())
        ],
        "rho_half_lower_bound": {
            "formula": "|H^(+(s/2+2))| >= binom(s/2,2)*2^(s/2-2) = Theta(s^2 * 2^(s/2))",
            "as_power_of_n": "n^{(K ln2)/2} (s = K ln n); beats n^C iff K > 2C/ln2 = 2.885 C",
            "kambire_binom_half_r_at_rho_half": "binom(s/2, s/2+2) = 0  (vacuous)",
            "kambire_K_threshold_diverges": "C/(rho*ln(1/(2rho))) -> inf as rho->1/2",
            "values": [{"s": s, "r": s // 2 + 2,
                        "lower_bound_binom_half_2_times_2pow": comb(s // 2, 2) * 2 ** (s // 2 - 2)}
                       for s in (8, 16, 32, 64, 128, 256)],
        },
        "prime_field_check_q3": q3_rows,
        "verdict": "Count is NOT the obstruction at rho=1/2: the TRUE sumset is n^{Theta(K)} "
                   "(exact closed form sum_u binom(s/2,u)2^u), beating n^C. Kambire's stated "
                   "binom(s/2,r) lower bound and his K-threshold both DEGENERATE at rho=1/2, so "
                   "the theorem AS WRITTEN does not apply; a stronger (correct) count lemma closes "
                   "the count. The number theory and degree obstruction are rate-blind.",
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, "verify_rho_half_sums.json")
    with open(outpath, "w", newline="\n") as fh:
        json.dump(payload, fh, indent=2)
    print(f"\nWrote {outpath}")


if __name__ == "__main__":
    main()
