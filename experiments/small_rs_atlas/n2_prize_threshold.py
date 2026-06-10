"""
n2_prize_threshold.py -- IS the (field-agnostic) extension counterexample PRIZE-LEVEL?

We established (n2_count_laws.py, exact) that the number of DISTINCT bad scalars for the
Kambire/KK25 line f=X^{rm}, g=X^{(r-1)m} on a smooth subgroup H=<xi> of order s is a
CHARACTERISTIC-ZERO cyclotomic quantity -- identical for prime fields and genuine
odd-char extension fields GF(p^2), modulo finite-field saturation.  Two counts, per RS
dimension convention:

  N0_sum(s,r) = |H^{(+r)}|  (distinct r-fold subset sums)   -- Kambire-literal (deg<=k);
                this is the count Kambire's Theorem 1 uses, with the RIGOROUS LOWER BOUND
                (Lemma 9 / KK25):  N0_sum(s,r) >= binom(s/2, r)  for s=2^t (phi(s)=s/2),
                provided the characteristic is large enough that the subset sums do NOT
                collide.  Measured N0_sum grows ~2^{c s} (c~0.75 at rho=1/4), EXCEEDING
                the binom(s/2,r) lower bound -- so binom(s/2,r) is a safe under-estimate.

  N0_fire(s,r) = #distinct e_1 over e_2=0 r-subsets  -- STRICT (deg<k); = (s/2-1)^2 at
                r=4, and 0 unless r==0 or 1 (mod 4).  This is the count for the strict
                RS[F,L,k] (dimension k) convention.

THE PRIZE CONDITION (brief): a genuine epsilon_mca violation needs the number of
DISTINCT bad scalars to EXCEED  2^{-128} * |F|.  Equivalently the bad-scalar FRACTION
N0/|F| > 2^{-128}.  We ask: does an (s,r) with rho=(r-2)/s in {1/2,1/4,1/8,1/16},
s <= 2^40 (prize domain cap), and delta = capacity - 2/s, achieve this over a 256-bit
and a 128-bit field -- FOR GENUINE EXTENSION FIELDS?

The field-agnosticism means the answer is the SAME for prime and genuine extension,
EXCEPT we must additionally certify the field is large enough that the sums stay
distinct (no saturation) -- the Kambire/KK25 distinctness condition p > phi(s)^{phi(s)}
(prime) and its genuine-extension analogue (resultant not divisible by p).  We track
BOTH the threshold (count > 2^{-128}|F|) and the distinctness feasibility.
"""

from __future__ import annotations

from math import comb, log2, lgamma


def log2_binom(n: int, k: int) -> float:
    """log2 of binom(n,k) via lgamma (handles huge n,k)."""
    if k < 0 or k > n:
        return float("-inf")
    if k == 0 or k == n:
        return 0.0
    return (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / 0.6931471805599453


def Hbin(x: float) -> float:
    if x <= 0 or x >= 1:
        return 0.0
    return -x * log2(x) - (1 - x) * log2(1 - x)


def log2_factorial(n: int) -> float:
    return lgamma(n + 1) / 0.6931471805599453


def distinctness_floor_bits(s: int, mode: str = "linnik") -> float:
    """log2 of the FIELD SIZE needed for the subset sums to stay distinct (so the realized
    count equals the char-0 count, no saturation).  Three calibrations:

      mode='linnik'   (DEFAULT, what Kambire actually uses): a good characteristic exists
                       in the Linnik window p in [4^s, 8^s], i.e. log2 p in [2s, 3s].  We
                       use the UPPER end 3s (8^s) as the rigorous floor: a field of >= 3s
                       bits admits a distinctness-good characteristic == 1 (mod s).  This
                       is the operative finite-field constraint (kambire.md: 8^s <= 2^b
                       => s <= b/3).  For a genuine extension GF(p^2) the same window
                       applies to p (so |F|=p^2 of >= 2*(... )); we conservatively keep 3s
                       on |F| as the binding feasibility bound (a SUPERSET of what's
                       strictly needed -- if anything the extension is EASIER, p^2 buys
                       headroom -- so a 'no witness' under 3s is robust).
      mode='linnik_lo' the optimistic lower end 2s (4^s).
      mode='lemma9'    the very loose KK25 Lemma 9 sufficient condition p > phi(s)^{phi(s)}
                       (phi(s)=s/2), i.e. (s/2) log2(s/2) bits.  Way larger than needed;
                       shown only to bound the pessimistic extreme.
    """
    if mode == "lemma9":
        phi = s // 2
        return phi * log2(phi) if phi >= 2 else 1.0
    if mode == "linnik_lo":
        return 2.0 * s
    return 3.0 * s  # linnik (default)


def best_witness_for_rate(rho_name: str, rho: float, field_bits: int,
                          s_cap_log2: int = 40, count_kind: str = "sum",
                          floor_mode: str = "linnik"):
    """Sweep s=2^t (t<=s_cap_log2) with r=rho*s+2 integer in [3,s).  A WITNESS needs BOTH:
      (i) count > 2^{field_bits-128}  (the threshold), using
          count_kind='sum' : rigorous lower bound binom(s/2,r) on N0_sum (vacuous, =0, at
                              rho>=1/2 where r>s/2 -- Kambire's bound does not cover 1/2);
          count_kind='sum_entropy' : the measured-growth estimate 2^{H(2rho)*s/2} (an
                              ASYMPTOTIC, not a proven lower bound -- reported as an upper
                              sanity check, since measured N0_sum exceeds binom(s/2,r));
      (ii) distinctness feasible: field_bits >= distinctness_floor_bits(s, floor_mode)
           (the subset sums actually stay distinct in the field).
    Returns (witness_or_None, diagnostic) where diagnostic records the BINDING constraint
    and the feasibility window [s_lo (count ok), s_hi (distinctness ok)]."""
    target = field_bits - 128
    # s_hi: largest s with distinctness OK.  floor=3s (linnik) => s <= field_bits/3.
    # s_lo: smallest s with count > target.
    witness = None
    s_count_ok = []      # t where count exceeds target
    s_distinct_ok = []   # t where distinctness OK
    for t in range(2, s_cap_log2 + 1):
        s = 1 << t
        rs = rho * s
        if abs(rs - round(rs)) > 1e-9:
            continue
        r = int(round(rs)) + 2
        if not (3 <= r < s):
            continue
        if count_kind == "sum":
            lc = log2_binom(s // 2, r) if r <= s // 2 else float("-inf")
        elif count_kind == "sum_entropy":
            lc = Hbin(2 * rho) * (s // 2)
        else:
            raise ValueError(count_kind)
        floor = distinctness_floor_bits(s, floor_mode)
        count_ok = lc > target
        dist_ok = field_bits >= floor
        if count_ok:
            s_count_ok.append(t)
        if dist_ok:
            s_distinct_ok.append(t)
        cand = dict(rho=rho_name, field_bits=field_bits, s=s, t=t, r=r,
                    rho_val=(r - 2) / s, capacity=1 - (r - 2) / s, eta=2 / s,
                    log2_count=lc, target_bits=target, floor_bits=floor,
                    count_ok=count_ok, dist_ok=dist_ok)
        if count_ok and dist_ok and witness is None:
            witness = cand
    diag = dict(rho=rho_name, field_bits=field_bits, count_kind=count_kind,
                floor_mode=floor_mode,
                count_ok_t=(min(s_count_ok) if s_count_ok else None),
                distinct_ok_t_max=(max(s_distinct_ok) if s_distinct_ok else None),
                feasible=(witness is not None))
    return witness, diag


def main():
    print("=" * 108)
    print("PRIZE-LEVEL THRESHOLD: distinct bad scalars > 2^{-128}*|F|  (the eps_mca violation)")
    print("  Count is FIELD-AGNOSTIC => the SAME (s,r) works for PRIME and GENUINE EXTENSION.")
    print("  A witness needs BOTH (i) count > 2^{b-128} AND (ii) field big enough for the")
    print("  subset sums to stay distinct: |F| >= 8^s (Kambire's Linnik window), i.e. b>=3s.")
    print("=" * 108)
    for field_bits in [128, 256]:
        tb = field_bits - 128
        print(f"\n#### {field_bits}-bit field  (threshold count > 2^{tb}; distinctness needs "
              f"s <= b/3 = {field_bits//3}, i.e. t <= {(field_bits//3).bit_length()-1 if field_bits//3>0 else 0}+) ####")
        print(f"  N0_sum >= binom(s/2,r)  (Kambire-literal deg<=k; rigorous lower bound):")
        for rho_name, rho in [("1/16", 1/16), ("1/8", 1/8), ("1/4", 1/4), ("1/2", 1/2)]:
            w, diag = best_witness_for_rate(rho_name, rho, field_bits, count_kind="sum")
            if w is not None:
                print(f"    rho={rho_name}: *** WITNESS *** s=2^{w['t']}={w['s']}, r={w['r']}, "
                      f"rho=(r-2)/s={w['rho_val']:.4f}, delta=capacity-2/s="
                      f"{w['capacity']:.6f}-{w['eta']:.2e}={w['capacity']-w['eta']:.6f}, "
                      f"count>=2^{w['log2_count']:.1f}>2^{tb}, distinctness floor 2^{w['floor_bits']:.0f}<=2^{field_bits}")
            else:
                # explain the binding constraint
                ct = diag["count_ok_t"]
                dt = diag["distinct_ok_t_max"]
                if rho == 0.5:
                    print(f"    rho={rho_name}: NO WITNESS -- binom(s/2,r)=0 for rho>=1/2 (r>s/2); "
                          f"Kambire's count bound is VACUOUS at 1/2 (needs the footnote-18 "
                          f"list-size variant, not this proximity-gap count).")
                elif ct is None:
                    print(f"    rho={rho_name}: NO WITNESS -- count never exceeds 2^{tb} within s<=2^40.")
                elif dt is None or (ct is not None and dt is not None and ct > dt):
                    s_lo = 1 << ct
                    s_hi = (1 << dt) if dt is not None else 0
                    print(f"    rho={rho_name}: NO WITNESS -- INCOMPATIBLE: count>2^{tb} needs "
                          f"s>=2^{ct}={s_lo}, but distinctness needs s<=2^{dt}={s_hi} "
                          f"(8^s<=2^{field_bits}).  The window is EMPTY (s_lo>s_hi): the field "
                          f"large enough for the sums to be distinct is too small to hold 2^{tb} of them.")
                else:
                    print(f"    rho={rho_name}: NO WITNESS (count_ok t>={ct}, distinct_ok t<={dt}).")
    print("\n" + "-" * 108)
    print("  CROSS-CHECK with measured-growth estimate 2^{H(2rho)*s/2} (asymptotic; measured")
    print("  N0_sum EXCEEDS this, so it is a conservative stand-in for the true count):")
    for field_bits in [128, 256]:
        tb = field_bits - 128
        msgs = []
        for rho_name, rho in [("1/4", 1/4), ("1/8", 1/8), ("1/16", 1/16)]:
            w, diag = best_witness_for_rate(rho_name, rho, field_bits,
                                            count_kind="sum_entropy")
            msgs.append(f"rho={rho_name}:{'WIT s=2^%d' % w['t'] if w else 'none'}")
        print(f"    {field_bits}-bit: " + "  ".join(msgs))
    print("\n" + "-" * 108)
    print("  Strict deg<k firing count N0_fire(s,4)=(s/2-1)^2 (rate rho=2/s, r=4 only):")
    for field_bits in [128, 256]:
        tb = field_bits - 128
        wit = None
        for t in range(2, 41):
            s = 1 << t
            cnt = (s // 2 - 1) ** 2
            if cnt > 0 and log2(cnt) > tb and field_bits >= distinctness_floor_bits(s):
                wit = (s, t, log2(cnt))
                break
        if wit:
            s, t, lc = wit
            print(f"    {field_bits}-bit: witness s=2^{t}={s} (rho=2/s={2/s:.2e}, NOT a target "
                  f"rate), count=2^{lc:.1f}>2^{tb}, floor 2^{distinctness_floor_bits(s):.0f}")
        else:
            mx_t = 40
            print(f"    {field_bits}-bit: NO witness with distinctness -- at the max distinct-"
                  f"feasible s=2^{field_bits//3 if False else int(field_bits/3)} the count "
                  f"(s/2-1)^2 ~ 2^{2*(log2(field_bits/6)):.1f} << 2^{tb} (firing count grows "
                  f"only ~s^2, far too slow).")
    # ---- field-bit boundary b* (largest b admitting a witness, rho=1/4 binom LB) ----
    exists = [b for b in range(80, 300)
              if best_witness_for_rate("1/4", 0.25, b, count_kind="sum")[0] is not None]
    b_star = max(exists) if exists else None
    print(f"\n  FIELD-BIT BOUNDARY (rho=1/4, binom LB, 8^s distinctness): witness exists iff")
    print(f"  b <= b* = {b_star} bits.  So 128-bit IS prize-level, 256-bit (and any b>={b_star+1}) is NOT.")
    print("=" * 108)
    return b_star


if __name__ == "__main__":
    import json
    import os
    b_star = main()
    # save a compact artifact
    out = {"meta": {"experiment": "N2_prize_threshold",
                    "condition": "distinct bad scalars > 2^{-128}*|F|",
                    "count_field_agnostic": True,
                    "distinctness_floor": "|F| >= 8^s (Kambire Linnik window), b>=3s"},
           "field_bit_boundary_b_star": b_star,
           "witnesses": {}}
    for fb in [128, 256]:
        out["witnesses"][fb] = {}
        for rn, rv in [("1/16", 1/16), ("1/8", 1/8), ("1/4", 1/4), ("1/2", 1/2)]:
            w, diag = best_witness_for_rate(rn, rv, fb, count_kind="sum")
            out["witnesses"][fb][rn] = {"witness": w, "feasible": w is not None,
                                        "diagnostic": diag}
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "n2_prize_threshold.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    print(f"Wrote {path}")
