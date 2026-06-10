"""
counterexample_extension.py -- DOES the Kambire/BCHKS near-capacity proximity-gap
counterexample EXTEND from PRIME fields to ODD-CHARACTERISTIC EXTENSION fields
GF(p^e)?  (Proximity-Prize sub-lemma N2 -- the genuine prize opening.)

Why this is the decisive experiment
------------------------------------
Every unconditional negative result (Kambire arXiv:2604.09724; BCHKS ePrint
2025/2055) lives over PRIME fields.  CGHLL26 (eprint 2026/532, App. A.5 Conj 2)
ASSERTS line-decodability -- hence SAFETY from this counterexample -- for
extension-field alphabets, with NO evidence.  The construction is MULTIPLICATIVE
(it lives on a smooth subgroup D=<omega> and an inner subgroup H=<xi>), so the
question is NOT about additive subgroups; it is whether the subgroup's
distinct-element subset-sums stay DISTINCT in F_{p^e} as they do in F_p.  Wave-2
tested only char-2 GF(2^m) -- which has ODD group order and thus NO smooth subgroup,
so it cannot even host the construction.  Odd-char GF(p^e) (which DOES have smooth
subgroups: |GF(9)*|=8, |GF(49)*|=48, |GF(289)*|=288, ...) is the untested case.

The construction, field-agnostically (faithful to counterexample_kambire.py)
----------------------------------------------------------------------------
On D=<omega> of order n with inner subgroup H=<xi> of order s, m=n/s cosets:
   r = rho*s + 2,   k = (r-2)m,   rho = (r-2)/s,
   f = X^{rm},  g = X^{(r-1)m},   delta = 1 - r/s,   capacity = 1 - rho,  eta = 2/s.
For r DISTINCT xi_1..xi_r in H,
   prod_{j=1}^r (X^m - xi_j) = X^{rm} - lambda X^{(r-1)m} + R(X),  lambda = sum xi_j,
with deg R <= (r-2)m = k.  When deg R < k, R is a codeword, so f - lambda g = -R on
the size-rm coset union  ==>  dist(f - lambda g, C) <= delta*n.  There are |H^{(+r)}|
such certified scalars (the distinct-element r-fold sumset of H).  Yet g=X^{(r-1)m}
has degree (r-1)m > k so S* <= (r-1)m < rm = CA threshold: NO correlated agreement.

THE prime signal we are testing for: on the SMOOTH prime subgroup,
   close_count(gamma over all F)  ==  predicted_sumset_size  >>  close_count(random).
We ask whether this same signal appears on odd-char EXTENSION smooth subgroups.

What is field-specific (and what we changed)
--------------------------------------------
counterexample_kambire.py computes the residual-degree / sumset with Python-int
`% p` arithmetic (prime-only).  Here `_residual_degree_F` and `sumset_H_F` do the
SAME polynomial algebra but in the ACTUAL FIELD via F.add/F.mul/F.neg, so lambda and
the R-coefficients are genuine GF(p^e) elements.  Everything else (RS, exact
distance, exact S* via branch-and-bound, the f=X^{rm}/g=X^{(r-1)m} line, the
delta-grid near capacity) is reused unchanged from the prime code.

The DECISIVE COMPARISON, matched on (q, n, k, delta):
  (a) smooth <omega> subgroup of GF(p^e)         -- the construction's native domain
  (b) a random size-n subset of GF(p^e)*         -- structure-free baseline
  (c) a PRIME field F_{p'} of comparable size with its own smooth subgroup
                                                  -- positive control (must fire)

Plus the ALGEBRAIC CRUX, independent of RS: subset-sum distinctness of <xi> in
GF(p^e) vs in F_p (the r-fold distinct sumset AND the full distinct-subset-sum count).

Everything EXACT (no decoder).  Results -> results/counterexample_extension.{json,csv}.
"""

from __future__ import annotations

import csv
import json
import itertools
import math
import os
import time
from dataclasses import dataclass, asdict

import numpy as np

from ff import PrimeField, FiniteField, _divisors, is_prime
from ff_ext import PrimePowerField
from rs import (build_codeword_book, domain_subgroup, domain_random, dist_to_code,
                min_distance)
from search_bad_lines import _max_common_agreement, _agreement_bits
from counterexample_kambire import (monomial_eval, measure_line_on_domain,
                                     delta_window, KParams, _ncr)


RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
QK_CAP = 3_000_000


# ===========================================================================
# Field-agnostic residual-degree / sumset (the prime code's logic, in-field).
# ===========================================================================
def _poly_mulmod_F(F: FiniteField, a: list[int], b: list[int]) -> list[int]:
    """Multiply two polynomials with FIELD-ELEMENT coefficients (no reduction)."""
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            out[i + j] = F.add(out[i + j], F.mul(ai, bj))
    return out


def _residual_degree_F(F: FiniteField, combo, m: int, r: int):
    """deg of R(X) = prod_j (X^m - xi_j) - (X^{rm} - lambda X^{(r-1)m}), IN-FIELD.

    Faithful in-field port of counterexample_kambire._residual_degree.  The Kambire
    identity cancels the X^{rm} and X^{(r-1)m} terms, leaving R of degree <= (r-2)m.
    A scalar lambda = sum xi_j is a CERTIFIED-close point iff some representation has
    deg R < k = (r-2)m (then f - lambda g = -R on the size-rm coset union, so it is
    delta-close).  Returns (deg R, lambda) with lambda a genuine field element.
    """
    poly = [1]                       # the constant polynomial 1
    lam = 0                          # field element accumulator for sum xi_j
    for xi in combo:
        lam = F.add(lam, int(xi))
        fac = [0] * (m + 1)
        fac[0] = F.neg(int(xi))      # (X^m - xi): constant term -xi, leading 1
        fac[m] = 1
        poly = _poly_mulmod_F(F, poly, fac)
    R = poly[:]
    # subtract (X^{rm} - lambda X^{(r-1)m})
    R[r * m] = F.sub(R[r * m], 1)
    R[(r - 1) * m] = F.add(R[(r - 1) * m], lam)
    deg = max((i for i, c in enumerate(R) if c != 0), default=-1)
    return deg, lam


def sumset_H_F(F: FiniteField, s: int, r: int, m: int, k: int) -> dict:
    """In-field analysis of the distinct-element r-fold sumset H^{(+r)} of H=<xi>.

    Mirrors counterexample_kambire.sumset_H but with FIELD arithmetic, so collisions
    are measured in GF(p^e), which is exactly the algebraic crux.  Returns dict with
    distinct lambda set, certified-close subset (deg R < k), and counts.
    """
    H = domain_subgroup(F, s).tolist()
    all_vals = set()
    predicted = set()
    ntuples = 0
    for combo in itertools.combinations(H, r):
        ntuples += 1
        deg, lam = _residual_degree_F(F, combo, m, r)
        all_vals.add(lam)
        if deg < k:
            predicted.add(lam)
    return {
        "all_lambda": sorted(all_vals),
        "predicted_lambda": sorted(predicted),
        "num_tuples": ntuples,
        "num_distinct": len(all_vals),
        "num_predicted": len(predicted),
    }


# ===========================================================================
# THE ALGEBRAIC CRUX (independent of the RS machinery).
# ===========================================================================
def subset_sum_distinctness(F: FiniteField, s: int,
                            full_cap_subsets: int = 2_000_000) -> dict:
    """Are the subset-sums of the order-s subgroup H=<xi> DISTINCT in F?

    This is the heart of N2.  We measure TWO distinctness statistics for H=<xi>:

      (1) r-fold distinct-element sumsets, for each r=2..min(s, R_MAX): how many of
          the C(s,r) sums  sum_{i in T,|T|=r} xi^i  are DISTINCT field values.  The
          construction's certified-close scalars are EXACTLY these (for the design r),
          so a collapse here is a collapse of the counterexample's "many close points".

      (2) FULL distinct-subset-sum property: over ALL 2^s subsets T (incl. empty),
          how many distinct values sum_{i in T} xi^i take.  2^s = 2^s distinct iff H
          has the "distinct subset sums" (Sidon-like) property -- the strongest form.
          (Capped: only computed when 2^s <= full_cap_subsets.)

    For PRIME fields, Kambire's argument (Z-independence of roots via a Linnik prime
    p == 1 mod n) makes these sums distinct asymptotically; the question is whether
    GF(p^e) -- which contains the subfield GF(p) and many F_p-linear relations --
    forces COLLISIONS.  We return raw counts so prime-vs-extension is directly
    comparable.

    NOTE on what "distinct" can possibly mean at small s: the sums live in F (size q),
    so at most q distinct values are possible; we report both the combinatorial count
    and the min(., q) ceiling so the comparison is honest about field-size saturation.
    """
    H = domain_subgroup(F, s).tolist()
    q = F.q

    # (1) r-fold distinct-element sumsets
    R_MAX = min(s, 5)
    rfold = {}
    for r in range(2, R_MAX + 1):
        vals = set()
        ntup = 0
        for combo in itertools.combinations(H, r):
            ntup += 1
            acc = 0
            for xi in combo:
                acc = F.add(acc, int(xi))
            vals.add(acc)
        rfold[r] = {
            "num_tuples": ntup,
            "num_distinct": len(vals),
            "max_possible_combinatorial": ntup,
            "field_ceiling": q,
            "all_distinct": len(vals) == ntup,
        }

    # (2) full distinct-subset-sum property (capped)
    full = None
    if (1 << s) <= full_cap_subsets:
        vals = set()
        # iterate subsets via bitmask; accumulate sums
        for mask in range(1 << s):
            acc = 0
            mm = mask
            idx = 0
            while mm:
                if mm & 1:
                    acc = F.add(acc, int(H[idx]))
                mm >>= 1
                idx += 1
            vals.add(acc)
        full = {
            "num_subsets": 1 << s,
            "num_distinct_sums": len(vals),
            "field_ceiling": q,
            "all_distinct": len(vals) == (1 << s),
            "saturates_field": len(vals) == q,
        }

    return {"field": F.name, "q": q, "s": s, "r_fold": rfold, "full_subset": full}


# ===========================================================================
# Subfield classification + firing-condition detector (the heart of N2).
# ===========================================================================
def subfield_orders(p: int, e: int) -> set:
    """All multiplicative-subgroup orders that live inside a PROPER subfield of
    GF(p^e), i.e. {divisors of p^d - 1 : d | e, d < e}.

    A smooth subgroup of order n is "subfield-contained" (n in this set) iff
    <omega> lies in a proper subfield GF(p^d) <= GF(p^e).  In that case the
    construction on <omega> is SECRETLY the GF(p^d) construction -- it carries NO
    genuinely-extension information.  We use this to separate honest extension tests
    from subfield artifacts.
    """
    out = set()
    for d in range(1, e):
        if e % d == 0:
            for nn in _divisors(p ** d - 1):
                out.add(nn)
    return out


def count_certified_triples(F: FiniteField, s: int, m: int, r: int, k: int) -> int:
    """Number of r-subsets of H=<xi> (order s) whose residual degree deg R < k.

    These are EXACTLY the construction's certified-close r-tuples.  For r=3, k=2m
    needs deg R < 2m; with m=k this is the e2-vanishing condition
    sum_{i<j} xi_i xi_j = 0 (the firing condition).  This counts the firing tuples
    via the SAME in-field residual-degree computation used everywhere else.
    """
    H = domain_subgroup(F, s).tolist()
    cnt = 0
    for combo in itertools.combinations(H, r):
        deg, _lam = _residual_degree_F(F, combo, m, r)
        if deg < k:
            cnt += 1
    return cnt


# ===========================================================================
# Prime control finder: a prime field of comparable size with a smooth subgroup.
# ===========================================================================
def find_prime_control(q_target: int, n: int, s: int, m: int, r: int, k: int,
                       qk_cap=QK_CAP, require_firing=True):
    """Find a PRIME control p' with an order-n subgroup (so p' == 1 mod n), p'^k <=
    qk_cap, size comparable to q_target.  Returns p' or None.

    POSITIVE CONTROL intent: the prime control exists to demonstrate the construction
    DOES fire at this (n,s,r,k) scale.  Whether it fires is prime-dependent (the
    'deg R < k' cancellation is an arithmetic coincidence -- e.g. at n=64,s=32,r=3
    p=193/257/449 fire but p=769/577 do NOT).  So with require_firing we prefer the
    SIZE-CLOSEST prime that actually fires (count_certified_triples > 0); if none in
    the window fires we fall back to the closest prime and flag it.
    """
    lo = max(n + 1, 5)
    hi = max(q_target * 4, n * 16, 1200)
    cands = []
    p = lo
    while p <= hi:
        if is_prime(p) and (p - 1) % n == 0 and (p ** k) <= qk_cap:
            cands.append(p)
        p += 1
    if not cands:
        return None
    if require_firing:
        firing = [pp for pp in cands
                  if count_certified_triples(PrimeField(pp), s, m, r, k) > 0]
        if firing:
            firing.sort(key=lambda pp: (abs(pp - q_target), pp))
            return firing[0]
    cands.sort(key=lambda pp: (abs(pp - q_target), pp))
    return cands[0]


# ===========================================================================
# Build the feasible (field, n, s, m, r, k) cases for an EXTENSION field.
# ===========================================================================
@dataclass
class ExtCase:
    p: int
    e: int
    q: int
    n: int
    s: int
    m: int
    r: int
    k: int
    rho: float
    delta: float
    capacity: float
    eta: float
    a_exp: int   # r*m


def enumerate_ext_cases(p: int, e: int, qk_cap=QK_CAP, smooth_only=True,
                        rho_below_half=True) -> list[ExtCase]:
    """All feasible (n,s,m,r,k) for GF(p^e) with q^k <= qk_cap, smooth n,s.

    Mirrors counterexample_kambire.enumerate_kambire_params' rules:
      n | q-1, s | n, m=n/s, r in [3, s), k=(r-2)m, rho=(r-2)/s < 1/2 (paper range).
    smooth_only requires n AND s to be powers of two (the deployed FFT-domain case,
    and the only regime where the prime construction's smooth signal is defined).
    """
    q = p ** e
    pm1 = q - 1
    out: list[ExtCase] = []
    seen = set()
    for n in _divisors(pm1):
        if n < 4:
            continue
        if smooth_only and (n & (n - 1)) != 0:
            continue
        for s in _divisors(n):
            if s < 2 or s >= n:
                continue
            if smooth_only and (s & (s - 1)) != 0:
                continue
            m = n // s
            for r in range(3, s):
                if rho_below_half and (r - 2) >= s / 2:
                    continue
                k = (r - 2) * m
                if k < 1:
                    continue
                if q ** k > qk_cap:
                    continue
                a_exp = r * m
                if a_exp > n:
                    continue
                rho = (r - 2) / s
                delta = 1.0 - r / s
                cap = 1.0 - rho
                eta = 2.0 / s
                key = (n, s, m, r, k)
                if key in seen:
                    continue
                seen.add(key)
                out.append(ExtCase(p=p, e=e, q=q, n=n, s=s, m=m, r=r, k=k, rho=rho,
                                   delta=delta, capacity=cap, eta=eta, a_exp=a_exp))
    return out


def _kparams_from_extcase(c: ExtCase) -> KParams:
    """Adapt an ExtCase to a KParams so we can reuse delta_window()."""
    return KParams(p=c.q, n=c.n, s=c.s, m=c.m, r=c.r, k=c.k, rho=c.rho,
                   delta=c.delta, capacity=c.capacity, eta=c.eta, a_exp=c.a_exp)


# ===========================================================================
# Run the N2 experiment for one extension field.
# ===========================================================================
def run_field(p: int, e: int, seed=0xC0FFEE, qk_cap=QK_CAP, verbose=True) -> list[dict]:
    rng = np.random.default_rng(seed + p * 1000 + e)
    Fext = PrimePowerField(p, e)
    cases = enumerate_ext_cases(p, e, qk_cap=qk_cap)
    # prefer larger n, then the construction-bearing s (largest sumset)
    cases.sort(key=lambda c: (-c.n, -c.s, c.r))
    sf_orders = subfield_orders(p, e)

    results = []
    for c in cases:
        kp = _kparams_from_extcase(c)
        grid = delta_window(kp)
        exp_f = c.r * c.m
        exp_g = (c.r - 1) * c.m

        D = domain_subgroup(Fext, c.n)

        # Is the INNER subgroup H=<xi> (order s) genuinely-extension, or does it lie
        # in a proper subfield (then the construction is a subfield artifact)?
        s_genuine = c.s not in sf_orders
        n_genuine = c.n not in sf_orders

        # ---- in-field predicted sumset on the EXTENSION smooth subgroup ----
        ss_ext = sumset_H_F(Fext, c.s, c.r, c.m, c.k)
        lambdas_ext = ss_ext["predicted_lambda"]
        ext_certified_triples = count_certified_triples(Fext, c.s, c.m, c.r, c.k)

        # ---- (a) smooth extension subgroup ----
        sub = measure_line_on_domain(Fext, D, c.k, exp_f, exp_g, grid,
                                     predicted_lambdas=lambdas_ext, neg_lambda=True)

        # ---- (b) random size-n subset of GF(p^e)* ----
        Lrand = domain_random(Fext, c.n, rng)
        rnd = measure_line_on_domain(Fext, Lrand, c.k, exp_f, exp_g, grid,
                                     predicted_lambdas=None)

        # ---- (c) prime control: comparable size, own smooth subgroup ----
        prime_summary = None
        pprime = find_prime_control(c.q, c.n, c.s, c.m, c.r, c.k, qk_cap=qk_cap)
        if pprime is not None and (pprime ** c.k) <= qk_cap:
            Fp = PrimeField(pprime)
            Dp = domain_subgroup(Fp, c.n)
            ss_p = sumset_H_F(Fp, c.s, c.r, c.m, c.k)
            lambdas_p = ss_p["predicted_lambda"]
            prm = measure_line_on_domain(Fp, Dp, c.k, exp_f, exp_g, grid,
                                         predicted_lambdas=lambdas_p, neg_lambda=True)
            prime_certified = count_certified_triples(Fp, c.s, c.m, c.r, c.k)
            prime_summary = {
                "prime": pprime, "field": Fp.name,
                "sumset_num_tuples": ss_p["num_tuples"],
                "sumset_num_distinct": ss_p["num_distinct"],
                "predicted_sumset_size": ss_p["num_predicted"],
                "certified_triples": prime_certified,
                "fired": prime_certified > 0,
                **prm,
            }

        # ---- algebraic crux: subset-sum distinctness ext vs prime ----
        crux_ext = subset_sum_distinctness(Fext, c.s)
        crux_prime = (subset_sum_distinctness(PrimeField(pprime), c.s)
                      if pprime is not None else None)

        # design-delta extraction (match prime code's at_design on each capacity)
        def at_design(meas, capacity, target_gap=c.eta):
            target_delta = capacity - target_gap
            return min(meas["per_delta"],
                       key=lambda rc: abs(rc["delta"] - target_delta))

        sub_design = at_design(sub, c.capacity)
        rnd_design = at_design(rnd, c.capacity)
        prime_design = (at_design(prime_summary, c.capacity) if prime_summary else None)

        rec = {
            "field_ext": Fext.name, "p": p, "e": e, "q": c.q,
            "params": asdict(c),
            "exp_f": exp_f, "exp_g": exp_g,
            "design_delta": round(c.delta, 6),
            # subfield classification (the crux): is H=<xi> genuinely extension?
            "s_genuine_extension": bool(s_genuine),
            "n_genuine_extension": bool(n_genuine),
            "subfield_orders": sorted(sf_orders),
            # extension sumset
            "ext_sumset_num_tuples": ss_ext["num_tuples"],
            "ext_sumset_num_distinct": ss_ext["num_distinct"],
            "ext_predicted_sumset_size": ss_ext["num_predicted"],
            "ext_certified_triples": ext_certified_triples,
            "predicted_sumset_max_possible": _ncr(c.s, c.r),
            "subgroup_ext": sub,
            "random_ext": rnd,
            "prime_control": prime_summary,
            "crux_ext": crux_ext,
            "crux_prime": crux_prime,
            "headline": {
                "s_genuine_extension": bool(s_genuine),
                # (a) smooth extension
                "ext_close_at_design": sub_design["close_count"],
                "ext_frac_close": sub_design["frac_close"],
                "ext_Sstar": sub["S_star"],
                "ext_ca_threshold": sub_design["ca_threshold"],
                "ext_is_bad": sub_design["is_bad_line"],
                "ext_meaningful": sub_design["meaningful_regime"],
                "ext_predicted_sumset": ss_ext["num_predicted"],
                "ext_certified_triples": ext_certified_triples,
                "ext_close_eq_pred": (sub_design["close_count"]
                                      == ss_ext["num_predicted"]),
                # (b) random extension
                "rand_close_at_design": rnd_design["close_count"],
                "rand_frac_close": rnd_design["frac_close"],
                "rand_Sstar": rnd["S_star"],
                "rand_is_bad": rnd_design["is_bad_line"],
                # (c) prime control
                "prime": pprime,
                "prime_fired": (prime_summary["fired"] if prime_summary else None),
                "prime_certified_triples": (prime_summary["certified_triples"]
                                            if prime_summary else None),
                "prime_close_at_design": (prime_design["close_count"]
                                          if prime_design else None),
                "prime_frac_close": (prime_design["frac_close"]
                                     if prime_design else None),
                "prime_predicted_sumset": (prime_summary["predicted_sumset_size"]
                                           if prime_summary else None),
                "prime_close_eq_pred": (
                    (prime_design["close_count"]
                     == prime_summary["predicted_sumset_size"])
                    if prime_summary else None),
                "prime_is_bad": (prime_design["is_bad_line"] if prime_design else None),
                # crux headline numbers (design r)
                "ext_rfold_distinct": crux_ext["r_fold"].get(c.r, {}).get("num_distinct"),
                "ext_rfold_tuples": crux_ext["r_fold"].get(c.r, {}).get("num_tuples"),
                "prime_rfold_distinct": (
                    crux_prime["r_fold"].get(c.r, {}).get("num_distinct")
                    if crux_prime else None),
            },
        }
        results.append(rec)
        if verbose:
            h = rec["headline"]
            pc = h["prime_close_at_design"]
            pcs = f"{pc}" if pc is not None else "n/a"
            genu = "GENU" if h["s_genuine_extension"] else "subf"
            print(f"{Fext.name:>8} n={c.n:>3} s={c.s:>3} m={c.m} r={c.r} k={c.k} "
                  f"rho={c.rho:.3f} d={c.delta:.4f} | Hinner={genu} "
                  f"cert3={h.get('ext_certified_triples', 0):>3} | "
                  f"EXT cl={h['ext_close_at_design']:>3}(pred {h['ext_predicted_sumset']:>3},"
                  f"={int(h['ext_close_eq_pred'])}) bad={int(h['ext_is_bad'])} "
                  f"S*={h['ext_Sstar']}/thr{h['ext_ca_threshold']} | "
                  f"RAND cl={h['rand_close_at_design']:>3} | "
                  f"PRIME({h['prime']}) fire={h.get('prime_fired')} cl={pcs}"
                  f"(pred {h['prime_predicted_sumset']})", flush=True)
    return results


# ===========================================================================
# Self-test: the construction identity holds IN-FIELD over an extension field.
# ===========================================================================
def _self_test():
    print("counterexample_extension self-test")
    print("=" * 70)

    # We need an odd-char extension field with a smooth subgroup big enough to host
    # the construction with a meaningful regime and p^k small.  GF(289)=17^2 has
    # |F*|=288=2^5*3^2, smooth subgroup of order n=32, s=16, m=2, r=3, k=2 (rho=1/16):
    # q^k = 289^2 = 83521 <= cap.  This MIRRORS the prime cases p=193/257 n=64? No;
    # we choose what's feasible: the r=3,k=2 line is the clean prime-signal analogue.
    F = PrimePowerField(17, 2)            # GF(289)
    n, s, m, r, k = 32, 16, 2, 3, 2
    assert (F.q - 1) % n == 0
    D = domain_subgroup(F, n)
    H = domain_subgroup(F, s)
    assert set(H.tolist()).issubset(set(D.tolist())), "H must be subgroup of D"

    exp_f, exp_g = r * m, (r - 1) * m
    f = monomial_eval(F, D, exp_f)
    g = monomial_eval(F, D, exp_g)
    book = build_codeword_book(F, D, k)
    assert min_distance(book) == n - k + 1, "RS over GF(289) must be MDS"
    print(f"  GF(289) n={n} k={k}: RS is MDS (d={n-k+1}) OK")

    # coset-union sanity (purely multiplicative; field-agnostic):
    H_list = H.tolist()
    combo0 = tuple(H_list[:r])
    union = []
    for xi in combo0:
        union.extend([int(a) for a in D.tolist() if F.pow(int(a), m) == xi])
    assert len(union) == r * m, f"coset union {len(union)} != r*m {r*m}"
    print(f"  coset structure: {r} cosets x m={m} -> union {len(union)} = r*m OK")

    # In-field identity: find a certified-close lambda (deg R < k) and verify the
    # word f - lambda*g is delta-close EXACTLY.
    ss = sumset_H_F(F, s, r, m, k)
    print(f"  in-field sumset: {ss['num_tuples']} tuples -> {ss['num_distinct']} "
          f"distinct lambda -> {ss['num_predicted']} certified-close (deg R < k)")
    delta_n = (1.0 - r / s) * n
    if ss["num_predicted"] >= 1:
        witness = ss["predicted_lambda"][0]
        gamma = F.neg(witness)
        w = F.add_vec(f, F.mul_scalar_vec(gamma, g))
        d = dist_to_code(book, w)
        assert d <= delta_n + 1e-9, \
            f"certified lambda={witness} dist {d} should be <= delta*n {delta_n}"
        print(f"  identity check: certified lambda={witness} gives dist {d} "
              f"<= delta*n={delta_n:.0f} OK")
    else:
        print("  (no certified-close lambda at this small case -- identity vacuous "
              "here; the full run uses larger s)")

    # No-CA: S* < CA threshold and bounded by deg(g).
    S_star, info = _max_common_agreement(book, f, g)
    assert info["exact"], "S* must be exact"
    pop_g, _ = _agreement_bits(book, g)
    assert int(pop_g.max()) <= exp_g, \
        f"g=X^{exp_g} cannot agree on more than deg(g)={exp_g} pts"
    ca_threshold = r * m
    assert S_star < ca_threshold, f"S*={S_star} must be < CA threshold {ca_threshold}"
    print(f"  no-CA check: S*={S_star} < CA threshold={ca_threshold}; "
          f"g max single-agreement={int(pop_g.max())} <= deg(g)={exp_g} OK")

    # Subset-sum crux runs and returns sane structure.
    crux = subset_sum_distinctness(F, s)
    rf = crux["r_fold"][r]
    print(f"  crux: GF(289) H=<xi> order {s}: {r}-fold sumset "
          f"{rf['num_distinct']} distinct / {rf['num_tuples']} tuples "
          f"(field ceiling {crux['q']})")
    assert rf["num_distinct"] >= 1
    if crux["full_subset"] is not None:
        fs = crux["full_subset"]
        print(f"        full subset-sums: {fs['num_distinct_sums']} distinct / "
              f"{fs['num_subsets']} subsets; all_distinct={fs['all_distinct']} "
              f"saturates_field={fs['saturates_field']}")

    # Prime control finder returns a usable prime (prefers one that FIRES).
    pc = find_prime_control(F.q, n, s, m, r, k)
    assert pc is not None and (pc - 1) % n == 0, "prime control must have order-n subgrp"
    fired = count_certified_triples(PrimeField(pc), s, m, r, k) > 0
    print(f"  prime control for q~{F.q}, n={n}: p'={pc} (p'-1 divisible by n; "
          f"fires={fired}) OK")

    print("=" * 70)
    print("ALL counterexample_extension SELF-TESTS PASSED")


# ===========================================================================
# Output writers + cross-field summary.
# ===========================================================================
def write_outputs(all_results, elapsed, out_dir=RESULTS_DIR):
    os.makedirs(out_dir, exist_ok=True)
    payload = {
        "meta": {
            "experiment": "N2_extension_field_counterexample",
            "question": ("Does the smooth-domain Kambire/BCHKS near-capacity "
                         "proximity-gap counterexample extend from PRIME fields to "
                         "ODD-CHARACTERISTIC EXTENSION fields GF(p^e)?"),
            "timestamp": time.strftime("%Y%m%d_%H%M%S"),
            "elapsed_sec": round(elapsed, 2),
            "qk_cap": QK_CAP,
            "convention": ("delta = radius. Line f + gamma*g, gamma over ALL of F. "
                           "EXACT distance via full codeword enumeration; S* exact "
                           "branch-and-bound. delta_design = capacity - 2/s = 1 - r/s."),
            "construction": ("f=X^{rm}, g=X^{(r-1)m} on D=<omega> (order n); H=<xi> "
                             "order s; k=(r-2)m; predicted close scalars = H^{(+r)} "
                             "(distinct-element r-fold sumset of H), computed IN-FIELD."),
            "signal": ("PRIME smooth subgroup shows close_count == predicted_sumset_size "
                       ">> random. We test whether odd-char EXTENSION smooth subgroups "
                       "reproduce this."),
        },
        "results": all_results,
    }
    jpath = os.path.join(out_dir, "counterexample_extension.json")
    with open(jpath, "w", newline="\n") as fh:
        json.dump(payload, fh, indent=2, default=str)

    # Flat CSV: one row per (ext field, case, domain, delta).
    cpath = os.path.join(out_dir, "counterexample_extension.csv")
    rows = []
    for rec in all_results:
        c = rec["params"]
        domains = [("subgroup_ext", "smooth_ext_subgroup", rec.get("subgroup_ext")),
                   ("random_ext", "random_ext_subset", rec.get("random_ext")),
                   ("prime_control", "prime_control_subgroup", rec.get("prime_control"))]
        for _key, label, meas in domains:
            if meas is None:
                continue
            # capacity is the same (matched n,k) for ext domains; for prime control too
            cap = c["capacity"]
            pred = (rec["ext_predicted_sumset_size"] if label != "prime_control_subgroup"
                    else meas.get("predicted_sumset_size"))
            for pd in meas["per_delta"]:
                rows.append({
                    "ext_field": rec["field_ext"], "p": rec["p"], "e": rec["e"],
                    "q": rec["q"],
                    "domain": label, "domain_field": meas["field"],
                    "s_genuine_extension": rec["s_genuine_extension"],
                    "n": meas["n"], "k": meas["k"],
                    "s": c["s"], "m": c["m"], "r": c["r"],
                    "rho": round(c["rho"], 4),
                    "capacity": round(cap, 4),
                    "design_delta": round(c["delta"], 4),
                    "eta_2_over_s": round(c["eta"], 4),
                    "exp_f": meas["exp_f"], "exp_g": meas["exp_g"],
                    "delta": pd["delta"],
                    "delta_minus_capacity": round(pd["delta"] - cap, 4),
                    "close_count": pd["close_count"],
                    "frac_close": pd["frac_close"],
                    "S_star": pd["S_star"],
                    "ca_threshold": pd["ca_threshold"],
                    "g_max_single_agreement": meas["g_max_single_agreement"],
                    "meaningful_regime": pd["meaningful_regime"],
                    "is_bad_line": pd["is_bad_line"],
                    "predicted_sumset_size": pred,
                })
    if rows:
        keys = list(rows[0].keys())
        with open(cpath, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, lineterminator="\n")
            w.writeheader()
            for row in rows:
                w.writerow(row)
    return jpath, cpath, len(rows)


def summarize(all_results):
    """The decisive comparison, SPLIT by whether the inner subgroup H=<xi> is
    genuinely-extension or a subfield artifact.

    A smooth subgroup of order s lying in a proper subfield GF(p^d) makes the
    construction LITERALLY the prime-power-subfield construction (no new info).  The
    honest N2 question is only about GENUINELY-EXTENSION inner subgroups.  We tabulate
    both, and the verdict counts firing ONLY among genuine-extension cases.
    """
    print("\n" + "=" * 104)
    print("DECISIVE COMPARISON (meaningful regime, design delta = capacity - 2/s = 1 - r/s)")
    print("  Q: does the smooth signal (close_count == predicted_sumset > 0, > random) "
          "fire on GENUINE-extension subgroups?")
    print("=" * 104)
    hdr = (f"  {'extfield':>8} {'n':>3} {'s':>3} {'r':>2} {'k':>2} {'Hgenu':>5} "
           f"{'cert3':>5} | {'EXTcl':>5} {'pred':>4} {'EXTbad':>6} | "
           f"{'RANDcl':>6} | {'prime':>5} {'Pfire':>5} {'Pcl':>4} {'Ppred':>5} {'Pbad':>4}")
    print(hdr)
    n_meaningful = 0
    genuine_total = 0
    genuine_fires = 0
    subfield_total = 0
    subfield_fires = 0
    prime_fires = 0
    prime_available = 0
    collapse_rows = []      # genuine-ext cases where prime fired but ext did not
    for rec in all_results:
        h = rec["headline"]
        c = rec["params"]
        if not h["ext_meaningful"]:
            continue
        n_meaningful += 1
        genu = h["s_genuine_extension"]
        ext_signal = (h["ext_close_eq_pred"] and h["ext_predicted_sumset"] > 0
                      and h["ext_close_at_design"] > h["rand_close_at_design"])
        prime_signal = bool(h.get("prime_fired")) and (h["prime_close_at_design"] or 0) > 0
        if genu:
            genuine_total += 1
            genuine_fires += int(ext_signal)
        else:
            subfield_total += 1
            subfield_fires += int(ext_signal)
        if h["prime"] is not None:
            prime_available += 1
            prime_fires += int(prime_signal)
        if genu and prime_signal and not ext_signal:
            collapse_rows.append(rec)
        print(f"  {rec['field_ext']:>8} {c['n']:>3} {c['s']:>3} {c['r']:>2} {c['k']:>2} "
              f"{str(genu):>5} {h.get('ext_certified_triples', 0):>5} | "
              f"{h['ext_close_at_design']:>5} {h['ext_predicted_sumset']:>4} "
              f"{str(h['ext_is_bad']):>6} | {h['rand_close_at_design']:>6} | "
              f"{str(h['prime']):>5} {str(h.get('prime_fired')):>5} "
              f"{str(h['prime_close_at_design']):>4} {str(h['prime_predicted_sumset']):>5} "
              f"{str(h['prime_is_bad']):>4}")
    print("-" * 104)
    print(f"  meaningful-regime cases: {n_meaningful}  "
          f"(genuine-ext inner H: {genuine_total}, subfield-contained H: {subfield_total})")
    print(f"  GENUINE-extension smooth signal fired (close==pred>0 AND > random): "
          f"{genuine_fires}/{genuine_total}")
    print(f"  SUBFIELD-artifact smooth signal fired (= prime/subfield in disguise): "
          f"{subfield_fires}/{subfield_total}")
    print(f"  PRIME control fired (positive control):                             "
          f"{prime_fires}/{prime_available}")
    print("-" * 104)
    if genuine_total == 0:
        print("  VERDICT: no genuine-extension MEANINGFUL case at these scales "
              "(all firing cases are subfield artifacts).")
    elif genuine_fires == 0:
        print("  VERDICT: the counterexample DID NOT FIRE on any genuine-extension "
              "smooth subgroup at these scales.")
        print("           Where the size-matched PRIME control DID fire, the genuine "
              "extension showed close_count = predicted = 0:")
        for rec in collapse_rows:
            h = rec["headline"]; c = rec["params"]
            print(f"             {rec['field_ext']} n={c['n']} s={c['s']} r={c['r']} "
                  f"k={c['k']}: EXT cert-triples={h.get('ext_certified_triples')} "
                  f"close={h['ext_close_at_design']} | PRIME({h['prime']}) fired, "
                  f"cert-triples={h.get('prime_certified_triples')} "
                  f"close={h['prime_close_at_design']}")
    else:
        print(f"  VERDICT: the counterexample FIRED on {genuine_fires} genuine-extension "
              f"smooth subgroup(s) -- the negative EXTENDS (opening closes).")
    return {"n_meaningful": n_meaningful, "genuine_total": genuine_total,
            "genuine_fires": genuine_fires, "subfield_total": subfield_total,
            "subfield_fires": subfield_fires, "prime_fires": prime_fires,
            "prime_available": prime_available, "collapse_rows": collapse_rows}


def run(fields=((3, 2), (5, 2), (7, 2), (3, 4), (11, 2), (13, 2), (17, 2),
                (19, 2), (23, 2), (31, 2)),
        seed=0xC0FFEE, qk_cap=QK_CAP, verbose=True):
    all_results = []
    t0 = time.time()
    for (p, e) in fields:
        if verbose:
            print(f"\n--- GF({p}^{e}) = GF({p**e}), |F*|={p**e-1} ---")
        all_results.extend(run_field(p, e, seed=seed, qk_cap=qk_cap, verbose=verbose))
    elapsed = time.time() - t0
    return all_results, elapsed


if __name__ == "__main__":
    import sys
    args = [a for a in sys.argv[1:] if not a.startswith("--out=")]
    out_dir = next((a.split("=", 1)[1] for a in sys.argv[1:]
                    if a.startswith("--out=")),
                   os.environ.get("EXT_OUT", RESULTS_DIR))
    _self_test()
    if args and args[0] == "test":
        sys.exit(0)
    print()
    print("=" * 70)
    print("RUNNING N2 EXTENSION EXPERIMENT")
    print("=" * 70)
    results, elapsed = run()
    summ = summarize(results)
    jpath, cpath, nrows = write_outputs(results, elapsed, out_dir=out_dir)
    print("-" * 70)
    print(f"Done in {elapsed:.1f}s. Wrote:\n  {jpath}\n  {cpath} ({nrows} rows)")
