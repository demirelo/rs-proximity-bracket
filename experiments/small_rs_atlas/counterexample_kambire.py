"""
counterexample_kambire.py -- INSTANTIATE the Kambiré / BCHKS near-capacity
proximity-gap counterexample at small scale, EXACTLY.

Background (Kambiré arXiv:2604.09724; BCHKS ePrint 2025/2055 Thm 1.13)
---------------------------------------------------------------------
Over a PRIME field F_p with a smooth multiplicative-subgroup domain D = <omega>
of order n, take a subgroup H = <xi> of order s | n, with m = n/s cosets.  Set

    r = rho*s + 2,     k = (r-2)*m,     rho = (r-2)/s,
    f = X^{r*m},       g = X^{(r-1)*m}        (monomial evaluations on D),
    delta = 1 - r/s,   capacity = 1 - rho,    eta = capacity - delta = 2/s.

The ENGINE is the polynomial identity (on the subgroup, X^a depends only on a mod n):
    prod_{j=1..r} (X^m - xi_j)  =  X^{rm} - (xi_1+...+xi_r) X^{(r-1)m} + R(X),
with deg R <= (r-2)m = k.  The r cosets  H_j = { a in D : a^m = xi_j }  union to a
set of size r*m on which the left side vanishes.  Hence for every
    lambda = xi_1 + ... + xi_r   (xi_j distinct in H; the "distinct-element sumset" H^{(+r)})
the word  f + (-lambda)*g  = X^{rm} - lambda X^{(r-1)m}  equals  -R(X)  on those r*m
points, i.e. it agrees with a degree-<=k polynomial on (1-delta)*n = r*m coordinates,
so  Delta(f + (-lambda) g, C) <= delta.  There are |H^{(+r)}| such scalars
(>= (s/2r)^r asymptotically; n^C for large s).  THIS IS MANY CLOSE POINTS.

Yet there is NO correlated agreement: g = X^{(r-1)m} has degree (r-1)m > k, so it
agrees with ANY single degree-<k codeword on at most (r-1)m points (a nonzero
polynomial of degree (r-1)m has at most that many roots), and (r-1)m < r*m = the
CA threshold.  So the joint common-agreement set S*(f,g) <= (r-1)m < (1-delta)n:
no single common set of density 1-delta explains the line.  ==> CA / proximity-gap
FAILS at delta = capacity - 2/s.

What this script does (EXACTLY -- full codeword enumeration, no decoder)
------------------------------------------------------------------------
For each feasible (p, n, s, m, r, k) with p^k <= a cap:
  * Build C = RS[F_p, D, k] on the smooth subgroup D = <omega>.
  * Construct f = X^{rm}, g = X^{(r-1)m}.
  * Enumerate the predicted sumset H^{(+r)} and check how many predicted lambda are
    actually delta-close (small fields lose some to collisions mod p).
  * Compute the GLOBAL close-count over ALL gamma in F_p (the honest line statistic).
  * Compute the EXACT S* = max common agreement, and classify the line: bad iff
    many close points AND S* < ceil((1-delta) n) in the meaningful regime.
  * Sweep delta around capacity (capacity - {a few}/s ... capacity) to locate the
    onset/width of the failure window and report close-count(delta), S*, CA-threshold.

Cross-checks for PRIME-SPECIFICITY:
  (a) RANDOM domain: same field, same n, same monomial exponents f=X^{rm}, g=X^{(r-1)m}
      evaluated on a RANDOM size-n subset of F_p^* (no subgroup structure).  The
      coset/sumset identity breaks, so the prediction is the bad line should
      LARGELY VANISH (few close points) -- the badness is structure-specific.
  (b) BINARY-EXTENSION subgroup of comparable size: same monomial line on a
      GF(2^m) multiplicative subgroup of order n' ~ n (where one exists).  Tests
      whether extension fields reproduce the prime-smooth badness.

All distances and S* are EXACT.  Results -> results/counterexample_kambire.{json,csv}.
"""

from __future__ import annotations

import csv
import json
import itertools
import os
import time
from dataclasses import dataclass, asdict

import numpy as np

from ff import PrimeField, BinaryExtensionField, FiniteField, _divisors, is_prime
from rs import (build_codeword_book, domain_subgroup, domain_random, dist_to_code,
                min_distance, CodewordBook)
from search_bad_lines import _max_common_agreement, _agreement_bits


RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")

# Cap on q^k for exact full enumeration of the codeword book (per the brief, ~3e6).
QK_CAP = 3_000_000


# ---------------------------------------------------------------------------
# Monomial evaluation on a domain (X^a on the elements of L).
# ---------------------------------------------------------------------------
def monomial_eval(F: FiniteField, L: np.ndarray, a: int) -> np.ndarray:
    """Evaluate X^a on every element of L (field powers)."""
    return np.array([F.pow(int(x), a) for x in L], dtype=np.int64)


# ---------------------------------------------------------------------------
# Parameter enumeration.
# ---------------------------------------------------------------------------
@dataclass
class KParams:
    p: int
    n: int          # |D| = order of the domain subgroup
    s: int          # |H| = order of the inner subgroup (s | n)
    m: int          # number of cosets = n / s
    r: int          # r = rho*s + 2
    k: int          # code dimension = (r-2)*m
    rho: float
    delta: float    # 1 - r/s  (the predicted failure radius)
    capacity: float
    eta: float      # capacity - delta = 2/s
    a_exp: int      # r*m = number of agreement coordinates = (1-delta)*n

    @property
    def qk(self) -> int:
        return self.p ** self.k


def enumerate_kambire_params(primes, smooth_only=True, qk_cap=QK_CAP) -> list[KParams]:
    """All feasible (p,n,s,m,r,k) for the construction with p^k <= qk_cap.

    smooth_only: require n AND s to be powers of two (deployed FFT-domain case).
    """
    out: list[KParams] = []
    seen = set()
    for p in primes:
        if not is_prime(p):
            continue
        pm1 = p - 1
        for n in _divisors(pm1):
            if n < 4:
                continue
            if smooth_only and (n & (n - 1)) != 0:
                continue                      # n not a power of two
            for s in _divisors(n):
                if s < 2 or s >= n:
                    continue
                if smooth_only and (s & (s - 1)) != 0:
                    continue                  # s not a power of two
                m = n // s
                for r in range(3, s):
                    if (r - 2) >= s / 2:       # rho < 1/2 strict (paper's range)
                        continue
                    k = (r - 2) * m
                    if k < 1:
                        continue
                    if p ** k > qk_cap:
                        continue
                    rho = (r - 2) / s
                    delta = 1.0 - r / s
                    cap = 1.0 - rho
                    eta = 2.0 / s
                    a_exp = r * m
                    if a_exp > n:
                        continue
                    key = (p, n, s, m, r, k)
                    if key in seen:
                        continue
                    seen.add(key)
                    out.append(KParams(p=p, n=n, s=s, m=m, r=r, k=k, rho=rho,
                                       delta=delta, capacity=cap, eta=eta,
                                       a_exp=a_exp))
    return out


# ---------------------------------------------------------------------------
# Core exact measurement on a single domain (subgroup / random / extension).
# ---------------------------------------------------------------------------
def measure_line_on_domain(F: FiniteField, L: np.ndarray, k: int,
                           exp_f: int, exp_g: int,
                           delta_grid: list[float],
                           predicted_lambdas: list[int] | None = None,
                           neg_lambda: bool = True) -> dict:
    """Build C=RS[F,L,k], the monomial line f=X^{exp_f}, g=X^{exp_g}, and measure
    EXACTLY over the whole field of scalars gamma:

      * per-delta global close-count #{gamma in F : dist(f+gamma g, C) <= delta n}
      * exact S* = max common agreement of (f,g)
      * (optional) how many of the PREDICTED lambda values are actually close
        (a predicted scalar contributes at gamma = -lambda when neg_lambda, matching
        the identity f - lambda g; we test that exact gamma).

    Returns a dict with all the numbers and the per-delta classification.
    """
    n = len(L)
    book = build_codeword_book(F, L, k)
    f = monomial_eval(F, L, exp_f)
    g = monomial_eval(F, L, exp_g)

    # Exact distance for EVERY gamma in F (q distance computations).
    dists = np.empty(F.q, dtype=np.int64)
    for gg in range(F.q):
        w = F.add_vec(f, F.mul_scalar_vec(gg, g))
        dists[gg] = dist_to_code(book, w)

    # Exact common agreement S* and the per-codeword max agreement of g alone
    # (the quantity that upper-bounds S* and certifies no-CA).
    S_star, info = _max_common_agreement(book, f, g)
    pop_g, _ = _agreement_bits(book, g)
    pop_f, _ = _agreement_bits(book, f)
    g_max_single = int(pop_g.max())
    f_max_single = int(pop_f.max())

    # Per-delta classification.
    per_delta = []
    for delta in delta_grid:
        thresh = delta * n + 1e-9
        close = int((dists <= thresh).sum())
        ca_threshold = int(np.ceil((1.0 - delta) * n - 1e-9))
        meaningful = ca_threshold > k
        ca_explained = S_star >= ca_threshold
        is_bad = meaningful and (not ca_explained) and (close > 1)
        rec = {
            "delta": round(float(delta), 6),
            "close_count": close,
            "frac_close": round(close / F.q, 6),
            "ca_threshold": ca_threshold,
            "S_star": int(S_star),
            "meaningful_regime": bool(meaningful),
            "ca_explained": bool(ca_explained),
            "is_bad_line": bool(is_bad),
        }
        per_delta.append(rec)

    # Predicted-lambda hit rate (only meaningful for the subgroup construction).
    pred = None
    if predicted_lambdas is not None:
        npred = len(predicted_lambdas)
        # use the construction's design delta = the first grid point that equals
        # 1 - exp_f-related? We just report distances at the predicted gammas.
        pred_dists = []
        for lam in predicted_lambdas:
            gamma = F.neg(lam) if neg_lambda else lam
            w = F.add_vec(f, F.mul_scalar_vec(int(gamma), g))
            pred_dists.append(int(dist_to_code(book, w)))
        pred_dists = np.array(pred_dists, dtype=np.int64)
        pred = {
            "num_predicted_lambda": npred,
            "pred_dist_min": int(pred_dists.min()) if npred else -1,
            "pred_dist_max": int(pred_dists.max()) if npred else -1,
            "pred_dist_mean": float(pred_dists.mean()) if npred else -1.0,
            "pred_dists": pred_dists.tolist(),
        }

    return {
        "q": F.q, "field": F.name, "n": n, "k": k,
        "exp_f": exp_f, "exp_g": exp_g,
        "num_codewords": book.num_codewords,
        "min_dist": int(min_distance(book)),
        "mds_ok": int(min_distance(book)) == (n - k + 1),
        "S_star": int(S_star), "S_exact": bool(info["exact"]),
        "g_max_single_agreement": g_max_single,
        "f_max_single_agreement": f_max_single,
        "dist_min_over_line": int(dists.min()),
        "dist_max_over_line": int(dists.max()),
        "per_delta": per_delta,
        "predicted": pred,
    }


def _polymul_mod(a, b, p):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] = (out[i + j] + ai * bj) % p
    return out


def _residual_degree(F: PrimeField, combo, m: int, r: int) -> int:
    """Degree of R(X) = prod_j (X^m - xi_j) - (X^{rm} - lambda X^{(r-1)m}).

    The Kambiré identity guarantees the X^{rm} and X^{(r-1)m} terms cancel, leaving
    R of degree <= (r-2)m.  A scalar lambda is a *theoretically-predicted close*
    point iff there is a representation with deg R < k = (r-2)m (then f - lambda g
    equals the degree-<k polynomial -R on the size-rm coset union, so it is
    delta-close).  At small m, deg R can land exactly at (r-2)m = k (not < k), in
    which case THAT representation does not certify closeness.
    """
    p = F.q
    poly = [1]
    lam = 0
    for xi in combo:
        lam = F.add(lam, int(xi))
        fac = [0] * (m + 1)
        fac[0] = (-int(xi)) % p
        fac[m] = 1
        poly = _polymul_mod(poly, fac, p)
    R = poly[:]
    R[r * m] = (R[r * m] - 1) % p
    R[(r - 1) * m] = (R[(r - 1) * m] + lam) % p
    deg = max((i for i, c in enumerate(R) if c % p != 0), default=-1)
    return deg, lam


def sumset_H(F: PrimeField, s: int, r: int, m: int, k: int) -> dict:
    """Analyse the distinct-element r-fold sumset H^{(+r)} of H = <xi> (order s).

    Returns dict with:
      * 'all_lambda'        : sorted DISTINCT field values lambda = xi_1+...+xi_r,
      * 'predicted_lambda'  : the subset of lambda admitting a representation with
                              deg R < k (the construction's certified-close scalars),
      * 'num_tuples', 'num_distinct', 'num_predicted'.

    Small fields collapse many tuples to the same lambda (collisions mod p) -- the
    phenomenon the asymptotic construction sidesteps via a Linnik prime choice.  At
    small m, some lambda fail the deg R < k test (residual degree lands at k); we
    separate the certified-close subset from the raw sumset honestly.
    """
    H = domain_subgroup(F, s).tolist()
    all_vals = set()
    predicted = set()
    ntuples = 0
    for combo in itertools.combinations(H, r):
        ntuples += 1
        deg, lam = _residual_degree(F, combo, m, r)
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


def delta_window(kp: KParams, n_below=4, n_above=1) -> list[float]:
    """Delta grid straddling the predicted failure radius and capacity.

    The construction predicts failure AT delta = capacity - 2/s = kp.delta.
    We also probe a few multiples of (1/s) below the design delta and up to
    capacity, all snapped to the 1/n distance lattice so close-counts are
    integer-clean.
    """
    n, s = kp.n, kp.s
    cap = kp.capacity
    pts = set()
    # the exact design radius and capacity
    pts.add(round(kp.delta, 6))
    pts.add(round(cap, 6))
    # a few steps of 1/n around the design delta toward capacity and below it
    for j in range(-n_below, n_above + 1):
        d = kp.delta + j / n
        if 0.02 < d < 0.999:
            pts.add(round(d, 6))
    # also half-way between Johnson and design delta, as a "well-below" anchor
    johnson = 1.0 - np.sqrt(kp.rho)
    mid = 0.5 * (johnson + kp.delta)
    if 0.02 < mid < 0.999:
        pts.add(round(mid, 6))
    return sorted(pts)


# ---------------------------------------------------------------------------
# Companion extension-field subgroup of comparable size.
# ---------------------------------------------------------------------------
def find_extension_companion(n_target: int, k: int, qk_cap=QK_CAP):
    """Find a GF(2^m) with a multiplicative subgroup of order n' close to n_target
    such that (2^m)^k <= qk_cap.  Returns (F, n_prime) or (None, None).

    GF(2^m)* has order 2^m - 1 (ODD), so its subgroup orders are all odd -- there is
    NO power-of-two-order (smooth) subgroup in a binary field (the very reason the
    deployed smooth case lives over prime / Mersenne-prime fields).  We therefore
    take the closest available subgroup order as an *odd-order* comparison domain.
    """
    best = None
    for m in (6, 7, 8):
        F = BinaryExtensionField(m)
        if F.q ** k > qk_cap:
            continue
        for nprime in _divisors(F.q - 1):
            if nprime < max(4, k + 1):
                continue
            if nprime > F.q - 1:
                continue
            cand = (abs(nprime - n_target), F, nprime, m)
            if best is None or cand[0] < best[0]:
                best = cand
    if best is None:
        return None, None
    return best[1], best[2]


# ---------------------------------------------------------------------------
# Driver: run the full instantiation experiment.
# ---------------------------------------------------------------------------
def run(primes=(17, 97, 193, 257), seed=0xCAFE, qk_cap=QK_CAP, verbose=True):
    rng = np.random.default_rng(seed)
    params = enumerate_kambire_params(primes, smooth_only=True, qk_cap=qk_cap)
    # Prefer the cases that actually have a meaningful regime and several cosets /
    # a non-trivial sumset.  Keep all, but sort so the headline (large n, smallest
    # eta) come first.
    params.sort(key=lambda kp: (kp.p, -kp.n, kp.s))

    results = []
    t0 = time.time()
    for kp in params:
        F = PrimeField(kp.p)
        D = domain_subgroup(F, kp.n)
        exp_f = kp.r * kp.m
        exp_g = (kp.r - 1) * kp.m

        # predicted sumset (only for the prime subgroup construction)
        ss = sumset_H(F, kp.s, kp.r, kp.m, kp.k)
        lambdas = ss["predicted_lambda"]
        grid = delta_window(kp)

        # (1) the smooth prime subgroup -- the construction's native domain
        sub = measure_line_on_domain(F, D, kp.k, exp_f, exp_g, grid,
                                      predicted_lambdas=lambdas, neg_lambda=True)

        # (2) RANDOM domain of the same size n over the same field (no structure)
        Lrand = domain_random(F, kp.n, rng)
        rnd = measure_line_on_domain(F, Lrand, kp.k, exp_f, exp_g, grid,
                                     predicted_lambdas=None)

        # (3) extension-field companion subgroup of comparable size.
        # GF(2^m)* has ODD order (2^m - 1), so it has NO power-of-two-order subgroup
        # and the literal coset/sumset construction cannot run.  The honest probe is:
        # take an extension subgroup of comparable n', a code of the SAME rate rho
        # (k' = round(rho*n')), and the analogous high-degree-monomial line
        # f=X^{a'}, g=X^{b'} with b' = (r-1)/r * a' scaled to n' (so g is a monomial
        # of degree > k', the no-CA engine), then measure at the SAME RELATIVE delta
        # window anchored on the EXTENSION's own capacity (1 - k'/n').  If badness is
        # prime-smooth-specific, the extension should NOT show the inflated close-count
        # that the smooth prime subgroup does at matched (rho, delta-from-capacity).
        ext_summary = None
        Fext, nprime = find_extension_companion(kp.n, kp.k, qk_cap=qk_cap)
        if Fext is not None and nprime is not None:
            kprime = max(1, int(round(kp.rho * nprime)))
            if (Fext.q ** kprime) <= qk_cap and kprime < nprime:
                Dext = domain_subgroup(Fext, nprime)
                # analogous monomials: f deg ~ r/(r-2) * k' (just-out high monomial),
                # g deg ~ (r-1)/(r-2) * k' -- both > k', g the no-CA driver.  Scale the
                # prime's exponents by n'/n and reduce mod n'.
                ef = (exp_f * nprime // kp.n) % nprime
                eg = (exp_g * nprime // kp.n) % nprime
                # ensure g is a genuine high-degree monomial (deg in (k', n')); if the
                # scaling collapses it, fall back to the smallest monomial above k'.
                if not (kprime < eg < nprime):
                    eg = (kprime + 1) % nprime if (kprime + 1) < nprime else kprime
                if not (kprime < ef < nprime):
                    ef = min(nprime - 1, kprime + 2)
                cap_ext = 1.0 - kprime / nprime
                # relative-delta window matched to the prime's (capacity - delta) gaps
                gaps = sorted({round(kp.capacity - d, 6) for d in grid})
                ext_grid = sorted({round(cap_ext - gp, 6) for gp in gaps
                                   if 0.02 < cap_ext - gp < 0.999})
                if not ext_grid:
                    ext_grid = [round(cap_ext - 1.0 / nprime, 6)]
                ext = measure_line_on_domain(Fext, Dext, kprime, ef, eg, ext_grid,
                                             predicted_lambdas=None)
                ext_summary = {"field": Fext.name, "n_prime": nprime,
                               "k_prime": kprime, "rho_ext": kprime / nprime,
                               "capacity_ext": cap_ext,
                               "exp_f_ext": ef, "exp_g_ext": eg, **ext}

        # find the close-count / S* at the DESIGN delta for the headline.  For prime
        # domains (subgroup, random) the design delta is kp.delta; for the extension
        # we match on the GAP BELOW CAPACITY (eta = 2/s) so the comparison is at the
        # same relative radius on each domain's own capacity.
        def at_design(meas, capacity=kp.capacity, target_gap=kp.eta):
            target_delta = capacity - target_gap
            best = min(meas["per_delta"],
                       key=lambda rc: abs(rc["delta"] - target_delta))
            return best

        sub_design = at_design(sub)
        rnd_design = at_design(rnd)

        rec = {
            "params": asdict(kp),
            "exp_f": exp_f, "exp_g": exp_g,
            "design_delta": round(kp.delta, 6),
            "sumset_num_tuples": ss["num_tuples"],
            "sumset_num_distinct": ss["num_distinct"],
            "predicted_sumset_size": ss["num_predicted"],
            "predicted_sumset_max_possible": _ncr(kp.s, kp.r),
            "subgroup": sub,
            "random": rnd,
            "extension": ext_summary,
            "headline": {
                "smooth_close_at_design": sub_design["close_count"],
                "smooth_Sstar": sub["S_star"],
                "smooth_ca_threshold": sub_design["ca_threshold"],
                "smooth_is_bad": sub_design["is_bad_line"],
                "smooth_meaningful": sub_design["meaningful_regime"],
                "random_close_at_design": rnd_design["close_count"],
                "random_Sstar": rnd["S_star"],
                "random_is_bad": rnd_design["is_bad_line"],
                "extension_close_at_design": (
                    at_design(ext_summary, ext_summary["capacity_ext"], kp.eta)["close_count"]
                    if ext_summary else None),
                "extension_frac_close": (
                    at_design(ext_summary, ext_summary["capacity_ext"], kp.eta)["frac_close"]
                    if ext_summary else None),
                "extension_is_bad": (
                    at_design(ext_summary, ext_summary["capacity_ext"], kp.eta)["is_bad_line"]
                    if ext_summary else None),
                "smooth_frac_close": sub_design["frac_close"],
                "random_frac_close": rnd_design["frac_close"],
            },
        }
        results.append(rec)
        if verbose:
            h = rec["headline"]
            extc = h['extension_close_at_design']
            extf = h['extension_frac_close']
            extstr = f"{extc}({extf:.3f})" if extc is not None else "n/a"
            print(f"p={kp.p:>4} n={kp.n:>3} s={kp.s:>3} m={kp.m} r={kp.r} k={kp.k} "
                  f"rho={kp.rho:.3f} d={kp.delta:.4f} cap={kp.capacity:.3f} "
                  f"eta={kp.eta:.4f} | "
                  f"SMOOTH cl={h['smooth_close_at_design']:>3}({h['smooth_frac_close']:.3f}) "
                  f"S*={h['smooth_Sstar']:>2}/thr{h['smooth_ca_threshold']:>2} "
                  f"bad={int(h['smooth_is_bad'])} "
                  f"| RAND cl={h['random_close_at_design']:>3}({h['random_frac_close']:.3f}) "
                  f"bad={int(h['random_is_bad'])} "
                  f"| EXT {Fext.name if Fext else '-'} cl={extstr} "
                  f"bad={h['extension_is_bad']}", flush=True)

    elapsed = time.time() - t0
    return results, elapsed


def summarize(results):
    """Aggregate the smooth-vs-random differential across all parameter sets.

    Reports, restricted to the MEANINGFUL regime at the design radius, the
    frac_close differential (smooth - random) per parameter set and overall, plus
    a verdict on whether the construction reproduces a SMOOTH-SPECIFIC bad line.
    """
    print("\n" + "=" * 78)
    print("SUMMARY: smooth-vs-random close-fraction differential at design delta")
    print("(meaningful regime; positive => smooth has MORE close points than random)")
    print("=" * 78)
    diffs = []
    smooth_bad = 0
    rand_bad = 0
    meaningful = 0
    smooth_strictly_more = 0
    print(f"  {'p':>4} {'n':>3} {'eta':>6} | {'smooth_frac':>11} {'rand_frac':>10} "
          f"{'diff':>7} {'sm_bad':>6} {'rn_bad':>6} {'ext_bad':>7}")
    for rec in results:
        h = rec["headline"]
        kp = rec["params"]
        if not h["smooth_meaningful"]:
            continue
        meaningful += 1
        sf, rf = h["smooth_frac_close"], h["random_frac_close"]
        diff = sf - rf
        diffs.append(diff)
        smooth_bad += int(h["smooth_is_bad"])
        rand_bad += int(h["random_is_bad"])
        if sf > rf + 1e-9:
            smooth_strictly_more += 1
        print(f"  {kp['p']:>4} {kp['n']:>3} {kp['eta']:>6.4f} | "
              f"{sf:>11.4f} {rf:>10.4f} {diff:>+7.4f} "
              f"{int(h['smooth_is_bad']):>6} {int(h['random_is_bad']):>6} "
              f"{str(h['extension_is_bad']):>7}")
    if diffs:
        import numpy as _np
        d = _np.array(diffs)
        print("-" * 78)
        print(f"  meaningful-regime parameter sets: {meaningful}")
        print(f"  close-fraction diff (smooth - random): mean={d.mean():+.4f} "
              f"max={d.max():+.4f} min={d.min():+.4f}")
        print(f"  smooth had STRICTLY more close points than random in "
              f"{smooth_strictly_more}/{meaningful} sets")
        print(f"  bad-line flagged: smooth {smooth_bad}/{meaningful}, "
              f"random {rand_bad}/{meaningful}")
        if d.max() > 0.05:
            print("  => At some scales the SMOOTH prime subgroup shows a SHARPLY larger "
                  "close-fraction\n     than random at matched (rho, delta) -- the "
                  "Kambire/BCHKS mechanism BITES.")
        else:
            print("  => No smooth-specific inflation at these scales.")
    return diffs


def _ncr(n, r):
    from math import comb
    return comb(n, r)


# ---------------------------------------------------------------------------
# Output writers.
# ---------------------------------------------------------------------------
def write_outputs(results, elapsed, out_dir=RESULTS_DIR):
    os.makedirs(out_dir, exist_ok=True)
    payload = {
        "meta": {
            "experiment": "kambire_bchks_counterexample_instantiation",
            "timestamp": time.strftime("%Y%m%d_%H%M%S"),
            "elapsed_sec": round(elapsed, 2),
            "qk_cap": QK_CAP,
            "convention": "delta = radius (ABF). Line f + gamma*g, gamma over all of F. "
                          "EXACT distance via full codeword enumeration; S* exact "
                          "branch-and-bound. delta_design = capacity - 2/s.",
            "construction": "f=X^{rm}, g=X^{(r-1)m} on D=<omega> (order n); H=<xi> "
                            "order s; k=(r-2)m; predicted close scalars = H^{(+r)} "
                            "(distinct-element r-fold sumset of H).",
        },
        "results": results,
    }
    jpath = os.path.join(out_dir, "counterexample_kambire.json")
    with open(jpath, "w", newline="\n") as f:
        json.dump(payload, f, indent=2, default=str)

    # Flat CSV: one row per (param set, delta) for subgroup vs random vs extension.
    cpath = os.path.join(out_dir, "counterexample_kambire.csv")
    rows = []
    for rec in results:
        kp = rec["params"]
        for domain_key, label in (("subgroup", "smooth_subgroup"),
                                   ("random", "random_subset"),
                                   ("extension", "ext_field_subgroup")):
            meas = rec[domain_key]
            if meas is None:
                continue
            # each domain's own capacity (extension may differ in n', k')
            dom_cap = meas.get("capacity_ext", kp["capacity"])
            for pd in meas["per_delta"]:
                rows.append({
                    "p": kp["p"], "domain_field": meas["field"],
                    "domain": label, "n": meas["n"], "k": meas["k"],
                    "s": kp["s"], "m": kp["m"], "r": kp["r"],
                    "rho_prime_construction": round(kp["rho"], 4),
                    "domain_capacity": round(dom_cap, 4),
                    "prime_design_delta": round(kp["delta"], 4),
                    "eta_2_over_s": round(kp["eta"], 4),
                    "exp_f": meas["exp_f"], "exp_g": meas["exp_g"],
                    "delta": pd["delta"],
                    "delta_minus_capacity": round(pd["delta"] - dom_cap, 4),
                    "close_count": pd["close_count"],
                    "frac_close": pd["frac_close"],
                    "S_star": pd["S_star"],
                    "ca_threshold": pd["ca_threshold"],
                    "g_max_single_agreement": meas["g_max_single_agreement"],
                    "meaningful_regime": pd["meaningful_regime"],
                    "is_bad_line": pd["is_bad_line"],
                    "predicted_sumset_size": rec["predicted_sumset_size"],
                })
    if rows:
        keys = list(rows[0].keys())
        with open(cpath, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, lineterminator="\n")
            w.writeheader()
            for row in rows:
                w.writerow(row)
    return jpath, cpath, len(rows)


# ===========================================================================
# Self-test: verify the construction identity exactly on the smallest case.
# ===========================================================================
def _self_test():
    print("counterexample_kambire self-test")
    print("=" * 64)
    # Smallest clean case p=17, n=16, s=8, m=2, r=4, k=4.
    p, n, s, m, r, k = 17, 16, 8, 2, 4, 4
    F = PrimeField(p)
    D = domain_subgroup(F, n)
    H = domain_subgroup(F, s)
    assert set(H.tolist()).issubset(set(D.tolist())), "H must be a subgroup of D"

    exp_f, exp_g = r * m, (r - 1) * m
    f = monomial_eval(F, D, exp_f)
    g = monomial_eval(F, D, exp_g)
    book = build_codeword_book(F, D, k)
    assert min_distance(book) == n - k + 1, "RS must be MDS"

    # The polynomial identity: for distinct xi_1..xi_r in H, the union of cosets
    # H_j = {a in D : a^m = xi_j} has size r*m, and  prod_j(X^m - xi_j) =
    # X^{rm} - lambda X^{(r-1)m} + R(X)  with deg R <= (r-2)m = k.  When deg R < k,
    # R is a codeword (deg < k), so  f - lambda g = -R  on the size-rm union, giving
    # agreement on rm = (1-delta)n points, i.e. dist <= delta*n.  Coset-union sanity:
    H_list = H.tolist()
    combo0 = tuple(H_list[:r])
    union = []
    for xi in combo0:
        union.extend([int(a) for a in D.tolist() if F.pow(int(a), m) == xi])
    assert len(union) == r * m, f"coset union size {len(union)} != r*m {r*m}"
    print(f"  coset structure: {r} cosets x size m={m} -> union {len(union)} = "
          f"(1-delta)n = r*m OK")

    # Find a CERTIFIED-close witness: a combo whose residual degree deg R < k, then
    # verify it is delta-close exactly (dist <= delta*n).
    delta_n = (1.0 - r / s) * n
    ss = sumset_H(F, s, r, m, k)
    assert ss["num_predicted"] >= 1, "expected at least one certified-close lambda"
    witness_lam = ss["predicted_lambda"][0]
    gamma = F.neg(witness_lam)
    w = F.add_vec(f, F.mul_scalar_vec(gamma, g))
    d = dist_to_code(book, w)
    assert d <= delta_n + 1e-9, \
        f"certified lambda={witness_lam} dist {d} should be <= delta*n {delta_n}"
    print(f"  identity check: certified lambda={witness_lam} (deg R < k) gives "
          f"dist {d} <= delta*n={delta_n:.0f} OK")
    print(f"  sumset: {ss['num_tuples']} tuples -> {ss['num_distinct']} distinct "
          f"lambda -> {ss['num_predicted']} certified-close (deg R < k)")

    # S* must be < CA threshold (no correlated agreement), and bounded by deg(g).
    S_star, info = _max_common_agreement(book, f, g)
    assert info["exact"], "S* must be exact for this small code"
    pop_g, _ = _agreement_bits(book, g)
    assert int(pop_g.max()) <= exp_g, \
        f"g=X^{exp_g} cannot agree with a codeword on more than deg(g)={exp_g} points"
    ca_threshold = int(np.ceil((1.0 - (1.0 - r / s)) * n - 1e-9))   # = r*m
    assert ca_threshold == r * m
    assert S_star < ca_threshold, \
        f"S*={S_star} must be < CA threshold {ca_threshold} (no CA)"
    print(f"  no-CA check: S*={S_star} < CA threshold={ca_threshold}; "
          f"g max single-agreement={int(pop_g.max())} <= deg(g)={exp_g} OK")

    # Global close-count must be >= 2 (many close points) => genuine bad line.
    dists = np.array([dist_to_code(book, F.add_vec(f, F.mul_scalar_vec(gg, g)))
                      for gg in range(F.q)])
    close = int((dists <= delta_n + 1e-9).sum())
    assert close >= 2, f"expected many close points, got {close}"
    print(f"  bad-line check: global close-count={close} (>=2), meaningful regime, "
          f"no CA => GENUINE BAD LINE on smooth prime subgroup OK")

    # Sanity: a RANDOM domain with the same monomials should (almost surely) have a
    # much smaller close-count -- the structure is what creates the badness.
    rng = np.random.default_rng(7)
    Lr = domain_random(F, n, rng)
    fr = monomial_eval(F, Lr, exp_f)
    gr = monomial_eval(F, Lr, exp_g)
    bookr = build_codeword_book(F, Lr, k)
    distsr = np.array([dist_to_code(bookr, F.add_vec(fr, F.mul_scalar_vec(gg, gr)))
                       for gg in range(F.q)])
    close_r = int((distsr <= delta_n + 1e-9).sum())
    print(f"  random-domain control: same monomials close-count={close_r} "
          f"(vs smooth {close}); structure-specificity {'OK' if close_r < close else 'WEAK'}")

    print("=" * 64)
    print("ALL counterexample_kambire SELF-TESTS PASSED")


if __name__ == "__main__":
    import sys
    args = [a for a in sys.argv[1:] if not a.startswith("--out=")]
    out_dir = next((a.split("=", 1)[1] for a in sys.argv[1:]
                    if a.startswith("--out=")),
                   os.environ.get("KAMBIRE_OUT", RESULTS_DIR))
    if args and args[0] == "test":
        _self_test()
    else:
        _self_test()
        print()
        print("=" * 64)
        print("RUNNING FULL INSTANTIATION (primes 17, 97, 193, 257)")
        print("=" * 64)
        results, elapsed = run()
        summarize(results)
        jpath, cpath, nrows = write_outputs(results, elapsed, out_dir=out_dir)
        print("-" * 64)
        print(f"Done in {elapsed:.1f}s. Wrote:\n  {jpath}\n  {cpath} ({nrows} rows)")
