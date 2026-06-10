"""
listsize_resolution.py
=======================

Self-contained resolver for **Proximity Prize sub-problem 2** — the *interleaved
list-size* grand challenge (ABF survey, eprint 2026/680, §1, boxed p.5):

    C := RS[F, L, k], L a smooth (power-of-two multiplicative-subgroup) domain,
    rho = k/n in {1/2, 1/4, 1/8, 1/16}, constant interleaving m, eps* = 2^-128.
    Determine the largest delta*_C in [0,1] such that

        | Lambda( C^{equiv m}, delta*_C ) |  <=  eps* * |F|,

    with a proof that for all delta > delta*_C the bound fails.

Here Lambda(C^{equiv m}, delta) is the list of m-tuples of RS codewords within
relative radius delta of a target word under the **columnwise / mixed metric**
d(W, c) = #{ i : exists j, w_j[i] != c_j[i] } (a column counts as an error if any
of the m rows disagrees there).  This is ABF Definition 2.9 (interleaved code
C^{equiv m}) with the list Lambda(C, delta) := max_f |{c in C : Delta(c,f) <= delta}|.

KEY STRUCTURAL FACT (the field-size lever).
-------------------------------------------
The constraint is `|Lambda| <= eps* * |F| = 2^-128 * |F|`.  Write B := 2^-128*|F|.

    * |F| <  2^128 :  B < 1.  Since any target that *is* an interleaved codeword
                      has |Lambda| >= 1, the constraint is violated for EVERY
                      delta in the support of a codeword.  delta*_C DOES NOT
                      EXIST.  This is exactly ABF's "assuming |F| sufficiently
                      large so that such a delta*_C exists" caveat.  -> DEGENERATE.
    * |F| == 2^128 :  B == 1.  Forces |Lambda| <= 1, i.e. *unique decoding of the
                      interleaved object*.  Rigorously delta*_C = (1-rho)/2 (the
                      unique-decoding radius): below it every Hamming ball of
                      radius delta < delta_min/2 contains <= 1 codeword (MDS), so
                      the interleaved list is also <= 1; at any larger delta a
                      worst-case target already has list >= 2.  -> BINDING.
    * |F| >  2^128 :  B = 2^(logF-128) >= 2^64, a large budget.  PROVEN bracket
                      [J - o(1), r_E]:  the Johnson list bound (a CONSTANT) gives a
                      PROVEN LOWER reach delta* >= J - eta_min with
                      eta_min = 1/(2 rho B^{1/m}) (Johnson proves J - eta_min at
                      finite budget, NOT J itself; J - o(1) is the explicit
                      lower-endpoint convention); the Elias/CS25 volume bound gives
                      a PROVEN UPPER ceiling delta* <= r_E = 1 - H_q(rho)
                      ~ (1-rho) - 1/log2 q (r_E is the FORMULA convention; the
                      mathematically exact proven object is the inverse-entropy
                      crossing H_q^{-1}(1-rho), which is NOT equal to the formula
                      at any deployed rate -- above it at rho=1/2, below at
                      rho<=1/4, gap <= 0.0017).  The value lies in [J - o(1), r_E],
                      NOT at r_E.  Reaching ~r_E from below is CONJECTURAL (needs a
                      worst-case large-list smooth-domain RS
                      list-decoding-beyond-Johnson theorem; ABF Sec 7.9 /
                      sub-lemma P'/(D2)); the volume-bound crossing this module
                      computes is that conjectural upper reach, NOT a proven
                      delta*.  -> PROVEN BRACKET [J - o(1), r_E] for large fields.

LIST-SIZE MODELS (each with its exact source and validity window).
------------------------------------------------------------------
We never invent constants.  Three rigorously-sourced single-RS list bounds, then
the interleaving relation |Lambda(C^{equiv m},delta)| <= |Lambda(C,delta)|^m
(ABF Definition 2.9; the upper inequality).

  (UD)  unique-decoding.  For delta < delta_min(C)/2 = (1-rho+1/n)/2, every ball
        of radius delta contains at most ONE codeword (classical MDS / Singleton).
        Hence |Lambda(C,delta)| = 1 and |Lambda(C^{equiv m},delta)| = 1.
        RIGOROUS, all fields, all m.  Source: minimum-distance / Singleton bound.

  (J)   Johnson (MDS form).  ABF Corollary 3.3 (from Johnson 1962, Theorem 3.2):
        for delta = 1 - sqrt(rho) - eta with eta > 0,
            |Lambda(C, delta)|  <=  1 / (2 * eta * rho).
        A CONSTANT independent of n.  Interleaved:
            |Lambda(C^{equiv m}, delta)|  <=  (1/(2*eta*rho))^m.
        RIGOROUS, all RS codes incl. smooth domains, all fields, all m.

  (CAP) capacity / Elias CEILING (NOT a worst-case upper bound).  Crites-Stewart
        2025/2046 Thm 1 + Lemma 1 (Elias 1957 list-decoding-capacity volume bound).
        IMPORTANT: this is an AVERAGE + LOWER bound, not a proven worst-case UPPER
        bound on this code's list.  CS25 Lemma 1 sandwiches the close-count for a
        uniformly RANDOM center (a first-moment / AVERAGE quantity); CS25 Thm 1 and
        Thm 7.4.1(ii) are LOWER bounds (some word / some code has a LARGE list).  An
        average count and a lower bound CANNOT upper-bound the worst case.  What it
        rigorously gives is the UPPER CEILING delta*_C^{(2)} <= r_E = 1 - H_q(rho):
        once delta exceeds r_E the worst-case list is q^{Omega(n)} >> B for EVERY
        code (deep holes x^k, 1/(x-a) exist on every domain), so no delta* can
        exceed r_E.  Below r_E it gives NO positive (lower) reach: it does not
        bound the worst-case list of this specific smooth code in (J, r_E).
        We model:  log_q |Lambda(C,delta)| = max(0, (H_q(delta) - (1-rho)) * n),
        interleaved x m, and use the CROSSING delta with budget B as the
        CONJECTURAL UPPER REACH (~ r_E) -- NOT a proven delta*.  Source: CS25 Thm 1.

The resolver reports, for delta*_C^{(2)} in the large-field regime:
  (a) the PROVEN LOWER reach (Johnson floor: list <= 1/(2 eta rho) <= B => delta* >= J - o(1)),
  (b) the PROVEN UPPER ceiling r_E (Elias/CS25), and
  (c) the volume-bound crossing, clearly labeled CONJECTURAL upper reach (assumes the
      open smooth-domain large-list bound; NOT a proven worst-case bound).
It labels the regime (DEGENERATE / BINDING / BRACKET[J - o(1), r_E]).

Import policy: this module imports ONLY proximity_parameters (johnson_radius,
qary_entropy, capacity, unique_decoding_radius and friends).  It does NOT touch
bounds.py.  All numbers below are derived from the three sourced bounds above.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Optional

from proximity_parameters import (
    as_fraction,
    capacity,
    johnson_radius,
    qary_entropy,
    unique_decoding_radius,
)

EPS_STAR_LOG2 = -128.0  # eps* = 2^-128 (the prize target)

# Cap the alphabet passed to H_q at a value strictly above every field size we
# resolve (256-bit + headroom), so no deployed field is ever capped.  WARNING
# (R29 referee item M2): the capacity column r_E = 1 - H_q(rho) depends on q at
# the ~1/log2(q) scale, so capping BELOW the actual field size silently shifts
# r_E — the old cap of 200 reported the 200-bit value 0.49500 instead of the
# true 256-bit value 0.49609 (rho = 1/2).  The cap exists only to keep
# math.log() off absurdly large ints; it must exceed max logF (here 256).
_QCAP_LOG2 = 320


def _qcap(logF: float) -> int:
    return 2 ** int(min(round(logF), _QCAP_LOG2))


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------


@dataclass
class Resolution:
    """delta*_C^{(2)} under one list-size model, with provenance."""

    rho: Fraction
    n: int
    m: int
    logF: float
    model: str                 # "UD" | "Johnson" | "capacity"
    regime: str                # "DEGENERATE" | "BINDING" | "BRACKET" (proven [J - o(1), r_E])
    delta_star: Optional[float]  # None when no feasible delta exists
    budget_log2: float         # log2(eps* * |F|) = logF - 128
    note: str = ""

    def __str__(self) -> str:  # pragma: no cover - display only
        d = "n/a" if self.delta_star is None else f"{self.delta_star:.5f}"
        return (f"rho={self.rho} n=2^{int(round(math.log2(self.n)))} m={self.m} "
                f"logF={self.logF:.0f} [{self.model}] {self.regime}: "
                f"delta*={d}  (budget 2^{self.budget_log2:+.0f})  {self.note}")


# ---------------------------------------------------------------------------
# The three list-size models  (log2 of the INTERLEAVED list size)
# ---------------------------------------------------------------------------


def log2_list_ud(rho: Fraction, n: int, m: int, delta: float) -> Optional[float]:
    """log2 |Lambda(C^{equiv m}, delta)| under the unique-decoding model.

    For delta < delta_min/2 = (1-rho+1/n)/2 the list is exactly 1 (log2 = 0).
    Outside that window this model says nothing (returns None).  Rigorous (MDS).
    """
    rho_f = float(rho)
    ud = (1.0 - rho_f + 1.0 / n) / 2.0  # delta_min(C)/2, MDS
    if delta < ud:
        return 0.0  # list size 1
    return None


def log2_list_johnson(rho: Fraction, n: int, m: int, delta: float) -> Optional[float]:
    """log2 |Lambda(C^{equiv m}, delta)| under the Johnson model (ABF Cor 3.3).

    delta = 1 - sqrt(rho) - eta, eta > 0  =>  |Lambda(C,delta)| <= 1/(2 eta rho),
    interleaved <= (1/(2 eta rho))^m.  Returns None at/above the Johnson radius.
    """
    rho_f = float(rho)
    J = johnson_radius(rho)
    eta = J - delta
    if eta <= 0:
        return None  # Johnson bound only valid strictly below J
    single = 1.0 / (2.0 * eta * rho_f)
    if single < 1.0:
        single = 1.0  # a list is at least 1 when the target is a codeword
    return m * math.log2(single)


def log2_list_capacity(rho: Fraction, n: int, m: int, delta: float,
                       logF: float) -> float:
    """log2 |Lambda(C^{equiv m}, delta)| under the Elias/CS25 worst-case model.

    log_q |Lambda(C,delta)| = max(0, (H_q(delta) - (1-rho)) * n); interleaved
    multiplies by m; convert to log2 via * log2 q.  Below the list-decoding-
    capacity radius 1 - H_q(rho) this is 0 (a constant list, dominated by the
    Johnson model); at/above it the list grows super-polynomially.  Source:
    Crites-Stewart 2025/2046 Thm 1 / Lemma 1 (Elias volume bound).
    """
    rho_f = float(rho)
    qcap = _qcap(logF)
    Hq = qary_entropy(Fraction(delta).limit_denominator(10 ** 12), qcap)
    expo_q = max(0.0, (Hq - (1.0 - rho_f)) * n)  # log_q of single-RS list
    return m * expo_q * logF  # log2 = log_q * log2(q) = expo * logF


# ---------------------------------------------------------------------------
# Resolver:  largest delta with model list <= eps* |F|, by bisection
# ---------------------------------------------------------------------------


def _budget_log2(logF: float) -> float:
    return logF + EPS_STAR_LOG2  # log2(2^-128 * |F|) = logF - 128


def resolve_ud(rho, n, m, logF) -> Resolution:
    """delta*_C under the unique-decoding model: rigorous, list <= 1.

    Feasible iff the budget admits list size 1, i.e. budget_log2 >= 0
    (|F| >= 2^128).  Then delta*_C = delta_min/2 = (1-rho+1/n)/2 (the largest
    delta with list provably 1).  Otherwise DEGENERATE.
    """
    rho = as_fraction(rho)
    n = int(n)
    m = int(m)
    logF = float(logF)
    B = _budget_log2(logF)
    ud = (1.0 - float(rho) + 1.0 / n) / 2.0
    if B < 0:
        return Resolution(rho, n, m, logF, "UD", "DEGENERATE", None, B,
                          note="2^-128|F| < 1: even list=1 is infeasible.")
    regime = "BINDING" if B == 0 else "BRACKET"
    return Resolution(rho, n, m, logF, "UD", regime, ud, B,
                      note="rigorous list<=1 up to delta_min/2 (MDS).")


def _bisect(rho, n, m, logF, lo, hi, log2_list_fn, iters=300) -> float:
    """Largest delta in (lo, hi) with log2_list_fn(delta) <= budget_log2(logF)."""
    B = _budget_log2(logF)
    f_lo = log2_list_fn(lo)
    if f_lo is None or f_lo > B:
        return lo  # not even feasible at the bottom of the window
    for _ in range(iters):
        mid = (lo + hi) / 2.0
        v = log2_list_fn(mid)
        if v is not None and v <= B:
            lo = mid
        else:
            hi = mid
    return lo


def resolve_johnson(rho, n, m, logF) -> Resolution:
    """delta*_C under the Johnson model (ABF Cor 3.3), rigorous for smooth RS.

    Largest delta = J - eta with (1/(2 eta rho))^m <= 2^-128|F|.  Closed form:
    eta_min = 1 / (2 rho * (2^-128|F|)^{1/m});  delta* = J - eta_min, clamped to
    the unique-decoding radius from below (the Johnson statement is only
    informative above UD; below UD the UD model already gives list 1).
    """
    rho = as_fraction(rho)
    n = int(n)
    m = int(m)
    logF = float(logF)
    rho_f = float(rho)
    J = johnson_radius(rho)
    B = _budget_log2(logF)
    ud = float(unique_decoding_radius(rho))
    # per-instance budget on the single-code list: 2^(B/m)
    per_instance_log2 = B / m
    if per_instance_log2 < 0:
        # budget < 1 per instance: 1/(2 eta rho) <= 2^(B/m) < 1 impossible
        return Resolution(rho, n, m, logF, "Johnson", "DEGENERATE", None, B,
                          note="per-instance budget < 1: Johnson list (a constant "
                               ">=1) cannot fit; only UD region survives.")
    eta_min = 1.0 / (2.0 * rho_f * (2.0 ** per_instance_log2))
    delta = J - eta_min
    if delta <= ud:
        # Johnson model gives nothing beyond unique decoding at this budget.
        regime = "BINDING" if B == 0 else "BRACKET"
        return Resolution(rho, n, m, logF, "Johnson", regime, ud, B,
                          note=f"Johnson eta_min={eta_min:.4f} pushes below UD; "
                               f"binding delta* = UD radius {ud:.5f}.")
    regime = "BRACKET"
    return Resolution(rho, n, m, logF, "Johnson", regime, delta, B,
                      note=f"delta* = J - eta_min, eta_min={eta_min:.6f}.")


def resolve_capacity(rho, n, m, logF) -> Resolution:
    """delta*_C under the Elias/CS25 worst-case capacity model.

    Largest delta with m*(H_q(delta)-(1-rho))*n*log2(q) <= 2^-128|F| in log2.
    Below the list-decoding-capacity radius the worst-case list is a constant, so
    this model is only informative as a NEGATIVE ceiling: it returns the radius at
    which the worst-case list crosses the budget, which for any sane n sits right
    at the list-decoding-capacity radius 1 - H_q(rho) ~ (1-rho) - 1/log2 q.
    """
    rho = as_fraction(rho)
    n = int(n)
    m = int(m)
    logF = float(logF)
    B = _budget_log2(logF)
    if B < 0:
        return Resolution(rho, n, m, logF, "capacity", "DEGENERATE", None, B,
                          note="2^-128|F| < 1: list-size constraint unsatisfiable.")
    lo = float(unique_decoding_radius(rho))
    hi = float(capacity(rho)) - 1e-12
    delta = _bisect(rho, n, m, logF, lo, hi,
                    lambda d: log2_list_capacity(rho, n, m, d, logF))
    regime = "BINDING" if B == 0 else "BRACKET"
    return Resolution(rho, n, m, logF, "capacity", regime, delta, B,
                      note="worst-case (adversarial-target) ceiling ~ "
                           "list-decoding-capacity radius 1 - H_q(rho).")


def resolve_all(rho, n, m, logF) -> dict:
    """All three model resolutions, framed as a PROVEN BRACKET (large fields).

    Honest framing of delta*_C^{(2)}:
      * DEGENERATE (|F| < 2^128): delta*_C does not exist.
      * BINDING (|F| == 2^128): delta*_C = UD radius (the PROVEN list<=1 reach).
      * BRACKET (|F| > 2^128): delta*_C lies in the PROVEN two-sided bracket
        [J - o(1), r_E].  Returned keys:
          - `proven_floor`      = Johnson-model value: PROVEN delta* >= J - eta_min
            with eta_min = 1/(2 rho B^{1/m}) (Johnson proves J - eta_min at finite
            budget, NOT J itself; J - o(1) is the explicit lower-endpoint
            convention; list <= 1/(2 eta rho) <= B; field-robust, smooth domains).
          - `proven_ceiling`    = r_E = 1 - H_q(rho): PROVEN delta* <= r_E
            (Elias/CS25; above the ceiling the worst-case list is q^{Omega(n)} >> B
            for every code).  Computed as the FORMULA 1 - H_q(rho) at the field's
            true q (R29/M2).  Convention (load-bearing): the mathematically exact
            proven ceiling object is the inverse-entropy CROSSING H_q^{-1}(1-rho),
            which is NOT equal to the formula at any deployed rate -- the crossing
            sits ABOVE the formula at rho=1/2 (+9.6e-5 at 31-bit, +1.7e-7 at
            256-bit) and BELOW it at rho<=1/4 (~2e-5 at 256-bit); the gap is
            <= 0.0017 everywhere deployed.  All quoted r_E values are values of
            the formula 1 - H_q(rho), not of the crossing.
          - `conjectural_upper` = the SAME volume-bound crossing (~ r_E), but read as
            an ACHIEVED reach.  This is NOT a proven delta*: the volume bound is an
            average + lower bound and does not upper-bound the worst-case list of
            this specific smooth code in (J, r_E).  Reaching ~r_E needs the open
            smooth-domain large-list bound (ABF Sec 7.9 / sub-lemma P'/(D2)).
      `binding` is kept for back-compat and is the model-appropriate honest value:
      DEGENERATE -> None; BINDING -> UD; BRACKET -> the PROVEN FLOOR (J).  The
      conjectural upper reach is reported separately and NEVER as `binding`.
    """
    ud = resolve_ud(rho, n, m, logF)
    joh = resolve_johnson(rho, n, m, logF)
    cap = resolve_capacity(rho, n, m, logF)
    B = _budget_log2(float(logF))
    proven_floor = None
    proven_ceiling = None
    conjectural_upper = None
    if B < 0:
        binding = None
        regime = "DEGENERATE"
    elif B == 0:
        binding = ud.delta_star
        regime = "BINDING"
    else:
        # PROVEN bracket [J - o(1), r_E].  Floor = Johnson (proven, = J - eta_min);
        # ceiling = the FORMULA r_E = 1 - H_q(rho) at the field's true q (R29/M2);
        # the exact proven object is the crossing H_q^{-1}(1-rho), whose offset
        # from the formula is rate-dependent in sign (above at rho=1/2, below at
        # rho<=1/4); the finite-n volume-bound crossing read as an ACHIEVED reach
        # is CONJECTURAL.  `binding` = the proven floor, never the conjectural reach.
        proven_floor = joh.delta_star          # PROVEN delta* >= J - o(1)
        proven_ceiling = 1.0 - qary_entropy(   # PROVEN delta* <= r_E (exact)
            as_fraction(rho), _qcap(float(logF)))
        conjectural_upper = cap.delta_star      # volume crossing ~r_E, read as reach -> CONJECTURAL
        binding = proven_floor
        regime = "BRACKET"
    return {"ud": ud, "johnson": joh, "capacity": cap,
            "binding": binding, "regime": regime, "budget_log2": B,
            "proven_floor": proven_floor, "proven_ceiling": proven_ceiling,
            "conjectural_upper": conjectural_upper}


# ---------------------------------------------------------------------------
# Table
# ---------------------------------------------------------------------------

FIELDS = [
    ("M31", 31),                 # Mersenne-31 prime field, 2^31 - 1
    ("Goldilocks", 64),          # 2^64 - 2^32 + 1
    ("128-bit", 128),            # e.g. M31^4 (~124) / BabyBear^4 / a 128-bit prime
    ("256-bit", 256),            # large extension / 256-bit prime
]
RATES = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 16)]
NS = [2 ** 16, 2 ** 20, 2 ** 24]


def print_table(m: int = 2) -> None:
    """Print the prize table: rho x field x n, with the delta*_C^{(2)} bracket."""
    print("=" * 108)
    print(f"Proximity Prize sub-problem 2 — interleaved list-size delta*_C^(2),  "
          f"m={m},  eps*=2^-128")
    print("Columns: J=Johnson radius 1-sqrt(rho), cap=Singleton capacity 1-rho, "
          "UD=(1-rho)/2.")
    print("delta* models: [UD] rigorous list<=1; [Joh] Johnson Cor 3.3 (PROVEN "
          "floor, smooth RS);")
    print("  [Cap*] Elias/CS25 volume-bound crossing ~ r_E: the CONJECTURAL upper "
          "reach if read as")
    print("  achieved (the volume bound is an average+lower bound, NOT a proven "
          "worst-case bound;")
    print("  reaching it needs the open smooth-domain large-list result).  The "
          "PROVEN UPPER CEILING")
    print("  delta* <= r_E (above it the worst-case list is q^Omega(n) >> B) is "
          "quoted in the bracket")
    print("  column as the FORMULA r_E = 1-H_q(rho) at the field's true q (R29/M2); "
          "the exact proven")
    print("  object is the crossing H_q^{-1}(1-rho) -- ABOVE the formula at rho=1/2 "
          "(+1.7e-7 at 256-bit),")
    print("  ~2e-5 BELOW it at rho<=1/4 (rate-dependent sign).")
    print("PROVEN delta*_C^(2) BRACKET (|F|>2^128) = [ d*[Joh] (floor, = J - "
          "eta_min) , r_E (ceiling) ].")
    print("Regime: DEGEN (|F|<2^128, no delta*) / BIND (|F|=2^128, list<=1, exact) / "
          "BRACKET (|F|>2^128).")
    print("=" * 108)
    header = (f"{'rho':>5} {'field':>11} {'logF':>5} {'n':>6} | "
              f"{'UD':>6} {'J':>6} {'cap':>6} | "
              f"{'d*[UD]':>7} {'d*[Joh]':>8} {'d*[Cap*]':>8} | "
              f"{'PROVEN d* (bracket / value)':>30} {'regime':>9}")
    print(header)
    print(f"{'':>5} {'':>11} {'':>5} {'':>6} | {'':>6} {'':>6} {'':>6} | "
          f"{'(rig.)':>7} {'(floor)':>8} {'(cross)':>8} | "
          f"{'floor=PROVEN; ceiling=exact r_E':>30} {'':>9}")
    print("-" * 108)
    for rho in RATES:
        J = johnson_radius(rho)
        cap_v = float(capacity(rho))
        ud_v = float(unique_decoding_radius(rho))
        for fname, logF in FIELDS:
            for n in NS:
                r = resolve_all(rho, n, m, logF)
                def fmt(x):
                    return "  -  " if (x is None) else f"{x:.5f}"
                # Honest headline: DEGEN -> none; BIND -> exact UD; BRACKET -> [J - o(1), r_E].
                if r["regime"] == "DEGENERATE":
                    headline = "does not exist"
                elif r["regime"] == "BINDING":
                    headline = f"{r['binding']:.5f} (exact, UD)"
                else:  # BRACKET
                    headline = (f"[{r['proven_floor']:.5f}, "
                                f"{r['proven_ceiling']:.5f}]")
                line = (f"{str(rho):>5} {fname:>11} {logF:>5} "
                        f"2^{int(round(math.log2(n))):<3} | "
                        f"{ud_v:>6.4f} {J:>6.4f} {cap_v:>6.4f} | "
                        f"{fmt(r['ud'].delta_star):>7} "
                        f"{fmt(r['johnson'].delta_star):>8} "
                        f"{fmt(r['capacity'].delta_star):>8} | "
                        f"{headline:>30} {r['regime']:>9}")
                print(line)
            print("-" * 108)


def print_lever_summary() -> None:
    """Print the field-size lever: what eps*|F| forces, per field."""
    print()
    print("=" * 100)
    print("THE FIELD-SIZE LEVER:  eps* * |F| = 2^-128 * |F|  is the entire list-size budget.")
    print("=" * 100)
    print(f"{'field':>12} {'logF':>5} {'budget = 2^-128|F|':>22} {'forces':>40}")
    print("-" * 100)
    for fname, logF in [("M31", 31), ("M31^2", 62), ("Goldilocks", 64),
                        ("M31^4", 124), ("128-bit", 128), ("192-bit", 192),
                        ("256-bit", 256)]:
        B = logF - 128
        if B < 0:
            forces = f"2^{B} < 1  => delta*_C DOES NOT EXIST (degenerate)"
        elif B == 0:
            forces = "= 1  => |Lambda|<=1 (unique decoding of interleaved obj)"
        else:
            forces = f"= 2^{B}  => |Lambda|<=2^{B} (loose)"
        print(f"{fname:>12} {logF:>5} {('2^%+d' % B):>22} {forces:>40}")
    print("-" * 100)
    print("Crossover: |F| = 2^128 is the knife-edge.  Below it sub-problem 2 is")
    print("vacuous (no delta*); at it, delta* = (1-rho)/2 (exact); above it, delta*")
    print("lies in the PROVEN BRACKET [J - o(1), r_E] with r_E = 1 - H_q(rho) ~")
    print("(1-rho) - 1/log2|F| (Johnson proves J - eta_min at finite budget, not J;")
    print("reaching ~r_E from the Johnson floor is CONJECTURAL, not proven).")


def print_rho_half_focus(m: int = 2) -> None:
    """The sharp rho=1/2 verdict the brief asks for, across fields, n=2^20."""
    print()
    print("=" * 100)
    print(f"SHARP VERDICT FOCUS:  rho = 1/2, n = 2^20, m = {m}")
    print(f"  Johnson J = {johnson_radius(Fraction(1,2)):.5f}, "
          f"capacity = {float(capacity(Fraction(1,2))):.5f}, "
          f"UD = {float(unique_decoding_radius(Fraction(1,2))):.5f}")
    print("=" * 100)
    for fname, logF in [("M31", 31), ("Goldilocks", 64), ("128-bit", 128),
                        ("192-bit", 192), ("256-bit", 256)]:
        r = resolve_all(Fraction(1, 2), 2 ** 20, m, logF)
        if r["regime"] == "DEGENERATE":
            verdict = "delta*_C^(2) = DOES NOT EXIST"
        elif r["regime"] == "BINDING":
            verdict = f"delta*_C^(2) = {r['binding']:.5f}  (exact, UD radius)"
        else:  # BRACKET
            verdict = (f"delta*_C^(2) in PROVEN BRACKET "
                       f"[{r['proven_floor']:.5f}, {r['proven_ceiling']:.5f}]  "
                       f"(upper reach ~{r['conjectural_upper']:.5f} CONJECTURAL)")
        print(f"  {fname:>9} (2^{logF:>3}): regime={r['regime']:>10}  {verdict}")
    print("-" * 100)
    print("  Reading: M31/Goldilocks -> DEGENERATE (2^-128|F|<1, no delta*).")
    print("           128-bit        -> BINDING, delta* = 0.25 (= UD radius, list<=1, exact).")
    print("           192/256-bit    -> PROVEN BRACKET [J - o(1), r_E]: floor J=0.29289 (Johnson,")
    print("                             PROVEN to J - eta_min), ceiling r_E ~ 0.4961 (256-bit) / ~0.4948")
    print("                             (192-bit) (Elias/CS25, PROVEN).  Reaching ~r_E from")
    print("                             the J floor is CONJECTURAL (open smooth-domain")
    print("                             large-list bound; ABF 7.9 / P'/(D2)), NOT proven.")


if __name__ == "__main__":
    print_lever_summary()
    print_table(m=2)
    print("\n\nSame table for m=3 (interleaving degree 3):")
    print_table(m=3)
    print_rho_half_focus(m=2)
