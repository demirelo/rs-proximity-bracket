"""
delta_star.py
=============

Produce the concrete **delta*_C brackets** — the project's headline numeric
deliverable — for the Proximity Prize parameter family.

For each (rho, field, n) we report:

* **provable-safe delta** (``delta_safe``): the LARGEST radius ``delta`` for
  which the best *VERIFIED* positive MCA bound certifies
  ``eps_mca(C, delta) <= 2^-128`` (single code, m = 1 interleaving for the
  list/soundness; the Bordage-Chiesa trade parameter is swept to push delta
  toward Johnson while keeping the error <= target).  We also report the query
  count ``t`` of the full composite soundness expression at that delta, and the
  resulting bits of security.  Because the best proven positive bound runs out
  at the Johnson radius ``1 - sqrt(rho)``, ``delta_safe`` never exceeds Johnson.

* **delta_unsafe** (``delta_unsafe``, FIELD-AGNOSTIC) = the per-field
  **assembled-lemma** unsafe ceiling (R13):

      ``delta_unsafe = (1 - rho) - 2/s_max(b)``

  where the Kambiré quotient parameter ``s`` must be a **power of two** (R13
  s-integrality; the smooth domain has order ``n = 2^t``) and ``s_max(b)`` is
  the largest such ``s`` passing (i) the KK25 distinctness calibration
  ``p > phi(s)^{phi(s)}``, (ii) the prize-threshold count ``N(s, rho) >
  2^{b-128}``, and (iv) the above-Johnson radius condition
  ``s > 2/(sqrt(rho)-rho)`` — see ``bounds.kambire_smax``.  ``s_max =
  16/16/32`` at ``b = 31/64/128``; at ``b = 256`` **no valid s exists**
  (count ``3^32 = 2^51 << 2^128``), so the mechanism establishes **no
  threshold ceiling** there (rendered as an unclosed bracket; the generic
  CS25/Elias ceiling at ~ r_E still stands).  Lemma status: ASSEMBLED from
  cited components (KK25 cited not re-proved; rho = 1/2 N1-conditional).  The
  pre-R13 continuum ``(1 - rho) - 6/log2|F|`` (Linnik-window arithmetic) is an
  asymptotic reference only — it is NOT the tabled value.  This applies to
  PRIME **and** genuine odd-characteristic EXTENSION fields identically: the
  near-capacity counterexample's bad-scalar count is a characteristic-zero
  cyclotomic invariant (``n2-verdict.md``), so ``GF(p^e)`` inherits the same
  mechanism as ``GF(p)`` — the earlier "extension is OPEN" reading is
  SUPERSEDED.

The open gap interval per family is ``[delta_safe, delta_unsafe]`` — the same
bracket for prime and extension (at 256-bit the upper end is not
threshold-established: ``[delta_safe, —)``).

We use only VERIFIED positive bounds for ``delta_safe`` so the headline number is
defensible; the BCHKS Thm 1.5 bound (verified=False, hidden C_rho) is reported
separately as an *optimistic* comparison.

IMPORTANT (do not overclaim): ``delta_safe`` and ``delta_unsafe`` are
**best-known brackets: proven floor + inferred ceiling**, NOT a resolution of
``delta*_C``.  The true ``delta*_C`` lies somewhere in
``[delta_safe, delta_unsafe]`` and pinning it is the open prize question.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from fractions import Fraction
from typing import List, Optional

import bounds as _bounds
from proximity_parameters import (
    bits,
    capacity,
    johnson_radius,
    list_decoding_capacity_radius,
    unique_decoding_radius,
)
from soundness import (
    DEFAULT_COST,
    CostModel,
)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "out")

TARGET_BITS = 128.0

# Field family (name, log2|F|, field_type).  Small exact fields computed
# cancellation-safely.  field_type is retained for labelling, but the capacity
# no-go is now FIELD-AGNOSTIC (see below).
from proximity_parameters import log2_pow2_minus_c  # noqa: E402

_LOG2_M31 = log2_pow2_minus_c(31, 1)            # 2^31 - 1
_LOG2_BABYBEAR = log2_pow2_minus_c(31, 2 ** 27)  # 15*2^27 + 1
_LOG2_GOLDILOCKS = log2_pow2_minus_c(64, 2 ** 32 - 1)  # 2^64 - 2^32 + 1

# (name, log2|F|, field_type).  The base SNARK fields are PRIME; deployments do
# the SOUNDNESS argument over an EXTENSION of the base field (degree-4/5 over
# M31/BabyBear, degree-2 over Goldilocks) to get a large enough |F|.  The
# capacity-window no-go is now FIELD-AGNOSTIC: it applies IDENTICALLY to prime
# and genuine odd-characteristic extension fields, because the near-capacity
# counterexample's bad-scalar count is a characteristic-zero cyclotomic
# invariant (n2-verdict.md).  So extension rows report the SAME delta_unsafe and
# the SAME bracket as the matching-size prime row (no "OPEN" upper end).  We
# include both prime and extension to show that the result is identical.
FIELDS = [
    ("Mersenne31",      _LOG2_M31, "prime"),
    ("BabyBear",        _LOG2_BABYBEAR, "prime"),
    ("Goldilocks",      _LOG2_GOLDILOCKS, "prime"),
    ("prime128",        128.0, "prime"),
    ("prime256",        256.0, "prime"),
    # extension fields used for the soundness argument (field_type="extension"):
    ("M31^4 (ext)",     4 * _LOG2_M31, "extension"),         # ~124 bits
    ("Goldilocks^2(ext)", 2 * _LOG2_GOLDILOCKS, "extension"),  # ~128 bits
    ("ext256",          256.0, "extension"),                 # generic 256-bit ext
]

RATES = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 16)]
LOG2NS = [16, 20, 24, 30]


# ---------------------------------------------------------------------------
# provable-safe delta via the best VERIFIED positive MCA bound, sweeping m
# ---------------------------------------------------------------------------

def _verified_mca_eps(rho: Fraction, n: int, logF: float, delta: float,
                      eta: Optional[float], m: int) -> Optional[float]:
    """Smallest eps_mca from any *VERIFIED* positive single-code MCA bound at
    these params, or ``None`` if no verified positive bound applies.

    Excludes the two capacity-region bounds (``unknown-beyond-johnson``, value
    1.0/verified=False, and ``proven-near-capacity-nogo``, value 1.0, regime
    "capacity") — neither is a positive certificate — and excludes the meta
    interleaving bound (we evaluate single-code here).
    """
    best = None
    for b in _bounds.applicable(rho, n, logF, delta, eta=eta, m=m, kind="mca",
                                include_meta=False):
        if not b.verified:             # drops unknown-beyond-johnson (OPEN band)
            continue
        if b.regime == "capacity":     # drops proven-near-capacity-nogo
            continue
        v = b.value(rho, n, logF, delta, eta=eta, m=m)
        if v >= 1.0:                   # vacuous
            continue
        if best is None or v < best:
            best = v
    return best


def provable_safe_delta(rho: Fraction, n: int, logF: float,
                        target_bits: float = TARGET_BITS,
                        m_max: int = 2_000_000) -> dict:
    """Largest delta with a VERIFIED positive bound giving eps_mca <= 2^-target.

    Strategy.  The only VERIFIED Johnson-regime positive bound is Bordage-Chiesa,
    whose validity ceiling ``dmax(m) = 1 - (1 + 1/(2m)) sqrt(rho)`` increases
    with the trade parameter ``m`` toward the Johnson radius, while its error
    ``(m+1/2)^7 n^2 d / (3 rho^{3/2} |F|)`` increases with ``m``.  So:

      * find the LARGEST ``m`` with ``eps_BC(m) <= 2^-target`` (binary search on
        the monotone-increasing error);
      * ``delta_safe`` = ``dmax(that m)``, capped at the Johnson radius.

    Below the unique-decoding radius the (verified) UD bound ``n/|F|`` may give a
    smaller error; we also test ``delta`` up to the UD radius via that bound.
    Returns the better (larger) achievable ``delta_safe`` and which bound/m it
    used.  ``feasible=False`` if NO verified positive bound reaches the target at
    ANY delta > 0 (i.e. the field is too small even at unique decoding).
    """
    target_eps_log2 = -target_bits     # want log2(eps) <= -target_bits
    J = johnson_radius(rho)
    udr = float(unique_decoding_radius(rho))

    result = {
        "feasible": False, "delta_safe": 0.0, "bound": None, "m": None,
        "eps_bits": float("-inf"), "regime": None,
    }

    # --- candidate 1: unique-decoding bound n/|F| (m-independent) ---
    # Valid up to udr; error is constant = n/|F|.  If it meets target, delta can
    # be the full udr.
    ud = _bounds.get_bound("unique-decoding")
    if ud.verified and ud.validity(rho, n, logF, udr):
        eps_ud = ud.value(rho, n, logF, udr)
        if eps_ud > 0 and bits(eps_ud) >= target_bits:
            result.update(feasible=True, delta_safe=udr, bound=ud.name,
                          m=None, eps_bits=bits(eps_ud),
                          regime="unique-decoding")

    # --- candidate 2: Bordage-Chiesa Johnson bound, sweep m ---
    bc = _bounds.get_bound("bordage-chiesa-johnson")

    def bc_eps_log2(m: int) -> float:
        # log2 of the BC error at trade parameter m (delta within window).
        d = max(1, (rho.numerator * n) // rho.denominator - 1)
        num_log2 = (math.log2(1.0 / 3.0)
                    + 7.0 * math.log2(m + 0.5)
                    + 2.0 * math.log2(n)
                    + math.log2(d)
                    - 1.5 * math.log2(float(rho)))
        return num_log2 - logF

    m_min = 3
    # Is even m=3 feasible?
    if bc_eps_log2(m_min) <= target_eps_log2:
        # binary search for the largest m with eps(m) <= target.
        lo, hi = m_min, m_max
        if bc_eps_log2(hi) <= target_eps_log2:
            best_m = hi
        else:
            while lo + 1 < hi:
                mid = (lo + hi) // 2
                if bc_eps_log2(mid) <= target_eps_log2:
                    lo = mid
                else:
                    hi = mid
            best_m = lo
        dmax = 1.0 - (1.0 + 1.0 / (2.0 * best_m)) * math.sqrt(float(rho))
        delta_bc = min(dmax, J)
        eps_bits_bc = -bc_eps_log2(best_m)
        if delta_bc > result["delta_safe"]:
            result.update(feasible=True, delta_safe=delta_bc,
                          bound=bc.name, m=best_m, eps_bits=eps_bits_bc,
                          regime="johnson")

    return result


# ---------------------------------------------------------------------------
# provable-unsafe delta (prime fields) and the full row
# ---------------------------------------------------------------------------

@dataclass
class DeltaStarRow:
    rho: str
    field: str
    field_type: str
    logF: float
    log2n: int
    n: int
    johnson: float
    capacity: float
    listdec_cap: float
    safe_feasible: bool
    delta_safe: float
    safe_bound: Optional[str]
    safe_m: Optional[int]
    safe_eps_bits: float
    t_at_safe: Optional[int]
    sec_bits_at_safe: float
    proof_kb_at_safe: Optional[float]
    delta_unsafe: Optional[float]   # field-agnostic Kambiré radius (prime & ext)
    gap_lo: float
    gap_hi: Optional[float]         # field-agnostic bracket upper end


def delta_star_row(rho: Fraction, field_name: str, logF: float,
                   field_type: str, log2n: int,
                   cost: CostModel = DEFAULT_COST,
                   target_bits: float = TARGET_BITS) -> DeltaStarRow:
    n = 1 << log2n
    J = johnson_radius(rho)
    cap = float(capacity(rho))
    # list-decoding capacity radius at q = 2^logF (use a cap so huge fields are
    # fine); this is the generic q-ary entropy radius, ~1/log2 q below Singleton.
    q = 2 ** min(int(round(logF)), 4096)
    ldc = list_decoding_capacity_radius(rho, q)

    safe = provable_safe_delta(rho, n, logF, target_bits=target_bits)

    t_at_safe = None
    sec_bits = float("-inf")
    proof_kb = None
    if safe["feasible"]:
        d = safe["delta_safe"]
        # Build the composite soundness floor from the VERIFIED certificate at
        # delta_safe.  The generic solver at m=1 would pick whichever registry
        # bound is cheapest, which can be the *unverified* BCHKS bound; for the
        # headline we instead use the proven ε_mca = 2^-safe["eps_bits"] (from
        # Bordage-Chiesa at its best trade m, or the verified UD bound), add the
        # (constant-VERIFY) interleaved list term at m=1, and bisect for the
        # minimum query count t with  floor + (1-δ)^t ≤ 2^-128.
        eps_mca = 2.0 ** (-safe["eps_bits"])
        eta = (johnson_radius(rho) - d) if d < johnson_radius(rho) else None
        lb = _bounds.best_listsize_bound(rho, n, logF, d, eta=eta, m=1,
                                         field_type=field_type)
        list_over_field = 0.0
        if lb is not None:
            ls = lb.value(rho, n, logF, d, eta=eta, m=1, field_type=field_type)
            list_over_field = min(1.0, 2.0 ** (math.log2(ls) - logF))
        floor = min(1.0, eps_mca + list_over_field)
        target = 2.0 ** (-target_bits)
        one_minus = 1.0 - d
        if floor <= target and one_minus > 0:
            slack = target - floor
            if slack <= 0:
                t_at_safe = 0
            else:
                t_at_safe = max(0, math.ceil(math.log(slack)
                                             / math.log(one_minus)))
            total = min(1.0, floor + one_minus ** t_at_safe)
            sec_bits = bits(total) if total > 0 else float("inf")
            proof_kb = cost.proof_size_bits(t_at_safe, n, logF) / 8 / 1024
        else:
            # ε_mca ≤ target by construction, so a floor above target can only
            # come from the list term |Λ|/|F| (or one_minus<=0). Report floor
            # bits; INF query count (no finite t helps).
            sec_bits = bits(floor) if floor > 0 else float("inf")

    # provable-unsafe delta: FIELD-AGNOSTIC Kambiré finite-field radius for ALL
    # field types (prime AND genuine odd-char extension) — the near-capacity
    # counterexample is a characteristic-zero cyclotomic invariant
    # (n2-verdict.md), so extensions inherit the same no-go as primes.
    # R13: (1-rho) - 2/s_max(b), s a power of two; None at 256-bit (no valid s
    # — the mechanism establishes no threshold ceiling there).
    delta_unsafe: Optional[float] = _bounds.kambire_unsafe_delta(rho, logF)

    return DeltaStarRow(
        rho=str(rho), field=field_name, field_type=field_type, logF=float(logF),
        log2n=log2n, n=n, johnson=J, capacity=cap, listdec_cap=ldc,
        safe_feasible=safe["feasible"], delta_safe=safe["delta_safe"],
        safe_bound=safe["bound"], safe_m=safe["m"],
        safe_eps_bits=safe["eps_bits"],
        t_at_safe=t_at_safe, sec_bits_at_safe=sec_bits, proof_kb_at_safe=proof_kb,
        delta_unsafe=delta_unsafe,
        gap_lo=safe["delta_safe"],
        gap_hi=delta_unsafe,
    )


def all_rows(rates=RATES, fields=FIELDS, log2ns=LOG2NS) -> List[DeltaStarRow]:
    rows = []
    for r in rates:
        for (fname, lf, ft) in fields:
            for ln in log2ns:
                rows.append(delta_star_row(r, fname, lf, ft, ln))
    return rows


# ---------------------------------------------------------------------------
# Markdown emitter
# ---------------------------------------------------------------------------

def _fmt(x, nd=5):
    if x is None:
        return "—"
    if isinstance(x, float) and (math.isinf(x) or math.isnan(x)):
        return "∞" if x > 0 else "−∞"
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def render_markdown(rows: List[DeltaStarRow]) -> str:
    L: List[str] = []
    L.append("# delta*_C tables — best-known brackets: "
             "proven floor + inferred ceiling")
    L.append("")
    L.append("Generated by `calculator/delta_star.py` "
             "(`python3 cli.py delta-star`). Target ε* = 2⁻¹²⁸. "
             "ABF δ-radius convention throughout.")
    L.append("")
    L.append("**Reading guide.** `delta_safe` = largest radius where the best "
             "*VERIFIED* positive MCA bound (Bordage-Chiesa Thm 9.2, or the "
             "unique-decoding bound n/|F| below the UD radius) certifies "
             "ε_mca ≤ 2⁻¹²⁸ (capped at the Johnson radius, the proven positive "
             "ceiling). `delta_unsafe` = the per-field **assembled-lemma** "
             "unsafe ceiling (R13): `(1−ρ) − 2/s_max(b)`, where the Kambiré "
             "quotient parameter `s` must be a **power of two** (R13 "
             "s-integrality) and `s_max(b)` is the largest such `s` passing "
             "(i) the KK25 distinctness calibration `p > φ(s)^{φ(s)}`, (ii) "
             "the prize-threshold count `N(s,ρ) > 2^{b−128}`, and (iv) the "
             "above-Johnson condition `s > 2/(√ρ−ρ)` — `s_max = 16/16/32` at "
             "`b = 31/64/128`; at `b = 256` **no valid s exists** (count "
             "`3³² = 2⁵¹ ≪ 2¹²⁸`), so **no threshold-established Kambiré-type "
             "ceiling below r_E exists** there (rendered `—`; the generic "
             "CS25/Elias ceiling at ≈ r_E still stands). Lemma status: "
             "**ASSEMBLED** (KK25 cited not re-proved; ρ = 1/2 "
             "N1-conditional); the pre-R13 continuum (1−ρ)−6/log₂|F| is an "
             "asymptotic reference only. This is **FIELD-AGNOSTIC**: it applies "
             "IDENTICALLY to **prime** and genuine odd-characteristic "
             "**extension** fields, because the near-capacity counterexample's "
             "bad-scalar count is a characteristic-zero cyclotomic invariant "
             "(`n2-verdict.md`) — so an extension row has the SAME δ_unsafe and "
             "the SAME bracket as the matching-size prime row. The open prize "
             "gap is `[delta_safe, delta_unsafe]` (at 256-bit: "
             "`[delta_safe, —)`, upper end not threshold-established). These "
             "are best-known brackets: **proven floor + assembled-lemma "
             "ceiling**, not a resolution of δ*_C.")
    L.append("")
    L.append("`t@safe` = query count of the full composite soundness "
             "(ε_mca + |Λ|/|F| + (1−δ)^t ≤ 2⁻¹²⁸) at `delta_safe`; "
             "`bits@safe` = achieved security bits there; INF = the "
             "t-independent floor already exceeds the target (no finite t "
             "helps — typically because the list term |Λ|/|F| or the certified "
             "ε_mca floor is too large for that field).")
    L.append("")

    # group by rho
    for r in sorted(set(row.rho for row in rows),
                    key=lambda s: Fraction(s), reverse=True):
        rr = Fraction(r)
        sub = [row for row in rows if row.rho == r]
        J = sub[0].johnson
        cap = sub[0].capacity
        L.append(f"## ρ = {r}  (Johnson 1−√ρ = {J:.5f}, "
                 f"Singleton capacity 1−ρ = {cap:.5f})")
        L.append("")
        L.append("| field | type | log₂\\|F\\| | n | δ_safe | safe via | "
                 "(BC m) | ε_mca bits | t@safe | bits@safe | proof KB | "
                 "δ_unsafe (field-agnostic, R13) | open gap |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for row in sub:
            via = row.safe_bound or "—"
            via = {"unique-decoding": "UD n/\\|F\\|",
                   "bordage-chiesa-johnson": "B-C Thm9.2"}.get(via, via)
            # delta_unsafe is field-agnostic (prime & extension); None means no
            # threshold-established ceiling (R13: no valid power-of-two s at
            # this field size — NOT the superseded "extension is OPEN" reading).
            unsafe = (_fmt(row.delta_unsafe) if row.delta_unsafe is not None
                      else "— (none established, R13)")
            if row.gap_hi is not None:
                gap = f"[{row.gap_lo:.4f}, {row.gap_hi:.4f}]"
            else:
                gap = f"[{row.gap_lo:.4f}, —)"
            t_s = (str(row.t_at_safe) if row.t_at_safe is not None
                   else ("INF" if row.safe_feasible else "—"))
            bsafe = (f"{row.sec_bits_at_safe:.1f}"
                     if row.safe_feasible and math.isfinite(row.sec_bits_at_safe)
                     else ("∞" if row.safe_feasible else "—"))
            pkb = (f"{row.proof_kb_at_safe:.1f}"
                   if row.proof_kb_at_safe is not None else "—")
            dsafe = (_fmt(row.delta_safe) if row.safe_feasible else "infeasible")
            ebits = (f"{row.safe_eps_bits:.1f}" if row.safe_feasible else "—")
            L.append(f"| {row.field} | {row.field_type} | {row.logF:.2f} | "
                     f"2^{row.log2n} | {dsafe} | {via} | "
                     f"{row.safe_m if row.safe_m else '—'} | {ebits} | "
                     f"{t_s} | {bsafe} | {pkb} | {unsafe} | {gap} |")
        L.append("")

    # headline block
    L.append("## Headline — ρ = 1/2 (the 64→128-bit prize question)")
    L.append("")
    L.append("At ρ = 1/2 the proven positive ceiling is the Johnson radius "
             "δ = 1−√(1/2) ≈ 0.29289. The single-query catch probability there "
             "is (1−δ) = √(1/2), so the per-query soundness factor is "
             "(1−δ)^t = (1/√2)^t and **128 queries give exactly "
             "(1/√2)¹²⁸ = 2⁻⁶⁴ — only 64 bits.** Reaching Singleton capacity "
             "δ → 1/2 would give (1/2)^t, i.e. (1/2)¹²⁸ = 2⁻¹²⁸ — the full "
             "128 bits. Closing Johnson→capacity at ρ = 1/2 is the headline "
             "64→128-bit jump.")
    L.append("")
    # show the (1-delta)^t exponent table for the safe delta at a couple fields
    L.append("| field | n | δ_safe | (1−δ_safe) | bits/query −log₂(1−δ_safe) | "
             "t for 128 bits via query term alone |")
    L.append("|---|---|---|---|---|---|")
    half = Fraction(1, 2)
    for row in [x for x in rows if x.rho == "1/2"
                and x.field in ("Mersenne31", "prime256")
                and x.log2n in (20, 30)]:
        d = row.delta_safe if row.safe_feasible else row.johnson
        one_minus = 1.0 - d
        bpq = -math.log2(one_minus)
        t_q = math.ceil(128.0 / bpq) if bpq > 0 else float("inf")
        note = "" if row.safe_feasible else " (δ_safe infeasible; using Johnson)"
        L.append(f"| {row.field} | 2^{row.log2n} | {d:.5f}{note} | "
                 f"{one_minus:.5f} | {bpq:.4f} | {t_q} |")
    L.append("")
    L.append("> The query-term column is the *best case* (it ignores the "
             "ε_mca + |Λ|/|F| floor). Over small fields (Mersenne31) the "
             "certified ε_mca floor exceeds 2⁻¹²⁸, so the composite is "
             "infeasible regardless of t — see the per-ρ tables. Only ~256-bit "
             "fields certify the Johnson radius at 2⁻¹²⁸ from the single-code "
             "MCA bound.")
    L.append("")
    L.append("### Why small fields are infeasible (the honest punchline)")
    L.append("")
    L.append("The verified Bordage-Chiesa MCA error is "
             "(m+½)⁷·n²·d / (3ρ^{3/2}·|F|). The n²·d numerator alone costs "
             "≈ 2·log₂n + log₂(ρn) bits (≈ 60 bits at n = 2²⁰, ≈ 90 bits at "
             "n = 2³⁰). To get ε_mca ≤ 2⁻¹²⁸ you therefore need "
             "log₂|F| ≳ 128 + 2·log₂n + log₂(ρn) + 7·log₂(m+½) − (3/2)log₂ρ, "
             "i.e. **roughly a ≥ 256-bit field** for n up to 2³⁰. Mersenne31 "
             "(31-bit), BabyBear (31-bit) and Goldilocks (64-bit) cannot "
             "certify ε_mca ≤ 2⁻¹²⁸ from this single-code bound at all — "
             "deployments over those fields reach 128-bit soundness via the "
             "(1−δ)^t query term plus protocol repetition / proof-of-work and "
             "by working in a large extension for the soundness argument, NOT "
             "from a single-code 2⁻¹²⁸ MCA certificate.")
    L.append("")
    return "\n".join(L)


def write_tables(path: Optional[str] = None) -> str:
    os.makedirs(OUT_DIR, exist_ok=True)
    if path is None:
        path = os.path.join(OUT_DIR, "delta_star_tables.md")
    rows = all_rows()
    md = render_markdown(rows)
    with open(path, "w") as fh:
        fh.write(md)
    return os.path.abspath(path)


__all__ = [
    "DeltaStarRow", "provable_safe_delta", "delta_star_row", "all_rows",
    "render_markdown", "write_tables", "FIELDS", "RATES", "LOG2NS",
]


if __name__ == "__main__":  # pragma: no cover
    p = write_tables()
    print(f"wrote {p}")
