"""
test_calculator.py
==================

Tests for the Proximity Prize parameter calculator.

Runs under pytest::

    python3 -m pytest tests/test_calculator.py -q

or as a plain script (a minimal __main__ harness collects + runs every
``test_*`` function and prints a single PASS/FAIL summary line)::

    python3 tests/test_calculator.py

The tests assert KNOWN ANCHORS (the radii are mathematical facts, independent
of any literature constant) plus structural properties of the soundness solver
(monotonicity, the bisection anchor) that hold for any sensible bound model.
"""

from __future__ import annotations

import math
import os
import sys
from fractions import Fraction

# Make the calculator package importable whether run via pytest from repo root
# or directly as a script from the tests/ dir.
_HERE = os.path.dirname(os.path.abspath(__file__))
_CALC = os.path.dirname(_HERE)
for p in (_CALC, _HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import bounds as _bounds  # noqa: E402
from proximity_parameters import (  # noqa: E402
    bits,
    capacity,
    johnson_radius,
    johnson_radius_exact,
    johnson_radius_rational,
    log2,
    log2_pow2_minus_c,
    qary_entropy,
    list_decoding_capacity_radius,
    unique_decoding_radius,
)
from soundness import (  # noqa: E402
    TARGET_ERROR_DEFAULT,
    evaluate_point,
    min_query_count,
    soundness_error,
    soundness_terms,
)

RATES = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 16)]
TOL = 1e-12


# ===========================================================================
# Radii anchors
# ===========================================================================

def test_unique_decoding_anchors():
    assert unique_decoding_radius(Fraction(1, 2)) == Fraction(1, 4)
    assert unique_decoding_radius(Fraction(1, 4)) == Fraction(3, 8)
    assert unique_decoding_radius(Fraction(1, 8)) == Fraction(7, 16)
    assert unique_decoding_radius(Fraction(1, 16)) == Fraction(15, 32)


def test_capacity_anchors():
    assert capacity(Fraction(1, 2)) == Fraction(1, 2)
    assert capacity(Fraction(1, 4)) == Fraction(3, 4)
    assert capacity(Fraction(1, 8)) == Fraction(7, 8)
    assert capacity(Fraction(1, 16)) == Fraction(15, 16)


def test_johnson_anchors_numeric():
    # 1 - sqrt(rho) for the four rates, within 1e-12.
    expected = {
        Fraction(1, 2): 1 - math.sqrt(0.5),
        Fraction(1, 4): 0.5,            # exact: 1 - 1/2
        Fraction(1, 8): 1 - math.sqrt(0.125),
        Fraction(1, 16): 0.75,          # exact: 1 - 1/4
    }
    for r, exp in expected.items():
        assert abs(johnson_radius(r) - exp) < TOL
        assert abs(float(johnson_radius_rational(r)) - exp) < TOL


def test_johnson_exact_perfect_squares():
    # rho a perfect-square reciprocal => Johnson is exactly rational.
    assert johnson_radius_exact(Fraction(1, 4)) == Fraction(1, 2)
    assert johnson_radius_exact(Fraction(1, 16)) == Fraction(3, 4)
    assert johnson_radius_rational(Fraction(1, 4)) == Fraction(1, 2)
    assert johnson_radius_rational(Fraction(1, 16)) == Fraction(3, 4)


def test_radii_ordering():
    # (1-rho)/2 < 1 - sqrt(rho) < 1 - rho  for every rho in (0,1).
    for r in RATES:
        ud = float(unique_decoding_radius(r))
        jr = johnson_radius(r)
        cap = float(capacity(r))
        assert ud < jr < cap, f"ordering failed at rho={r}: {ud} {jr} {cap}"


def test_radii_ordering_random_rates():
    # Strict ordering should hold for arbitrary rates, not just dyadic ones.
    for num, den in [(1, 3), (2, 5), (3, 7), (5, 11), (9, 10), (1, 100)]:
        r = Fraction(num, den)
        ud = float(unique_decoding_radius(r))
        jr = johnson_radius(r)
        cap = float(capacity(r))
        assert ud < jr < cap, f"ordering failed at rho={r}"


# ===========================================================================
# q-ary entropy
# ===========================================================================

def test_qary_entropy_zero():
    for q in (2, 3, 16, 256, 2 ** 31):
        assert qary_entropy(0, q) == 0.0


def test_qary_entropy_binary_peak():
    # H_2(1/2) = 1 exactly (max of binary entropy).
    assert abs(qary_entropy(Fraction(1, 2), 2) - 1.0) < TOL


def test_qary_entropy_peak_at_one_minus_inv_q():
    # H_q is maximized at x = 1 - 1/q with value 1.
    for q in (2, 4, 8):
        x = Fraction(q - 1, q)
        assert abs(qary_entropy(x, q) - 1.0) < 1e-9


def test_qary_entropy_monotone_increasing_on_left_branch():
    # H_q strictly increasing on [0, 1 - 1/q].
    for q in (2, 4, 16):
        xmax = 1.0 - 1.0 / q
        xs = [i / 200 * xmax for i in range(201)]
        vals = [qary_entropy(Fraction(x).limit_denominator(10 ** 9), q) for x in xs]
        for a, b in zip(vals, vals[1:]):
            assert b >= a - 1e-12, f"H_{q} not monotone increasing on left branch"


def test_qary_capacity_radius_tends_to_capacity_for_large_q():
    # As q -> infinity, the q-ary list-dec capacity radius -> 1 - rho, but the
    # convergence is only O(1/log2 q): the residual H_q(delta) carries a
    # delta*log_q(q-1) term that decays like 1/log2(q).  We therefore (a) check
    # the radius is below capacity and within a logq-matched tolerance, and
    # (b) check the gap *halves* when log2(q) doubles (the 1/log2 q rate).
    for r in RATES:
        cap = float(capacity(r))
        logq = 60
        radius = list_decoding_capacity_radius(r, 2 ** logq)
        assert radius < cap + 1e-9
        # gap ~ C / log2(q); at logq=60 expect gap < ~1.5/60 ~ 0.025 (loose).
        assert abs(radius - cap) < 2.0 / logq, (
            f"capacity radius {radius} too far from {cap} at q=2^{logq}")


def test_qary_capacity_radius_convergence_rate():
    # Gap to capacity halves as log2(q) doubles (the O(1/log2 q) rate).
    r = Fraction(1, 2)
    cap = float(capacity(r))
    g1 = abs(list_decoding_capacity_radius(r, 2 ** 40) - cap)
    g2 = abs(list_decoding_capacity_radius(r, 2 ** 80) - cap)
    # g2 should be about half of g1.
    assert 0.4 < g2 / g1 < 0.6, f"convergence rate off: g1={g1} g2={g2}"


# ===========================================================================
# log2 / bits helpers
# ===========================================================================

def test_log2_fraction_exact():
    assert abs(log2(Fraction(1, 8)) - (-3.0)) < TOL
    assert abs(log2(Fraction(1024)) - 10.0) < TOL


def test_bits_of_target():
    assert abs(bits(Fraction(1, 2 ** 128)) - 128.0) < 1e-9


def test_log2_pow2_minus_c_cancellation_safe():
    # log2(2^31 - 1) should be just below 31 and *not* exactly 31.
    v = log2_pow2_minus_c(31, 1)
    assert 30.999999 < v < 31.0
    # Goldilocks 2^64 - 2^32 + 1: log2 just below 64.
    g = log2_pow2_minus_c(64, 2 ** 32 - 1)
    assert 63.99 < g < 64.0
    # c = 0 gives exactly the exponent.
    assert log2_pow2_minus_c(31, 0) == 31.0


# ===========================================================================
# bounds registry sanity
# ===========================================================================

def test_registry_has_core_bounds():
    names = {b.name for b in _bounds.all_bounds()}
    for expected in ("unique-decoding", "bordage-chiesa-johnson",
                     "bchks-johnson",
                     "unknown-beyond-johnson", "proven-near-capacity-nogo",
                     "interleave-mca-union",
                     "interleave-listsize"):
        assert expected in names
    # The capacity region is now modelled by TWO distinct bounds (the OPEN band
    # and the PROVEN near-capacity no-go). The earlier SINGLE broad
    # "capacity-nogo" — which mislabelled the whole window as a proven no-go —
    # is superseded and must no longer be registered. (Nor the even-earlier
    # prime/extension split.)
    assert "capacity-nogo" not in names
    assert "capacity-prime-nogo" not in names
    assert "capacity-extension-open" not in names


def test_ud_bound_applies_below_udr_only():
    rho = Fraction(1, 4)
    udr = float(unique_decoding_radius(rho))
    b = _bounds.get_bound("unique-decoding")
    assert b.validity(rho, 1 << 20, 128, udr - 0.01)       # inside
    assert not b.validity(rho, 1 << 20, 128, udr + 0.01)   # outside


def test_johnson_bound_applies_below_johnson_with_slack():
    rho = Fraction(1, 4)
    J = johnson_radius(rho)
    # BCHKS Thm 1.5 is valid right up to J - eta.
    b = _bounds.get_bound("bchks-johnson")
    # delta = J - 0.05, slack eta = 0.05 -> applies.
    assert b.validity(rho, 1 << 20, 128, J - 0.05, eta=0.05)
    # at/above Johnson radius -> does not apply.
    assert not b.validity(rho, 1 << 20, 128, J + 0.001, eta=0.05)


def test_bordage_chiesa_window_widens_with_m():
    # Bordage-Chiesa validity ceiling 1 - (1 + 1/(2m)) sqrt(rho) increases with
    # m toward the Johnson radius; a delta just below Johnson needs large m.
    rho = Fraction(1, 4)
    J = johnson_radius(rho)
    b = _bounds.get_bound("bordage-chiesa-johnson")
    delta = J - 0.005  # very close to Johnson
    # small m (3): window is well below Johnson -> does NOT reach this delta.
    assert not b.validity(rho, 1 << 20, 256, delta, eta=0.005, m=3)
    # large m: window approaches Johnson -> reaches this delta.
    assert b.validity(rho, 1 << 20, 256, delta, eta=0.005, m=5000)
    # error is monotone increasing in m (the (m+1/2)^7 inflation).
    e_small = b.value(rho, 1 << 20, 256, J - 0.1, eta=0.1, m=3)
    e_large = b.value(rho, 1 << 20, 256, J - 0.1, eta=0.1, m=300)
    assert e_large > e_small


def test_open_band_is_unknown_not_a_proven_nogo():
    # KEY SEMANTIC FIX. In the OPEN band (J, r_E) the applicable capacity-region
    # MCA bound is `unknown-beyond-johnson`: vacuous 1.0 (cannot certify) but
    # NOT a proven no-go (verified=False). The proven no-go must NOT fire here.
    rho = Fraction(1, 4)
    J = johnson_radius(rho)
    logF = 256.0
    r_E = _bounds.nogo_split_radius(rho, logF)
    assert J < r_E  # there is a genuine open band to test in
    mid_open = (J + r_E) / 2
    ub = _bounds.get_bound("unknown-beyond-johnson")
    pn = _bounds.get_bound("proven-near-capacity-nogo")
    for ft in ("prime", "extension", None):
        # the OPEN bound applies and is vacuous (1.0)
        assert ub.validity(rho, 1 << 20, logF, mid_open, field_type=ft)
        assert ub.value(rho, 1 << 20, logF, mid_open, field_type=ft) == 1.0
        # the PROVEN no-go does NOT apply in the open band
        assert not pn.validity(rho, 1 << 20, logF, mid_open, field_type=ft)
    # It is NOT a proven no-go: verified is False (an honest "don't know").
    assert ub.verified is False
    # The single broad "capacity-nogo" is gone (replaced by this open bound).
    names = {b.name for b in _bounds.all_bounds()}
    assert "capacity-nogo" not in names
    # And the applicable-mca selector in the open band returns the OPEN bound,
    # not a proven no-go.
    app = [b.name for b in _bounds.applicable(rho, 1 << 20, logF, mid_open,
                                              kind="mca")
           if b.regime == "capacity"]
    assert app == ["unknown-beyond-johnson"]


def test_near_capacity_is_a_proven_nogo():
    # At delta >= r_E the applicable capacity-region MCA bound is
    # `proven-near-capacity-nogo`: vacuous 1.0 AND verified=True (a genuine
    # proven failure). The OPEN bound must NOT fire here.
    rho = Fraction(1, 4)
    cap = float(capacity(rho))
    logF = 256.0
    r_E = _bounds.nogo_split_radius(rho, logF)
    assert r_E < cap
    mid_proven = (r_E + cap) / 2
    ub = _bounds.get_bound("unknown-beyond-johnson")
    pn = _bounds.get_bound("proven-near-capacity-nogo")
    for ft in ("prime", "extension", None):
        assert pn.validity(rho, 1 << 20, logF, mid_proven, field_type=ft)
        assert pn.value(rho, 1 << 20, logF, mid_proven, field_type=ft) == 1.0
        # the OPEN bound does NOT apply at/above the ceiling
        assert not ub.validity(rho, 1 << 20, logF, mid_proven, field_type=ft)
    # It IS a proven no-go.
    assert pn.verified is True
    # right AT the ceiling r_E the proven no-go takes over (boundary belongs to
    # delta >= r_E), and the open band excludes it -> exactly one applies.
    assert pn.validity(rho, 1 << 20, logF, r_E)
    assert not ub.validity(rho, 1 << 20, logF, r_E)
    app = [b.name for b in _bounds.applicable(rho, 1 << 20, logF, mid_proven,
                                              kind="mca")
           if b.regime == "capacity"]
    assert app == ["proven-near-capacity-nogo"]


def test_capacity_region_split_is_field_agnostic():
    # Both new bounds are FIELD-AGNOSTIC (n2-verdict.md): prime, extension and
    # unknown field types give IDENTICAL validity for each, in both bands.
    rho = Fraction(1, 4)
    cap = float(capacity(rho))
    J = johnson_radius(rho)
    logF = 256.0
    r_E = _bounds.nogo_split_radius(rho, logF)
    mid_open = (J + r_E) / 2
    mid_proven = (r_E + cap) / 2
    ub = _bounds.get_bound("unknown-beyond-johnson")
    pn = _bounds.get_bound("proven-near-capacity-nogo")
    for bnd, d in ((ub, mid_open), (pn, mid_proven)):
        v_prime = bnd.validity(rho, 1 << 20, logF, d, field_type="prime")
        v_ext = bnd.validity(rho, 1 << 20, logF, d, field_type="extension")
        v_none = bnd.validity(rho, 1 << 20, logF, d)
        assert v_prime == v_ext == v_none == True
    # And neither superseded bound (the single broad no-go, or the old
    # extension-OPEN) exists.
    names = {b.name for b in _bounds.all_bounds()}
    assert "capacity-nogo" not in names
    assert "capacity-extension-open" not in names


def test_capacity_region_partitions_window_exactly():
    # Across the whole (Johnson, capacity) window, EXACTLY one capacity-region
    # MCA bound applies at every delta — no overlap, no gap, including at the
    # r_E seam. (The split must tile the window the old single bound covered.)
    for rho in (Fraction(1, 2), Fraction(1, 4), Fraction(1, 8)):
        J = johnson_radius(rho)
        cap = float(capacity(rho))
        logF = 256.0
        r_E = _bounds.nogo_split_radius(rho, logF)
        pts = [J + (cap - J) * i / 500 for i in range(501)]
        pts += [r_E - 1e-9, r_E, r_E + 1e-9]
        for d in pts:
            if not (J <= d < cap):
                continue
            app = [b for b in _bounds.applicable(rho, 1 << 20, logF, d,
                                                 kind="mca")
                   if b.regime == "capacity"]
            assert len(app) == 1, (
                f"expected exactly 1 capacity bound at rho={rho} "
                f"delta={d}, got {[b.name for b in app]}")


def test_kambire_unsafe_delta_below_capacity():
    # R13: the per-field established-unsafe radius is (1-rho) - 2/s_max(b) with
    # s_max a POWER OF TWO (s_max = 16/16/32 at b = 31/64/128), strictly below
    # capacity and strictly above Johnson (hypothesis (iv)).  At b = 256 no
    # valid s exists (count 3^32 = 2^51 << 2^128) -> None (no threshold-
    # established ceiling; the generic CS25/Elias ceiling takes over as the
    # band split, see nogo_split_radius).
    rho = Fraction(1, 2)
    cap = float(capacity(rho))
    J = johnson_radius(rho)
    for logF, s_expected in ((31.0, 16), (64.0, 16), (128.0, 32)):
        assert _bounds.kambire_smax(rho, logF) == s_expected
        d = _bounds.kambire_unsafe_delta(rho, logF)
        assert J < d < cap
        assert abs(d - (cap - 2.0 / s_expected)) < 1e-12
    assert _bounds.kambire_smax(rho, 256.0) is None
    assert _bounds.kambire_unsafe_delta(rho, 256.0) is None
    # The band split falls back to the generic CS25/Elias ceiling at 256-bit
    # and equals the assembled-lemma value where the lemma fires.
    split_256 = _bounds.nogo_split_radius(rho, 256.0)
    assert J < split_256 < cap
    assert abs(split_256 - 0.49609) < 5e-4  # ~ H_q^{-1}(1-rho) at q = 2^256
    assert abs(_bounds.nogo_split_radius(rho, 128.0) - 0.4375) < 1e-12


def test_no_infinite_recursion_in_meta_bound():
    # Regression: the interleaving meta-bound must not recurse into itself.
    rho = Fraction(1, 4)
    J = johnson_radius(rho)
    sel = _bounds.interleaved_mca(rho, 1 << 20, 256, J - 0.05, eta=0.05, m=8)
    assert sel is not None
    # m>1 should select the composition meta-bound.
    assert sel.name == "interleave-mca-union"
    v = sel.value(rho, 1 << 20, 256, J - 0.05, eta=0.05, m=8)
    assert 0.0 < v <= 1.0


def test_interleave_scales_error_with_m():
    rho = Fraction(1, 4)
    J = johnson_radius(rho)
    delta = J - 0.05
    e1 = soundness_terms(rho, 256, 1 << 20, delta, m=1, t=0).eps_mca
    e4 = soundness_terms(rho, 256, 1 << 20, delta, m=4, t=0).eps_mca
    # union bound -> linear in m (within float tolerance).
    assert abs(e4 - 4 * e1) < 1e-9 * max(e4, 4 * e1) + 1e-300


# ===========================================================================
# soundness solver: monotonicity + bisection anchor
# ===========================================================================

def test_t_monotone_nonincreasing_in_delta():
    # More slack (larger delta) => fewer queries needed, over a field big
    # enough that the floor stays below target throughout.
    rho = Fraction(1, 4)
    logF = 256
    ts = []
    deltas = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.49]
    for d in deltas:
        sol = min_query_count(rho, logF, 1 << 20, d, m=1)
        assert sol.feasible, f"expected feasible at delta={d}"
        ts.append(sol.t)
    for a, b in zip(ts, ts[1:]):
        assert b <= a, f"t not non-increasing in delta: {list(zip(deltas, ts))}"


def test_t_increases_as_field_shrinks():
    # Tighter field (smaller logF) => larger floor & weaker query budget =>
    # t should not decrease (and typically increases) for fixed delta, as long
    # as both remain feasible.
    rho = Fraction(1, 8)
    delta = johnson_radius(rho) - 0.05
    sol_big = min_query_count(rho, 400, 1 << 20, delta, m=1)
    sol_small = min_query_count(rho, 300, 1 << 20, delta, m=1)
    assert sol_big.feasible and sol_small.feasible
    assert sol_small.t >= sol_big.t


def test_t_increases_as_target_tightens():
    # Tighter target (more bits) => more queries, for a fixed feasible point.
    rho = Fraction(1, 8)
    delta = johnson_radius(rho) - 0.05
    logF = 400
    sol_128 = min_query_count(rho, logF, 1 << 20, delta, m=1,
                              target_error=Fraction(1, 2 ** 128))
    sol_200 = min_query_count(rho, logF, 1 << 20, delta, m=1,
                              target_error=Fraction(1, 2 ** 200))
    assert sol_128.feasible and sol_200.feasible
    assert sol_200.t >= sol_128.t


def test_bisection_anchor():
    # The solver must return t with error(t) <= target < error(t-1).
    rho = Fraction(1, 4)
    logF = 256
    delta = 0.40
    target = TARGET_ERROR_DEFAULT
    sol = min_query_count(rho, logF, 1 << 20, delta, m=1, target_error=target)
    assert sol.feasible and sol.t is not None and sol.t >= 1
    t = sol.t
    err_t = soundness_error(rho, logF, 1 << 20, delta, 1, t)
    err_tm1 = soundness_error(rho, logF, 1 << 20, delta, 1, t - 1)
    assert err_t <= float(target), f"error(t)={err_t} > target"
    assert err_tm1 > float(target), f"error(t-1)={err_tm1} <= target (not minimal)"


def test_infeasible_when_floor_exceeds_target():
    # Small field where the t-independent floor already exceeds the target:
    # the solver must report infeasible (no finite t helps).
    rho = Fraction(1, 4)
    delta = johnson_radius(rho) - 0.05
    sol = min_query_count(rho, 64, 1 << 20, delta, m=1)   # Goldilocks-ish
    assert not sol.feasible
    assert sol.t is None
    assert sol.floor_bits < 128.0


def test_capacity_window_is_infeasible():
    # In the Johnson->capacity no-go window, the MCA floor is vacuous (1.0),
    # so no query count can reach the target regardless of field size.
    rho = Fraction(1, 4)
    J = johnson_radius(rho)
    cap = float(capacity(rho))
    delta = (J + cap) / 2
    sol = min_query_count(rho, 1024, 1 << 20, delta, m=1)
    assert not sol.feasible


def test_query_term_decreases_with_t():
    # The composite total is strictly decreasing in t down to the floor.
    rho = Fraction(1, 4)
    logF = 256
    delta = 0.40
    prev = None
    for t in (0, 10, 50, 100, 200):
        tot = soundness_error(rho, logF, 1 << 20, delta, 1, t)
        if prev is not None:
            assert tot <= prev + 1e-18
        prev = tot


def test_proof_size_increases_with_t_and_n():
    # Proof size grows with query count and with block length (deeper Merkle).
    from soundness import CostModel
    cost = CostModel()
    p_small_t = cost.proof_size_bits(50, 1 << 20, 128)
    p_large_t = cost.proof_size_bits(100, 1 << 20, 128)
    assert p_large_t > p_small_t
    p_small_n = cost.proof_size_bits(50, 1 << 20, 128)
    p_large_n = cost.proof_size_bits(50, 1 << 28, 128)
    assert p_large_n > p_small_n


def test_evaluate_point_feasible_has_proof_size():
    rho = Fraction(1, 4)
    row = evaluate_point(rho, logF=256, log2n=20, delta=0.40, m=1)
    assert row.feasible
    assert row.t is not None and row.t >= 1
    assert row.proof_size_kb is not None and row.proof_size_kb > 0


# ===========================================================================
# delta*_C provable brackets (the headline deliverable)
# ===========================================================================

def test_delta_star_small_field_infeasible_big_field_reaches_johnson():
    import delta_star as ds
    rho = Fraction(1, 2)
    J = johnson_radius(rho)
    # Mersenne31 (31-bit): single-code MCA certificate cannot reach 2^-128.
    safe_m31 = ds.provable_safe_delta(rho, 1 << 20, log2_pow2_minus_c(31, 1))
    assert not safe_m31["feasible"]
    # 256-bit field: reaches essentially the Johnson radius via Bordage-Chiesa.
    safe_256 = ds.provable_safe_delta(rho, 1 << 20, 256.0)
    assert safe_256["feasible"]
    assert safe_256["bound"] == "bordage-chiesa-johnson"
    assert safe_256["delta_safe"] <= J + 1e-9          # never exceeds Johnson
    assert safe_256["delta_safe"] > J - 0.01           # but reaches close to it
    assert safe_256["eps_bits"] >= 128.0 - 1e-6        # certifies the target


def test_delta_star_prime_and_extension_have_identical_unsafe():
    # FIELD-AGNOSTIC (n2-verdict.md): extension delta_unsafe == prime
    # delta_unsafe at matched field size. The earlier "extension is OPEN"
    # behavior is SUPERSEDED. R13: at 256-bit NO threshold-established ceiling
    # exists (no valid power-of-two s) — both prime and extension report None;
    # at 128-bit both report the assembled-lemma (1-rho) - 2/32 = 0.4375.
    import delta_star as ds
    rho = Fraction(1, 2)
    cap = float(capacity(rho))
    prime_row = ds.delta_star_row(rho, "prime256", 256.0, "prime", 20)
    ext_row = ds.delta_star_row(rho, "ext256", 256.0, "extension", 20)
    # R13: 256-bit has no threshold-established ceiling — identically for both.
    assert prime_row.delta_unsafe is None
    assert ext_row.delta_unsafe is None
    # Same provable-safe delta too (the positive bound is field-type agnostic).
    assert abs(prime_row.delta_safe - ext_row.delta_safe) < 1e-12
    # At 128-bit the assembled lemma fires (s_max = 32) — identically for both.
    prime_128 = ds.delta_star_row(rho, "prime128", 128.0, "prime", 20)
    ext_128 = ds.delta_star_row(rho, "Goldilocks^2(ext)", 128.0, "extension", 20)
    assert prime_128.delta_unsafe is not None
    assert ext_128.delta_unsafe is not None
    assert abs(prime_128.delta_unsafe - (cap - 2.0 / 32.0)) < 1e-12
    assert abs(ext_128.delta_unsafe - prime_128.delta_unsafe) < 1e-12
    assert ext_128.delta_unsafe < cap


def test_delta_star_prime_and_extension_have_identical_bracket():
    # Explicit: prime and extension produce the IDENTICAL [delta_safe,
    # delta_unsafe] bracket at every target rate, over a matched-size field.
    # R13: at 256-bit the upper end is None (no threshold-established ceiling)
    # for BOTH; at 128-bit both report the assembled-lemma (1-rho) - 2/32.
    import delta_star as ds
    for rho in RATES:
        cap = float(capacity(rho))
        prime_row = ds.delta_star_row(rho, "prime256", 256.0, "prime", 20)
        ext_row = ds.delta_star_row(rho, "ext256", 256.0, "extension", 20)
        # R13: identical (absent) unsafe radius at 256-bit.
        assert prime_row.delta_unsafe is None
        assert ext_row.delta_unsafe is None
        assert prime_row.gap_hi is None and ext_row.gap_hi is None
        # Identical lower bracket endpoint.
        assert abs(prime_row.gap_lo - ext_row.gap_lo) < 1e-12
        # Identical (present) unsafe radius at 128-bit.
        prime_128 = ds.delta_star_row(rho, "prime128", 128.0, "prime", 20)
        ext_128 = ds.delta_star_row(rho, "Goldilocks^2(ext)", 128.0,
                                    "extension", 20)
        assert prime_128.delta_unsafe is not None
        assert ext_128.delta_unsafe is not None
        assert abs(prime_128.delta_unsafe - (cap - 2.0 / 32.0)) < 1e-12
        assert abs(ext_128.delta_unsafe - prime_128.delta_unsafe) < 1e-12


def test_delta_star_unsafe_ground_truth_values():
    # Ground-truth delta_unsafe = (1-rho) - 2/s_max(b) (R13 assembled lemma,
    # s a POWER OF TWO; s_max = 16/16/32 at b = 31/64/128, None at 256-bit),
    # applied to BOTH prime and extension. Spot-check the headline values from
    # problem-ledger.md / technical-note.md (R13).
    import delta_star as ds
    # rho=1/2 across field sizes (prime values, now also extension):
    cases_half = {
        128.0: 0.4375,    # 128-bit (s_max = 32; old continuum 0.45313)
        log2_pow2_minus_c(64, 2 ** 32 - 1): 0.375,  # Goldilocks (s_max = 16)
        log2_pow2_minus_c(31, 1): 0.375,            # M31 31-bit (s_max = 16)
        log2_pow2_minus_c(31, 2 ** 27): 0.375,      # BabyBear (s_max = 16)
    }
    for logF, expected in cases_half.items():
        d = _bounds.kambire_unsafe_delta(Fraction(1, 2), logF)
        assert abs(d - expected) < 1e-12, f"rho=1/2 logF={logF}: {d} != {expected}"
    # R13: 256-bit — no threshold-established ceiling at ANY rate.
    for rho in RATES:
        assert _bounds.kambire_unsafe_delta(rho, 256.0) is None
    # rho in {1/4,1/8,1/16} at 128-bit (s_max = 32):
    cases_128 = {
        Fraction(1, 4): 0.6875,
        Fraction(1, 8): 0.8125,
        Fraction(1, 16): 0.875,
    }
    for rho, expected in cases_128.items():
        d = _bounds.kambire_unsafe_delta(rho, 128.0)
        assert abs(d - expected) < 1e-12, f"rho={rho} @128: {d} != {expected}"


def test_delta_star_johnson_query_term_is_64_bits_at_128_queries():
    # The crux of the headline: at rho=1/2, delta=Johnson, (1-delta)^128 = 2^-64.
    rho = Fraction(1, 2)
    J = johnson_radius(rho)
    one_minus = 1.0 - J               # = sqrt(1/2)
    assert abs(one_minus - math.sqrt(0.5)) < 1e-12
    assert abs(-128.0 * math.log2(one_minus) - 64.0) < 1e-9
    # At capacity delta=1/2 it would be 128 bits.
    assert abs(-128.0 * math.log2(1.0 - 0.5) - 128.0) < 1e-9


# ===========================================================================
# plain-script fallback runner
# ===========================================================================

def _run_all() -> int:
    """Collect and run every test_* in this module; print one summary line."""
    g = globals()
    tests = sorted(name for name in g
                   if name.startswith("test_") and callable(g[name]))
    failures = []
    for name in tests:
        try:
            g[name]()
        except Exception as exc:  # noqa: BLE001
            failures.append((name, repr(exc)))
    n = len(tests)
    if failures:
        for name, err in failures:
            print(f"FAIL  {name}: {err}")
        print(f"\n{len(failures)}/{n} FAILED, {n - len(failures)}/{n} passed")
        return 1
    print(f"PASSED {n}/{n} tests")
    return 0


if __name__ == "__main__":
    raise SystemExit(_run_all())
