"""
bounds.py
=========

A **registry** of known Reed-Solomon proximity-gap / mutual-correlated-agreement
(MCA) error bounds and interleaved list-size bounds, with the literature
constants now **filled from the Wave-1 extraction** (see ``../literature/notes``)
and each one tagged ``verified=True`` (traceable to a source note) or
``verified=False`` (a clearly-flagged ``VERIFY`` placeholder for a hidden /
not-yet-pinned constant).

Convention
----------
We carry the **ABF δ-radius convention everywhere**: ``delta`` is the relative
proximity *radius* (fraction of disagreeing coordinates), ``rho = k/n`` is the
rate, the Johnson radius is ``J = 1 - sqrt(rho)`` and Singleton capacity is
``1 - rho``.  BCHKS (eprint 2025/2055) uses the OPPOSITE convention (their
``delta`` = code minimum distance ``1 - rho``, their ``gamma`` = radius); we
translate its statements into the ABF convention here and document the clash in
``../problem-ledger.md``.

Design goals
------------
* Every bound is a self-describing :class:`Bound` object carrying metadata
  (name, source tag, the formula as a string, a ``validity`` predicate, a
  ``value`` callable, a ``verified`` flag and free-text notes).
* Bounds are looked up by *regime* (unique-decoding / Johnson / capacity) and by
  *kind* (``"mca"`` for a single code, ``"list"`` for an interleaved list size).
  ``soundness.py`` consumes them purely through this interface.
* A ledger-driven script can patch the module-level constants in one place
  (``bounds.CONSTANTS`` or the module globals) **without changing any call
  site**; the validity windows and value callables read the globals live.

Parameter glossary (shared signature ``(rho, n, logF, delta, eta, m)`` plus an
optional ``field_type``)
----------------------------------------------------------------------
* ``rho``        : code rate ``k/n`` (Fraction or float), in (0, 1).
* ``n``          : block length ``|L|`` (int).
* ``logF``       : ``log2|F|`` of the field (float).  We always carry the *log*
                   of the field size, never ``|F|`` itself, so 256-bit fields
                   are fine.
* ``delta``      : relative proximity radius being tested (float in (0, 1)).
* ``eta``        : Johnson slack, ``delta = (1 - sqrt(rho)) - eta`` with
                   ``eta > 0``.  Only meaningful for the Johnson-regime bounds;
                   pass ``None`` else.
* ``m``          : interleaving / batch parameter.  For the *positive* RS bounds
                   it doubles as the Bordage-Chiesa integer trade parameter
                   ``m >= 3`` (larger ``m`` => closer to Johnson, worse error);
                   for the interleaving composition it is the number of
                   interleaved instances.  (These two roles coincide in the
                   protocol use: a single integer parameter is swept.)
* ``field_type`` : ``"prime"`` | ``"extension"`` | ``None`` (unknown).  Retained
                   in the signatures for back-compat, but the capacity-region
                   bounds are now **FIELD-AGNOSTIC**: the smooth-domain
                   near-capacity counterexample is a characteristic-zero
                   cyclotomic invariant (``n2-verdict.md``), so genuine
                   odd-characteristic extension fields ``GF(p^e)`` inherit the
                   SAME proven near-capacity no-go as prime fields ``GF(p)``.
                   Both capacity-region bounds therefore fire identically for
                   ``"prime"``, ``"extension"`` and ``None``; the parameter no
                   longer gates them.

The capacity-region split (the honest semantics)
------------------------------------------------
Above the Johnson radius there is NOT one homogeneous "no-go". The window
``[J, capacity)`` is split into TWO regions that the calculator now models with
two distinct bounds, so the semantics match the math (and the project's own
thesis that the band ``(J, r_E)`` is OPEN):

* ``unknown-beyond-johnson`` for ``J < delta < nogo_split_radius`` — the
  **OPEN** band between the Johnson radius and the established-unsafe ceiling
  (R13: ``(1-rho) - 2/s_max(b)`` per the assembled power-of-two lemma where it
  fires; at 256-bit, the generic CS25/Elias ceiling). NO
  certified *positive* MCA theorem is known here, but neither is a proven
  impossibility: the region is genuinely open. The bound returns a *vacuous*
  ``1.0`` meaning "cannot certify security from the current bounds", which is
  NOT a proven no-go. ``verified=False`` (there is no positive theorem to
  verify).
* ``proven-near-capacity-nogo`` for ``delta >= nogo_split_radius`` — the
  thin near-capacity band where the proximity gap / CA provably FAILS
  (Kambiré / BCHKS / Crites-Stewart). The bound returns ``1.0`` as a genuine
  no-go, ``verified=True`` (with the per-field ceiling carrying the assembled
  lemma's soft spots: KK25 cited not re-proved, rho=1/2 N1-conditional).

``value(...)`` returns a probability ``in (0, 1]`` for ``kind == "mca"`` and a
**list size** (a count ``>= 1``) for ``kind == "list"``.  The composite model in
``soundness.py`` converts a list size to an error via ``list_size / |F|``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Callable, Dict, List, Optional

from proximity_parameters import (
    as_fraction,
    johnson_radius,
    list_decoding_capacity_radius,
    unique_decoding_radius,
)

# ---------------------------------------------------------------------------
# Bound container
# ---------------------------------------------------------------------------

# A validity predicate and a value function share this signature.
ValidityFn = Callable[..., bool]
ValueFn = Callable[..., float]


@dataclass
class Bound:
    """One proximity-gap / MCA / list-size bound with full provenance."""

    name: str
    regime: str                      # "unique-decoding" | "johnson" | "capacity"
    kind: str                        # "mca" | "list"
    source_tag: str                  # short citation handle
    formula_str: str                 # human-readable formula
    _validity: ValidityFn
    _value: ValueFn
    verified: bool = False           # True only once checked against the paper
    notes: str = ""
    is_meta: bool = False            # True for *composition* bounds (e.g. the
    #                                  interleaving union bound) whose value
    #                                  defers to other registered single-code
    #                                  bounds.  Excluded from single-code
    #                                  selection to avoid self-recursion.

    def validity(self, rho, n, logF, delta, eta=None, m=1,
                 field_type: Optional[str] = None) -> bool:
        """Whether this bound is *claimed to apply* at these parameters."""
        try:
            return bool(self._validity(rho=as_fraction(rho), n=int(n),
                                       logF=float(logF), delta=float(delta),
                                       eta=(None if eta is None else float(eta)),
                                       m=int(m), field_type=field_type))
        except Exception:
            return False

    def value(self, rho, n, logF, delta, eta=None, m=1,
              field_type: Optional[str] = None) -> float:
        """Evaluate the bound (error probability for ``mca``; list size for
        ``list``).  Caller is responsible for first checking :meth:`validity`.
        """
        return float(self._value(rho=as_fraction(rho), n=int(n),
                                  logF=float(logF), delta=float(delta),
                                  eta=(None if eta is None else float(eta)),
                                  m=int(m), field_type=field_type))

    # convenience -----------------------------------------------------------
    def describe(self) -> str:
        flag = "VERIFIED" if self.verified else "VERIFY"
        return (f"[{flag}] {self.name} ({self.source_tag}) "
                f"regime={self.regime} kind={self.kind}\n"
                f"    formula: {self.formula_str}\n"
                f"    notes  : {self.notes}")


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

_REGISTRY: Dict[str, Bound] = {}


def register_bound(bound: Bound, *, overwrite: bool = False) -> Bound:
    """Add a :class:`Bound` to the global registry, keyed by ``bound.name``."""
    if bound.name in _REGISTRY and not overwrite:
        raise KeyError(f"bound {bound.name!r} already registered "
                       f"(pass overwrite=True to replace)")
    _REGISTRY[bound.name] = bound
    return bound


def get_bound(name: str) -> Bound:
    return _REGISTRY[name]


def all_bounds() -> List[Bound]:
    return list(_REGISTRY.values())


def lookup(regime: str, kind: str) -> List[Bound]:
    """All registered bounds matching a ``regime`` and ``kind``."""
    return [b for b in _REGISTRY.values()
            if b.regime == regime and b.kind == kind]


def applicable(rho, n, logF, delta, eta=None, m=1, kind: str = "mca",
               include_meta: bool = True,
               field_type: Optional[str] = None) -> List[Bound]:
    """All registered bounds of the given ``kind`` whose ``validity`` holds.

    Set ``include_meta=False`` to exclude composition meta-bounds (used by the
    meta-bounds themselves to select an underlying single-code bound without
    recursing into themselves).
    """
    return [b for b in _REGISTRY.values()
            if b.kind == kind
            and (include_meta or not b.is_meta)
            and b.validity(rho, n, logF, delta, eta=eta, m=m,
                           field_type=field_type)]


def best_mca_bound(rho, n, logF, delta, eta=None, m=1,
                   include_meta: bool = True,
                   field_type: Optional[str] = None) -> Optional[Bound]:
    """Among applicable ``mca`` bounds, the one giving the smallest error.

    Returns ``None`` if no registered MCA bound claims to apply at ``delta``
    (e.g. ``delta`` is beyond every regime's validity window).  Pass
    ``include_meta=False`` to consider only single-code (non-composition)
    bounds.

    NOTE: the capacity-region bounds (``unknown-beyond-johnson`` and
    ``proven-near-capacity-nogo``) both report a *vacuous* error of 1.0 (no
    security claimable), so a positive bound — when one applies — is always
    preferred by the ``min``.  The two differ only in *meaning*:
    ``unknown-beyond-johnson`` is "cannot certify (OPEN region)", whereas
    ``proven-near-capacity-nogo`` is a proven impossibility.
    """
    cands = applicable(rho, n, logF, delta, eta=eta, m=m, kind="mca",
                       include_meta=include_meta, field_type=field_type)
    if not cands:
        return None
    return min(cands, key=lambda b: b.value(rho, n, logF, delta, eta=eta, m=m,
                                            field_type=field_type))


# ===========================================================================
# Concrete bounds.  Constants now FILLED from ../literature/notes; each carries
# a `verified` flag (True = traceable to a source note; False = flagged VERIFY
# placeholder for a hidden / not-yet-pinned constant).
# ===========================================================================
#
# We work with E := |F| = 2**logF only through logs to stay exact for 256-bit
# fields.  A term "A / |F|" therefore becomes "A * 2**(-logF)"; we compute it in
# log space when A itself is a polynomial in n to avoid overflow, then exp back.


def _over_field(numerator_log2: float, logF: float) -> float:
    """Return ``2**(numerator_log2) / 2**logF`` clamped to (0, 1].

    Both arguments are *log2* quantities; the ratio is ``2**(num - logF)``.
    Clamped at 1.0 because a "probability" bound exceeding 1 is vacuous (the
    regime simply gives no information there).
    """
    return min(1.0, 2.0 ** (numerator_log2 - logF))


def _degree(rho: Fraction, n: int) -> int:
    """RS message degree ``d = k - 1 = floor(rho*n) - 1`` (>= 1)."""
    k = (rho.numerator * n) // rho.denominator
    return max(1, k - 1)


# --- 1. Unique-decoding regime --------------------------------------------
#
# For delta <= (1 - rho)/2 the correlated-agreement statement is classical and
# tight up to the leading constant: the (M)CA error is O(n / |F|).  We use the
# explicit form  C_UD * n / |F|  with C_UD = 1.
#
# This matches the BCIKS20 unique-decoding error eps_U = n/q (transcribed in
# crites-stewart.md, "the true [BCIKS Thm 1.2/1.4] error ... unique-decoding
# delta in (0,(1-rho)/2] -> eps_U = n/q").  Numerator n (not the degree d) and
# leading constant 1 are the standard statement.

_C_UD = 1.0   # leading constant in the unique-decoding (M)CA error  (= 1).


def _ud_validity(*, rho, n, logF, delta, eta, m, field_type):
    return 0 < delta <= float(unique_decoding_radius(rho))


def _ud_value(*, rho, n, logF, delta, eta, m, field_type):
    # error <= C_UD * n / |F|   (single code).  log2 numerator = log2(C_UD * n).
    num_log2 = math.log2(_C_UD * n)
    return _over_field(num_log2, logF)


register_bound(Bound(
    name="unique-decoding",
    regime="unique-decoding",
    kind="mca",
    source_tag="BCIKS20-UDR (via crites-stewart.md)",
    formula_str="eps_mca(C, delta) <= n / |F|   for delta <= (1-rho)/2",
    _validity=_ud_validity,
    _value=_ud_value,
    verified=True,
    notes=("Classical correlated-agreement up to the unique-decoding radius; "
           "BCIKS20 error eps_U = n/|F| (transcribed in crites-stewart.md). "
           "Leading constant C_UD = 1, numerator n."),
))


# --- 2. Johnson regime, Bordage-Chiesa (the headline positive MCA bound) ----
#
# Bordage-Chiesa-Guan-Manzur (eprint 2025/2051, Thm 2 / Thm 9.2), transcribed in
# bordage-chiesa.md.  For C = RS[F, D, k], rho = k/n, the univariate-powers /
# affine-line generator (the FRI/STIR/WHIR base case, d = degree = k - 1), and
# any integer m >= 3:
#
#   eps_MCA(delta) <= (m + 1/2)^7 * n^2 * d / (3 * rho^{3/2} * |F|)
#       valid for   delta <= 1 - (1 + 1/(2m)) * sqrt(rho).
#
# As m -> infinity the validity window -> the Johnson radius 1 - sqrt(rho).
# Larger m moves delta closer to Johnson but inflates the error (the (m+1/2)^7
# factor): the m-vs-delta trade is the lever the soundness solver sweeps.
# d = k - 1 = floor(rho*n) - 1.  This matches BCIKS20 correlated agreement and
# extends it to MCA (loss-free, = ABF's eps_mca).

_C_BC = 1.0 / 3.0     # leading constant 1/(3 rho^{3/2}) -> the "1/3" part.
_EXP_M_BC = 7.0       # exponent of (m + 1/2).
_EXP_N_BC = 2.0       # exponent of n.
_BC_M_MIN = 3         # smallest admissible Bordage-Chiesa trade parameter.


def _bc_dmax(rho: Fraction, m: int) -> float:
    """Bordage-Chiesa validity ceiling 1 - (1 + 1/(2m)) sqrt(rho) for trade m."""
    return 1.0 - (1.0 + 1.0 / (2.0 * m)) * math.sqrt(float(rho))


def _bc_validity(*, rho, n, logF, delta, eta, m, field_type):
    mm = max(int(m), _BC_M_MIN)
    if delta <= 0:
        return False
    # Applies up to 1 - (1 + 1/(2m)) sqrt(rho); strictly below the Johnson
    # radius for every finite m, reaching it only in the m -> inf limit.
    return delta <= _bc_dmax(rho, mm) + 1e-12


def _bc_value(*, rho, n, logF, delta, eta, m, field_type):
    mm = max(int(m), _BC_M_MIN)
    d = _degree(rho, n)
    # log2 numerator = log2(1/3) + 7 log2(m+1/2) + 2 log2(n) + log2(d)
    #                  - (3/2) log2(rho).
    num_log2 = (math.log2(_C_BC)
                + _EXP_M_BC * math.log2(mm + 0.5)
                + _EXP_N_BC * math.log2(n)
                + math.log2(d)
                - 1.5 * math.log2(float(rho)))
    return _over_field(num_log2, logF)


register_bound(Bound(
    name="bordage-chiesa-johnson",
    regime="johnson",
    kind="mca",
    source_tag="Bordage-Chiesa-2025/2051-Thm9.2",
    formula_str=("eps_mca(C, delta) <= (m+1/2)^7 * n^2 * d / (3 * rho^{3/2} * |F|)"
                 "   for delta <= 1 - (1 + 1/(2m)) sqrt(rho),  m >= 3 integer, "
                 "d = k-1 = floor(rho*n)-1"),
    _validity=_bc_validity,
    _value=_bc_value,
    verified=True,
    notes=("Headline positive MCA bound for smooth-domain RS (any domain D, "
           "incl. multiplicative subgroups; any field). Loss-free (= ABF "
           "eps_mca). m>=3 is the integer trade parameter: larger m pushes "
           "delta toward Johnson 1-sqrt(rho) but inflates error via (m+1/2)^7. "
           "Constants 1/3, exp 7 (on m+1/2), exp 2 (on n), degree d in "
           "numerator all transcribed verbatim in bordage-chiesa.md "
           "(Thm 9.2 / Lemma 9.3, univariate-powers/affine-line case, d=1 per "
           "generator -> here the RS message degree d=k-1)."),
))


# --- 3. Johnson regime, BCHKS Thm 1.5 (better n-scaling, hidden rho-constant) -
#
# BCHKS (eprint 2025/2055, Thm 1.5; = ABF Thm 4.12 core), transcribed in
# bchks.md.  In the ABF convention: for delta in [0, 1 - sqrt(rho)), set
# eta = 1 - sqrt(rho) - delta; the correlated-agreement set size a = O_rho(n/η^5)
# suffices for ZERO proximity loss, giving an MCA error
#
#   eps_mca(delta) <= C_rho * n / (eta^5 * |F|)        up to delta = J - eta.
#
# The n-scaling is LINEAR (vs Bordage-Chiesa's n^2 * d ~ n^3), so for large n
# this is the stronger bound IF the hidden constant C_rho is benign.  But the
# constant is an O_rho(.) (bchks.md flags it; the leading constant in BCHKS's a
# is 1/(48(1-delta)^{3/2}) = 1/(48 rho^{3/2}) only in the eta << sqrt(rho)
# asymptotic).  We register C_rho = 1 as a PLACEHOLDER and keep verified=False
# so this bound is never silently trusted as the binding one.

_C_RHO_BCHKS = 1.0    # VERIFY: hidden O_rho(.) constant in BCHKS Thm 1.5 (=1 placeholder).
_EXP_ETA_BCHKS = 5.0  # exponent of 1/eta in BCHKS Thm 1.5 (a = O_rho(n/eta^5)).


def _bchks_validity(*, rho, n, logF, delta, eta, m, field_type):
    if eta is None or eta <= 0:
        return False
    J = johnson_radius(rho)
    if not (0 < delta < J):
        return False
    return delta <= J - eta + 1e-12


def _bchks_value(*, rho, n, logF, delta, eta, m, field_type):
    # eps <= C_rho * n / (eta^5 * |F|).
    num_log2 = (math.log2(_C_RHO_BCHKS) + math.log2(n)
                - _EXP_ETA_BCHKS * math.log2(eta))
    return _over_field(num_log2, logF)


register_bound(Bound(
    name="bchks-johnson",
    regime="johnson",
    kind="mca",
    source_tag="BCHKS-2025/2055-Thm1.5",
    formula_str=("eps_mca(C, delta) <= C_rho * n / (eta^5 * |F|)   "
                 "for delta <= 1 - sqrt(rho) - eta,  eta > 0   "
                 "[ABF convention; BCHKS states a = O_rho(n/eta^5), eps* = 0]"),
    _validity=_bchks_validity,
    _value=_bchks_value,
    verified=False,
    notes=("VERIFY the hidden O_rho(.) constant C_rho (set to 1.0 placeholder). "
           "BCHKS Thm 1.5 gives a = O_rho(n/eta^5) far-set size for ZERO "
           "proximity loss (eps* = 0), translating to MCA error "
           "C_rho * n / (eta^5 |F|); the eta-exponent 5 is firm (Table 1, "
           "'this work', <Johnson, a = Theta_delta(n/eta^5)), the leading "
           "constant is an unspecified O_rho(.) (asymptotic leading term "
           "1/(48 rho^{3/2}) only for eta << sqrt(rho)). LINEAR in n (vs "
           "Bordage-Chiesa's n^2 d), so asymptotically stronger, but kept "
           "unverified until C_rho is pinned. NOTE the BCHKS<->ABF convention "
           "flip: BCHKS delta = 1-rho (min-distance), gamma = radius; here "
           "translated to ABF radius delta and rho = k/n."),
))


# --- 4. Capacity region: SPLIT into an OPEN band + a PROVEN near-capacity no-go
#
# Above the Johnson radius the window [J, capacity) is NOT homogeneous, and the
# calculator must not paint it with a single "proven no-go" brush — that would
# contradict the project's own thesis that the band (J, r_E) is OPEN.  We split
# it at the finite-field established-unsafe ceiling ``nogo_split_radius(rho,
# logF)`` — the R13 assembled-lemma value ``(1-rho) - 2/s_max(b)`` where it
# exists (s_max = 16/16/32 at b = 31/64/128), else (b = 256: no threshold-
# established Kambiré-type ceiling) the generic CS25/Elias list-decoding-
# capacity ceiling — into two clearly-distinct bounds:
#
#   (4a) ``unknown-beyond-johnson`` for  J < delta < r_E   — the OPEN band.
#        There is NO certified positive MCA theorem here (Bordage-Chiesa / BCHKS
#        run out at J), but there is ALSO no proven impossibility: this is the
#        genuinely open region the prize asks about.  We register a VACUOUS bound
#        (value 1.0, verified=False) so the solver never extrapolates the Johnson
#        form past its proven range — the 1.0 means "cannot certify security from
#        the current bounds", NOT "proven insecure".
#
#   (4b) ``proven-near-capacity-nogo`` for  delta >= r_E   — the thin band just
#        below Singleton where the proximity gap / CA provably FAILS:
#          * Kambiré (arXiv 2604.09724, Thm 1): C = RS[F_p, <ω>, k], n = 2^t,
#            prime p ≡ 1 mod n, p < n^A, rho in (0, 1/2).  At
#            delta = (1-rho) - 2/(K log n) there are >= n^C scalars z with
#            Δ(f+zg, C) <= delta yet Δ([f,g], C^2) > delta — a direct
#            proximity-gap + CA failure.  rho = 1/2 excluded as written.
#          * Crites-Stewart (eprint 2025/2046, Cor 1): ε_ca = 1 once delta
#            reaches the list-decoding capacity 1 - H_q(rho) (~ 1/log2 q below
#            Singleton), q >= n.
#          * BCHKS Thm 1.13 (cond. on Conj 1.12; M31 UNCONDITIONAL at rho ≈ 1/2):
#            prime multiplicative-subgroup RS at gamma = delta_min - Θ(1/log n)
#            forces proximity loss Θ(1/log n).
#        This is a genuine no-go (value 1.0, verified=True).  R13: the split
#        radius is the ASSEMBLED-lemma ceiling ``(1-rho) - 2/s_max(b)`` (s a
#        power of two; s_max = 16/16/32 at b = 31/64/128) where it exists, else
#        the generic CS25/Elias list-decoding-capacity ceiling (b = 256: the
#        Kambiré mechanism is sub-threshold there).  The pre-R13 continuum
#        ``(1-rho) - 6/log2|F|`` is an asymptotic reference only — flagged in
#        the notes.
#
# FIELD-AGNOSTIC (the decisive correction; see n2-verdict.md).  The Kambiré-type
# near-capacity counterexample's distinct-bad-scalar count is a
# CHARACTERISTIC-ZERO CYCLOTOMIC INVARIANT (e_1, e_2 are reductions mod the
# characteristic p of fixed elements of Z[ξ_s], independent of the extension
# degree e).  A primitive s-th root of unity lives in F_p when p ≡ 1 (mod s) and
# genuinely in F_{p^2}\F_p when p ≡ -1 (mod s); in BOTH cases the reduction
# realizes the SAME char-0 count.  So genuine odd-characteristic extension fields
# GF(p^e) inherit the SAME near-capacity no-go as prime fields GF(p) (the
# obstruction is MULTIPLICATIVE — cyclotomic subset sums — not the char-2
# additive subspace mechanism, and it lives in F_{p^e} identically to F_p).
# BOTH new bounds therefore fire for prime, extension and unknown field types
# identically; field_type is ignored by both.

# Pre-R13 continuum Kambiré gap below capacity: eta_min ≈ 6 / log2|F| (from
# needing a prime p in [4^s, 8^s] with 8^s <= |F| => s <= log2|F|/3 =>
# eta = 2/s).  R13: superseded as the operative ceiling by the power-of-two
# s_max(b) lemma (see kambire_smax / kambire_unsafe_delta); retained ONLY as an
# explicitly-labeled asymptotic reference.
_KAMBIRE_GAP_CONST = 6.0   # pre-R13 continuum reference; NOT the operative ceiling (R13).


def _capacity_window(rho: Fraction, delta: float) -> bool:
    """The whole region above Johnson, up to (open) Singleton capacity."""
    return johnson_radius(rho) <= delta < 1.0 - float(rho)


def nogo_split_radius(rho, logF) -> float:
    """The radius splitting the OPEN band from the PROVEN near-capacity no-go.

    R13: the established-unsafe ceiling is the assembled-lemma value
    ``(1-rho) - 2/s_max(b)`` (:func:`kambire_unsafe_delta`) where it exists.
    Where no admissible power-of-two ``s`` exists (b = 256: the Kambiré
    mechanism is sub-threshold at the prize target), the **generic CS25/Elias
    list-decoding-capacity ceiling** ``H_q^{-1}(1 - rho)`` takes over — it
    stands at every field size (CS25 Cor 1: eps_ca = 1 there), so the no-go
    band never disappears, it just starts at the generic ceiling instead of a
    Kambiré-type constant-eta one.  Clamped to be >= Johnson.
    """
    r = as_fraction(rho)
    d = kambire_unsafe_delta(r, logF)
    if d is None:
        q = 2 ** min(int(round(float(logF))), 4096)
        d = list_decoding_capacity_radius(r, q)
    return max(johnson_radius(r), d)


# --- 4a. OPEN band:  J < delta < r_E  (no positive theorem; NOT a proven no-go).

def _unknown_beyond_johnson_validity(*, rho, n, logF, delta, eta, m, field_type):
    # Fires in the OPEN band between the Johnson radius (the proven positive
    # ceiling) and the finite-field provable-unsafe ceiling r_E.  FIELD-AGNOSTIC
    # (field_type ignored): the region is open identically for prime, extension
    # and unknown fields.  Lower edge is the Johnson radius (inclusive — the
    # positive Johnson-regime bounds reach up to but not past J); upper edge is
    # r_E (exclusive — at/above r_E the PROVEN no-go takes over).
    if not _capacity_window(rho, delta):
        return False
    split = nogo_split_radius(rho, logF)
    return delta < split - 1e-12


def _vacuous_value(*, rho, n, logF, delta, eta, m, field_type):
    # No certified POSITIVE bound here, so no finite query count certifies
    # security.  Vacuous 1.0 = "cannot certify", NOT a proven impossibility.
    return 1.0


register_bound(Bound(
    name="unknown-beyond-johnson",
    regime="capacity",
    kind="mca",
    source_tag="open-region (no positive theorem)",
    formula_str=("OPEN band (FIELD-AGNOSTIC — prime & extension): "
                 "1-sqrt(rho) <= delta < (1-rho)-2/s_max(b) (R13 assembled "
                 "lemma, s a power of two; s_max = 16/16/32 at b = 31/64/128; "
                 "at b = 256 the band runs to the generic CS25/Elias ceiling). "
                 "No certified positive MCA bound; value = 1.0 (VACUOUS, "
                 "cannot certify), NOT a proven no-go."),
    _validity=_unknown_beyond_johnson_validity,
    _value=_vacuous_value,
    verified=False,
    notes=("No certified *positive* MCA bound is available between Johnson and "
           "the established no-go split; the region is **OPEN**. The vacuous "
           "1.0 means 'cannot certify security from current bounds', NOT a "
           "proven impossibility. This is the band the Proximity Prize asks "
           "about: above the proven Johnson positive ceiling and below the "
           "near-capacity no-go split (R13: (1-rho) - 2/s_max(b) where the "
           "assembled lemma fires; at 256-bit no Kambiré-type ceiling is "
           "threshold-established and the split is the generic CS25/Elias "
           "ceiling ~ 1 - H_q(rho)). The positive bounds (Bordage-Chiesa "
           "Thm 9.2, BCHKS Thm 1.5) run out exactly at the Johnson radius, and "
           "no impossibility is known until the split — so neither security NOR "
           "insecurity is certified here. FIELD-AGNOSTIC: open identically for "
           "prime, genuine odd-characteristic extension, and unknown field "
           "types (the matching near-capacity no-go is itself field-agnostic, "
           "n2-verdict.md). verified=False because there is no positive theorem "
           "to verify (this is an honest 'don't know', not a checked bound)."),
))


# --- 4b. PROVEN near-capacity no-go:  delta >= r_E.

def _proven_near_capacity_nogo_validity(*, rho, n, logF, delta, eta, m, field_type):
    # Fires ONLY in the thin near-capacity band delta >= r_E where the proximity
    # gap / CA provably fails (Kambiré / CS / BCHKS).  FIELD-AGNOSTIC.
    if not _capacity_window(rho, delta):
        return False
    split = nogo_split_radius(rho, logF)
    return delta >= split - 1e-12


def _nogo_value(*, rho, n, logF, delta, eta, m, field_type):
    # Proven failure: the proximity gap fails outright here.
    return 1.0


register_bound(Bound(
    name="proven-near-capacity-nogo",
    regime="capacity",
    kind="mca",
    source_tag="Kambiré 2604.09724 / BCHKS / CS25 (R13 assembled per-field ceiling)",
    formula_str=("PROVEN no-go (FIELD-AGNOSTIC — prime & extension): eps_ca = 1 "
                 "/ proximity-gap fails for delta >= (1-rho)-2/s_max(b) (R13 "
                 "assembled lemma, s a power of two; s_max = 16/16/32 at "
                 "b = 31/64/128; at b = 256 no Kambiré-type ceiling is "
                 "threshold-established and the split is the generic CS25/Elias "
                 "ceiling ~ 1 - H_q(rho) ~ (1-rho) - 1/log2 q, which stands at "
                 "every field size); Kambiré failure at "
                 "delta = (1-rho) - 2/(K log n)"),
    _validity=_proven_near_capacity_nogo_validity,
    _value=_nogo_value,
    verified=True,
    notes=("PROVEN no-go in the thin near-capacity band delta >= split, with "
           "split = (1-rho) - 2/s_max(b) per the R13 assembled lemma (s a power "
           "of two; s_max = 16/16/32 at b = 31/64/128) where it fires, else "
           "(b = 256) the generic CS25/Elias list-decoding-capacity ceiling. "
           "Returns a genuine vacuous error 1.0 "
           "(the proximity gap fails outright). Backed by Kambiré Thm 1 (smooth "
           "D=<ω>, n=2^t, prime p, rho in (0,1/2): >= n^C bad scalars at "
           "delta=(1-rho)-2/(K log n)); Crites-Stewart Cor 1 (eps_ca=1 at the "
           "list-decoding capacity 1-H_q(rho), needs q>=n); BCHKS Thm 1.13 "
           "(prime mult.-subgroup, loss Θ(1/log n); M31 UNCONDITIONAL at "
           "rho≈1/2 with a CA failure at delta=1/2). FIELD-AGNOSTIC: prime "
           "= extension. EXTENSION INHERITANCE (n2-verdict.md): the "
           "near-capacity counterexample's bad-scalar count is a "
           "CHARACTERISTIC-ZERO CYCLOTOMIC INVARIANT, so genuine GF(p^e) "
           "realizes the same count as GF(p) (verified exactly; decoder-free "
           "certified bad lines on genuine GF(31^2), GF(127^2) at "
           "delta=1/2,3/4). Lemma status: ASSEMBLED (KK25 distinctness cited "
           "not re-proved; rho=1/2 N1-conditional; deg-convention tracked). "
           "The pre-R13 continuum (1-rho) - 6/log2|F| is an asymptotic "
           "reference only. The R13 gap 2/s_max(b) is a CONSTANT (= 0.125 at "
           "31/64-bit, 0.0625 at 128-bit), NOT o(1) at deployed sizes — so this "
           "no-go forecloses only a thin band just below capacity; the rho=1/2 "
           "boundary is excluded by Kambiré as written (covered instead by the "
           "unconditional M31 BCHKS witness)."),
))


# --- 5. Interleaving / batching: MCA composition (Bordage-Chiesa Lem 10.1) --
#
# Bordage-Chiesa Lemma 10.1 (transcribed in bordage-chiesa.md): if G has MCA for
# C with error eps_MCA, then G has MCA for the m-interleaving C^m with error
# m * eps_MCA.  This wraps an underlying single-code MCA bound (chosen at the
# same m, which for the positive bound also sets the Bordage-Chiesa trade
# parameter) and multiplies by m.  It is a *meta-bound*: its value defers to the
# best applicable single-code bound.

def _interleave_mca_validity(*, rho, n, logF, delta, eta, m, field_type):
    if m < 1:
        return False
    # Applicable exactly when some *single-code* MCA bound applies.
    # include_meta=False prevents recursing into this meta-bound itself.
    return best_mca_bound(rho, n, logF, delta, eta=eta, m=m,
                          include_meta=False, field_type=field_type) is not None


def _interleave_mca_value(*, rho, n, logF, delta, eta, m, field_type):
    base = best_mca_bound(rho, n, logF, delta, eta=eta, m=m,
                          include_meta=False, field_type=field_type)
    if base is None:
        return 1.0
    single = base.value(rho, n, logF, delta, eta=eta, m=m, field_type=field_type)
    return min(1.0, m * single)


register_bound(Bound(
    name="interleave-mca-union",
    regime="johnson",  # nominal; validity defers to underlying single-code bound
    kind="mca",
    source_tag="Bordage-Chiesa-2025/2051-Lem10.1",
    formula_str="eps_mca(C^{equiv m}, delta) <= m * eps_mca(C, delta)",
    _validity=_interleave_mca_validity,
    _value=_interleave_mca_value,
    verified=True,
    is_meta=True,
    notes=("Bordage-Chiesa Lemma 10.1 (verbatim in bordage-chiesa.md): MCA for "
           "the m-interleaving C^m has error <= m * eps_MCA(C). Linear-in-m "
           "union over the interleaved instances, wrapping the best applicable "
           "single-code MCA bound. A benign factor below the Johnson radius; "
           "NOTE m here also drives the underlying Bordage-Chiesa trade "
           "parameter."),
))


# --- 6. Interleaved list size |Lambda(C^{equiv m}, delta)| (Johnson bound) ---
#
# Within the Johnson regime the list size for a single RS code obeys the
# Johnson bound.  For delta = J(rho) - eta with eta > 0, the standard RS Johnson
# list-size bound is poly(n)/poly(eta); a defensible registry value is
#
#       |Lambda| <= n / eta
#
# (the classical Johnson-bound list size is O(1/eta) up to a factor n in the
# worst case; the exact constant / exponent is the subject of the deep
# list-size analysis being done separately).  We multiply by m for the
# interleaved code as a conservative union over instances.  Returns a COUNT.
#
# VERIFY: the exact constant and the precise eta / n / m exponents.  Agent 2.2
# does the deep list-size analysis; this is a defensible placeholder
# (verified=False) so soundness has a non-trivial list term.

_C_L = 1.0    # leading constant in the Johnson list-size bound.
_P_L = 1.0    # VERIFY exponent of n.
_Q_L = 1.0    # VERIFY exponent of 1/eta.
_S_L = 1.0    # VERIFY exponent of m (interleaving union).


def _listsize_validity(*, rho, n, logF, delta, eta, m, field_type):
    if eta is None or eta <= 0 or m < 1:
        return False
    J = johnson_radius(rho)
    return 0 < delta < J and delta <= J - eta + 1e-12


def _listsize_value(*, rho, n, logF, delta, eta, m, field_type):
    # Returns a list-size COUNT (>= 1), not a probability.
    log2_count = (math.log2(_C_L) + _P_L * math.log2(n)
                  - _Q_L * math.log2(eta) + _S_L * math.log2(m))
    return max(1.0, 2.0 ** log2_count)


register_bound(Bound(
    name="interleave-listsize",
    regime="johnson",
    kind="list",
    source_tag="Johnson-list-bound (constant VERIFY; deep analysis = Agent 2.2)",
    formula_str="|Lambda(C^{equiv m}, delta)| <= C_L * n^p * m^s / eta^q",
    _validity=_listsize_validity,
    _value=_listsize_value,
    verified=False,
    notes=("Johnson-bound RS list size within delta <= J(rho) - eta. Defensible "
           "registry value |Lambda| <= n * m / eta (C_L=1, p=1, q=1, s=1). "
           "VERIFY the exact constant and the n / (1/eta) / m exponents — the "
           "deep list-size analysis is done separately (Agent 2.2). "
           "soundness.py converts to error via |Lambda| / |F|."),
))


def best_listsize_bound(rho, n, logF, delta, eta=None, m=1,
                        field_type: Optional[str] = None) -> Optional[Bound]:
    """Among applicable ``list`` bounds, the one giving the smallest list size."""
    cands = applicable(rho, n, logF, delta, eta=eta, m=m, kind="list",
                       field_type=field_type)
    if not cands:
        return None
    return min(cands, key=lambda b: b.value(rho, n, logF, delta, eta=eta, m=m,
                                            field_type=field_type))


def interleaved_mca(rho, n, logF, delta, eta=None, m=1,
                    field_type: Optional[str] = None) -> Optional[Bound]:
    """The MCA bound to use for the *interleaved* code ``C^{equiv m}``.

    This is the protocol-relevant selector used by :mod:`soundness`:

    * For ``m == 1`` it returns the best applicable single-code MCA bound.
    * For ``m > 1`` it returns the best applicable *composition* (meta) bound,
      which carries the dependence on ``m`` (Bordage-Chiesa Lemma 10.1, the
      linear union bound).  If no meta-bound applies it falls back to the
      single-code bound.

    Returning the :class:`Bound` (not a number) keeps provenance: the caller
    reads both its ``name`` and its ``value(...)``.
    """
    if int(m) <= 1:
        return best_mca_bound(rho, n, logF, delta, eta=eta, m=1,
                              include_meta=False, field_type=field_type)
    # m > 1: prefer a composition meta-bound.
    metas = [b for b in applicable(rho, n, logF, delta, eta=eta, m=m, kind="mca",
                                   field_type=field_type)
             if b.is_meta]
    if metas:
        return min(metas, key=lambda b: b.value(rho, n, logF, delta, eta=eta,
                                                m=m, field_type=field_type))
    return best_mca_bound(rho, n, logF, delta, eta=eta, m=m,
                          include_meta=False, field_type=field_type)


# ---------------------------------------------------------------------------
# Kambiré finite-field negative radius (for the ledger's provable-unsafe delta)
# ---------------------------------------------------------------------------

def kambire_smax(rho, logF) -> Optional[int]:
    """Largest admissible power-of-two quotient parameter ``s_max(b)`` (R13).

    The assembled per-field lemma (technical-note.md §3.2, problem-ledger.md §4)
    requires the Kambiré quotient parameter ``s = n/m`` to be a **power of two**
    (R13 integrality — the smooth domain has order ``n = 2^t``) satisfying:

      (i)  distinctness (KK25; cited, not re-proved):  ``p > phi(s)^{phi(s)}``,
           with ``phi(s) = s/2`` for power-of-two ``s >= 4``;
      (ii) count:  ``N(s, rho) > 2^{b-128}`` at ``b = log2|F|`` bits.  We check
           the rho=1/2 antipodal count ``3^{s/2}`` (the largest of the per-rate
           counts; for rho < 1/2 Kambiré's N0_sum is smaller still).  At every
           deployed sub-256-bit size the threshold ``2^{b-128} <= 1`` is trivial,
           and at b = 256 even ``3^{s/2}`` fails for every distinctness-allowed
           ``s`` (``3^32 = 2^51 << 2^128``) — so the proxy reproduces the
           lemma's tabled ``s_max(b)`` at all deployed sizes and rates;
      (iv) radius/size:  the firing radius stays strictly above Johnson,
           ``(1-rho) - 2/s > J = 1 - sqrt(rho)``, i.e. ``s > 2/(sqrt(rho)-rho)``
           (threshold <= 10.7 over the deployed rates, so ``s >= 16`` suffices).

    Returns the largest such power of two, or ``None`` when no ``s`` satisfies
    (i)+(ii)+(iv) — in particular at ``b = 256``, where the mechanism is
    infeasible at the prize threshold and establishes **no threshold ceiling**.
    Tabled values: ``s_max = 16/16/32`` at ``b = 31/64/128``.
    """
    r = as_fraction(rho)
    rho_f = float(r)
    b = float(logF)
    iv_threshold = 2.0 / (math.sqrt(rho_f) - rho_f)  # hypothesis (iv)
    best: Optional[int] = None
    s = 4
    while s <= 1 << 20:
        phi = s // 2                                # phi(2^k) = 2^(k-1)
        distinct_ok = b > phi * math.log2(phi)      # (i)  p > phi(s)^phi(s)
        count_ok = (s / 2.0) * math.log2(3.0) > b - 128.0  # (ii) 3^{s/2} > 2^{b-128}
        radius_ok = s > iv_threshold                # (iv) above Johnson
        if distinct_ok and count_ok and radius_ok:
            best = s
        s <<= 1
    return best


def kambire_unsafe_delta(rho, logF) -> Optional[float]:
    """Smallest radius at which the no-go is established, finite-field (R13).

    FIELD-AGNOSTIC (n2-verdict.md): the same provable-unsafe radius holds for
    prime AND genuine odd-characteristic extension fields, because the
    near-capacity counterexample's bad-scalar count is a characteristic-zero
    cyclotomic invariant — extensions inherit the prime no-go identically.

    R13 (s-integrality) correction: the quotient parameter ``s`` must be a
    **power of two**, so the honest per-field ceiling is the assembled lemma's

        ``delta*_C <= (1 - rho) - 2/s_max(b)``

    with ``s_max(b)`` from :func:`kambire_smax` (= 16/16/32 at b = 31/64/128).
    Returns ``None`` when no admissible ``s`` exists — at ``b = 256`` the count
    ``3^32 = 2^51 << 2^128`` is sub-threshold, so **no threshold-established
    Kambiré-type (near-capacity, constant-eta) ceiling below r_E exists** (the
    generic CS25/Elias ceiling at ~ r_E = 1 - H_q(rho) still stands at every
    field size — see the capacity-region bounds' split fallback).

    Hypothesis (iv) guarantees the returned radius is strictly above Johnson.
    The pre-R13 continuum ``(1 - rho) - 6/log2|F|`` (from the ``[4^s, 8^s]``
    Linnik window arithmetic) is retained only as an explicitly-labeled
    asymptotic reference (``_KAMBIRE_GAP_CONST``), NOT returned here.

    Status: ASSEMBLED from cited components (Kambiré Thm 1 count + KK25
    distinctness + threshold arithmetic + R13 integrality); soft spots: KK25
    cited not re-proved; rho = 1/2 N1-conditional; deg-convention tracked.
    """
    r = as_fraction(rho)
    cap = 1.0 - float(r)
    s_max = kambire_smax(r, logF)
    if s_max is None:
        return None
    return cap - 2.0 / s_max


# ---------------------------------------------------------------------------
# VERIFY surface
# ---------------------------------------------------------------------------

def verify_flags() -> List[Dict[str, str]]:
    """Structured list of every unverified bound + its VERIFY note.

    Used by the CLI and the README example so the human can see exactly which
    literature constants still need pinning.
    """
    out: List[Dict[str, str]] = []
    for b in _REGISTRY.values():
        if not b.verified:
            out.append({
                "name": b.name,
                "source_tag": b.source_tag,
                "formula": b.formula_str,
                "notes": b.notes,
            })
    return out


# Module-level dict of the tunable constants, so a ledger-driven script can
# patch them in one place (e.g. bounds.CONSTANTS['_C_RHO_BCHKS'] = ...).
CONSTANTS = {
    # unique decoding
    "_C_UD": _C_UD,
    # Bordage-Chiesa Johnson MCA (VERIFIED constants)
    "_C_BC": _C_BC, "_EXP_M_BC": _EXP_M_BC, "_EXP_N_BC": _EXP_N_BC,
    "_BC_M_MIN": _BC_M_MIN,
    # BCHKS Thm 1.5 (hidden C_rho -> VERIFY)
    "_C_RHO_BCHKS": _C_RHO_BCHKS, "_EXP_ETA_BCHKS": _EXP_ETA_BCHKS,
    # Pre-R13 continuum Kambiré gap constant (asymptotic reference ONLY; the
    # operative R13 ceiling is (1-rho) - 2/s_max(b), see kambire_smax)
    "_KAMBIRE_GAP_CONST": _KAMBIRE_GAP_CONST,
    # list size (constants VERIFY)
    "_C_L": _C_L, "_P_L": _P_L, "_Q_L": _Q_L, "_S_L": _S_L,
}


__all__ = [
    "Bound",
    "register_bound",
    "get_bound",
    "all_bounds",
    "lookup",
    "applicable",
    "best_mca_bound",
    "best_listsize_bound",
    "interleaved_mca",
    "kambire_smax",
    "kambire_unsafe_delta",
    "nogo_split_radius",
    "verify_flags",
    "CONSTANTS",
]


if __name__ == "__main__":  # pragma: no cover
    for b in all_bounds():
        print(b.describe())
        print()
