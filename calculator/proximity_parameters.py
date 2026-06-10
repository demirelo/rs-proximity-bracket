"""
proximity_parameters.py
=======================

Core radii / entropy / log helpers for Reed-Solomon proximity-gap analysis.

A Reed-Solomon code ``C = RS[F, L, k]`` has evaluation domain ``L`` with
``|L| = n`` and message dimension ``k``.  Its *rate* is ``rho = k / n`` and we
work throughout with the **relative** Hamming distance ``delta in [0, 1]``
(fraction of disagreeing coordinates).

The three radii that organize the whole subject (relative, as fractions of
``n``) are, for ``rho in (0, 1)``::

    unique decoding   : (1 - rho) / 2
    Johnson           : 1 - sqrt(rho)
    list-dec capacity : 1 - rho

and they satisfy the strict ordering

    (1 - rho)/2  <  1 - sqrt(rho)  <  1 - rho     for rho in (0, 1).

Arithmetic conventions
----------------------
* The unique-decoding radius and capacity are *rational* in ``rho`` so we keep
  them exact with :class:`fractions.Fraction` when the input is rational.
* The Johnson radius involves ``sqrt(rho)``; we expose both an exact symbolic
  value (via :mod:`sympy`) and a high-precision ``Fraction`` rationalization.
* Security-bit counts are computed as ``-log2(error)``.  We deliberately route
  these through :func:`math.log2` / :mod:`mpmath` on a single positive quantity
  rather than subtracting two large floats, to avoid catastrophic cancellation.

Nothing in this module hard-codes a literature constant; the only "physics"
here is the definition of the radii and of the q-ary entropy function.  All
proximity-gap *constants* live in ``bounds.py`` and are flagged ``VERIFY``.
"""

from __future__ import annotations

import math
from fractions import Fraction
from numbers import Rational
from typing import Union

import sympy as sp

Number = Union[int, float, Fraction]

# ---------------------------------------------------------------------------
# Rate normalization
# ---------------------------------------------------------------------------


def as_fraction(x: Number) -> Fraction:
    """Coerce ``x`` to an exact :class:`Fraction`.

    Floats are converted via :meth:`Fraction.from_float` then limited to a
    sane denominator so that, e.g., ``0.5`` becomes ``Fraction(1, 2)`` rather
    than a 53-bit dyadic monster.  Pass integers / ``Fraction`` / strings for
    exactness.
    """
    if isinstance(x, Fraction):
        return x
    if isinstance(x, int):
        return Fraction(x)
    if isinstance(x, str):
        return Fraction(x)
    if isinstance(x, float):
        # Limit denominator: assume the user meant a "nice" rate like 1/8.
        return Fraction(x).limit_denominator(10**9)
    raise TypeError(f"cannot interpret {x!r} as a rational rate")


def _check_rate(rho: Fraction) -> None:
    if not (0 < rho < 1):
        raise ValueError(f"rate rho must lie strictly in (0, 1); got {rho}")


# ---------------------------------------------------------------------------
# The three radii
# ---------------------------------------------------------------------------


def unique_decoding_radius(rho: Number) -> Fraction:
    """Relative unique-decoding radius ``(1 - rho) / 2`` (exact)."""
    r = as_fraction(rho)
    _check_rate(r)
    return (1 - r) / 2


def capacity(rho: Number) -> Fraction:
    """Relative list-decoding *capacity* radius ``1 - rho`` (exact).

    This is the information-theoretic limit: beyond ``1 - rho`` the list of
    nearby codewords can be super-polynomial / the ball can contain more than
    ``poly`` codewords for generic codes.
    """
    r = as_fraction(rho)
    _check_rate(r)
    return 1 - r


def johnson_radius_exact(rho: Number) -> sp.Expr:
    """Johnson radius ``1 - sqrt(rho)`` as an exact :mod:`sympy` expression.

    Use this when you want to keep the surd symbolic (e.g. ``rho = 1/4`` gives
    exactly ``1/2``).  For numeric work use :func:`johnson_radius` (float) or
    :func:`johnson_radius_rational` (high-precision Fraction).
    """
    r = as_fraction(rho)
    _check_rate(r)
    return 1 - sp.sqrt(sp.Rational(r.numerator, r.denominator))


def johnson_radius(rho: Number) -> float:
    """Johnson radius ``1 - sqrt(rho)`` as a float (double precision)."""
    r = as_fraction(rho)
    _check_rate(r)
    return 1.0 - math.sqrt(r.numerator / r.denominator)


def johnson_radius_rational(rho: Number, max_denominator: int = 10**15) -> Fraction:
    """Johnson radius as a rational approximation good to ``max_denominator``.

    Computed at high precision through :mod:`sympy` (50 digits) so the result
    is faithful well past double precision before the final
    :meth:`Fraction.limit_denominator`.
    """
    expr = johnson_radius_exact(rho)
    # Evaluate the (possibly still-symbolic) surd to 50 significant digits,
    # then rationalize.  evalf handles ``1 - sqrt(2)/2`` whereas sp.Float does
    # not accept an unevaluated symbolic expression.
    high = expr.evalf(50)
    val = sp.Rational(high)
    return Fraction(int(val.p), int(val.q)).limit_denominator(max_denominator)


# ---------------------------------------------------------------------------
# q-ary entropy and the capacity/entropy radius
# ---------------------------------------------------------------------------


def qary_entropy(x: Number, q: Number) -> float:
    """The q-ary entropy function ``H_q(x)`` as a float.

    Convention (standard, e.g. Guruswami-Rudra-Sudan lecture notes)::

        H_q(x) = x * log_q(q - 1) - x * log_q(x) - (1 - x) * log_q(1 - x)

    Defined for ``x in [0, 1]`` and integer alphabet size ``q >= 2``.  With this
    convention ``H_q`` increases on ``[0, 1 - 1/q]`` from ``H_q(0) = 0`` to its
    maximum ``H_q(1 - 1/q) = 1`` and is the right object for the
    *list-decoding capacity* statement: a q-ary code of rate ``R`` is list
    decodable (with polynomial lists) up to radius ``delta`` essentially when
    ``R <= 1 - H_q(delta)``.

    Edge cases: ``H_q(0) = 0`` and ``H_q(1) = log_q(q - 1)`` are returned
    exactly without hitting ``log(0)``.
    """
    xf = as_fraction(x)
    if not (0 <= xf <= 1):
        raise ValueError(f"H_q argument x must lie in [0, 1]; got {xf}")
    qi = int(q)
    if qi < 2:
        raise ValueError(f"alphabet size q must be >= 2; got {q}")

    logq = math.log(qi)

    def term(p: Fraction) -> float:
        # -p * log_q(p), with the standard limit 0*log(0) = 0.
        if p == 0:
            return 0.0
        pf = p.numerator / p.denominator
        return -pf * (math.log(pf) / logq)

    # base term is x * log_q(q-1); for q == 2 this is x * log_q(1) = 0.
    base = (xf.numerator / xf.denominator) * (math.log(qi - 1) / logq)
    return base + term(xf) + term(1 - xf)


def qary_entropy_inverse(target: float, q: Number, lo: float = 0.0,
                         hi: float | None = None, iters: int = 200) -> float:
    """Smallest ``x in [0, 1-1/q]`` with ``H_q(x) >= target`` (bisection).

    Used by :func:`list_decoding_capacity_radius`.  Because ``H_q`` is strictly
    increasing on ``[0, 1 - 1/q]`` this inverse is well defined for
    ``target in [0, 1]``.
    """
    qi = int(q)
    if hi is None:
        hi = 1.0 - 1.0 / qi
    target = max(0.0, min(1.0, target))
    for _ in range(iters):
        mid = (lo + hi) / 2
        if qary_entropy(Fraction(mid).limit_denominator(10**12), qi) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def list_decoding_capacity_radius(rho: Number, q: Number) -> float:
    """The q-ary list-decoding *capacity* radius for rate ``rho``.

    Convention documented here precisely:  we return the ``delta`` solving

        1 - H_q(delta) = rho        <=>        H_q(delta) = 1 - rho,

    i.e. the relative radius at which a rate-``rho`` q-ary code sits exactly on
    the list-decoding capacity curve.  As ``q -> infinity`` this radius tends to
    the RS capacity ``1 - rho`` (since ``H_q(delta) -> delta`` in that limit).
    The convergence is, however, only ``O(1/log2 q)``: the entropy carries a
    ``delta * log_q(q - 1)`` term whose deviation from ``delta`` decays like
    ``1/log2(q)``.  So even at ``q = 2**40`` the radius is still ~0.025 below
    capacity for ``rho = 1/2`` — do NOT read this as a tight RS bound; it is the
    generic q-ary entropy curve (tested in test_qary_capacity_radius_*).

    NOTE: this is the generic *q-ary* capacity radius, NOT the RS-specific
    proximity-gap threshold.  The RS radii above (unique / Johnson / capacity)
    are the operationally relevant ones for FRI/STIR/WHIR; this function is
    provided for the entropy-volume / counterexample discussion.
    """
    r = as_fraction(rho)
    _check_rate(r)
    return qary_entropy_inverse(1.0 - (r.numerator / r.denominator), q)


# ---------------------------------------------------------------------------
# log2 / bit helpers  (cancellation-safe)
# ---------------------------------------------------------------------------


def log2(x: Number) -> float:
    """``log2(x)`` for a positive number, exact-friendly for Fractions.

    For a :class:`Fraction` we compute ``log2(num) - log2(den)`` via
    :func:`math.log2` on the *integer* parts; this is a difference of two
    well-conditioned logs of exact integers (no cancellation of nearly-equal
    large floats), and Python's ``math.log2`` handles big ints exactly enough
    for our bit-budget purposes.
    """
    if isinstance(x, Fraction):
        if x <= 0:
            raise ValueError("log2 of non-positive Fraction")
        return math.log2(x.numerator) - math.log2(x.denominator)
    xf = float(x)
    if xf <= 0:
        raise ValueError("log2 of non-positive value")
    return math.log2(xf)


def bits(error: Number) -> float:
    """Security bits ``= -log2(error)`` for a positive error probability.

    This is the safe way to talk about security level: we take ``log2`` of the
    single positive quantity ``error`` and negate, never subtracting two large
    bit-counts.
    """
    return -log2(error)


def log2_pow2_minus_c(exp: int, c: int) -> float:
    """High-precision ``log2(2**exp - c)`` for field sizes like 2**31 - 1.

    Catastrophic-cancellation-safe: ``log2(2**exp - c) = exp + log2(1 - c/2**exp)``
    and the second term is computed with :func:`math.log1p` on the small ratio
    ``-c / 2**exp`` rather than by logging the huge integer ``2**exp - c`` (which
    would already have rounded).  Returns a value within a few ULPs of the truth.
    """
    if c == 0:
        return float(exp)
    ratio = -c / (2.0 ** exp)
    return exp + math.log1p(ratio) / math.log(2.0)


__all__ = [
    "as_fraction",
    "unique_decoding_radius",
    "capacity",
    "johnson_radius",
    "johnson_radius_exact",
    "johnson_radius_rational",
    "qary_entropy",
    "qary_entropy_inverse",
    "list_decoding_capacity_radius",
    "log2",
    "bits",
    "log2_pow2_minus_c",
]


if __name__ == "__main__":  # pragma: no cover - quick manual smoke
    for rho in (Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 16)):
        print(f"rho={rho}: UD={unique_decoding_radius(rho)} "
              f"J~={johnson_radius(rho):.6f} cap={capacity(rho)}")
