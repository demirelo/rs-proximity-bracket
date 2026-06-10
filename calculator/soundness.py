"""
soundness.py
============

Compose proximity-gap bounds into a FRI/STIR/WHIR-style **soundness error** as a
function of the query count ``t``, solve for the minimum integer ``t`` reaching a
security target (default ``2**-128``), and report **proof-size proxies**.

Soundness model (simplified, DOCUMENTED, and intentionally conservative)
------------------------------------------------------------------------
We use the composite shape from the research plan:

    soundness_error(t)  ~  eps_mca(C^{equiv m}, delta)              [batching / proximity]
                        +  |Lambda(C^{equiv m}, delta)| / |F|        [interleaved list, folding randomness]
                        +  (1 - delta)**t                            [t-fold query / sampling phase]

* ``eps_mca`` comes from the best applicable MCA bound (``bounds.best_mca_bound``,
  which already accounts for interleaving ``m`` via the union meta-bound).
* the list term comes from the best applicable list-size bound
  (``bounds.best_listsize_bound``), divided by ``|F| = 2**logF``.
* the query term ``(1 - delta)**t`` is the probability that all ``t`` independent
  queries miss the >= ``delta`` fraction of corrupted positions; this is the only
  ``t``-dependent term, so increasing ``t`` drives the total down to the
  ``t``-independent floor ``eps_mca + list/|F|``.

CAVEATS (read before trusting a number):
* This is a *composite proxy*, not any specific protocol's exact soundness.
  Real FRI/STIR/WHIR analyses have additional lower-order terms (per-round
  folding error, proof-of-work grinding, repetition structure, batching slack)
  and constants that differ.  Use this to compare *parameter regimes*, not to
  certify a deployed system.
* If the ``t``-independent floor already exceeds the target, NO finite ``t``
  helps; the solver reports infeasibility.  This is exactly what happens
  anywhere above the Johnson radius (floor = 1.0): in the OPEN band
  ``(J, r_E)`` the floor is 1.0 because no certified *positive* MCA bound
  applies (``unknown-beyond-johnson`` — "cannot certify", not a proven no-go),
  and at/above ``r_E`` it is 1.0 because the proximity gap provably fails
  (``proven-near-capacity-nogo``).
* All bound constants are provisional (see ``bounds.verify_flags()``).

Everything here is computed in log space where it matters; the only place we add
probabilities is the final composite (three terms of comparable-or-smaller
magnitude), which we do in linear space after confirming none individually
overflows.
"""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass, asdict
from fractions import Fraction
from typing import Callable, List, Optional, Sequence

import bounds as _bounds
from proximity_parameters import (
    as_fraction,
    bits,
    johnson_radius,
    unique_decoding_radius,
    capacity,
)

TARGET_BITS_DEFAULT = 128
TARGET_ERROR_DEFAULT = Fraction(1, 2 ** 128)


# ---------------------------------------------------------------------------
# Cost model parameters
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CostModel:
    """Proof-size cost model (a proxy, in BITS).

    A single query in a FRI/STIR/WHIR-style protocol opens a Merkle leaf and
    supplies its authentication path.  We model the per-query cost as

        path:    log_arity(n)  inner nodes  x  hash_size_bits          (Merkle path)
        leaf:    arity * field_elem_bits     (the opened coset / leaf alphabet)

    and the total proof size proxy as

        t * (path_bits + leaf_bits)  +  commitment_overhead_bits.

    All knobs are explicit so the caller can match a concrete instantiation.

    Defaults: 256-bit hash (e.g. a 2-to-1 compression with 256-bit digests),
    arity = 2 (binary Merkle tree / rate-1/2 folding), field element size taken
    from ``logF`` at call time, modest fixed commitment overhead.
    """

    hash_size_bits: int = 256
    arity: int = 2
    commitment_overhead_bits: int = 256          # roots + misc; a proxy
    field_elem_bits_override: Optional[int] = None  # else use ceil(logF)

    def merkle_depth(self, n: int) -> int:
        """Number of hashes on an authentication path = ceil(log_arity(n))."""
        if self.arity < 2:
            raise ValueError("arity must be >= 2")
        return math.ceil(math.log(n, self.arity))

    def field_elem_bits(self, logF: float) -> int:
        if self.field_elem_bits_override is not None:
            return self.field_elem_bits_override
        return math.ceil(logF)

    def per_query_bits(self, n: int, logF: float) -> float:
        path = self.merkle_depth(n) * self.hash_size_bits
        leaf = self.arity * self.field_elem_bits(logF)
        return float(path + leaf)

    def proof_size_bits(self, t: int, n: int, logF: float) -> float:
        return t * self.per_query_bits(n, logF) + self.commitment_overhead_bits


DEFAULT_COST = CostModel()


# ---------------------------------------------------------------------------
# Composite soundness error
# ---------------------------------------------------------------------------


@dataclass
class SoundnessTerms:
    """The three composite terms at a given ``t`` plus their sum and bits."""

    eps_mca: float
    list_over_field: float
    query_term: float
    total: float
    total_bits: float
    mca_bound_name: Optional[str]
    list_bound_name: Optional[str]
    floor: float            # eps_mca + list_over_field  (t-independent)
    floor_bits: float


def _query_term(delta: float, t: int) -> float:
    """``(1 - delta)**t`` computed stably for large ``t``.

    For large ``t`` this underflows to 0.0 in linear space; we therefore compute
    it as ``2**(t * log2(1 - delta))`` and only materialize the linear value,
    which is exactly what the composite sum needs (a true 0.0 floor is fine).
    """
    if not (0.0 < delta < 1.0):
        raise ValueError(f"delta must be in (0,1); got {delta}")
    log2_term = t * math.log2(1.0 - delta)
    if log2_term < -1100:          # below double's exponent range -> effectively 0
        return 0.0
    return 2.0 ** log2_term


def soundness_terms(rho, logF, n, delta, m=1, t: int = 0,
                    mca_bound: Optional[_bounds.Bound] = None,
                    list_bound: Optional[_bounds.Bound] = None,
                    eta: Optional[float] = None,
                    field_type: Optional[str] = None) -> SoundnessTerms:
    """Evaluate the composite soundness error and its parts at query count ``t``.

    ``mca_bound`` / ``list_bound`` may be supplied explicitly; otherwise the best
    applicable registered bound is selected automatically.  ``eta`` (Johnson
    slack) is derived from ``delta`` as ``J(rho) - delta`` when not given and
    when ``delta`` lies below the Johnson radius (so the Johnson-regime bounds
    can fire); pass it explicitly to override.

    ``field_type`` (``"prime"`` | ``"extension"`` | ``None``) is forwarded to the
    bounds registry for back-compat, but the two capacity-region bounds are now
    **FIELD-AGNOSTIC**: the proven near-capacity no-go applies identically over
    prime AND genuine odd-characteristic extension fields (the counterexample is
    a characteristic-zero cyclotomic invariant — ``n2-verdict.md``), so the
    whole ``(Johnson, capacity)`` window is vacuous (no security claimable)
    regardless of ``field_type``.  The window is split: in the OPEN band
    ``(J, r_E)`` the applicable MCA bound is ``unknown-beyond-johnson`` (vacuous
    1.0 = "cannot certify", NOT a proven no-go), and at/above ``r_E`` it is
    ``proven-near-capacity-nogo`` (a genuine proven no-go).  Either way the
    ``eps_mca`` floor is 1.0, so the *behavior* (and every number) is unchanged;
    only the reported ``mca_bound_name`` distinguishes the two.
    """
    rho_f = as_fraction(rho)
    delta = float(delta)
    J = johnson_radius(rho_f)
    if eta is None and delta < J:
        eta = J - delta
    # else leave eta as given/None (UD regime doesn't need it).

    # --- MCA term ---
    # Use the *interleaved* selector so the batching factor m is reflected
    # (best_mca_bound would otherwise pick the m-independent single-code bound
    # whenever m*single > single, silently dropping the batching penalty).
    if mca_bound is None:
        mca_bound = _bounds.interleaved_mca(rho_f, n, logF, delta, eta=eta, m=m,
                                            field_type=field_type)
    if mca_bound is None:
        eps_mca = 1.0
        mca_name = None
    else:
        eps_mca = mca_bound.value(rho_f, n, logF, delta, eta=eta, m=m,
                                  field_type=field_type)
        mca_name = mca_bound.name

    # --- list term ---
    if list_bound is None:
        list_bound = _bounds.best_listsize_bound(rho_f, n, logF, delta, eta=eta,
                                                 m=m, field_type=field_type)
    if list_bound is None:
        list_over_field = 0.0      # no list term contributed (e.g. UD-only model)
        list_name = None
    else:
        list_size = list_bound.value(rho_f, n, logF, delta, eta=eta, m=m,
                                     field_type=field_type)
        # |Lambda| / |F| = 2**(log2(list_size) - logF)
        list_over_field = min(1.0, 2.0 ** (math.log2(list_size) - float(logF)))
        list_name = list_bound.name

    # --- query term ---
    qterm = _query_term(delta, t)

    floor = min(1.0, eps_mca + list_over_field)
    total = min(1.0, floor + qterm)
    return SoundnessTerms(
        eps_mca=eps_mca,
        list_over_field=list_over_field,
        query_term=qterm,
        total=total,
        total_bits=(bits(total) if total > 0 else math.inf),
        mca_bound_name=mca_name,
        list_bound_name=list_name,
        floor=floor,
        floor_bits=(bits(floor) if floor > 0 else math.inf),
    )


def soundness_error(rho, logF, n, delta, m, t, **kw) -> float:
    """Just the composite soundness error (a probability) at query count ``t``."""
    return soundness_terms(rho, logF, n, delta, m=m, t=t, **kw).total


# ---------------------------------------------------------------------------
# Solve for minimum t  (bisection)
# ---------------------------------------------------------------------------


@dataclass
class QuerySolution:
    feasible: bool
    t: Optional[int]                 # minimum t achieving the target (if feasible)
    achieved_bits: float             # security bits at t (inf if total==0)
    floor_bits: float                # t-independent floor in bits
    target_bits: float
    note: str = ""


def min_query_count(rho, logF, n, delta, m=1,
                    target_error: Fraction = TARGET_ERROR_DEFAULT,
                    t_max: int = 1 << 20, **kw) -> QuerySolution:
    """Minimum integer ``t`` with ``soundness_error(t) <= target_error``.

    Strategy:
    1. Compute the ``t``-independent floor.  If ``floor > target`` the problem is
       infeasible for *any* ``t`` (return ``feasible=False``).
    2. Otherwise the total ``= floor + (1 - delta)**t`` is strictly decreasing in
       ``t`` toward ``floor``, so bisect on ``t in [0, t_max]`` for the smallest
       ``t`` meeting the target.  The solver guarantees the standard anchor
       ``error(t) <= target < error(t - 1)``.

    Monotonicity note: because only the query term depends on ``t`` and it is
    monotone non-increasing, ``t`` itself is monotone non-increasing in ``delta``
    (larger ``delta`` => each query is more likely to catch corruption) over any
    range where the floor stays below target — this is asserted in the tests.
    """
    target = float(target_error)
    target_b = bits(target_error)

    terms0 = soundness_terms(rho, logF, n, delta, m=m, t=0, **kw)
    floor = terms0.floor
    floor_b = terms0.floor_bits

    if floor > target:
        return QuerySolution(
            feasible=False, t=None,
            achieved_bits=floor_b, floor_bits=floor_b, target_bits=target_b,
            note=(f"infeasible: t-independent floor = {floor_b:.2f} bits "
                  f"< target {target_b:.2f} bits (eps_mca + list/|F| too large "
                  f"at delta={float(delta):.6f}). No finite t helps."),
        )

    # floor <= target: a sufficiently large t works.  Check upper bound first.
    def err(t: int) -> float:
        return soundness_error(rho, logF, n, delta, m, t, **kw)

    if err(t_max) > target:
        return QuerySolution(
            feasible=False, t=None,
            achieved_bits=bits(err(t_max)), floor_bits=floor_b,
            target_bits=target_b,
            note=(f"target not reached even at t_max={t_max}; "
                  f"increase t_max (delta={float(delta):.6f} very small)."),
        )

    # Bisection for the smallest t in [0, t_max] meeting the target.
    lo, hi = 0, t_max
    if err(0) <= target:
        t_star = 0
    else:
        while lo + 1 < hi:
            mid = (lo + hi) // 2
            if err(mid) <= target:
                hi = mid
            else:
                lo = mid
        t_star = hi

    achieved = err(t_star)
    return QuerySolution(
        feasible=True, t=t_star,
        achieved_bits=(bits(achieved) if achieved > 0 else math.inf),
        floor_bits=floor_b, target_bits=target_b,
        note=("ok" if t_star > 0 else "floor alone already meets target (t=0)"),
    )


# ---------------------------------------------------------------------------
# Combined row: t and proof-size at one (rho, logF, n, delta, m)
# ---------------------------------------------------------------------------


@dataclass
class ParamRow:
    rho: str
    logF: float
    log2n: int
    n: int
    m: int
    delta: float
    eta: Optional[float]
    feasible: bool
    t: Optional[int]
    achieved_bits: float
    floor_bits: float
    proof_size_bits: Optional[float]
    proof_size_kb: Optional[float]
    mca_bound: Optional[str]
    list_bound: Optional[str]
    note: str


def evaluate_point(rho, logF, log2n, delta, m=1,
                   cost: CostModel = DEFAULT_COST,
                   target_error: Fraction = TARGET_ERROR_DEFAULT,
                   **kw) -> ParamRow:
    """Full evaluation at one parameter point: solve t, then size the proof."""
    n = 1 << int(log2n)
    rho_f = as_fraction(rho)
    J = johnson_radius(rho_f)
    eta = (J - float(delta)) if float(delta) < J else None

    sol = min_query_count(rho_f, logF, n, delta, m=m,
                          target_error=target_error, **kw)
    terms = soundness_terms(rho_f, logF, n, delta, m=m,
                            t=(sol.t or 0), **kw)

    if sol.feasible and sol.t is not None:
        ps_bits = cost.proof_size_bits(sol.t, n, logF)
        ps_kb = ps_bits / 8 / 1024
    else:
        ps_bits = None
        ps_kb = None

    return ParamRow(
        rho=str(rho_f), logF=float(logF), log2n=int(log2n), n=n, m=int(m),
        delta=float(delta), eta=eta,
        feasible=sol.feasible, t=sol.t,
        achieved_bits=sol.achieved_bits, floor_bits=sol.floor_bits,
        proof_size_bits=ps_bits, proof_size_kb=ps_kb,
        mca_bound=terms.mca_bound_name, list_bound=terms.list_bound_name,
        note=sol.note,
    )


# ---------------------------------------------------------------------------
# Sweep delta from small up to capacity
# ---------------------------------------------------------------------------


def sweep_delta(rho, logF, log2n, m=1,
                delta_min: float = 0.01,
                n_points: int = 60,
                cost: CostModel = DEFAULT_COST,
                target_error: Fraction = TARGET_ERROR_DEFAULT,
                include_capacity: bool = True,
                **kw) -> List[ParamRow]:
    """Sweep ``delta`` from ``delta_min`` up toward capacity ``1 - rho``.

    Returns one :class:`ParamRow` per sampled ``delta``.  The upper end is just
    below capacity (``capacity - tiny``) since at/above capacity no positive
    bound exists.  Useful for t(delta) and proof-size(delta) curves / heatmaps.
    """
    rho_f = as_fraction(rho)
    cap = float(capacity(rho_f))
    hi = cap - 1e-6 if include_capacity else johnson_radius(rho_f) - 1e-6
    hi = max(hi, delta_min)
    rows: List[ParamRow] = []
    if n_points < 2:
        n_points = 2
    for i in range(n_points):
        delta = delta_min + (hi - delta_min) * i / (n_points - 1)
        rows.append(evaluate_point(rho_f, logF, log2n, delta, m=m,
                                   cost=cost, target_error=target_error, **kw))
    return rows


# ---------------------------------------------------------------------------
# CSV writer
# ---------------------------------------------------------------------------

CSV_FIELDS = list(ParamRow.__annotations__.keys())


def write_rows_csv(rows: Sequence[ParamRow], path: str) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))
    return os.path.abspath(path)


__all__ = [
    "CostModel", "DEFAULT_COST",
    "SoundnessTerms", "soundness_terms", "soundness_error",
    "QuerySolution", "min_query_count",
    "ParamRow", "evaluate_point", "sweep_delta",
    "write_rows_csv", "CSV_FIELDS",
    "TARGET_ERROR_DEFAULT", "TARGET_BITS_DEFAULT",
]


if __name__ == "__main__":  # pragma: no cover
    # quick smoke: one point a hair below Johnson for rho=1/4 over a 128-bit field
    rho = Fraction(1, 4)
    row = evaluate_point(rho, logF=128, log2n=20,
                         delta=johnson_radius(rho) - 0.02, m=1)
    import pprint
    pprint.pprint(asdict(row))
