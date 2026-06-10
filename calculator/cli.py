#!/usr/bin/env python3
"""
cli.py
======

Command-line front-end for the Proximity Prize parameter calculator.

Subcommands
-----------
* ``tables``  — print, for each rate ``rho in {1/2,1/4,1/8,1/16}`` and a couple
                of field sizes, the three radii and the required query count /
                proof-size proxy at a chosen ``delta`` (default ``delta = Johnson
                - eta`` with a small ``eta``).  Also prints the bare radii table
                and the VERIFY surface.
* ``sweep``   — dump a CSV over ``(rho, delta, logF, log2n, m)`` of ``t`` and
                proof size to ``calculator/out/`` for later heatmaps.
* ``bounds``  — list the registered proximity-gap bounds with provenance.
* ``delta-star`` — write the delta*_C provable-bracket tables (the headline
                numeric deliverable) to ``out/delta_star_tables.md`` and print a
                short ρ=1/2 summary.

Examples
--------
    python3 cli.py tables
    python3 cli.py tables --eta 0.05 --log2n 24
    python3 cli.py sweep --rhos 1/2,1/4 --logFs 31,64,128 --log2ns 20,24 --m 1
    python3 cli.py bounds
    python3 cli.py delta-star
"""

from __future__ import annotations

import argparse
import os
import sys
from fractions import Fraction
from typing import List

import bounds as _bounds
from proximity_parameters import (
    as_fraction,
    capacity,
    johnson_radius,
    johnson_radius_rational,
    unique_decoding_radius,
    log2_pow2_minus_c,
)
from soundness import (
    CostModel,
    DEFAULT_COST,
    TARGET_ERROR_DEFAULT,
    TARGET_BITS_DEFAULT,
    evaluate_point,
    sweep_delta,
    write_rows_csv,
)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "out")

# Named fields the user is likely to care about, expressed as log2|F|.
# The three "exact" small fields are computed cancellation-safely.
NAMED_FIELDS = {
    "Mersenne31": log2_pow2_minus_c(31, 1),               # 2^31 - 1
    "BabyBear": log2_pow2_minus_c(31, 2 ** 27) ,          # 2^31 - 2^27 + 1 ~ 15*2^27+1
    "Goldilocks": log2_pow2_minus_c(64, 2 ** 32 - 1),     # 2^64 - 2^32 + 1
    "prime128": 128.0,
    "prime256": 256.0,
}

DEFAULT_RATES = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 16)]


# ---------------------------------------------------------------------------
# parsing helpers
# ---------------------------------------------------------------------------

def _parse_fraction_list(s: str) -> List[Fraction]:
    return [as_fraction(tok.strip()) for tok in s.split(",") if tok.strip()]


def _parse_float_list(s: str) -> List[float]:
    return [float(tok.strip()) for tok in s.split(",") if tok.strip()]


def _parse_int_list(s: str) -> List[int]:
    return [int(tok.strip()) for tok in s.split(",") if tok.strip()]


# ---------------------------------------------------------------------------
# radii table
# ---------------------------------------------------------------------------

def print_radii_table(rates: List[Fraction]) -> None:
    print("=" * 72)
    print("Reed-Solomon radii (relative Hamming distance)")
    print("=" * 72)
    print(f"{'rho':>6} | {'unique-dec (1-r)/2':>20} | {'Johnson 1-sqrt(r)':>20} "
          f"| {'capacity 1-r':>14}")
    print("-" * 72)
    for r in rates:
        ud = unique_decoding_radius(r)
        jr = johnson_radius_rational(r)
        cap = capacity(r)
        print(f"{str(r):>6} | {str(ud):>9} = {float(ud):.6f} "
              f"| {float(jr):>20.6f} | {str(cap):>6} = {float(cap):.4f}")
    print("-" * 72)
    print("ordering check: (1-r)/2 < 1-sqrt(r) < 1-r  for all rho in (0,1)")
    print()


# ---------------------------------------------------------------------------
# tables subcommand
# ---------------------------------------------------------------------------

def cmd_tables(args: argparse.Namespace) -> int:
    rates = _parse_fraction_list(args.rhos) if args.rhos else DEFAULT_RATES

    print_radii_table(rates)

    # Field set for the t / proof-size table.
    if args.logFs:
        fields = [("logF=%g" % x, x) for x in _parse_float_list(args.logFs)]
    else:
        fields = [(name, lf) for name, lf in NAMED_FIELDS.items()
                  if name in ("Goldilocks", "prime128", "prime256")]

    cost = CostModel(hash_size_bits=args.hash_bits, arity=args.arity)

    print("=" * 100)
    print(f"Query count t and proof-size proxy at delta = Johnson(rho) - eta,  "
          f"eta = {args.eta}")
    print(f"target = 2^-{TARGET_BITS_DEFAULT}, m = {args.m}, log2(n) = {args.log2n}, "
          f"hash = {args.hash_bits} bits, arity = {args.arity}")
    print("=" * 100)
    header = (f"{'rho':>6} | {'field':>11} | {'logF':>7} | {'delta':>8} | "
              f"{'t':>7} | {'sec_bits':>8} | {'floor_b':>8} | "
              f"{'proof_KB':>9} | {'mca bound':>16}")
    print(header)
    print("-" * len(header))
    for r in rates:
        J = float(johnson_radius(r))
        delta = J - args.eta
        if delta <= 0:
            delta = J / 2  # guard tiny rates / large eta
        for fname, lf in fields:
            row = evaluate_point(r, logF=lf, log2n=args.log2n, delta=delta,
                                 m=args.m, cost=cost,
                                 target_error=TARGET_ERROR_DEFAULT)
            t_str = str(row.t) if row.feasible else "INF"
            sb = f"{row.achieved_bits:.1f}" if row.feasible else "-"
            pk = f"{row.proof_size_kb:.2f}" if row.proof_size_kb is not None else "-"
            print(f"{str(r):>6} | {fname:>11} | {lf:>7.2f} | {delta:>8.5f} | "
                  f"{t_str:>7} | {sb:>8} | {row.floor_bits:>8.1f} | "
                  f"{pk:>9} | {str(row.mca_bound):>16}")
    print("-" * len(header))
    print("INF in t => the t-independent floor (eps_mca + |Lambda|/|F|) already "
          "exceeds the target;\nno query count can rescue it (see 'floor_b').")
    print()

    if not args.no_verify:
        _print_verify_surface()
    return 0


def _print_verify_surface() -> None:
    print("=" * 72)
    print("VERIFY surface — literature constants/exponents still to be filled")
    print("=" * 72)
    for i, vf in enumerate(_bounds.verify_flags(), 1):
        print(f"[{i}] {vf['name']}  (source: {vf['source_tag']})")
        print(f"    formula: {vf['formula']}")
        print(f"    {vf['notes']}")
    print()
    print("Tunable constants (patch via bounds.CONSTANTS or module globals):")
    print("   ", _bounds.CONSTANTS)
    print()


# ---------------------------------------------------------------------------
# sweep subcommand
# ---------------------------------------------------------------------------

def cmd_sweep(args: argparse.Namespace) -> int:
    rates = _parse_fraction_list(args.rhos) if args.rhos else DEFAULT_RATES
    logFs = _parse_float_list(args.logFs) if args.logFs else [64.0, 128.0]
    log2ns = _parse_int_list(args.log2ns) if args.log2ns else [20]
    ms = _parse_int_list(args.ms) if args.ms else [1]

    cost = CostModel(hash_size_bits=args.hash_bits, arity=args.arity)
    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = (args.out if args.out
                else os.path.join(OUT_DIR, "sweep.csv"))

    all_rows = []
    for r in rates:
        for lf in logFs:
            for ln in log2ns:
                for mm in ms:
                    rows = sweep_delta(r, logF=lf, log2n=ln, m=mm,
                                       delta_min=args.delta_min,
                                       n_points=args.points, cost=cost,
                                       target_error=TARGET_ERROR_DEFAULT,
                                       include_capacity=not args.up_to_johnson)
                    all_rows.extend(rows)

    path = write_rows_csv(all_rows, out_path)
    feas = sum(1 for r in all_rows if r.feasible)
    print(f"wrote {len(all_rows)} rows ({feas} feasible) to:\n  {path}")
    print(f"grid: rhos={[str(x) for x in rates]} logFs={logFs} "
          f"log2ns={log2ns} ms={ms} points/curve={args.points}")
    print("columns:", ", ".join(
        __import__("soundness").CSV_FIELDS))
    return 0


# ---------------------------------------------------------------------------
# bounds subcommand
# ---------------------------------------------------------------------------

def cmd_bounds(_args: argparse.Namespace) -> int:
    print("Registered proximity-gap / MCA / list-size bounds")
    print("=" * 72)
    for b in _bounds.all_bounds():
        print(b.describe())
        print()
    return 0


# ---------------------------------------------------------------------------
# delta-star subcommand
# ---------------------------------------------------------------------------

def cmd_delta_star(args: argparse.Namespace) -> int:
    import delta_star as ds
    path = args.out if args.out else None
    written = ds.write_tables(path)
    rows = ds.all_rows()
    # Quick console summary: the rho=1/2 brackets for M31 and a 256-bit field.
    print(f"wrote delta*_C tables to:\n  {written}\n")
    print("ρ=1/2 provable brackets [delta_safe, delta_unsafe] (n = 2^20; "
          "δ_unsafe field-agnostic — prime = extension):")
    for row in rows:
        if row.rho == "1/2" and row.log2n == 20 and row.field in (
                "Mersenne31", "prime256", "ext256"):
            safe = (f"{row.delta_safe:.5f}" if row.safe_feasible
                    else "infeasible")
            unsafe = (f"{row.delta_unsafe:.5f}"
                      if row.delta_unsafe is not None else "OPEN")
            t_s = (str(row.t_at_safe) if row.t_at_safe is not None
                   else ("INF" if row.safe_feasible else "—"))
            print(f"  {row.field:>12} ({row.field_type:>9}): "
                  f"delta_safe={safe:>10}  t@safe={t_s:>5}  "
                  f"delta_unsafe={unsafe:>10}")
    print()
    print("Johnson(1/2) = 0.29289; capacity(1/2) = 0.5 -> at Johnson "
          "(1-delta)^t = (1/sqrt2)^t, so t=128 gives only 2^-64 (64 bits). "
          "See the markdown for the full family + headline.")
    return 0


# ---------------------------------------------------------------------------
# argparse wiring
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="proximity-calc",
        description="Reproducible parameter calculator for RS proximity gaps.")
    sub = p.add_subparsers(dest="cmd", required=True)

    # tables
    pt = sub.add_parser("tables", help="print radii + t/proof-size tables")
    pt.add_argument("--rhos", type=str, default="",
                    help="comma rates, e.g. '1/2,1/4,1/8,1/16' (default these)")
    pt.add_argument("--logFs", type=str, default="",
                    help="comma log2|F| values; default Goldilocks/128/256")
    pt.add_argument("--log2n", type=int, default=20, help="log2 of block length n")
    pt.add_argument("--eta", type=float, default=0.05,
                    help="Johnson slack: delta = Johnson(rho) - eta (default .05)")
    pt.add_argument("--m", type=int, default=1, help="interleaving parameter")
    pt.add_argument("--hash-bits", type=int, default=256, help="Merkle hash size")
    pt.add_argument("--arity", type=int, default=2, help="Merkle/folding arity")
    pt.add_argument("--no-verify", action="store_true",
                    help="suppress the VERIFY surface print")
    pt.set_defaults(func=cmd_tables)

    # sweep
    ps = sub.add_parser("sweep", help="dump CSV over (rho,delta,logF,log2n,m)")
    ps.add_argument("--rhos", type=str, default="")
    ps.add_argument("--logFs", type=str, default="")
    ps.add_argument("--log2ns", type=str, default="")
    ps.add_argument("--ms", type=str, default="")
    ps.add_argument("--points", type=int, default=60,
                    help="delta samples per curve (default 60)")
    ps.add_argument("--delta-min", type=float, default=0.01)
    ps.add_argument("--up-to-johnson", action="store_true",
                    help="cap the sweep at the Johnson radius instead of capacity")
    ps.add_argument("--hash-bits", type=int, default=256)
    ps.add_argument("--arity", type=int, default=2)
    ps.add_argument("--out", type=str, default="",
                    help="output CSV path (default calculator/out/sweep.csv)")
    ps.set_defaults(func=cmd_sweep)

    # bounds
    pb = sub.add_parser("bounds", help="list registered bounds + provenance")
    pb.set_defaults(func=cmd_bounds)

    # delta-star
    pd = sub.add_parser(
        "delta-star",
        help="write the delta*_C provable-bracket tables (out/delta_star_tables.md)")
    pd.add_argument("--out", type=str, default="",
                    help="output markdown path (default out/delta_star_tables.md)")
    pd.set_defaults(func=cmd_delta_star)

    return p


def main(argv: List[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
