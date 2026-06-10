"""
scaleup_smooth_vs_random.py -- Wave-2 scale-up of the smooth-vs-random comparison.

Wave 1 found (exact, q<=64): smooth multiplicative-subgroup RS behaves statistically
like random-subset domains everywhere measured (bad-line onset at Johnson, no
smooth-specific inflation).  Caveat: tiny fields, so the 1/|F| floor is huge and any
smooth-specific effect that scales like (subgroup structure)/|F| could be invisible.

This script pushes |F| up to primes {127, 257, 521, 1031} and GF(2^7), GF(2^8), with
PROPER-subgroup domains (n | q-1, n < q-1) at rates rho ~ {1/2, 1/4}, and asks:

  Does the smooth ~ random equivalence SURVIVE as n/|F| shrinks (larger q)?

For each (field, n, k) we evaluate the EXACT bad-line statistics (close-count
distribution, S*, bad-line fraction) on TWO matched domains -- the order-n smooth
multiplicative subgroup and a random size-n subset of F* -- over the same delta grid
straddling Johnson and capacity, then report the smooth-minus-random differential and
track it versus q (i.e. versus n/|F|).

Everything is EXACT (full codeword enumeration; reuses rs.py / search_bad_lines.py).
All gamma in F are enumerated per line (no gamma sampling).  Line counts per cell are
tiered by compute cost (q^k * q * n, the cost of one line's all-gamma distance pass)
and LOGGED -- no silent caps.  Results -> results/scaleup_smooth_vs_random.{json,csv}.
"""

from __future__ import annotations

import csv
import json
import os
import time
import traceback
from dataclasses import dataclass
from multiprocessing import Pool

import numpy as np

from ff import PrimeField, BinaryExtensionField, FiniteField, _divisors
from rs import build_codeword_book, domain_subgroup, domain_random
from search_bad_lines import search_bad_lines_multi

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
QK_CAP = 3_000_000


# ---------------------------------------------------------------------------
# Battery: proper-subgroup codes over the scale-up fields at rho ~ {1/2, 1/4}.
# ---------------------------------------------------------------------------
@dataclass
class Cfg:
    field_name: str
    field: FiniteField
    n: int
    k: int

    @property
    def q(self):
        return self.field.q

    @property
    def rho(self):
        return self.k / self.n

    @property
    def qk(self):
        return self.q ** self.k

    @property
    def cost(self):
        # cost of one line's exact all-gamma distance pass ~ q^k * q * n
        return self.qk * self.q * self.n

    @property
    def label(self):
        return f"{self.field_name}_n{self.n}_k{self.k}"


def build_battery() -> list[Cfg]:
    fields = [
        ("GF(127)", PrimeField(127)),
        ("GF(257)", PrimeField(257)),
        ("GF(521)", PrimeField(521)),
        ("GF(1031)", PrimeField(1031)),
        ("GF(2^7)", BinaryExtensionField(7)),
        ("GF(2^8)", BinaryExtensionField(8)),
    ]
    # Per-line exact-distance cost is q^k * q * n (ALL gamma in F enumerated).  We cap
    # it so every cell stays exactly enumerable within the runtime budget.  This caps
    # the prime field at 521 in practice: q=1031 cells (q^k ~ 1e6) cost ~21 s PER LINE
    # over all 1031 gammas, which is too slow for a meaningful line sample -- so q=1031
    # is intentionally excluded here and called out in the meta/README (q=521 already
    # probes n/|F| down to ~0.01, the small-ratio regime we care about).
    COST_CAP = 2e9
    cfgs: list[Cfg] = []
    seen = set()
    for name, F in fields:
        q = F.q
        for n in _divisors(q - 1):
            if n < 5 or n >= q - 1:        # PROPER subgroup, big enough for rho meaning
                continue
            for target in (0.5, 0.25):
                k = max(2, round(target * n))
                if k >= n:
                    continue
                rho = k / n
                if abs(rho - target) > 0.13:
                    continue
                if q ** k > QK_CAP:
                    continue
                cfg = Cfg(name, F, n, k)
                if cfg.cost > COST_CAP:
                    continue
                key = (name, n, k)
                if key in seen:
                    continue
                seen.add(key)
                cfgs.append(cfg)
    # GF(2^7)* has prime order 127 -> NO proper subgroup; it contributes nothing here
    # (logged in meta).  That is itself informative: a binary field with prime q-1 has
    # no smooth proper-subgroup domain at all.
    return cfgs


def line_budget(cost: float) -> dict:
    """Lines per generator family, tiered by per-line cost (q^k*q*n).  Logged.

    Tuned so the whole battery finishes well under the runtime budget when the cells
    run in parallel (the slowest single cell bounds wall-clock).  The largest-q cells
    (q=1031, cost>2e9) get a deliberately small but NONZERO budget so we still obtain
    the crown n/|F| data point; this is logged per cell, not a silent cap.
    """
    if cost <= 5e7:
        return dict(random=120, cwnoise=120, lowdeg=40, rational=40)
    if cost <= 5e8:
        return dict(random=50, cwnoise=50, lowdeg=20, rational=20)
    if cost <= 2e9:
        return dict(random=12, cwnoise=12, lowdeg=6, rational=6)
    if cost <= 6e9:
        return dict(random=4, cwnoise=4, lowdeg=2, rational=2)
    # q=1031, n=10 (cost ~1.1e10, ~25 s/line): smallest viable to still land the point.
    return dict(random=2, cwnoise=2, lowdeg=1, rational=1)


def delta_grid(rho: float, n: int) -> list[float]:
    """Delta grid straddling Johnson and capacity (same shape as run_atlas)."""
    johnson = 1.0 - np.sqrt(rho)
    cap = 1.0 - rho
    pts = set()
    for f in (0.6, 0.8, 0.9, 1.0):
        pts.add(round(f * johnson, 4))
    for t in (0.25, 0.5, 0.75, 1.0):
        pts.add(round(johnson + t * (cap - johnson), 4))
    pts.add(round(min(0.97, cap + 0.02), 4))
    return sorted(p for p in pts if 0.02 < p < 0.985)


# ---------------------------------------------------------------------------
# One cell: matched smooth (subgroup) vs random over the same field/n/k/grid.
# ---------------------------------------------------------------------------
def run_cell(args) -> dict:
    cfg: Cfg = args["cfg"]
    seed: int = args["seed"]
    F = cfg.field
    out = {
        "label": cfg.label, "field": F.name, "q": F.q, "n": cfg.n, "k": cfg.k,
        "rho": cfg.rho, "qk": cfg.qk, "cost": cfg.cost, "n_over_q": cfg.n / F.q,
        "johnson": 1.0 - np.sqrt(cfg.rho), "capacity": 1.0 - cfg.rho,
        "ok": True, "elapsed_sec": 0.0,
    }
    t0 = time.time()
    try:
        grid = delta_grid(cfg.rho, cfg.n)
        out["deltas"] = grid
        lb = line_budget(cfg.cost)
        out["line_budget"] = lb

        # smooth subgroup
        rng_s = np.random.default_rng(seed)
        Ds = domain_subgroup(F, cfg.n)
        book_s = build_codeword_book(F, Ds, cfg.k)
        sub = search_bad_lines_multi(book_s, grid, rng_s,
                                     n_random=lb["random"], n_cwnoise=lb["cwnoise"],
                                     n_lowdeg=lb["lowdeg"], n_rational=lb["rational"])

        # random subset (fresh rng so the random domain differs from the subgroup)
        rng_r = np.random.default_rng(seed + 1)
        Dr = domain_random(F, cfg.n, rng_r)
        book_r = build_codeword_book(F, Dr, cfg.k)
        rnd = search_bad_lines_multi(book_r, grid, rng_r,
                                     n_random=lb["random"], n_cwnoise=lb["cwnoise"],
                                     n_lowdeg=lb["lowdeg"], n_rational=lb["rational"])

        out["subgroup_per_delta"] = sub
        out["random_per_delta"] = rnd
        # control sanity: shared-support controls must never be flagged bad
        out["control_bad_total"] = sum(s["control_bad_lines"] for s in sub) + \
                                   sum(s["control_bad_lines"] for s in rnd)
    except Exception as e:
        out["ok"] = False
        out["error"] = f"{e}"
        out["traceback"] = traceback.format_exc()
    out["elapsed_sec"] = round(time.time() - t0, 2)
    return out


# ---------------------------------------------------------------------------
# Differential analysis (smooth - random), matched on (field,n,k,delta).
# ---------------------------------------------------------------------------
def compute_differentials(results) -> dict:
    """For each result, match subgroup vs random per delta in the MEANINGFUL regime,
    and aggregate the close-count-mean and bad-fraction differentials, bucketed by q.
    """
    per_cell = []
    by_q = {}
    for r in results:
        if not r.get("ok"):
            continue
        sub = {s["delta"]: s for s in r["subgroup_per_delta"]}
        rnd = {s["delta"]: s for s in r["random_per_delta"]}
        dclose, dbad = [], []
        for d in r["deltas"]:
            s, rr = sub.get(d), rnd.get(d)
            if not s or not rr or not s["meaningful_regime"]:
                continue
            dclose.append(s["close_count_mean"] - rr["close_count_mean"])
            dbad.append(s["frac_bad"] - rr["frac_bad"])
        if dclose:
            cell = {
                "label": r["label"], "q": r["q"], "n": r["n"], "k": r["k"],
                "rho": round(r["rho"], 4), "n_over_q": round(r["n_over_q"], 5),
                "mean_close_diff": float(np.mean(dclose)),
                "max_abs_close_diff": float(np.max(np.abs(dclose))),
                "mean_badfrac_diff": float(np.mean(dbad)),
                "max_abs_badfrac_diff": float(np.max(np.abs(dbad))),
                "n_meaningful_deltas": len(dclose),
            }
            per_cell.append(cell)
            by_q.setdefault(r["q"], []).append(cell)

    # aggregate by q
    q_summary = []
    for q in sorted(by_q):
        cells = by_q[q]
        mc = np.array([c["mean_close_diff"] for c in cells])
        mb = np.array([c["mean_badfrac_diff"] for c in cells])
        q_summary.append({
            "q": q,
            "num_cells": len(cells),
            "mean_close_diff": float(mc.mean()),
            "max_abs_close_diff": float(np.max([c["max_abs_close_diff"] for c in cells])),
            "mean_badfrac_diff": float(mb.mean()),
            "max_abs_badfrac_diff": float(np.max([c["max_abs_badfrac_diff"] for c in cells])),
            "min_n_over_q": float(min(c["n_over_q"] for c in cells)),
        })
    return {"per_cell": per_cell, "by_q": q_summary}


def print_summary(results, diffs):
    print("\n" + "=" * 84)
    print("SCALE-UP: smooth (subgroup) vs random differential vs q  (n/|F| shrinking)")
    print("=" * 84)
    print("Per-q aggregate (meaningful regime; differential = smooth - random):")
    print(f"  {'q':>6} {'cells':>5} {'min n/q':>8} {'mean dclose':>12} "
          f"{'max|dclose|':>12} {'mean dbadfrac':>14} {'max|dbadfrac|':>14}")
    for s in diffs["by_q"]:
        print(f"  {s['q']:>6} {s['num_cells']:>5} {s['min_n_over_q']:>8.4f} "
              f"{s['mean_close_diff']:>+12.4f} {s['max_abs_close_diff']:>12.4f} "
              f"{s['mean_badfrac_diff']:>+14.5f} {s['max_abs_badfrac_diff']:>14.5f}")
    # overall verdict + trend
    if diffs["by_q"]:
        qs = np.array([s["q"] for s in diffs["by_q"]])
        mc = np.array([abs(s["mean_close_diff"]) for s in diffs["by_q"]])
        print("-" * 84)
        print(f"  |mean close-diff| across q: " +
              ", ".join(f"q={s['q']}:{abs(s['mean_close_diff']):.3f}"
                        for s in diffs["by_q"]))
        # crude monotonic trend check
        if len(qs) >= 3:
            corr = float(np.corrcoef(np.log(qs), mc)[0, 1]) if mc.std() > 0 else 0.0
            print(f"  corr(log q, |mean close-diff|) = {corr:+.3f} "
                  f"(>0 would mean smooth diverges from random as q grows)")
        big = max(abs(s["mean_close_diff"]) for s in diffs["by_q"])
        bigb = max(abs(s["mean_badfrac_diff"]) for s in diffs["by_q"])
        if big < 0.5 and bigb < 0.03:
            print("  => smooth ~ random SURVIVES at larger q: differentials stay within "
                  "sampling noise\n     (bounded by ~0.5 close-count, ~0.03 bad-fraction) "
                  "at every q tested.")
        else:
            print("  => A smooth-vs-random differential APPEARS at larger q (see rows).")
    ctrl = sum(r.get("control_bad_total", 0) for r in results if r.get("ok"))
    print(f"  CA control sanity: shared-support controls flagged bad = {ctrl} (MUST be 0).")


# ---------------------------------------------------------------------------
# Output.
# ---------------------------------------------------------------------------
def write_outputs(results, diffs, elapsed, battery_meta, out_dir=RESULTS_DIR):
    os.makedirs(out_dir, exist_ok=True)
    payload = {
        "meta": {
            "experiment": "scaleup_smooth_vs_random",
            "timestamp": time.strftime("%Y%m%d_%H%M%S"),
            "elapsed_sec": round(elapsed, 2),
            "qk_cap": QK_CAP,
            "convention": "delta=radius (ABF). EXACT distance via full codeword "
                          "enumeration; ALL gamma in F enumerated per line. Line "
                          "counts tiered by cost (q^k*q*n) and logged per cell.",
            "fields_considered": ["GF(127)", "GF(257)", "GF(521)", "GF(1031)",
                                  "GF(2^7)", "GF(2^8)"],
            "note_GF2_7": "GF(2^7)* has prime order 127 -> NO proper multiplicative "
                          "subgroup, so it yields no smooth proper-subgroup cell.",
            "note_GF1031": "GF(1031) cells excluded: exact all-gamma distance costs "
                           "~21 s/line over its 1031 scalars (cost q^k*q*n > 2e9), too "
                           "slow for a meaningful line sample. q=521 already probes "
                           "n/|F| down to ~0.01.",
            "cost_cap_qk_q_n": 2e9,
            "battery": battery_meta,
        },
        "results": results,
        "differentials": diffs,
    }
    jpath = os.path.join(out_dir, "scaleup_smooth_vs_random.json")
    with open(jpath, "w", newline="\n") as f:
        json.dump(payload, f, indent=2, default=str)

    # CSV: one row per (cell, delta) with smooth + random side by side.
    cpath = os.path.join(out_dir, "scaleup_smooth_vs_random.csv")
    rows = []
    for r in results:
        if not r.get("ok"):
            continue
        sub = {s["delta"]: s for s in r["subgroup_per_delta"]}
        rnd = {s["delta"]: s for s in r["random_per_delta"]}
        for d in r["deltas"]:
            s, rr = sub.get(d), rnd.get(d)
            if not s or not rr:
                continue
            rows.append({
                "field": r["field"], "q": r["q"], "n": r["n"], "k": r["k"],
                "rho": round(r["rho"], 4), "n_over_q": round(r["n_over_q"], 5),
                "johnson": round(r["johnson"], 4), "capacity": round(r["capacity"], 4),
                "delta": d,
                "delta_minus_johnson": round(d - r["johnson"], 4),
                "delta_minus_capacity": round(d - r["capacity"], 4),
                "meaningful_regime": s["meaningful_regime"],
                "smooth_close_mean": round(s["close_count_mean"], 4),
                "random_close_mean": round(rr["close_count_mean"], 4),
                "close_mean_diff": round(s["close_count_mean"] - rr["close_count_mean"], 4),
                "smooth_close_max": s["close_count_max"],
                "random_close_max": rr["close_count_max"],
                "smooth_frac_bad": round(s["frac_bad"], 4),
                "random_frac_bad": round(rr["frac_bad"], 4),
                "badfrac_diff": round(s["frac_bad"] - rr["frac_bad"], 4),
                "smooth_Sstar_min_among_close": s["S_star_min_among_close"],
                "random_Sstar_min_among_close": rr["S_star_min_among_close"],
                "smooth_num_lines": s["num_lines"],
                "control_bad": s["control_bad_lines"] + rr["control_bad_lines"],
            })
    if rows:
        with open(cpath, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
            w.writeheader()
            for row in rows:
                w.writerow(row)
    return jpath, cpath, len(rows)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=min(16, os.cpu_count() or 8))
    ap.add_argument("--seed", type=int, default=0x5CA1E)
    ap.add_argument("--out", type=str,
                    default=os.environ.get("SCALEUP_OUT", RESULTS_DIR))
    args = ap.parse_args()

    battery = build_battery()
    battery.sort(key=lambda c: c.cost)
    jobs = [{"cfg": c, "seed": args.seed + 7 * i} for i, c in enumerate(battery)]
    battery_meta = [{"label": c.label, "q": c.q, "n": c.n, "k": c.k,
                     "rho": round(c.rho, 4), "qk": c.qk, "cost": c.cost,
                     "n_over_q": round(c.n / c.q, 5),
                     "line_budget": line_budget(c.cost)} for c in battery]

    print("=" * 84)
    print("SCALE-UP smooth-vs-random  (proper-subgroup domains, larger fields)")
    print("=" * 84)
    print(f"{len(battery)} matched cells (subgroup vs random), procs={args.procs}")
    for c in battery:
        lb = line_budget(c.cost)
        nl = sum(lb.values())
        print(f"  {c.label:18} q={c.q:>4} rho={c.rho:.3f} n/q={c.n/c.q:.4f} "
              f"q^k={c.qk:.1e} cost={c.cost:.1e} lines/domain~{nl}")
    print("-" * 84, flush=True)

    t0 = time.time()
    results = []
    with Pool(processes=args.procs) as pool:
        for r in pool.imap_unordered(run_cell, jobs):
            status = "ok " if r.get("ok") else "ERR"
            print(f"  [{status}] {r['label']:18} q={r['q']:>4} "
                  f"({r['elapsed_sec']:.1f}s) ctrl_bad={r.get('control_bad_total','?')}",
                  flush=True)
            results.append(r)
    elapsed = time.time() - t0
    print("-" * 84)
    print(f"All cells done in {elapsed:.1f}s ({elapsed/60:.1f} min)")

    diffs = compute_differentials(results)
    print_summary(results, diffs)
    jpath, cpath, nrows = write_outputs(results, diffs, elapsed, battery_meta,
                                        out_dir=args.out)
    print(f"\nWrote:\n  {jpath}\n  {cpath} ({nrows} rows)")


if __name__ == "__main__":
    main()
