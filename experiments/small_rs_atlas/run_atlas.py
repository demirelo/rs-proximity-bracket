"""
run_atlas.py -- orchestrate the small-field Reed-Solomon proximity-gap atlas.

For a battery of (field, n, k) codes and several evaluation-domain TYPES
  {smooth subgroup, multiplicative coset, random subset, full F*}
at rates rho ~ {1/2, 1/4, ...} where feasible, and several proximity radii delta
straddling the Johnson radius (1-sqrt(rho)) and capacity (1-rho), we run:

  * the bad-line search  (search_bad_lines.search_bad_lines)
  * the interleaved list-size search for m in {2,3} (search_interleaved)

and write structured JSON + flat CSV to results/.

The CENTRAL comparison the prize cares about: at a FIXED (rho, delta, |F|), does
the SMOOTH (subgroup/coset) domain show MORE bad lines / LARGER interleaved lists
than RANDOM-subset / FULL domains?  We therefore evaluate all four domain types on
the SAME field and the SAME n,k so the only thing changing is domain structure.

Performance budget
------------------
Exact distance is O(q^k * n) per word, so we TIER the per-cell sample counts by
q^k (documented in BUDGET below) to keep total wall-clock under ~15 min on 16
cores.  Parallelism: one process per (config, domain-type) cell via a
multiprocessing Pool; each worker builds its own codeword book (uint8, cheap to
fork).  All sampling counts are logged explicitly -- NO silent caps.
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
import argparse
import traceback
from dataclasses import dataclass, asdict
from multiprocessing import Pool

import numpy as np

from ff import PrimeField, BinaryExtensionField, FiniteField, _divisors
from rs import (build_codeword_book, domain_subgroup, domain_coset,
                domain_full, domain_random)
from search_bad_lines import search_bad_lines, search_bad_lines_multi
from search_interleaved import search_interleaved


RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


# ---------------------------------------------------------------------------
# The battery of codes.  Each entry: (field_spec, n, k, label).
# Chosen so q^k stays <= ~5e6 (exact distance feasible) and we get rates near
# 1/2 and 1/4, over BOTH prime and binary-extension fields, with n having
# multiple divisors (so coset/subgroup domains of the right size exist).
# ---------------------------------------------------------------------------
@dataclass
class CodeConfig:
    field_spec: str
    field: FiniteField
    n: int
    k: int
    label: str

    @property
    def qk(self) -> int:
        return self.field.q ** self.k

    @property
    def rho(self) -> float:
        return self.k / self.n


def _mk(field: FiniteField, n: int, k: int) -> CodeConfig:
    return CodeConfig(field.name, field, n, k,
                      f"{field.name}_n{n}_k{k}")


def build_battery() -> list[CodeConfig]:
    """Construct the list of codes to test.

    For each field we pick n | (q-1) (so a smooth subgroup of order n exists) with
    enough structure, and k giving rho ~ 1/2 and ~1/4.
    """
    F31 = PrimeField(31)      # q-1 = 30 = 2*3*5 ; divisors incl 5,6,10,15,30
    F16 = BinaryExtensionField(4)  # q-1 = 15 = 3*5 ; divisors 3,5,15
    F32 = BinaryExtensionField(5)  # q-1 = 31 (prime) ; only n=31 subgroup
    F64 = BinaryExtensionField(6)  # q-1 = 63 = 7*9 ; divisors 7,9,21,63
    F13 = PrimeField(13)      # q-1 = 12 = 2^2*3 ; divisors 4,6,12
    F61 = PrimeField(61)      # q-1 = 60 ; lots of divisors (rich smooth structure)

    cfgs: list[CodeConfig] = []

    # --- GF(2^4), n=15 (full subgroup) ---
    cfgs.append(_mk(F16, 15, 7))   # rho ~ 0.467 (~1/2)  q^k=2.7e8 -> interleaved only
    cfgs.append(_mk(F16, 15, 4))   # rho ~ 0.267 (~1/4)  q^k=6.5e4
    cfgs.append(_mk(F16, 15, 3))   # rho = 0.20         q^k=4.1e3 (cheap)

    # --- GF(31), n=15 and n=10 ---
    cfgs.append(_mk(F31, 15, 4))   # rho ~ 0.267        q^k=9.2e5
    cfgs.append(_mk(F31, 10, 5))   # rho = 0.5          q^k=2.86e7 -> interleaved only
    cfgs.append(_mk(F31, 10, 3))   # rho = 0.3          q^k=2.97e4
    cfgs.append(_mk(F31, 6,  3))   # rho = 0.5          q^k=2.97e4

    # --- GF(13), n=12 ---
    cfgs.append(_mk(F13, 12, 6))   # rho = 0.5          q^k=4.83e6
    cfgs.append(_mk(F13, 12, 3))   # rho = 0.25         q^k=2197 (cheap)

    # --- GF(2^5), n=31 (full subgroup, prime order) ---
    cfgs.append(_mk(F32, 31, 4))   # rho ~ 0.129        q^k=1.05e6

    # --- GF(2^6), n=21 and n=9 ---
    cfgs.append(_mk(F64, 21, 3))   # rho ~ 0.143        q^k=2.6e5
    cfgs.append(_mk(F64, 9,  4))   # rho ~ 0.444        q^k=1.68e7 -> interleaved only
    cfgs.append(_mk(F64, 9,  3))   # rho = 0.333        q^k=2.6e5

    # --- GF(61), n=20,15,12,10 (rich smooth structure, prime field) ---
    cfgs.append(_mk(F61, 20, 5))   # rho = 0.25         q^k=8.4e8 -> interleaved only
    cfgs.append(_mk(F61, 15, 4))   # rho ~ 0.267        q^k=1.38e7 -> interleaved only
    cfgs.append(_mk(F61, 12, 3))   # rho = 0.25         q^k=2.27e5
    cfgs.append(_mk(F61, 10, 5))   # rho = 0.5          q^k=8.4e8 -> interleaved only
    cfgs.append(_mk(F61, 10, 3))   # rho = 0.3          q^k=2.27e5

    # --- Low-k (k=2) codes for m=2 AND m=3 interleaved coverage.  q^(3k)=q^6 is
    #     small, so m=3 list-size searches are feasible here (the higher-k codes
    #     above have q^(3k) too large for m=3).  These are also PROPER subgroups
    #     (n < q-1) where smooth and random domains genuinely differ. ---
    cfgs.append(_mk(F31, 15, 2))   # rho=0.133  q^k=961    proper subgroup, m=3 ok
    cfgs.append(_mk(F31, 10, 2))   # rho=0.2    q^k=961    proper subgroup, m=3 ok
    cfgs.append(_mk(F31, 6,  2))   # rho=0.333  q^k=961    proper subgroup, m=3 ok
    cfgs.append(_mk(F61, 12, 2))   # rho=0.167  q^k=3721   proper subgroup, m=3 ok
    cfgs.append(_mk(F61, 10, 2))   # rho=0.2    q^k=3721   proper subgroup, m=3 ok
    cfgs.append(_mk(F64, 21, 2))   # rho=0.095  q^k=4096   proper subgroup, m=3 ok
    cfgs.append(_mk(F64, 9,  2))   # rho=0.222  q^k=4096   proper subgroup, m=3 ok
    cfgs.append(_mk(F16, 15, 2))   # rho=0.133  q^k=256    n=q-1, m=3 ok (q^6=1.7e7)

    return cfgs


# ---------------------------------------------------------------------------
# Per-cell sampling budget, tiered by q^k (exact-distance cost ~ q^k * n).
# ---------------------------------------------------------------------------
def line_budget(qk: int) -> dict:
    """Number of lines to sample per generator for the bad-line search.

    Tiered by q^k because exact distance costs O(q^k * n) PER line (each line
    triggers q distance computations).  Costs measured on this machine:
      q^k ~ 4e3  -> ~0.1 ms/line       q^k ~ 3e4 -> ~22 ms/line
      q^k ~ 2e5  -> ~0.2 s/line        q^k ~ 1e6 -> ~1.0 s/line
      q^k ~ 5e6  -> ~3 s/line
    Returns a dict (random/cwnoise/lowdeg/rational), or None if q^k is too big
    for an exact line search (then we run interleaved only and log the skip).
    """
    if qk <= 1e4:
        return dict(random=500, cwnoise=500, lowdeg=250, rational=250)
    if qk <= 5e4:
        return dict(random=250, cwnoise=250, lowdeg=120, rational=120)
    if qk <= 3e5:
        return dict(random=80, cwnoise=80, lowdeg=40, rational=40)
    if qk <= 1.2e6:
        return dict(random=20, cwnoise=20, lowdeg=10, rational=10)
    if qk <= 5e6:
        return dict(random=8, cwnoise=8, lowdeg=4, rational=4)
    return None    # too big for exact line search


def interleaved_budget(qk_m: int) -> dict | None:
    """Sampling for interleaved search, tiered by q^(k*m) (list enumeration cost).

    qk_m = q^(k*m) = number of interleaved codewords (the list lives in this set).
    Branch-and-bound makes this fast when lists are small, but we bound targets.
    """
    if qk_m <= 1e6:
        return dict(n_random=120, n_cwnoise=120, n_shared=120)
    if qk_m <= 1e8:
        return dict(n_random=60, n_cwnoise=60, n_shared=60)
    if qk_m <= 1e10:
        return dict(n_random=30, n_cwnoise=30, n_shared=30)
    if qk_m <= 5e10:
        # Branch-and-bound stays fast even here (lists are small except near
        # capacity); a handful of targets suffices to find the structured worst
        # case (shared-support boundary) plus the random/cw+noise baseline.
        return dict(n_random=12, n_cwnoise=12, n_shared=12)
    return None


# ---------------------------------------------------------------------------
# Delta grid: anchored at Johnson and capacity.
# ---------------------------------------------------------------------------
def delta_grid(rho: float, n: int) -> list[float]:
    """A grid of delta values straddling Johnson (1-sqrt(rho)) and capacity (1-rho).

    We snap to multiples of 1/n where helpful (distances are integers /n) but keep
    the raw anchor values too.  Points: below Johnson, at Johnson, between Johnson
    and capacity, at capacity, and slightly beyond capacity.
    """
    johnson = 1.0 - np.sqrt(rho)
    cap = 1.0 - rho
    pts = set()
    # fractions of johnson
    for f in (0.6, 0.8, 0.9, 1.0):
        pts.add(round(f * johnson, 4))
    # between johnson and capacity
    for t in (0.25, 0.5, 0.75, 1.0):
        pts.add(round(johnson + t * (cap - johnson), 4))
    # slightly beyond capacity
    pts.add(round(min(0.95, cap + 0.5 * (1 - cap) * 0.2), 4))
    # ensure within (0,1)
    grid = sorted(p for p in pts if 0.02 < p < 0.98)
    return grid


# ---------------------------------------------------------------------------
# Domain builders for a given config + domain type.
# ---------------------------------------------------------------------------
def make_domain(cfg: CodeConfig, dtype: str, rng: np.random.Generator):
    """Return (L, ok, note).  ok=False if this domain type is infeasible for cfg."""
    F, n = cfg.field, cfg.n
    if dtype == "subgroup":
        if (F.q - 1) % n != 0:
            return None, False, f"n={n} does not divide q-1={F.q-1}"
        return domain_subgroup(F, n), True, "order-n multiplicative subgroup"
    if dtype == "coset":
        if (F.q - 1) % n != 0:
            return None, False, f"n={n} does not divide q-1={F.q-1}"
        # choose a shift outside the subgroup
        H = set(domain_subgroup(F, n).tolist())
        shift = next((s for s in range(2, F.q) if s not in H), 1)
        return domain_coset(F, n, shift), True, f"coset shift={shift}"
    if dtype == "random":
        if n > F.q - 1:
            return None, False, f"n={n} > |F*|={F.q-1}"
        return domain_random(F, n, rng), True, "random subset of F*"
    if dtype == "full":
        if n != F.q - 1:
            return None, False, f"full domain has size q-1={F.q-1} != n={n}"
        return domain_full(F), True, "full F*"
    raise ValueError(dtype)


# ---------------------------------------------------------------------------
# One worker: evaluate one (config, domain-type) cell.
# ---------------------------------------------------------------------------
def run_cell(args) -> dict:
    """Worker entry point.  args = (cfg, dtype, seed, do_lines, ms_for_inter)."""
    cfg: CodeConfig = args["cfg"]
    dtype: str = args["dtype"]
    seed: int = args["seed"]
    F = cfg.field
    rng = np.random.default_rng(seed)

    out = {
        "label": cfg.label, "field": F.name, "q": F.q, "n": cfg.n, "k": cfg.k,
        "rho": cfg.rho, "qk": cfg.qk, "domain_type": dtype,
        "johnson": 1.0 - np.sqrt(cfg.rho), "capacity": 1.0 - cfg.rho,
        "ok": True, "note": "", "bad_lines": None, "interleaved": None,
        "elapsed_sec": 0.0,
    }
    t0 = time.time()
    try:
        L, ok, note = make_domain(cfg, dtype, rng)
        out["note"] = note
        if not ok:
            out["ok"] = False
            out["elapsed_sec"] = time.time() - t0
            return out

        deltas = delta_grid(cfg.rho, cfg.n)
        out["deltas"] = deltas

        # ---------- bad-line search (if q^k small enough) ----------
        lb = line_budget(cfg.qk)
        if lb is not None:
            book = build_codeword_book(F, L, cfg.k)
            # Single sampling pass classified across the whole delta grid.
            line_results = search_bad_lines_multi(
                book, deltas, rng,
                n_random=lb["random"], n_cwnoise=lb["cwnoise"],
                n_lowdeg=lb["lowdeg"], n_rational=lb["rational"],
                verbose=False)
            out["bad_lines"] = {
                "budget": lb,
                "per_delta": line_results,
                "any_bad": any(s["num_bad_lines"] > 0 for s in line_results),
                "max_bad_over_deltas": max(s["num_bad_lines"] for s in line_results),
            }
        else:
            out["bad_lines"] = {"skipped": True,
                                "reason": f"q^k={cfg.qk:.2e} too big for exact line search"}

        # ---------- interleaved search for m=2,3 ----------
        inter = {}
        for m in (2, 3):
            qk_m = F.q ** (cfg.k * m)
            ib = interleaved_budget(qk_m)
            if ib is None:
                inter[f"m{m}"] = {"skipped": True,
                                  "reason": f"q^(k*m)={qk_m:.2e} too big"}
                continue
            # interleaved needs a codeword book too
            book_i = build_codeword_book(F, L, cfg.k)
            summ = search_interleaved(book_i, m, deltas, rng,
                                      n_random=ib["n_random"],
                                      n_cwnoise=ib["n_cwnoise"],
                                      n_shared=ib["n_shared"], verbose=False)
            inter[f"m{m}"] = summ
        out["interleaved"] = inter

    except Exception as e:  # never let one cell kill the whole run
        out["ok"] = False
        out["note"] = f"EXCEPTION: {e}"
        out["traceback"] = traceback.format_exc()
    out["elapsed_sec"] = round(time.time() - t0, 2)
    return out


# ---------------------------------------------------------------------------
# Orchestration.
# ---------------------------------------------------------------------------
def build_jobs(battery, domain_types, base_seed) -> list[dict]:
    jobs = []
    sid = base_seed
    for cfg in battery:
        for dtype in domain_types:
            jobs.append({"cfg": cfg, "dtype": dtype, "seed": sid})
            sid += 1
    return jobs


def flatten_for_csv(results: list[dict]) -> list[dict]:
    """One CSV row per (cell, delta): the headline numbers for analysis."""
    rows = []
    for r in results:
        if not r.get("ok"):
            continue
        base = {k: r[k] for k in
                ("label", "field", "q", "n", "k", "rho", "qk", "domain_type",
                 "johnson", "capacity")}
        deltas = r.get("deltas", [])
        bl = r.get("bad_lines") or {}
        bl_per = {s["delta"]: s for s in bl.get("per_delta", [])} if "per_delta" in bl else {}
        inter = r.get("interleaved") or {}
        for d in deltas:
            row = dict(base)
            row["delta"] = d
            row["delta_minus_johnson"] = round(d - base["johnson"], 4)
            row["delta_minus_capacity"] = round(d - base["capacity"], 4)
            s = bl_per.get(d)
            if s:
                row["line_meaningful"] = s["meaningful_regime"]
                row["line_num_bad"] = s["num_bad_lines"]
                row["line_frac_bad"] = round(s["frac_bad"], 4)
                row["line_control_bad"] = s["control_bad_lines"]
                row["line_close_max"] = s["close_count_max"]
                row["line_close_mean"] = round(s["close_count_mean"], 3)
                row["line_close_p90"] = round(s["close_count_p90"], 3)
                row["line_frac_2plus_close"] = round(s["frac_lines_with_2plus_close"], 4)
                row["line_Sstar_min_among_close"] = s["S_star_min_among_close"]
                row["line_Sstar_mean"] = round(s["S_star_mean"], 3)
                row["line_ca_threshold"] = s["ca_threshold"]
                row["line_num_lines"] = s["num_lines"]
            else:
                row["line_num_bad"] = ""
            for m in (2, 3):
                key = f"m{m}"
                im = inter.get(key)
                if im and not im.get("skipped"):
                    mb = im["max_list_by_delta"].get(str(d))
                    row[f"inter_m{m}_maxlist"] = mb["size"] if mb else ""
                    row[f"inter_m{m}_kind"] = mb["kind"] if mb else ""
                    meanb = im["mean_list_by_delta"].get(str(d))
                    row[f"inter_m{m}_meanlist"] = round(meanb, 3) if meanb is not None else ""
                else:
                    row[f"inter_m{m}_maxlist"] = ""
            rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description="RS proximity-gap small-field atlas")
    ap.add_argument("--procs", type=int, default=min(16, os.cpu_count() or 8))
    ap.add_argument("--seed", type=int, default=0xA71A5)
    ap.add_argument("--quick", action="store_true",
                    help="reduced battery for a fast smoke run")
    ap.add_argument("--out", type=str, default=RESULTS_DIR)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    battery = build_battery()
    if args.quick:
        battery = [c for c in battery if c.qk <= 1e5][:5]

    domain_types = ["subgroup", "coset", "random", "full"]
    jobs = build_jobs(battery, domain_types, args.seed)

    # Count feasible jobs for the progress log.
    print("=" * 70)
    print("Reed-Solomon proximity-gap SMALL-FIELD ATLAS")
    print("=" * 70)
    print(f"Battery: {len(battery)} codes x {len(domain_types)} domain types "
          f"= {len(jobs)} cells")
    print(f"Procs: {args.procs}   base seed: {hex(args.seed)}")
    print("Codes (q^k = exact-distance cost driver):")
    for c in battery:
        lb = line_budget(c.qk)
        print(f"  {c.label:22} rho={c.rho:.3f} q^k={c.qk:.2e} "
              f"johnson={1-np.sqrt(c.rho):.3f} cap={1-c.rho:.3f} "
              f"lines={'yes' if lb else 'NO(too big)'}")
    print("-" * 70, flush=True)

    t_start = time.time()
    results = []
    done = 0
    with Pool(processes=args.procs) as pool:
        for r in pool.imap_unordered(run_cell, jobs):
            done += 1
            status = "ok " if r.get("ok") else "SKIP/ERR"
            extra = ""
            if r.get("ok"):
                bl = r.get("bad_lines") or {}
                nb = bl.get("max_bad_over_deltas")
                if bl.get("skipped"):
                    extra = "lines=skip"
                elif nb is not None:
                    extra = f"maxbad={nb}"
            print(f"[{done:3}/{len(jobs)}] {status} {r['label']:22} "
                  f"{r['domain_type']:9} {extra:14} "
                  f"({r['elapsed_sec']:.1f}s)  {r.get('note','')[:40]}", flush=True)
            results.append(r)

    elapsed = time.time() - t_start
    print("-" * 70)
    print(f"All cells done in {elapsed:.1f}s ({elapsed/60:.1f} min)")

    # ---- write outputs ----
    ts = time.strftime("%Y%m%d_%H%M%S")
    json_path = os.path.join(args.out, "atlas_results.json")
    csv_path = os.path.join(args.out, "atlas_results.csv")

    # Strip non-serializable field objects from configs in the dump.
    def clean(r):
        return r  # run_cell already returns only plain types
    payload = {
        "meta": {
            "timestamp": ts, "elapsed_sec": round(elapsed, 1),
            "procs": args.procs, "seed": args.seed, "quick": args.quick,
            "num_codes": len(battery), "num_cells": len(jobs),
            "domain_types": domain_types,
            "numpy": np.__version__,
            "note": ("All distances/list-sizes are EXACT (full codeword "
                     "enumeration). Sampling is over WORDS/LINES/TARGETS only; "
                     "counts are recorded per cell. gamma is enumerated over all "
                     "of F for every line (no gamma sampling)."),
        },
        "results": [clean(r) for r in results],
    }
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2, default=str)
    print(f"Wrote {json_path}")

    rows = flatten_for_csv(results)
    if rows:
        # union of all keys for a stable header
        keys = []
        for row in rows:
            for kk in row:
                if kk not in keys:
                    keys.append(kk)
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for row in rows:
                w.writerow(row)
        print(f"Wrote {csv_path} ({len(rows)} rows)")

    # ---- quick top-level analysis printed to console ----
    _print_headline_analysis(results)


def _collect_bl(results):
    """Index bad-line per-delta summaries by (field,n,k,delta,domain)."""
    idx = {}
    for r in results:
        if not r.get("ok"):
            continue
        bl = r.get("bad_lines") or {}
        for s in bl.get("per_delta", []):
            key = (r["field"], r["n"], r["k"], s["delta"])
            idx.setdefault(key, {})[r["domain_type"]] = s
    return idx


def _print_headline_analysis(results):
    """Print the smooth-vs-random/full DIFFERENTIAL comparison and onset points.

    The research question is NOT "are there bad lines" (near capacity every domain
    has them -- list-decoding regime) but "do SMOOTH domains show MORE/larger
    badness than RANDOM/FULL at the SAME (rho,delta,|F|)".  So every comparison
    here is matched on (field,n,k,delta) and contrasts domain types.
    """
    print("=" * 78)
    print("HEADLINE ANALYSIS  (smooth=subgroup/coset  vs  random/full)")
    print("=" * 78)

    bl_idx = _collect_bl(results)

    # ---- 1. Differential close-count: smooth minus random, matched cells ----
    # For each (field,n,k,delta) with BOTH subgroup and random present and in the
    # MEANINGFUL regime, compute mean close-count difference and frac_bad diff.
    print("\n[1] Matched smooth-vs-random close-count & bad-fraction differential")
    print("    (meaningful regime only: ca_threshold>k; positive => smooth worse)")
    diffs_close = []   # (subgroup - random) mean close-count
    diffs_badfrac = []
    diffs_close_coset = []
    n_matched = 0
    for key, dmap in sorted(bl_idx.items()):
        if "subgroup" not in dmap or "random" not in dmap:
            continue
        sg, rd = dmap["subgroup"], dmap["random"]
        if not sg.get("meaningful_regime"):
            continue
        n_matched += 1
        diffs_close.append(sg["close_count_mean"] - rd["close_count_mean"])
        diffs_badfrac.append(sg["frac_bad"] - rd["frac_bad"])
        if "coset" in dmap:
            diffs_close_coset.append(dmap["coset"]["close_count_mean"]
                                     - rd["close_count_mean"])
    if n_matched:
        dc = np.array(diffs_close); db = np.array(diffs_badfrac)
        print(f"    matched (field,n,k,delta) cells: {n_matched}")
        print(f"    close-count mean diff (subgroup - random): "
              f"mean={dc.mean():+.3f}  max={dc.max():+.3f}  min={dc.min():+.3f}")
        print(f"    bad-fraction  diff (subgroup - random): "
              f"mean={db.mean():+.4f}  max={db.max():+.4f}  min={db.min():+.4f}")
        if diffs_close_coset:
            dcc = np.array(diffs_close_coset)
            print(f"    close-count mean diff (coset    - random): "
                  f"mean={dcc.mean():+.3f}  max={dcc.max():+.3f}")
        # Verdict
        if abs(dc.mean()) < 0.5 and abs(db.mean()) < 0.02:
            print("    => SMOOTH ~ RANDOM (no smooth-specific line badness): "
                  "supports the POSITIVE route.")
        elif dc.mean() > 0.5 or db.mean() > 0.02:
            print("    => SMOOTH shows MORE line badness than random: "
                  "possible counterexample seed.")
        else:
            print("    => SMOOTH shows LESS badness than random.")
    else:
        print("    (no matched meaningful-regime cells)")

    # ---- 2. Control sanity across all cells ----
    ctrl_bad_total = 0
    cells = 0
    for r in results:
        if not r.get("ok"):
            continue
        for s in (r.get("bad_lines") or {}).get("per_delta", []):
            cells += 1
            ctrl_bad_total += s["control_bad_lines"]
    print(f"\n[2] CA control sanity: shared-support controls flagged bad = "
          f"{ctrl_bad_total} (MUST be 0) across {cells} (cell,delta) entries.")

    # ---- 3. Badness onset vs Johnson/capacity, by domain type ----
    print("\n[3] BAD-LINE onset: smallest delta (meaningful regime) with frac_bad>0,")
    print("    reported relative to Johnson (1-sqrt rho) and capacity (1-rho).")
    onset_by_dt = {dt: [] for dt in ("subgroup", "coset", "random", "full")}
    for r in results:
        if not r.get("ok"):
            continue
        dt = r["domain_type"]
        johnson, cap = r["johnson"], r["capacity"]
        per = sorted((r.get("bad_lines") or {}).get("per_delta", []),
                     key=lambda s: s["delta"])
        onset = None
        for s in per:
            if s["meaningful_regime"] and s["num_bad_lines"] > 0:
                onset = s["delta"]
                break
        if onset is not None:
            onset_by_dt[dt].append((onset, onset - johnson, onset - cap))
    for dt in ("subgroup", "coset", "random", "full"):
        v = onset_by_dt[dt]
        if v:
            rj = np.mean([x[1] for x in v]); rc = np.mean([x[2] for x in v])
            print(f"    {dt:9}: {len(v):2} cells, mean onset (delta-Johnson)="
                  f"{rj:+.3f}, (delta-capacity)={rc:+.3f}")
        else:
            print(f"    {dt:9}: no bad lines in meaningful regime (good)")

    # ---- 4. Interleaved list onset, matched smooth-vs-random ----
    print("\n[4] Interleaved (m=2) list-size onset: smallest delta with max-list>1.")
    onset_il = {dt: [] for dt in ("subgroup", "coset", "random", "full")}
    for r in results:
        if not r.get("ok"):
            continue
        inter = (r.get("interleaved") or {}).get("m2")
        if not inter or inter.get("skipped"):
            continue
        dt = r["domain_type"]; johnson = r["johnson"]; cap = r["capacity"]
        onset = None
        for d in sorted(r.get("deltas", [])):
            mb = inter["max_list_by_delta"].get(str(d))
            if mb and mb["size"] > 1:
                onset = d
                break
        if onset is not None:
            onset_il[dt].append((onset, onset - johnson, onset - cap))
    for dt in ("subgroup", "coset", "random", "full"):
        v = onset_il[dt]
        if v:
            rj = np.mean([x[1] for x in v]); rc = np.mean([x[2] for x in v])
            print(f"    {dt:9}: {len(v):2} cells, mean onset (delta-Johnson)="
                  f"{rj:+.3f}, (delta-capacity)={rc:+.3f}")
        else:
            print(f"    {dt:9}: no list>1 observed in tested delta range")

    # ---- 5. Max interleaved list head-to-head, same (field,n,k) ----
    print("\n[5] Max m=2 interleaved list over all tested delta, by domain "
          "(same field,n,k):")
    by_code = {}
    for r in results:
        if not r.get("ok"):
            continue
        inter = (r.get("interleaved") or {}).get("m2")
        if not inter or inter.get("skipped"):
            continue
        key = (r["field"], r["n"], r["k"])
        mx = max((inter["max_list_by_delta"][str(d)]["size"]
                  for d in r.get("deltas", [])
                  if str(d) in inter["max_list_by_delta"]), default=0)
        by_code.setdefault(key, {})[r["domain_type"]] = mx
    print(f"    {'code':22} {'subgroup':>9} {'coset':>9} {'random':>9} {'full':>9}")
    smooth_worse = 0; total_cmp = 0
    for key, d in sorted(by_code.items()):
        lbl = f"{key[0]}_n{key[1]}_k{key[2]}"
        print(f"    {lbl:22} {str(d.get('subgroup','-')):>9} "
              f"{str(d.get('coset','-')):>9} "
              f"{str(d.get('random','-')):>9} {str(d.get('full','-')):>9}")
        if "subgroup" in d and "random" in d:
            total_cmp += 1
            if d["subgroup"] > d["random"]:
                smooth_worse += 1
    if total_cmp:
        print(f"    -> smooth subgroup had a STRICTLY larger max list than random "
              f"in {smooth_worse}/{total_cmp} matched codes.")
    print("=" * 78)


if __name__ == "__main__":
    main()
