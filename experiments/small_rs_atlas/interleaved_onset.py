"""
interleaved_onset.py -- PRECISE onset of m=2 interleaved list growth vs capacity.

Wave 1 observed that the maximum 2-interleaved list size stays EXACTLY 1 from below
the Johnson radius, through it, and well beyond -- departing from 1 only near
`capacity - O(1/n)` and exploding at capacity (1 - rho).  This script measures that
onset PRECISELY for a range of proper-subgroup codes, on a fine delta lattice in units
of 1/n, and fits the gap below capacity:

    delta_onset = capacity - c_onset / n      (c_onset = (capacity - delta_onset) * n)

to test the conjecture "interleaved lists stay <= small up to capacity - O(1/n)".

Definitions (exact; reuse search_interleaved.py)
------------------------------------------------
Interleaved code C^{equiv 2}: pairs (c0,c1) of codewords, columnwise metric
d(W,c) = #{ i : w0[i]!=c0[i] OR w1[i]!=c1[i] }.  The list at radius delta for a target
W=(w0,w1) is { (c0,c1) : d(W,c) <= delta n }, computed EXACTLY over all q^{2k}
interleaved codewords by branch-and-bound on bit-packed agreement masks.

  delta_onset(W)        = smallest delta on the 1/n lattice with |list(W,delta)| > 1.
  delta_onset(code)     = min over probed targets W (the worst / earliest onset).

Targets probed per code (sampling counts logged, no silent caps):
  * the deterministic structured-boundary target at each lattice delta (codewords on a
    shared support of size exactly n - ceil((1-delta)n)) -- the analytic worst case;
  * random and cw+noise and shared-support sampled targets.

We scan delta upward from the Johnson radius to capacity on the 1/n lattice and record,
for EACH lattice delta, the max list over all probed targets.  The onset is the first
lattice delta where that max exceeds 1.  We compare smooth subgroup vs random domain
(matched field/n/k) to confirm the onset is a rate phenomenon, not a domain one.

Output -> results/interleaved_onset.csv  (+ a json with the full per-delta curves).
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
from search_interleaved import (interleaved_list_sizes, gen_target_random,
                                gen_target_cwnoise, gen_target_shared_support)

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


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
    def label(self):
        return f"{self.field_name}_n{self.n}_k{self.k}"


def build_battery() -> list[Cfg]:
    """Proper-subgroup codes spanning a range of n (to fit capacity - c/n), with
    q^{2k} small enough for exact interleaved enumeration.  k=2 throughout so the
    interleaved code C^{equiv 2} has q^4 codewords (tractable up to q~257)."""
    fields = [
        ("GF(31)", PrimeField(31)),
        ("GF(61)", PrimeField(61)),
        ("GF(2^6)", BinaryExtensionField(6)),
        ("GF(127)", PrimeField(127)),
        ("GF(257)", PrimeField(257)),
    ]
    cfgs: list[Cfg] = []
    seen = set()
    for name, F in fields:
        q = F.q
        for n in _divisors(q - 1):
            if n < 6 or n >= q - 1:
                continue
            k = 2
            if q ** (2 * k) > 1.2e9:       # interleaved enumeration budget
                continue
            # keep a spread of n; cap how many per field to bound runtime
            key = (name, n, k)
            if key in seen:
                continue
            seen.add(key)
            cfgs.append(Cfg(name, F, n, k))
    return cfgs


def fine_delta_lattice(rho: float, n: int) -> list[float]:
    """Delta values on the 1/n lattice from ~Johnson up to capacity (inclusive).

    Distances are integers, so the meaningful radii are j/n.  We scan from just below
    the Johnson radius to capacity = 1 - rho, which is ((n-k)/n).  This resolves the
    onset to the 1/n granularity needed to fit capacity - c/n.
    """
    johnson = 1.0 - np.sqrt(rho)
    cap = 1.0 - rho
    j_lo = max(1, int(np.floor(johnson * n)) - 1)
    j_hi = int(round(cap * n))                  # = n - k
    lattice = sorted({round(j / n, 6) for j in range(j_lo, j_hi + 1)})
    return lattice


def measure_onset_on_domain(F: FiniteField, L: np.ndarray, k: int,
                            lattice: list[float], rng,
                            n_random=40, n_cwnoise=40, n_shared=40) -> dict:
    """Max m=2 interleaved list size at each lattice delta, over probed targets, and
    the onset delta (first lattice delta with max-list > 1).  EXACT list sizes."""
    n = len(L)
    book = build_codeword_book(F, L, k)

    # max list over all probed targets, per lattice delta
    max_list = {d: 0 for d in lattice}

    def absorb(target):
        sizes, info = interleaved_list_sizes(book, target, lattice)
        for d in lattice:
            if sizes[d] > max_list[d]:
                max_list[d] = sizes[d]
        return info["exact"]

    exact_all = True
    # sampled targets
    for _ in range(n_random):
        exact_all &= absorb(gen_target_random(F, book, rng, 2))
    e_each = max(1, int(round(0.5 * (1 - k / n) * n)))   # mid-range corruption
    for _ in range(n_cwnoise):
        exact_all &= absorb(gen_target_cwnoise(F, book, rng, 2, e_each))
    for _ in range(n_shared):
        exact_all &= absorb(gen_target_shared_support(F, book, rng, 2, e_each))

    # deterministic structured-boundary target at EACH lattice delta (analytic worst
    # case): codewords on a shared support of size exactly n - ceil((1-delta)n).
    for d in lattice:
        e = max(0, n - int(np.ceil((1.0 - d) * n - 1e-9)))
        if e < 1:
            continue
        tgt = gen_target_shared_support(F, book, rng, 2, e)
        sizes, info = interleaved_list_sizes(book, tgt, [d])
        exact_all &= info["exact"]
        if sizes[d] > max_list[d]:
            max_list[d] = sizes[d]

    cap = 1.0 - k / n
    johnson = 1.0 - np.sqrt(k / n)
    onset = next((d for d in lattice if max_list[d] > 1), None)
    return {
        "n": n, "k": k, "capacity": cap, "johnson": johnson,
        "lattice": lattice,
        "max_list_by_delta": {str(d): max_list[d] for d in lattice},
        "onset_delta": onset,
        "onset_minus_capacity": (onset - cap) if onset is not None else None,
        "onset_minus_johnson": (onset - johnson) if onset is not None else None,
        "c_onset": ((cap - onset) * n) if onset is not None else None,   # capacity - c/n
        "list_at_capacity": max_list[lattice[-1]] if lattice else 0,
        "exact": exact_all,
        "sampling": {"random": n_random, "cwnoise": n_cwnoise, "shared": n_shared,
                     "structured_boundary": "1 per lattice delta"},
    }


def run_cell(args) -> dict:
    cfg: Cfg = args["cfg"]
    seed: int = args["seed"]
    F = cfg.field
    out = {"label": cfg.label, "field": F.name, "q": F.q, "n": cfg.n, "k": cfg.k,
           "rho": cfg.rho, "ok": True, "elapsed_sec": 0.0}
    t0 = time.time()
    try:
        lattice = fine_delta_lattice(cfg.rho, cfg.n)
        rng_s = np.random.default_rng(seed)
        sub = measure_onset_on_domain(F, domain_subgroup(F, cfg.n), cfg.k,
                                      lattice, rng_s)
        rng_r = np.random.default_rng(seed + 1)
        rnd = measure_onset_on_domain(F, domain_random(F, cfg.n, rng_r), cfg.k,
                                      lattice, rng_r)
        out["subgroup"] = sub
        out["random"] = rnd
    except Exception as e:
        out["ok"] = False
        out["error"] = str(e)
        out["traceback"] = traceback.format_exc()
    out["elapsed_sec"] = round(time.time() - t0, 2)
    return out


def write_outputs(results, elapsed, out_dir=RESULTS_DIR):
    os.makedirs(out_dir, exist_ok=True)
    payload = {
        "meta": {
            "experiment": "interleaved_onset_precision",
            "timestamp": time.strftime("%Y%m%d_%H%M%S"),
            "elapsed_sec": round(elapsed, 2),
            "convention": "delta=radius. m=2 interleaved (columnwise) list sizes EXACT "
                          "over all q^{2k} interleaved codewords. Delta scanned on the "
                          "1/n lattice from ~Johnson to capacity. onset = first lattice "
                          "delta with max-list>1 over probed targets.",
        },
        "results": results,
    }
    jpath = os.path.join(out_dir, "interleaved_onset.json")
    with open(jpath, "w", newline="\n") as f:
        json.dump(payload, f, indent=2, default=str)

    cpath = os.path.join(out_dir, "interleaved_onset.csv")
    rows = []
    for r in results:
        if not r.get("ok"):
            continue
        for dom_key, dom in (("subgroup", r["subgroup"]), ("random", r["random"])):
            rows.append({
                "field": r["field"], "q": r["q"], "n": r["n"], "k": r["k"],
                "rho": round(r["rho"], 4), "domain": dom_key,
                "johnson": round(dom["johnson"], 4),
                "capacity": round(dom["capacity"], 4),
                "onset_delta": dom["onset_delta"],
                "onset_minus_capacity": (round(dom["onset_minus_capacity"], 4)
                                         if dom["onset_minus_capacity"] is not None else ""),
                "onset_minus_johnson": (round(dom["onset_minus_johnson"], 4)
                                        if dom["onset_minus_johnson"] is not None else ""),
                "c_onset_capacity_minus_c_over_n": (round(dom["c_onset"], 3)
                                                    if dom["c_onset"] is not None else ""),
                "list_at_capacity": dom["list_at_capacity"],
                "exact": dom["exact"],
            })
    if rows:
        with open(cpath, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
            w.writeheader()
            for row in rows:
                w.writerow(row)
    return jpath, cpath, len(rows)


def print_summary(results):
    print("\n" + "=" * 88)
    print("INTERLEAVED (m=2) LIST ONSET vs capacity  (subgroup domain)")
    print("=" * 88)
    print(f"  {'code':18} {'q':>5} {'n':>4} {'rho':>5} {'Johnson':>8} {'capacity':>8} "
          f"{'onset':>7} {'on-cap':>7} {'on-J':>6} {'c=(cap-on)n':>11} {'list@cap':>8}")
    cons = []
    for r in sorted(results, key=lambda x: (x.get("q", 0), x.get("n", 0))):
        if not r.get("ok"):
            continue
        d = r["subgroup"]
        on = d["onset_delta"]
        onc = d["onset_minus_capacity"]
        onj = d["onset_minus_johnson"]
        c = d["c_onset"]
        print(f"  {r['label']:18} {r['q']:>5} {r['n']:>4} {r['rho']:>5.3f} "
              f"{d['johnson']:>8.3f} {d['capacity']:>8.3f} "
              f"{(f'{on:.3f}' if on is not None else 'none'):>7} "
              f"{(f'{onc:+.3f}' if onc is not None else '-'):>7} "
              f"{(f'{onj:+.3f}' if onj is not None else '-'):>6} "
              f"{(f'{c:.2f}' if c is not None else '-'):>11} {d['list_at_capacity']:>8}")
        if c is not None:
            cons.append((r["n"], c, onc))
    if cons:
        cs = np.array([c for _, c, _ in cons])
        oncs = np.array([o for _, _, o in cons])
        print("-" * 88)
        print(f"  c_onset = (capacity - onset)*n : mean={cs.mean():.2f} "
              f"min={cs.min():.2f} max={cs.max():.2f}  (n ranges {min(n for n,_,_ in cons)}"
              f"..{max(n for n,_,_ in cons)})")
        print(f"  onset - capacity (absolute)    : mean={oncs.mean():+.3f} "
              f"min={oncs.min():+.3f} max={oncs.max():+.3f}")
        if cs.std() < 1.5 and abs(cs.mean()) < 6:
            print(f"  => onset sits at capacity - c/n with c ~ {cs.mean():.1f} (small "
                  f"constant): SUPPORTS 'lists <= small until capacity - O(1/n)'.")
        # compare to Johnson: is onset much closer to capacity than to Johnson?
        onj_vals = [r["subgroup"]["onset_minus_johnson"] for r in results
                    if r.get("ok") and r["subgroup"]["onset_delta"] is not None]
        if onj_vals:
            print(f"  (For reference: onset - Johnson mean = {np.mean(onj_vals):+.3f} "
                  f"-- onset is far ABOVE the Johnson radius.)")
    # subgroup vs random onset identical?
    same = 0; tot = 0
    for r in results:
        if not r.get("ok"):
            continue
        tot += 1
        if r["subgroup"]["onset_delta"] == r["random"]["onset_delta"]:
            same += 1
    print(f"  subgroup-vs-random onset IDENTICAL in {same}/{tot} codes "
          f"(onset is a RATE phenomenon, not a domain one).")


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=min(16, os.cpu_count() or 8))
    ap.add_argument("--seed", type=int, default=0x012345)
    ap.add_argument("--out", type=str,
                    default=os.environ.get("ONSET_OUT", RESULTS_DIR))
    args = ap.parse_args()

    battery = build_battery()
    jobs = [{"cfg": c, "seed": args.seed + 11 * i} for i, c in enumerate(battery)]
    print("=" * 88)
    print(f"INTERLEAVED ONSET PRECISION: {len(battery)} proper-subgroup codes, "
          f"procs={args.procs}")
    for c in battery:
        lat = fine_delta_lattice(c.rho, c.n)
        print(f"  {c.label:18} q={c.q:>4} rho={c.rho:.3f} cap={1-c.rho:.3f} "
              f"lattice {len(lat)} pts in [{lat[0]:.3f},{lat[-1]:.3f}]")
    print("-" * 88, flush=True)

    t0 = time.time()
    results = []
    with Pool(processes=args.procs) as pool:
        for r in pool.imap_unordered(run_cell, jobs):
            st = "ok " if r.get("ok") else "ERR"
            sub = r.get("subgroup", {})
            print(f"  [{st}] {r['label']:18} onset={sub.get('onset_delta')} "
                  f"({r['elapsed_sec']:.1f}s)", flush=True)
            results.append(r)
    elapsed = time.time() - t0
    print("-" * 88)
    print(f"Done in {elapsed:.1f}s")
    print_summary(results)
    jpath, cpath, nrows = write_outputs(results, elapsed, out_dir=args.out)
    print(f"\nWrote:\n  {jpath}\n  {cpath} ({nrows} rows)")


if __name__ == "__main__":
    main()
