# Proximity Prize — Reproducible Parameter Calculator

A small, dependency-light (numpy / sympy / stdlib only) calculator for
Reed–Solomon **proximity-gap** parameters as used in hash-based SNARKs
(FRI / STIR / WHIR-style). It computes the three governing radii, composes
proximity-gap / MCA / list-size bounds into a soundness error, solves for the
minimum query count `t` reaching a `2^-128` target, and reports proof-size
proxies.

> **Status of the numbers.** The *radii* are exact mathematical facts. The
> *proximity-gap constants* are now **filled from the Wave-1 literature
> extraction** (`../literature/notes`). Each registered bound carries a
> `verified` flag: **VERIFIED** = traceable to a source note (the
> Bordage-Chiesa Thm 9.2 MCA error, the unique-decoding `n/|F|`, the
> Bordage-Chiesa Lem 10.1 interleaving, and the **field-agnostic** Kambiré /
> Crites-Stewart / BCHKS **proven near-capacity no-go**), **VERIFY** = a
> clearly-flagged placeholder for a *hidden / not-yet-pinned* constant or an
> *open* region (the BCHKS Thm 1.5 hidden `C_rho`, the list-size constant, and
> the `unknown-beyond-johnson` OPEN band — see below). The API is designed
> so call sites never change when a constant is pinned. Run
> `python3 cli.py bounds` / `bounds.verify_flags()` to see what is still
> unverified. **The capacity region above Johnson is modelled by TWO distinct
> bounds**, so the calculator's semantics match the math (the band between
> Johnson and the established ceiling is **OPEN**, not a proven no-go):
> `unknown-beyond-johnson` for the open band (vacuous `1.0` = "cannot certify",
> `verified=False`), and `proven-near-capacity-nogo` only at/above the split
> `nogo_split_radius` (the established failure, `verified=True`). **R13
> (s-integrality):** the split is the assembled-lemma per-field ceiling
> `(1−ρ)−2/s_max(b)` with `s` a **power of two** (`s_max = 16/16/32` at
> `b = 31/64/128`); at `b = 256` **no valid `s` exists** (count
> `3³² = 2⁵¹ ≪ 2¹²⁸` — no threshold-established Kambiré-type ceiling) and the
> split falls back to the generic CS25/Elias list-decoding-capacity ceiling,
> which stands at every field size. Both bounds are **field-agnostic** (prime
> and genuine odd-char extension fields share the same `δ_unsafe` / same
> field-agnostic mechanism; `../n2-verdict.md`); the assembled lemma's soft
> spots (KK25 distinctness cited not re-proved; ρ = 1/2 N1-conditional) are
> carried in `proven-near-capacity-nogo`'s notes, and the pre-R13 continuum
> `(1−ρ)−6/log₂|F|` is retained only as an asymptotic reference. The δ\*\_C
> tables and headline are **best-known provable brackets**, not a resolution
> of δ\*\_C.

---

## Modules

| File | Purpose |
|------|---------|
| `proximity_parameters.py` | Core radii (`unique_decoding_radius`, `johnson_radius*`, `capacity`), `qary_entropy` + `list_decoding_capacity_radius`, and cancellation-safe `log2` / `bits` / `log2_pow2_minus_c` helpers. Exact `Fraction`/`sympy` arithmetic. |
| `bounds.py` | A **registry** of proximity-gap / MCA / list-size bounds. Each is a self-describing `Bound` (name, source tag, formula string, `validity(...)`, `value(...)`, `verified` flag, notes). `register_bound`, `lookup`, `best_mca_bound`, `interleaved_mca`, `best_listsize_bound`, `verify_flags`. |
| `soundness.py` | Composite soundness error as a function of `t`; `min_query_count` (bisection solver); `CostModel` proof-size proxy; `evaluate_point` and `sweep_delta`; CSV writer. Threads an optional `field_type` (`"prime"`/`"extension"`/`None`) for back-compat; the two capacity-region bounds are now **field-agnostic** (prime and extension behave identically). Above Johnson the `ε_mca` floor is `1.0` either way (open band → `unknown-beyond-johnson`, "cannot certify"; near capacity → `proven-near-capacity-nogo`), so the numbers are unchanged — only the reported `mca_bound` name distinguishes the two. |
| `delta_star.py` | The **δ\*\_C provable-bracket** generator: `provable_safe_delta` (sweeps the Bordage-Chiesa trade parameter `m` for the largest δ with a *verified* positive bound ≤ `2^-128`), `kambire_unsafe_delta` established-unsafe radius `(1−ρ)−2/s_max(b)` (R13: `s` a power of two, `s_max = 16/16/32` at `b = 31/64/128`; `None` at 256-bit — no threshold-established ceiling; field-agnostic; `bounds.nogo_split_radius` supplies the open-band/no-go split, falling back to the generic CS25/Elias ceiling where the lemma does not fire), `delta_star_row`/`all_rows`, and the markdown emitter `write_tables` → `out/delta_star_tables.md`. |
| `cli.py` | `tables`, `sweep`, `bounds`, `delta-star` subcommands. |
| `tests/test_calculator.py` | pytest **and** plain-`__main__` runnable; asserts known anchors + solver properties + the δ\*\_C brackets. |

---

## The three radii (relative Hamming distance)

For rate `rho = k/n in (0,1)`:

```
unique decoding   = (1 - rho) / 2          # exact rational
Johnson           = 1 - sqrt(rho)          # sympy-exact / high-precision rational
list-dec capacity = 1 - rho                # exact rational
```

with the strict ordering `(1-rho)/2 < 1 - sqrt(rho) < 1 - rho`.

`qary_entropy(x, q) = x*log_q(q-1) - x*log_q(x) - (1-x)*log_q(1-x)` (standard
convention; `H_q(0)=0`, increasing on `[0, 1-1/q]`, max `1` at `1-1/q`).
`list_decoding_capacity_radius(rho, q)` solves `H_q(delta) = 1 - rho`; it tends
to the RS capacity `1-rho` only as `O(1/log2 q)` (documented in the docstring).

---

## Soundness cost model (and its caveats)

The composite soundness error at query count `t` is the conservative shape from
the research plan:

```
soundness_error(t)  ~  eps_mca(C^{equiv m}, delta)                 [batching / proximity term]
                    +  |Lambda(C^{equiv m}, delta)| / |F|           [interleaved list term]
                    +  (1 - delta)^t                                [t-fold query phase]
```

* `eps_mca` is taken from the best applicable MCA bound via
  `bounds.interleaved_mca(...)`, which for `m>1` selects the interleaving
  composition (currently the linear union bound `m * eps_mca(C,delta)`).
* the list term uses the best applicable list-size bound, divided by `|F|`.
* `(1-delta)^t` is the only `t`-dependent term; it drives the total down to the
  **`t`-independent floor** `eps_mca + |Lambda|/|F|`.

`min_query_count` first computes that floor. If `floor > target` the problem is
**infeasible for any `t`** (reported explicitly, with the floor in bits);
otherwise it bisects on `t in [0, t_max]` for the smallest `t` with
`error(t) <= target`, guaranteeing the anchor `error(t) <= target < error(t-1)`.

**Proof-size proxy** (`CostModel`, in bits):

```
per_query_bits = merkle_depth(n) * hash_size_bits        # authentication path
               + arity * field_elem_bits                  # opened leaf / coset
proof_size_bits = t * per_query_bits + commitment_overhead_bits
```

Defaults: `hash_size_bits=256`, `arity=2`, `field_elem_bits=ceil(logF)`,
`commitment_overhead_bits=256`. All are constructor knobs.

### Caveats (read before trusting a number)

1. **Composite proxy, not a protocol.** Real FRI/STIR/WHIR soundness has extra
   lower-order terms (per-round folding error, grinding/PoW, repetition
   structure, batching slack) and different constants. Use this to compare
   *regimes*, not to certify a deployment.
2. **Constants now filled (with explicit VERIFY flags).** The headline positive
   Johnson bound is Bordage-Chiesa Thm 9.2,
   `eps_mca ≤ (m+1/2)^7 · n^2 · d / (3 ρ^{3/2} |F|)` valid for
   `delta ≤ 1 - (1+1/(2m))√ρ` (verified). Its `n^2 · d` numerator is large, so
   over small fields the floor exceeds `2^-128` and the solver reports `INF`
   (infeasible) — this is *real*, not a placeholder artefact: certifying a
   single-code `2^-128` MCA needs `log2|F| ≳ 256` for `n` up to `2^30`. The
   still-`VERIFY` items are the BCHKS Thm 1.5 hidden `C_rho` (set to 1) and the
   Johnson list-size constant/exponents.
3. **The capacity region is split: OPEN band + proven near-capacity no-go.**
   The window `[1 - sqrt(rho), 1 - rho)` is NOT a single homogeneous "no-go" —
   that would contradict the project's own thesis that the band below the
   established ceiling is **OPEN**. It is split at `nogo_split_radius` into two
   distinct bounds:
   * **`unknown-beyond-johnson`** (`verified=False`) for the OPEN band from
     `J` up to the split: no certified *positive* MCA theorem applies
     (Bordage-Chiesa / BCHKS run out at Johnson), but **no impossibility is
     known either**. Its vacuous `1.0` means *"cannot certify security from
     current bounds"*, **NOT a proven no-go**. This is the band the Proximity
     Prize asks about.
   * **`proven-near-capacity-nogo`** (`verified=True`) only at/above the split:
     the thin band just below Singleton where the proximity gap / CA
     **provably fails** — Kambiré Thm 1 / Crites-Stewart Cor 1 / BCHKS
     Thm 1.13. Returns a genuine `1.0`. **R13:** the split is the
     assembled-lemma per-field ceiling `(1−ρ)−2/s_max(b)` (`s` a **power of
     two**; `s_max = 16/16/32` at `b = 31/64/128`; soft spots: KK25 cited not
     re-proved, ρ = 1/2 N1-conditional); at `b = 256` no valid `s` exists and
     the split falls back to the generic CS25/Elias list-decoding-capacity
     ceiling. The pre-R13 continuum `(1−ρ)−6/log₂|F|` (from the `8^s ≤ |F|`
     Linnik window) is an asymptotic reference only.

   Both are **field-agnostic** (prime AND genuine odd-characteristic extension
   alike), with the extension case established by the characteristic-zero
   cyclotomic-invariant argument of `../n2-verdict.md` (the bad-scalar count is
   identical for `GF(p)` and genuine `GF(p^e)`). Both return `1.0`, so the
   solver never extrapolates the Johnson form past its proven range and every
   number is unchanged from the earlier single-bound model. The earlier single
   broad `capacity-nogo` (which mislabelled the whole window as a proven no-go),
   and the even-earlier `capacity-prime-nogo` / `capacity-extension-open` split,
   are both superseded and removed.
4. **`log_q` convergence.** `list_decoding_capacity_radius` is the generic q-ary
   entropy radius; it approaches `1-rho` only slowly (`O(1/log2 q)`).

---

## Plugging in verified constants from `../literature`

Two equivalent routes, neither of which touches any call site:

* **Patch the tunables.** All exponents/constants live in module globals and a
  mirror dict `bounds.CONSTANTS`. To pin the BCHKS Thm 1.5 hidden constant once
  it is read off the paper and promote the bound to verified:

  ```python
  import bounds
  bounds._C_RHO_BCHKS = 1.0 / 48.0          # ledger value (example)
  bounds.get_bound("bchks-johnson").verified = True
  ```

* **Register a new/replacement bound.** Build a `Bound(...)` with the exact
  `validity` window and `value` formula and `register_bound(b, overwrite=True)`.
  `soundness.py` will pick it up automatically through `best_mca_bound` /
  `interleaved_mca` / `best_listsize_bound`.

Run `python3 cli.py bounds` or `bounds.verify_flags()` any time to see what is
still unverified.

---

## Example commands

```bash
# radii table + t/proof-size at delta = Johnson(rho) - eta, plus VERIFY surface
python3 cli.py tables

# tighter slack, larger domain, explicit field sizes (as log2|F|)
python3 cli.py tables --eta 0.02 --log2n 24 --logFs 64,128,256

# CSV sweep over (rho, delta, logF, log2n, m) -> calculator/out/sweep.csv
python3 cli.py sweep --rhos 1/2,1/4 --logFs 64,128,256 --log2ns 20,24 --ms 1,4 --points 60

# list every registered bound with provenance (shows VERIFIED vs VERIFY)
python3 cli.py bounds

# write the delta*_C provable-bracket tables -> out/delta_star_tables.md
python3 cli.py delta-star

# tests (either works)
python3 -m pytest tests/test_calculator.py -q
python3 tests/test_calculator.py
```

### Sample radii output

```
   rho |   unique-dec (1-r)/2 |    Johnson 1-sqrt(r) |   capacity 1-r
   1/2 |       1/4 = 0.250000 |             0.292893 |    1/2 = 0.5000
   1/4 |       3/8 = 0.375000 |             0.500000 |    3/4 = 0.7500
   1/8 |      7/16 = 0.437500 |             0.646447 |    7/8 = 0.8750
  1/16 |     15/32 = 0.468750 |             0.750000 |  15/16 = 0.9375
```

### Sample `t(delta)` rows (n = 2^20, m = 1, delta = Johnson − 0.05)

```
   rho |       field |    logF |    delta |       t | sec_bits |  floor_b |  proof_KB |        mca bound
   1/4 |  Goldilocks |   64.00 |  0.45000 |     INF |        - |     22.4 |         - |    bchks-johnson
   1/4 |    prime256 |  256.00 |  0.45000 |     149 |    128.5 |    214.4 |    102.47 |    bchks-johnson
```

(`INF` in the `t` column means the `t`-independent floor already exceeds the
target — no query count helps; see the `floor_b` column. Goldilocks (64-bit)
is infeasible for a single-code `2^-128` MCA; the 256-bit field reaches the
target. Note: just below the unique-decoding radius — e.g. ρ=1/2 at
δ=0.243 < UDR=0.25 — the cheaper `unique-decoding` bound `n/|F|` is selected
instead of the Johnson bounds.)

---

## Named fields

`cli.py` knows these (carried as `log2|F|`, computed cancellation-safely for the
exact ones):

| Field | Definition | log2\|F\| |
|-------|------------|-----------|
| Mersenne31 | `2^31 - 1` | ≈ 31.0 |
| BabyBear | `15·2^27 + 1 = 2^31 - 2^27 + 1` | ≈ 30.91 |
| Goldilocks | `2^64 - 2^32 + 1` | ≈ 64.0 |
| prime128 | 128-bit prime | 128.0 |
| prime256 | 256-bit prime | 256.0 |

Extension fields: pass any `--logFs` value (the calculator only ever uses
`log2|F|`, so any field / extension is parameterized by its bit-length). The
δ\*\_C generator (`delta_star.py`) still tags each field `prime`/`extension` for
labelling, but the capacity-window split (the `unknown-beyond-johnson` OPEN band
and the `proven-near-capacity-nogo` ceiling) is now **field-agnostic** — both
bounds apply identically over prime **and** genuine odd-characteristic extension
fields (Kambiré / Crites-Stewart / BCHKS for primes; the characteristic-zero
cyclotomic-invariant argument of `../n2-verdict.md` for extensions). It includes
`M31^4 (ext)`, `Goldilocks^2 (ext)` and a generic `ext256` to show that the
extension rows reproduce the matching-size prime `δ_unsafe = (1−ρ)−2/s_max(b)`
exactly (same bracket; R13 — at 256-bit both report no threshold-established
ceiling). Note `δ_unsafe` is the established-no-go *ceiling*, not the top of an
open interval: the genuinely **open** band sits between Johnson and it.
