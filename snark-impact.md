# Proximity Prize — SNARK Protocol Impact (Workstream F)

Date: 2026-06-02 (values aligned 2026-06-10 to the R13 s-integrality correction;
canonical per-field values: `calculator/out/delta_star_tables.md`). Convention:
**ABF δ-radius** (δ = relative proximity *radius*,
ρ = k/n, Johnson J = 1−√ρ, Singleton capacity = 1−ρ, list-decoding capacity
1−H_q(ρ)). Target ε* = 2⁻¹²⁸, n = 2²⁰ headline (n = 2²⁴, 2³⁰ trends noted).

This note translates the established MCA / list-size brackets into concrete
hash-based-SNARK (FRI / STIR / WHIR) parameters: query count `t`, Merkle
authentication cost, proof-size proxy, and verifier work — and asks whether
switching code families buys anything after overheads. **Every headline number
is reproduced from a calculator run; the commands and raw output are in the
[Appendix](#appendix--calculator-runs-every-headline-number).** The calculator
was run, not edited.

> **TL;DR.** Provably today, deployed smooth-domain RS is stuck at the Johnson
> radius for MCA. At ρ=1/2 that pins the per-query soundness factor at
> (1/√2)ᵗ — so you need **t ≈ 257 queries** for 128 bits where the capacity
> ideal would need **t = 128**. That is a flat **~2× proof-size penalty** at
> every rate (≈ **89 KB** extra at ρ=1/2, n=2²⁰, 256-bit field). Switching to
> capacity-reaching **folded / subspace-design RS does NOT win**: the folding
> overhead s ≳ 1/η² symbols per query inflates the proof **3.7×–230×**, dwarfing
> the query saving. The query/list phase must run over a **≥256-bit field**
> (where the list term is *provably* negligible at the operating Johnson radius —
> sub-problem 2's achievable radius is a proven bracket [J − o(1), r_E] and is never the
> binding constraint, MCA is); sub-problem 1 stays stuck at Johnson ⇒ live with
> ~2× queries, or recover the 64-bit gap with grinding/PoW.

---

## 1. The soundness model (ABF Lemma 6.6 composite) and its cost-model caveats

### 1.1 The model

All deployed protocols' round-by-round knowledge-soundness error reduces (ABF
survey, eprint 2026/680, Lemma 6.6) to the composite

```
soundness_error  =  max( ε_mca(C, δ) + |Λ(C^{≡2}, δ)| / |F| ,  (1 − δ)^t ).
```

Three quantities:

* **ε_mca(C, δ)** — the mutual-correlated-agreement error of the (interleaved)
  RS code: the batching/folding-soundness term. This is sub-problem 1.
* **|Λ(C^{≡2}, δ)| / |F|** — the interleaved list-size term (here the pair code
  C^{≡2}): the cost of the random folding combination landing in a large list.
  This is sub-problem 2.
* **(1 − δ)ᵗ** — the *query phase*: the probability that all `t` independent
  spot-checks miss the ≥ δ-fraction of corrupted positions. This is the **only
  `t`-dependent term** and the lever the prover/verifier pay for in proof size.

The calculator (`soundness.py`) uses the slightly more conservative *sum* shape
`ε_mca + |Λ|/|F| + (1−δ)ᵗ` and solves for the smallest `t` reaching the target;
the first two terms form a `t`-independent **floor**. If the floor already
exceeds ε* no finite `t` helps (the solver reports `INF`).

The whole design tension: **only `δ` controls the per-query catch rate
(1−δ).** Raising the provable-safe δ toward capacity makes each query worth more
bits, so fewer queries (smaller proofs) suffice. The Johnson barrier on δ is
therefore *directly* a proof-size cost.

### 1.2 Cost-model caveats (read before trusting a number)

1. **Composite proxy, not a protocol.** Real FRI/STIR/WHIR soundness carries
   extra lower-order terms (per-round folding error, grinding/PoW contributions,
   the round-by-round vs many-round repetition structure, batching slack) and
   different leading constants. The numbers here compare **regimes**; they do not
   certify any specific deployed system.
2. **Proof-size proxy** (`CostModel`, in bits):
   ```
   per_query_bits   = merkle_depth(n) · hash_size_bits   (authentication path)
                    + arity · field_elem_bits            (opened leaf / coset)
   proof_size_bits  = t · per_query_bits + commitment_overhead_bits
   ```
   Defaults used throughout: `hash_size_bits = 256`, `arity = 2` (binary Merkle),
   `field_elem_bits = ⌈log₂|F|⌉`, `commitment_overhead_bits = 256`. This is a
   *single-commitment* per-query proxy: it omits the multi-round STIR/WHIR
   commitment chain, FRI's per-round roots, and any batching of openings. It is
   deliberately linear in `t` so the Johnson-vs-capacity comparison is apples to
   apples.
3. **VERIFY constants in the floor.** The certified ε_mca floor uses the
   *verified* Bordage-Chiesa Thm 9.2 bound; the asymptotically-better BCHKS
   Thm 1.5 bound is `VERIFY` (hidden `C_ρ` set to 1). The list term uses a
   `VERIFY` Johnson-list placeholder `|Λ| ≤ n·m/η`. These affect *which field is
   feasible at all* (§2, §5) but **not** the Johnson-vs-capacity query ratio
   (§3), which depends only on δ through (1−δ)ᵗ. See §7.
4. **The ε_mca floor `= 1.0` above Johnson splits into two bounds — and `1.0`
   does NOT always mean "proven unsafe".** Above the Johnson radius the calculator's
   `ε_mca` floor is `1.0`, but it is reported under **two different bound names** with
   genuinely different meanings. In the **OPEN band** above Johnson and below the
   no-go split the bound is
   **`unknown-beyond-johnson`** (`verified=False`): the `1.0` means **"cannot certify
   — the region is OPEN"**, *not* a proven impossibility (no positive theorem is
   available there, but no proven no-go either). Only at **`δ ≥ (1−ρ)−2/s_max(b)`** —
   the per-field **assembled-lemma ceiling (R13)**, `s` a power of two, `s_max = 16/16/32`
   at `b = 31/64/128`; at `b = 256` no Kambiré-type ceiling is threshold-established and
   the split is the generic CS25/Elias ceiling ≈ r_E — does the bound become
   **`proven-near-capacity-nogo`** — a genuine (near-capacity) no-go. So a
   "floor = 1.0 / infeasible" cell in the open band
   is an *open* region, not a refutation; this is the two-bound model (`n2-verdict.md`,
   `calculator/bounds.py`). The δ_unsafe column elsewhere in this doc is the
   `proven-near-capacity-nogo` threshold (the R13 assembled-lemma ceiling, §6).

---

## 2. Per-(ρ, field) parameter tables at n = 2²⁰

Fields: Mersenne31 (2³¹−1, ≈31.0 b, prime), BabyBear (≈30.91 b, prime),
Goldilocks (2⁶⁴−2³²+1, ≈64.0 b, prime), a 128-bit prime, a 256-bit prime, and a
256-bit **extension** (`ext256`, e.g. M31⁸ / BabyBear⁸ / Goldilocks⁴). There is
now **no prime-vs-extension distinction**: the near-capacity no-go's *mechanism* is
**field-agnostic** (the Kambiré-type counterexample is a characteristic-zero cyclotomic
invariant, so genuine odd-char extensions inherit the same no-go as primes —
`n2-verdict.md`, §5/§6), and the feasibility of a single-code 2⁻¹²⁸ MCA *certificate*
depends only on log₂|F| anyway. Prime and extension rows are therefore identical at
equal log₂|F|. (The region above Johnson is governed by **two** calculator bounds:
`unknown-beyond-johnson` in the OPEN band between Johnson and the no-go split —
`ε_mca=1.0` = "cannot certify, OPEN", not a proven no-go — and
`proven-near-capacity-nogo` only at `δ ≥` the R13 assembled-lemma ceiling
`(1−ρ)−2/s_max(b)` (at 256-bit: no Kambiré-type ceiling is threshold-established,
the split is the generic CS25/Elias ceiling ≈ r_E); see §1.2 caveat 4. The δ_unsafe
column below is the former assembled-lemma threshold where it exists.)

Two readings of "δ" are reported, because they answer different questions:

* **δ_safe** (provable-safe radius, from `cli.py delta-star` /
  `out/delta_star_tables.md`): the largest δ at which the best **verified**
  positive MCA bound certifies ε_mca ≤ 2⁻¹²⁸. Capped at Johnson (the proven
  positive ceiling). The companion `t@safe` is the full composite query count
  there.
* **δ = Johnson − η** (from `cli.py tables`, η=0.05): a fixed near-Johnson radius,
  to show the floor behaviour and Merkle/proof costs uniformly.

### 2.1 Provable-safe table (the headline deliverable), n = 2²⁰

Source: `out/delta_star_tables.md` (`python3 cli.py delta-star`). `t@safe` is the
composite query count; "proof KB" uses the cost model above at that field's
log₂|F|. INF / infeasible = the `t`-independent floor exceeds 2⁻¹²⁸ (the
single-code MCA certificate is impossible at that field size — *not* a
placeholder artefact; see §2.3).

| ρ | field | log₂\|F\| | δ_safe | safe via | BC m | ε_mca bits | **t@safe** | bits | **proof KB** | δ_unsafe (field-agnostic, R13) | open gap |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1/2 | Mersenne31 | 31.0 | infeasible | — | — | — | — | — | — | 0.37500 | [0.000, 0.375] |
| 1/2 | BabyBear | 30.91 | infeasible | — | — | — | — | — | — | 0.37500 | [0.000, 0.375] |
| 1/2 | Goldilocks | 64.0 | infeasible | — | — | — | — | — | — | 0.37500 | [0.000, 0.375] |
| 1/2 | prime128 | 128.0 | infeasible | — | — | — | — | — | — | 0.43750 | [0.000, 0.438] |
| 1/2 | **prime256** | 256.0 | **0.29251** | B-C Thm9.2 | 934 | 128.0 | **272** | 128.0 | **187.0** | — (none established, R13) | [0.293, —) |
| 1/2 | **ext256** | 256.0 | **0.29251** | B-C Thm9.2 | 934 | 128.0 | **272** | 128.0 | **187.0** | — (none established, R13) | [0.293, —) |
| 1/4 | Goldilocks | 64.0 | infeasible | — | — | — | — | — | — | 0.62500 | [0.000, 0.625] |
| 1/4 | prime128 | 128.0 | infeasible | — | — | — | — | — | — | 0.68750 | [0.000, 0.688] |
| 1/4 | **prime256** | 256.0 | **0.49972** | B-C Thm9.2 | 889 | 128.0 | **136** | 128.0 | **93.5** | — (none established, R13) | [0.500, —) |
| 1/4 | **ext256** | 256.0 | **0.49972** | B-C Thm9.2 | 889 | 128.0 | **136** | 128.0 | **93.5** | — (none established, R13) | [0.500, —) |
| 1/8 | **prime256** | 256.0 | **0.64624** | B-C Thm9.2 | 846 | 128.0 | **91** | 128.0 | **62.6** | — (none established, R13) | [0.646, —) |
| 1/8 | **ext256** | 256.0 | **0.64624** | B-C Thm9.2 | 846 | 128.0 | **91** | 128.0 | **62.6** | — (none established, R13) | [0.646, —) |
| 1/16 | **prime256** | 256.0 | **0.74984** | B-C Thm9.2 | 805 | 128.0 | **68** | 128.0 | **46.8** | — (none established, R13) | [0.750, —) |
| 1/16 | **ext256** | 256.0 | **0.74984** | B-C Thm9.2 | 805 | 128.0 | **68** | 128.0 | **46.8** | — (none established, R13) | [0.750, —) |

(Mersenne31 / BabyBear / Goldilocks / prime128 rows for ρ ∈ {1/4,1/8,1/16} are
all `infeasible` for the single-code certificate, identical in form to the ρ=1/2
rows; their δ_unsafe values are in `out/delta_star_tables.md`.)

**The δ_unsafe column is the per-field *assembled-lemma* ceiling (R13)**
`(1−ρ) − 2/s_max(b)`: the Kambiré quotient parameter `s` must be a **power of
two** (s-integrality), and `s_max(b)` is the largest such `s` passing the KK25
distinctness calibration `p > φ(s)^{φ(s)}`, the prize-threshold count
`N(s,ρ) > 2^{b−128}`, and the above-Johnson condition `s > 2/(√ρ−ρ)` — giving
`s_max = 16/16/32` at `b = 31/64/128`. At `b = 256` **no valid `s` exists** (the
count `3³² = 2⁵¹ ≪ 2¹²⁸` fails the threshold), so no Kambiré-type ceiling below
`r_E` is threshold-established there — only the generic CS25/Elias ceiling at
≈ `r_E` stands. Lemma status: ASSEMBLED (KK25 cited not re-proved; ρ=1/2
N1-conditional); the pre-R13 continuum `(1−ρ)−6/log₂|F|` is an asymptotic
reference only. δ_safe (the proven Johnson floor) and the unconditional
M31 δ=1/2 witness are the proven anchors; "field-agnostic" holds because the
counterexample's bad-scalar count is a characteristic-zero cyclotomic invariant
(`n2-verdict.md`), so prime = extension at equal log₂|F|.

**The honest punchline (verified bound).** A single-code 2⁻¹²⁸ MCA certificate
is a **field-size constraint**. The verified Bordage-Chiesa error
(m+½)⁷·n²·d / (3ρ^{3/2}·|F|) has an n²·d numerator costing ≈ 2·log₂n + log₂(ρn)
bits (≈ 60 b at n=2²⁰, ≈ 90 b at n=2³⁰). Hence you need roughly

```
log₂|F| ≳ 128 + 2·log₂n + log₂(ρn) + 7·log₂(m+½) − 1.5·log₂ρ   ≈  256-bit field
```

for n up to 2³⁰. The calculator confirms the threshold: at ρ=1/2, n=2²⁰ the
minimum log₂|F| is **201 bits** (B-C trade m=3) rising to **257 bits** (m≈900,
the value that pushes δ_safe right up to Johnson); at n=2³⁰ it is **231–287
bits**. So **M31, BabyBear (31 b) and Goldilocks (64 b) cannot certify
ε_mca ≤ 2⁻¹²⁸ from this bound at all** — only ~256-bit fields can.

### 2.2 Uniform near-Johnson table (δ = Johnson − 0.05), n = 2²⁰

Source: `python3 cli.py tables` (m=1, default cost model). Shows the floor and
the chosen MCA bound directly. `floor_b` = the `t`-independent floor in bits.

| ρ | field | log₂\|F\| | δ | t | sec_bits | floor_b | proof KB | mca bound |
|---|---|---|---|---|---|---|---|---|
| 1/2 | Goldilocks | 64.0 | 0.24289 | INF | — | 39.6 | — | unique-decoding |
| 1/2 | prime128 | 128.0 | 0.24289 | INF | — | 103.6 | — | unique-decoding |
| 1/2 | prime256 | 256.0 | 0.24289 | **319** | 128.1 | 231.6 | **219.34** | unique-decoding |
| 1/4 | Goldilocks | 64.0 | 0.45000 | INF | — | 22.4 | — | bchks-johnson |
| 1/4 | prime128 | 128.0 | 0.45000 | INF | — | 86.4 | — | bchks-johnson |
| 1/4 | prime256 | 256.0 | 0.45000 | **149** | 128.5 | 214.4 | **102.47** | bchks-johnson |
| 1/8 | prime256 | 256.0 | 0.59645 | **98** | 128.3 | 214.4 | **67.41** | bchks-johnson |
| 1/16 | prime256 | 256.0 | 0.70000 | **74** | 128.5 | 214.4 | **50.91** | bchks-johnson |

(Goldilocks/prime128 at ρ ∈ {1/8,1/16} are likewise INF, floor_b ≈ 86.4 / 22.4;
full grid in the Appendix.) Note at ρ=1/2 the chosen δ=0.243 sits just *below*
the unique-decoding radius 0.25, so the cheaper UD bound `n/|F|` is selected (and
its floor is 231.6 b at 256-bit — feasible). Because n=2²⁰ the `tables` run uses
the conservative *unverified* `bchks-johnson` bound where it applies, so its `t`
values differ slightly from the verified-only `delta-star` table; both agree that
only the 256-bit field is feasible.

### 2.3 n = 2²⁴ and n = 2³⁰ trend

`δ_safe` *decreases* slightly as n grows (the B-C trade m must shrink to keep
ε_mca ≤ 2⁻¹²⁸ against the larger n²·d, so the validity ceiling
1−(1+1/2m)√ρ drops a little below Johnson). The query count `t@safe` is nearly
n-independent (it tracks δ_safe ≈ Johnson), but the proof **size** grows because
`merkle_depth(n)` grows linearly in log₂n. Concretely, at ρ=1/2, prime256:

| n | δ_safe | t@safe | proof KB | BC m |
|---|---|---|---|---|
| 2¹⁶ | 0.29278 | 276 | 155.3 | 3068 |
| 2²⁰ | 0.29251 | 272 | 187.0 | 934 |
| 2²⁴ | 0.29165 | 270 | 219.4 | 284 |
| 2³⁰ | 0.28537 | 273 | 273.0 | 47 |

So the proof grows ~10 KB per extra log₂n (the Merkle path lengthens), while `t`
barely moves — exactly the smooth-domain trend you expect: query count is set by
the *radius*, proof size by the *tree depth*.

---

## 3. Proof-size DELTA: "stuck at Johnson" vs "reach capacity"

This is the central impact result. Isolate the query term: to reach `b` bits from
(1−δ)ᵗ alone needs `t = ⌈b / (−log₂(1−δ))⌉`. The per-query bit yield is
`−log₂(1−δ)`.

### 3.1 Query-count delta (query-term only, 128 bits)

Source: query-term computation in the Appendix (radii from
`proximity_parameters`).

| ρ | Johnson δ | bits/query | **t @ Johnson** | capacity δ | bits/query | **t @ capacity** | **extra queries (J − cap)** |
|---|---|---|---|---|---|---|---|
| **1/2** | 0.29289 | 0.5000 | **257** | 0.500 | 1.0000 | **128** | **+129 (≈2.0×)** |
| 1/4 | 0.50000 | 1.0000 | 128 | 0.750 | 2.0000 | 64 | +64 (2.0×) |
| 1/8 | 0.64645 | 1.5000 | 86 | 0.875 | 3.0000 | 43 | +43 (2.0×) |
| 1/16 | 0.75000 | 2.0000 | 64 | 0.938 | 4.0000 | 32 | +32 (2.0×) |

**The ρ=1/2 story, exactly.** At the Johnson radius the catch probability per
query is 1−δ = √(1/2), so the query factor is (1/√2)ᵗ and **128 queries give only
(1/√2)¹²⁸ = 2⁻⁶⁴ — 64 bits.** Reaching the Singleton capacity δ→1/2 gives (1/2)ᵗ,
so (1/2)¹²⁸ = 2⁻¹²⁸ — the full 128 bits at **half the queries**. The Johnson
limitation **doubles the query count** at ρ=1/2 (257 vs 128). The same 2× holds
at every rate (it is the ratio −log₂(1−cap) / −log₂(1−J) = 2 for these RS rates,
since 1−J = √ρ and 1−cap = ρ, so the ratio is exactly log(ρ)/log(√ρ) = 2).

### 3.2 Proof-size delta in KB (cost model, log₂|F|=256, binary Merkle)

Source: `CostModel` run in the Appendix.

| ρ | n | per-query bits | t@J | **KB @ Johnson** | t@cap | **KB @ capacity** | **ΔKB (J−cap)** | ratio |
|---|---|---|---|---|---|---|---|---|
| **1/2** | 2²⁰ | 5632 | 257 | **176.72** | 128 | **88.03** | **+88.69** | 2.01× |
| 1/2 | 2²⁴ | 6656 | 257 | 208.84 | 128 | 104.03 | +104.81 | 2.01× |
| 1/2 | 2³⁰ | 8192 | 257 | 257.03 | 128 | 128.03 | +129.00 | 2.01× |
| 1/4 | 2²⁰ | 5632 | 128 | 88.03 | 64 | 44.03 | +44.00 | 2.00× |
| 1/8 | 2²⁰ | 5632 | 86 | 59.16 | 43 | 29.59 | +29.56 | 2.00× |
| 1/16 | 2²⁰ | 5632 | 64 | 44.03 | 32 | 22.03 | +22.00 | 2.00× |

**Being stuck at Johnson costs ≈ 89 KB at ρ=1/2, n=2²⁰** (and grows with n: ~105
KB at 2²⁴, ~129 KB at 2³⁰), a flat ~2× of the proof. The cost shrinks at lower
rates in absolute terms (22 KB at ρ=1/16) but stays a 2× factor.

### 3.3 Proof-size delta in verifier hash / eval counts

Verifier work is dominated by re-hashing authentication paths:
`hashes ≈ t · merkle_depth(n)`, plus O(arity) field evals per query for the
folding-consistency check. Source: Appendix.

| ρ | n | merkle depth | hashes @ Johnson | hashes @ capacity | **extra verifier hashes** |
|---|---|---|---|---|---|
| **1/2** | 2²⁰ | 20 | 5140 | 2560 | **+2580** |
| 1/2 | 2³⁰ | 30 | 7710 | 3840 | +3870 |
| 1/4 | 2²⁰ | 20 | 2560 | 1280 | +1280 |
| 1/8 | 2²⁰ | 20 | 1720 | 860 | +860 |
| 1/16 | 2²⁰ | 20 | 1280 | 640 | +640 |

The verifier does **~2580 extra hash compressions** (and ~2× the
folding-consistency field evals, since those scale with `t`) per proof at ρ=1/2,
n=2²⁰, purely because of the Johnson limitation. At n=2³⁰ it is ~3870 extra.

---

## 4. Effect of interleaving `m` and folding

### 4.1 Interleaving `m` (batching many RS instances into C^{≡m})

Two places `m` enters the composite:

* **MCA union (Bordage-Chiesa Lemma 10.1, verified):**
  ε_mca(C^{≡m}, δ) ≤ m · ε_mca(C, δ). A *benign* linear factor: it costs only
  **+log₂ m bits** of field (m=2 → +1 b, m=16 → +4 b). Below Johnson this is
  never the bottleneck.
* **Interleaved list (the upper relation |Λ(C^{≡m},δ)| ≤ |Λ(C,δ)|ᵐ):** the list
  *log* is multiplied by `m`. So the list term contributes
  m · log₂|Λ(C,δ)| − log₂|F| bits to the floor. Since the single-code Johnson
  list is a constant 1/(2ηρ) independent of n, even mᵗʰ-power it stays tiny
  against a ≥256-bit field. **At the operating Johnson radius this is proven**
  (Johnson Cor. 3.3 gives a constant list, so |Λ(C^{≡m},J)|/|F| ≪ 2⁻¹²⁸); the
  largest δ at which the list term *alone* stays ≤2⁻¹²⁸ is a proven bracket
  [J − o(1), r_E] (§5.3, `listsize_resolution.py`) and is never the binding constraint.

In the calculator, `m` also drives the Bordage-Chiesa **trade parameter** (larger
m pushes the validity window toward Johnson but inflates the (m+½)⁷ error). That
(m+½)⁷ factor is what forces the ≥256-bit field at the *very top* of the Johnson
window: e.g. ρ=1/2, n=2²⁰ needs log₂|F| ≥ 201 b at m=3 but **257 b at m≈900** (to
reach δ_safe=0.2925 ≈ Johnson). So interleaving is cheap for *soundness* but the
B-C trade is what makes "Johnson exactly" expensive in field size.

**Net:** for the deployed regime (constant m, e.g. 2), interleaving adds ≤ a few
field bits and does not change the §3 query-count story.

### 4.2 Folding / alphabet overhead for capacity-reaching code families

The capacity-reaching positive results are for **other code classes**, and they
buy their larger δ at an alphabet/folding cost:

* **Folded / subspace-design RS (Goyal-Guruswami 2025/2054):** reach MCA radius
  δ = (1−R) − η ≈ capacity−η, but the folding parameter is **s ≳ 1/η²**, and the
  code is over the alphabet **F^s** — i.e. *each query opens s field symbols* (a
  folded coset), not 1. To approach capacity (η→0) the per-query leaf cost blows
  up quadratically in 1/η.
* **Random-punctured RS (Gao-Cai 2025/870):** reaches near capacity but needs a
  **random** evaluation domain — not the smooth (FFT) domain deployments require.
  Off the table for plain deployed RS.

The quantitative consequence is §5.

---

## 5. Three implementation choices compared — does any proof-size win survive?

We compare at the headline point **ρ=1/2, n=2²⁰, log₂|F|=256**, against the
baseline of deployed smooth RS at the Johnson radius (t=257, arity-2 leaf,
**176.72 KB**).

### 5.1 (A) Deployed smooth-domain RS — the baseline

Provable-safe δ = Johnson = 0.293. t = 257, proof **176.72 KB**, floor feasible
only at ≥256-bit field. **This is what you can prove today.** Everything else is
measured against it.

### 5.2 (B) Folded / subspace-design RS (Goyal-Guruswami) — reaches capacity, but…

Each query opens s ≳ 1/η² symbols over F^s. Compare the folded proof at
δ = capacity−η against the smooth-RS Johnson baseline. Source: Appendix.

| η | folded δ = cap−η | s ≈ 1/η² | t_fold | leaf bits/query | path bits/query | **proof KB** | **vs smooth-RS @ Johnson** |
|---|---|---|---|---|---|---|---|
| 0.100 | 0.4000 | 100 | 174 | 25 600 | 5 120 | 652.5 | **3.69×** |
| 0.050 | 0.4500 | 400 | 149 | 102 400 | 5 120 | 1 955.7 | **11.07×** |
| 0.025 | 0.4750 | 1 600 | 138 | 409 600 | 5 120 | 6 986.3 | **39.53×** |
| 0.010 | 0.4900 | 10 000 | 132 | 2 560 000 | 5 120 | 41 332.5 | **233.89×** |

**The proof-size win does NOT survive the folding overhead.** Reaching capacity
cuts the query count roughly in half (257 → ~130–174), but the s ≳ 1/η² symbols
opened *per query* dominate the proof: even at the loosest η=0.10 (which only
reaches δ=0.40, still short of capacity 0.5), the folded proof is **3.7× larger**
than smooth-RS at Johnson; pushing
δ toward capacity (η small) makes it catastrophically worse (40×–230×). The
folding-leaf cost is the killer: 25.6 Kb–2.56 Mb *per query* of opened coset vs
the 512 bits (2 symbols) at arity-2 smooth RS.

(Caveats: this uses the calculator's same per-query Merkle proxy with the leaf
enlarged to s symbols; it ignores any subspace-design-specific list/commitment
savings and any proof-of-work that folded constructions also use. The exact GG
constant in s ≳ 1/η² is a `VERIFY` item — but even s = 1/η² *exactly* gives the
3.7×–230× above, and any realistic constant > 1 makes it worse, so the qualitative
verdict is robust.)

### 5.3 (C) Large-field sub-problem-2 regime — list term provably negligible at Johnson, MCA still stuck

If you run the query/interleaving phase over a **large extension (≥256-bit)** —
which deployments already do — then at the **operating Johnson radius** the list
term is **provably negligible**: the Johnson list bound (ABF Cor. 3.3) gives a
constant single-code list, so |Λ(C^{≡2}, J)|/|F| ≪ 2⁻¹²⁸ over any ≥256-bit field.
The list term is therefore **never the binding constraint** — MCA is.

How far could the list term *alone* stay ≤2⁻¹²⁸? That is the sub-problem-2
quantity δ\*_C^{(2)}, governed by the budget B = 2⁻¹²⁸·|F|. The honest answer is a
**proven bracket** (`listsize_resolution.py`):

| ρ (256-bit, n=2²⁰, m=2) | δ\*_C^{(2)} ∈ [J − o(1), r_E] | regime |
|---|---|---|
| 1/2 | **[0.293, 0.496]** | BRACKET (upper reach conjectural) |
| 1/4 | [0.500, 0.747] | BRACKET (upper reach conjectural) |
| 1/8 | [0.646, 0.873] | BRACKET (upper reach conjectural) |
| 1/16 | [0.750, 0.936] | BRACKET (upper reach conjectural) |

The lower end (Johnson floor) is **proven** — to `J − η_min` with
`η_min = 1/(2ρB^{1/m})`, written `J − o(1)`; the upper end r_E (formula
convention `r_E = 1−H_q(ρ)`; the mathematically exact proven object is the
inverse-entropy crossing `H_q^{-1}(1−ρ)`, the two differing only beyond the
displayed precision here) is a
**proven ceiling** (above it the worst-case list is q^{Ω(n)} ≫ B); whether
δ\*_C^{(2)} actually *reaches* ≈r_E is **conjectural** (it needs the open
smooth-domain large-list bound — see technical-note §4). For proof size this gap is
**moot**: MCA caps the usable δ at Johnson regardless, and the list term is already
provably fine there. The field-size lever is sharp: **|F| = 2¹²⁸ is the
knife-edge.** Below it (M31 2⁻⁹⁷, Goldilocks 2⁻⁶⁴, even M31⁴ ≈ 2⁻⁴) the budget
B < 1 and **δ\*_C^{(2)} does not exist** — the list term *cannot* be driven to
2⁻¹²⁸ at all, so the phase must run in a ≥256-bit extension. At exactly 2¹²⁸,
δ\*_C^{(2)} = (1−ρ)/2 (unique decoding, exact).

**But this does not rescue MCA.** The list term is **decoupled** from sub-problem 1
(the BCHKS Thm 1.9 barrier ties MCA-beyond-Johnson to RS list-decoding at list
size q, a separate hard problem; sub-problem 2 only needs list ≤ 2⁻¹²⁸·q ≪ q, so it
is strictly easier and not gated by Thm 1.9). So in regime (C): the *list* term is
provably fine up to Johnson (and bracketed [J − o(1), r_E] beyond), but the *MCA* term
still caps the provable-safe δ at Johnson, and the query count stays at the §3
t≈257. **Big field neutralises the list term; it does not move the query count.**

### 5.4 Verdict

| choice | provable δ (ρ=1/2) | query count | proof @ n=2²⁰ | net vs baseline |
|---|---|---|---|---|
| **(A) smooth RS @ Johnson** | 0.293 | 257 | 176.7 KB | baseline |
| **(B) folded/subspace-design RS** | →0.49 (cap−η) | ~130–174 | 652 KB – 41 MB | **3.7×–230× WORSE** |
| **(C) smooth RS, ≥256-bit field** | 0.293 (MCA); list ∈ [J−o(1), r_E] | 257 | 176.7 KB | list provably negligible at J; **MCA unchanged** |

**No proof-size win survives.** (B)'s query saving is obliterated by the s ≳ 1/η²
folding alphabet. (C) is what deployments already do and is necessary (sub-problem
2 is degenerate below 2¹²⁸), but it leaves the MCA query count exactly where the
Johnson barrier puts it. The only thing that would genuinely shrink the proof is a
**positive MCA theorem for plain smooth-domain RS above Johnson at ρ=1/2** — the
open prize question.

---

## 6. Concrete recommendation for hash-based SNARK designers

**1. Safe radius to use today: the Johnson radius, δ ≤ 1−√ρ** (ρ=1/2 → 0.293,
ρ=1/4 → 0.5, ρ=1/8 → 0.646, ρ=1/16 → 0.75). This is the proven positive ceiling
for MCA on deployed smooth RS. Do **not** operate in the (Johnson, capacity)
window: the whole band is uncertifiable (OPEN), and its top is a near-capacity
no-go (Kambiré / Crites-Stewart / BCHKS). The
*proven* part is the **unconditional M31 witness, which fails right at δ=1/2**
(BCHKS, ρ≈1/2); the per-field no-go ceilings are the **assembled-lemma (R13)**
values `δ_unsafe = (1−ρ)−2/s_max(b)` — at ρ=1/2: **0.375** (31/64-bit) and
**0.4375** (128-bit) — with `s` a power of two and `s_max` set by the KK25
distinctness calibration plus the prize-threshold count (lemma ASSEMBLED: KK25
cited not re-proved; ρ=1/2 N1-conditional). At **256-bit no Kambiré-type ceiling
is threshold-established** (count `3³² ≪ 2¹²⁸`); the generic CS25/Elias ceiling
at ≈ r_E stands. The pre-R13 continuum `(1−ρ)−6/log₂|F|` (e.g. 0.477 at 256-bit)
is an asymptotic reference only, *not* a proven endpoint. The
*mechanism* of this no-go is **field-agnostic** — it holds **identically over prime
and genuine odd-characteristic extension base fields** (the counterexample is a
characteristic-zero cyclotomic invariant, `n2-verdict.md`). Extensions are **not**
a safe haven: same δ_unsafe, same band.

**2. Required query count (128-bit soundness, query term):** **t ≈ 257 at ρ=1/2**
(257/128/86/64 at ρ = 1/2 / 1/4 / 1/8 / 1/16). This is ~**2× the capacity ideal**
(128/64/43/32). Budget the proof accordingly: ≈ **177 KB at ρ=1/2, n=2²⁰** in the
arity-2 / 256-bit cost model, growing ~10 KB per extra log₂n.

**3. Run the query/interleaving phase over a large extension field
(≥256-bit / |F| > 2¹²⁸).** This is not optional: (i) the single-code MCA
*certificate* is infeasible below ~256 bits (the n²·d numerator); (ii) sub-problem
2 is **degenerate** below 2¹²⁸ (the list term can't reach 2⁻¹²⁸). Over a ≥256-bit
extension, **the list term is provably negligible at the operating Johnson radius**
(Johnson list bound ⇒ |Λ(C^{≡2},J)|/|F| ≪ 2⁻¹²⁸); the largest radius the list term
alone tolerates is a proven bracket [J − o(1), r_E] (≈[0.293, 0.496] at ρ=1/2), but that
is moot since MCA already binds at Johnson. This matches deployed FRI/STIR/WHIR
practice (small FFT base field, big extension for the query phase).

**4. Do NOT switch to folded / subspace-design RS for proof size.** They reach
capacity for MCA, but the s ≳ 1/η² folding alphabet makes the proof **3.7×–230×
larger** than smooth-RS-at-Johnson (§5.2). The query-count win is real but is
swamped by per-query coset openings. (They may still be worth it for *other*
reasons — linear-time prover, recursion structure — but not to shrink the proof.)

**5. To recover the ~64-bit Johnson shortfall without doubling queries, use
grinding / proof-of-work.** At ρ=1/2 the query-ideal t=128 yields only 64 bits at
Johnson; adding `w` PoW bits makes up the rest. Tradeoff (Appendix):

| queries t (ρ=1/2, Johnson) | query bits | PoW bits w to reach 128 | prover PoW cost |
|---|---|---|---|
| 128 | 64.0 | 64.0 | 2⁶⁴ hashes (heavy) |
| 160 | 80.0 | 48.0 | 2⁴⁸ |
| 192 | 96.0 | 32.0 | 2³² (cheap) |
| 257 | 128.5 | 0 | none |

A common sweet spot is t≈192 + 32 bits of grinding (2³² prover hashes, negligible)
to land the proof between the Johnson (257-query) and capacity (128-query) sizes.

**Practical upshot (one line):** *operate at the Johnson radius over a ≥256-bit
extension field — the list term is provably negligible there (sub-problem 2 is not
the bottleneck), but sub-problem 1 (MCA) is stuck at Johnson, so pay ~2× the
queries vs the capacity ideal (≈257 vs 128 at ρ=1/2, ≈89 KB / ~2580 verifier
hashes extra at n=2²⁰), or buy back the 64-bit gap with grinding/PoW; switching to
folded/subspace-design RS does NOT reduce proof size after the s ≳ 1/η² alphabet
overhead.*

---

## 7. Caveats and VERIFY items affecting the numbers

| Item | Effect on the numbers | Status |
|---|---|---|
| **ABF Lemma 6.6 composite is a proxy** | Real FRI/STIR/WHIR add per-round folding error, repetition structure, batching slack, different constants. Use for regime comparison, not certification. | model caveat |
| **Cost model is single-commitment, linear-in-t** | Omits STIR/WHIR multi-round commitment chains and FRI per-round roots; omits opening batching. Absolute KB are proxies; the **Johnson-vs-capacity ratio (§3) is robust** (depends only on (1−δ)ᵗ). | model caveat |
| **BCHKS Thm 1.5 hidden C_ρ** (`bchks-johnson`, VERIFY, set =1) | Its **linear-n** error would relax the field-size floor: if C_ρ≈1, ρ=1/2,n=2²⁰ needs only log₂\|F\| ≈ 170 b (vs 201–257 b for verified B-C). **Could make ~192-bit fields feasible for the certificate** — but unverified, so the δ_safe table uses only the verified B-C bound. Does **not** affect §3 (query ratio). | **VERIFY** C_ρ |
| **Johnson list-size placeholder** `|Λ| ≤ n·m/η` (`interleave-listsize`, VERIFY) | Sets the list term in the floor. Over ≥256-bit fields it is negligible (confirmed by the independent `listsize_resolution.py` giving ≈capacity); matters only for the exact floor at marginal field sizes. | **VERIFY** constant + exponents |
| **Goyal-Guruswami folding constant** in s ≳ 1/η² | §5.2's 3.7×–230× uses s = ⌈1/η²⌉ exactly. A larger real constant makes folded RS *worse*; the qualitative "does not win" verdict is robust to it. | **VERIFY** constant (verdict robust) |
| **Near-capacity no-go: field-agnostic mechanism, assembled-lemma ceiling (R13); OPEN band is OPEN** | The near-capacity no-go `δ_unsafe = (1−ρ)−2/s_max(b)` is the per-field **assembled-lemma ceiling (R13)** (`s` a power of two; `s_max = 16/16/32` at `b = 31/64/128`; at 256-bit **none threshold-established** — count `3³² ≪ 2¹²⁸` — the generic CS25/Elias ceiling at ≈ r_E stands; KK25 cited not re-proved, ρ=1/2 N1-conditional; the pre-R13 continuum `(1−ρ)−6/log₂\|F\|` is an asymptotic reference only) whose *mechanism* holds for **both** prime and genuine odd-char extension fields — the Kambiré-type counterexample is a characteristic-zero cyclotomic invariant (`n2-verdict.md`), sub-threshold at 256-bit and not a Conj-2 refutation. The region above Johnson splits into two bounds: in the **OPEN band** between Johnson and the no-go split the floor `1.0` is `unknown-beyond-johnson` = *cannot certify (OPEN)*, **not** a proven no-go; only at `δ ≥` the assembled-lemma ceiling is it `proven-near-capacity-nogo`. (The *positive* MCA theorem above Johnson remains open, for all field types alike — not an extension-specific opening.) δ_safe stays at the proven Johnson floor; the recommendation and the gap are identical for prime and extension. | mechanism field-agnostic; ceiling **ASSEMBLED (R13)**; OPEN band OPEN |
| **Kambiré ρ=1/2 boundary** | Kambiré states ρ ∈ (0,1/2) open; ρ=1/2 exactly is covered by the unconditional M31 BCHKS witness (prime). Affects the δ_unsafe column at ρ=1/2, not δ_safe or t. | verified (via BCHKS M31) |
| **List-decoding-capacity radius** (`list_decoding_capacity_radius`) | The generic q-ary entropy radius; approaches Singleton only as O(1/log₂q). §5.3's δ\*_C^{(2)} values carry this gap (0.496 not 0.5 at ρ=1/2, 256-bit). | exact (documented convergence) |

---

## Appendix — calculator runs (every headline number)

Run from `calculator/`. Calculator was executed, not modified.

### A.1 `python3 cli.py tables`  (radii + t/proof-size at δ=Johnson−0.05, m=1, n=2²⁰)

```
========================================================================
Reed-Solomon radii (relative Hamming distance)
========================================================================
   rho |   unique-dec (1-r)/2 |    Johnson 1-sqrt(r) |   capacity 1-r
------------------------------------------------------------------------
   1/2 |       1/4 = 0.250000 |             0.292893 |    1/2 = 0.5000
   1/4 |       3/8 = 0.375000 |             0.500000 |    3/4 = 0.7500
   1/8 |      7/16 = 0.437500 |             0.646447 |    7/8 = 0.8750
  1/16 |     15/32 = 0.468750 |             0.750000 |  15/16 = 0.9375
------------------------------------------------------------------------
ordering check: (1-r)/2 < 1-sqrt(r) < 1-r  for all rho in (0,1)

====================================================================================================
Query count t and proof-size proxy at delta = Johnson(rho) - eta,  eta = 0.05
target = 2^-128, m = 1, log2(n) = 20, hash = 256 bits, arity = 2
====================================================================================================
   rho |       field |    logF |    delta |       t | sec_bits |  floor_b |  proof_KB |        mca bound
--------------------------------------------------------------------------------------------------------
   1/2 |  Goldilocks |   64.00 |  0.24289 |     INF |        - |     39.6 |         - |  unique-decoding
   1/2 |    prime128 |  128.00 |  0.24289 |     INF |        - |    103.6 |         - |  unique-decoding
   1/2 |    prime256 |  256.00 |  0.24289 |     319 |    128.1 |    231.6 |    219.34 |  unique-decoding
   1/4 |  Goldilocks |   64.00 |  0.45000 |     INF |        - |     22.4 |         - |    bchks-johnson
   1/4 |    prime128 |  128.00 |  0.45000 |     INF |        - |     86.4 |         - |    bchks-johnson
   1/4 |    prime256 |  256.00 |  0.45000 |     149 |    128.5 |    214.4 |    102.47 |    bchks-johnson
   1/8 |  Goldilocks |   64.00 |  0.59645 |     INF |        - |     22.4 |         - |    bchks-johnson
   1/8 |    prime128 |  128.00 |  0.59645 |     INF |        - |     86.4 |         - |    bchks-johnson
   1/8 |    prime256 |  256.00 |  0.59645 |      98 |    128.3 |    214.4 |     67.41 |    bchks-johnson
  1/16 |  Goldilocks |   64.00 |  0.70000 |     INF |        - |     22.4 |         - |    bchks-johnson
  1/16 |    prime128 |  128.00 |  0.70000 |     INF |        - |     86.4 |         - |    bchks-johnson
  1/16 |    prime256 |  256.00 |  0.70000 |      74 |    128.5 |    214.4 |     50.91 |    bchks-johnson
--------------------------------------------------------------------------------------------------------
```

VERIFY surface (same run): `bchks-johnson` (hidden C_ρ=1), `interleave-listsize`
(constant/exponents), and `unknown-beyond-johnson` (the OPEN-band marker,
`verified=False`) are the unverified items. The capacity region is now **two bounds**:
`unknown-beyond-johnson` (fires in the OPEN band between Johnson and the no-go split,
`ε_mca=1.0` = "cannot certify — region OPEN", *not* a proven no-go) and
`proven-near-capacity-nogo` (fires only at `δ ≥ (1−ρ)−2/s_max(b)` — the R13
assembled-lemma ceiling, `s` a power of two, `s_max = 16/16/32` at `b = 31/64/128`;
at 256-bit no Kambiré-type ceiling is threshold-established and the split is the
generic CS25/Elias ceiling ≈ r_E). The latter's *mechanism* is **field-agnostic**
(`n2-verdict.md`): the Kambiré-type negative applies identically to prime and genuine
odd-char extension fields, so an extension's δ_unsafe equals the matching-size prime
value. (The old single `capacity-nogo` — and the earlier split bound that
left the extension upper end unresolved — have been **removed**; the whole window is no
longer mislabelled as a proven no-go.) The pre-R13 continuum 6/log₂|F| gap is retained
in the `proven-near-capacity-nogo` notes as an asymptotic reference only.
Tunable constants: `{_C_UD:1.0, _C_BC:0.333…, _EXP_M_BC:7.0, _EXP_N_BC:2.0,
_BC_M_MIN:3, _C_RHO_BCHKS:1.0, _EXP_ETA_BCHKS:5.0, _KAMBIRE_GAP_CONST:6.0,
_C_L:1.0, _P_L:1.0, _Q_L:1.0, _S_L:1.0}`.

### A.2 `python3 cli.py delta-star`  (provable-safe brackets; writes out/delta_star_tables.md)

```
ρ=1/2 provable brackets [delta_safe, delta_unsafe] (n = 2^20; δ_unsafe field-agnostic — prime = extension):
    Mersenne31 (    prime): delta_safe=infeasible  t@safe=    —  delta_unsafe=   0.37500
      prime256 (    prime): delta_safe=   0.29251  t@safe=  272  delta_unsafe=      OPEN
        ext256 (extension): delta_safe=   0.29251  t@safe=  272  delta_unsafe=      OPEN

Johnson(1/2) = 0.29289; capacity(1/2) = 0.5 -> at Johnson (1-delta)^t = (1/sqrt2)^t,
so t=128 gives only 2^-64 (64 bits).
```

The `ext256` and `prime256` rows are **identical**: per `n2-verdict.md` the
Kambiré-type no-go is **field-agnostic**, so at every field size the
genuine-extension `delta_unsafe` equals the matching-size prime value
(0.375 at 31/64-bit, 0.4375 at 128-bit). At 256-bit both print `OPEN`: no
Kambiré-type ceiling below r_E is threshold-established (R13), so the bracket
upper end is open there. (The earlier split bound that left the extension
upper end unresolved has been removed from the calculator.)

Full per-(ρ, field, n) table (incl. proof KB, BC m, ε_mca bits, n ∈ {2¹⁶,2²⁰,2²⁴,2³⁰})
is in **`calculator/out/delta_star_tables.md`**. Key ρ=1/2 prime256 rows:
2²⁰ → δ_safe 0.29251, t 272, **187.0 KB**, BC m 934; 2³⁰ → δ_safe 0.28537, t 273,
**273.0 KB**, BC m 47. ρ=1/4/1/8/1/16 prime256/ext256 at 2²⁰: t = 136 / 91 / 68,
proof = 93.5 / 62.6 / 46.8 KB.

### A.3 `python3 listsize_resolution.py`  (sub-problem 2: the field-size lever)

```
THE FIELD-SIZE LEVER:  eps* * |F| = 2^-128 * |F|  is the entire list-size budget.
         M31    31    2^-97 < 1  => delta*_C DOES NOT EXIST (degenerate)
       M31^2    62    2^-66 < 1  => delta*_C DOES NOT EXIST (degenerate)
  Goldilocks    64    2^-64 < 1  => delta*_C DOES NOT EXIST (degenerate)
       M31^4   124    2^-4  < 1  => delta*_C DOES NOT EXIST (degenerate)
     128-bit   128    = 1   => |Lambda|<=1 (unique decoding of interleaved obj)
     192-bit   192    = 2^64  => |Lambda|<=2^64 (loose)
     256-bit   256    = 2^128 => |Lambda|<=2^128 (loose)
Crossover: |F| = 2^128 is the knife-edge.  Below it sub-problem 2 is
vacuous (no delta*); at it, delta* = (1-rho)/2 (exact); above it, delta*
lies in the PROVEN BRACKET [J - o(1), r_E] with r_E = 1 - H_q(rho) ~
(1-rho) - 1/log2|F| (Johnson proves J - eta_min at finite budget, not J;
reaching ~r_E from the Johnson floor is CONJECTURAL, not proven).

SHARP VERDICT FOCUS:  rho = 1/2, n = 2^20, m = 2
  Johnson J = 0.29289, capacity = 0.50000, UD = 0.25000
        M31 (2^ 31): regime=DEGENERATE  delta*_C^(2) = DOES NOT EXIST
  Goldilocks (2^ 64): regime=DEGENERATE  delta*_C^(2) = DOES NOT EXIST
    128-bit (2^128): regime=   BINDING  delta*_C^(2) = 0.25000  (exact, UD radius)
    192-bit (2^192): regime=   BRACKET  delta*_C^(2) in PROVEN BRACKET [0.29289, 0.49479]  (upper reach ~0.49479 CONJECTURAL)
    256-bit (2^256): regime=   BRACKET  delta*_C^(2) in PROVEN BRACKET [0.29289, 0.49609]  (upper reach ~0.49609 CONJECTURAL)
```

(Proven floor = Johnson, in the `J − o(1)` convention — Johnson proves `J − η_min`
at finite budget; proven ceiling = r_E, quoted in the formula convention
`r_E = 1 − H_q(ρ)` at the field's true q (R29/M2; the exact proven object is the
inverse-entropy crossing `H_q^{-1}(1−ρ)`); the upper reach to ≈r_E is
conjectural. m=3 identical to resolved precision. 256-bit ceiling r_E =
0.49609 / 0.74683 / 0.87288 / 0.93618 for ρ = 1/2 / 1/4 / 1/8 / 1/16.)

### A.4 `python3 cli.py bounds`  (registered bounds + provenance)

The capacity-region no-go is now **split into two bounds** (this is the key honesty
fix): `unique-decoding` (VERIFIED, n/|F|), `bordage-chiesa-johnson`
(VERIFIED, (m+½)⁷n²d/(3ρ^{3/2}|F|)), `bchks-johnson` (**VERIFY** C_ρ, C_ρ·n/(η⁵|F|)),
**`unknown-beyond-johnson`** (`verified=False`; fires in the OPEN band between Johnson
and the no-go split,
value `ε_mca = 1.0` meaning **"cannot certify — the region is OPEN"**, *not* a proven
impossibility), **`proven-near-capacity-nogo`** (`verified=True`; fires **only** at
`δ ≥ (1−ρ)−2/s_max(b)` — the R13 assembled-lemma ceiling, `s` a power of two,
`s_max = 16/16/32` at `b = 31/64/128`; at 256-bit no Kambiré-type ceiling is
threshold-established and the split is the generic CS25/Elias ceiling ≈ r_E — and is
**field-agnostic**:
applies identically to prime and genuine odd-char extension fields, `n2-verdict.md`),
`interleave-mca-union` (VERIFIED, m·ε_mca), `interleave-listsize` (**VERIFY**,
C_L·nᵖ·mˢ/η^q). Full text in the run output. **The two capacity bounds both report
`ε_mca = 1.0` above Johnson, but they mean different things:** in the open band the `1.0`
is `unknown-beyond-johnson` = *cannot certify (OPEN)*, while only at `δ ≥` the
assembled-lemma
ceiling does the `1.0` become `proven-near-capacity-nogo` = a genuine (near-capacity)
no-go. (The old single `capacity-nogo` — and the earlier `capacity-prime-nogo` /
`capacity-extension-open` split — mislabelled the *whole* window as a proven no-go and
are gone. The pre-R13 continuum 6/log₂|F| gap is retained in the
`proven-near-capacity-nogo` notes as an asymptotic reference only.)

### A.5 Derived proof-size / query computations (using `proximity_parameters`, `soundness.CostModel`)

**Query-term-only t and bits/query (128-bit target), Johnson vs capacity vs UD:**
```
  rho        J    t@J   bpq@J |    cap  t@cap  bpq@cap |      UD   t@UD | extra t (J vs cap)
  1/2  0.29289    257  0.5000 |  0.500    128   1.0000 | 0.25000    309 |                129
  1/4  0.50000    128  1.0000 |  0.750     64   2.0000 | 0.37500    189 |                 64
  1/8  0.64645     86  1.5000 |  0.875     43   3.0000 | 0.43750    155 |                 43
 1/16  0.75000     64  2.0000 |  0.938     32   4.0000 | 0.46875    141 |                 32
```

**Proof-size KB (CostModel: hash=256, arity=2, logF=256, overhead=256), Johnson vs capacity:**
```
  rho  log2n  depth  pq_bits |   t@J     KB@J |  t@cap   KB@cap |  dKB(J-cap)  ratio
  1/2     20     20     5632 |   257   176.72 |    128    88.03 |       88.69   2.01
  1/2     24     24     6656 |   257   208.84 |    128   104.03 |      104.81   2.01
  1/2     30     30     8192 |   257   257.03 |    128   128.03 |      129.00   2.01
  1/4     20     20     5632 |   128    88.03 |     64    44.03 |       44.00   2.00
  1/8     20     20     5632 |    86    59.16 |     43    29.59 |       29.56   2.00
 1/16     20     20     5632 |    64    44.03 |     32    22.03 |       22.00   2.00
```

**Verifier hash count (t · merkle_depth), Johnson vs capacity:**
```
  rho  log2n  depth |   t@J  hashes@J |  t@cap  hashes@cap |  extra hashes
  1/2     20     20 |   257      5140 |    128        2560 |          2580
  1/2     30     30 |   257      7710 |    128        3840 |          3870
  1/4     20     20 |   128      2560 |     64        1280 |          1280
 1/16     20     20 |    64      1280 |     32         640 |           640
```

**Folded-RS (s ≈ 1/η²) vs smooth-RS @ Johnson (176.72 KB), ρ=1/2, n=2²⁰, logF=256:**
```
   eta delta_fold s~1/eta^2  t_fold  leaf_bits  path_bits  KB_fold  vs smoothJ
 0.100     0.4000       100     174      25600       5120   652.53       3.69x
 0.050     0.4500       400     149     102400       5120  1955.66      11.07x
 0.025     0.4750      1600     138     409600       5120  6986.28      39.53x
 0.010     0.4900     10000     132    2560000       5120 41332.53     233.89x
```

**Field-size floor (verified B-C vs hypothetical pinned BCHKS), ρ=1/2:**
```
 B-C rho=1/2 n=2^20 m=3: need logF >= 201.2 bits
 B-C rho=1/2 n=2^20 m=900: need logF >= 257.2 bits   (pushes delta_safe to ~Johnson)
 B-C rho=1/2 n=2^30 m=3: need logF >= 231.2 bits
 B-C rho=1/2 n=2^30 m=900: need logF >= 287.2 bits
 BCHKS (if C_rho=1) rho=1/2 n=2^20 eta=0.05: need logF >= 169.6 bits   [VERIFY]
 BCHKS (if C_rho=1) rho=1/2 n=2^30 eta=0.01: need logF >= 191.2 bits   [VERIFY]
```

**Grinding/PoW to recover the Johnson gap (ρ=1/2, 0.5 bits/query):**
```
  t= 128 queries @ Johnson -> 64.0 query bits; need w=64.0 PoW bits to reach 128
  t= 192 queries @ Johnson -> 96.0 query bits; need w=32.0 PoW bits to reach 128
  t= 257 queries @ Johnson -> 128.5 query bits; need w=0.0 PoW bits to reach 128
```
