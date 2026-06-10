# Proximity Prize — Formal Problem Ledger (Workstream A)

Date: 2026-06-02. Convention: **ABF δ-radius** throughout (δ = relative proximity
*radius*, ρ = k/n, Johnson J = 1−√ρ, Singleton capacity = 1−ρ; list-decoding
capacity radius `R_cap := H_q^{-1}(1−ρ)`, the exact inverse-entropy crossing —
the R41 naming of the proven Elias/CS25 ceiling object, previously written `r_E`;
the nearby formula `1−H_q(ρ)` is its labeled **large-q approximation**, the two
differing by ≤ 0.0017 at deployed rates with rate-dependent sign). Sources are the
Wave-1 extraction notes in `literature/notes/*.md`; every numeric constant below
is either traceable to a note or flagged `VERIFY`. Numbers in the tables are
produced by the calculator (`calculator/proximity_parameters.py`,
`calculator/delta_star.py`) using exact `Fraction`/`sympy`/cancellation-safe
`log2` arithmetic — not hand-rounded.

The prize (ABF survey, eprint 2026/680) poses two grand challenges for the
deployed code `C = RS[F, L, k]` with `L` a **smooth** domain (multiplicative
subgroup / coset of order a power of two), `ρ = k/|L| ∈ {1/2, 1/4, 1/8, 1/16}`,
`ε* = 2⁻¹²⁸`, in the regime `k ≤ 2⁴⁰`, `|F| < 2²⁵⁶`:

1. **Grand MCA challenge.** Determine the largest `δ*_C` with
   `ε_mca(C, δ*_C) ≤ ε*`, **with a proof that for all `δ > δ*_C`,
   `ε_mca(C, δ) > ε*`.**
2. **Grand list-decoding challenge.** For constant interleaving `m`, determine
   the largest `δ*_C` with `|Λ(C^{≡m}, δ*_C)| ≤ ε*·|F|`, with the analogous
   "fails for all `δ > δ*_C`" proof.

Both feed the deployed protocols' round-by-round knowledge-soundness error,
which has the exact composite form (ABF Lemma 6.6)

```
soundness_error  =  max( ε_mca(C, δ) + |Λ(C^{≡2}, δ)|/|F| ,  (1 − δ)^t ).
```

> **This ledger records best-known PROVABLE BRACKETS, not a resolution.** For
> every rate family we give `[δ_known_positive, δ_known_negative]`. The true
> `δ*_C` lies inside that bracket; pinning it is the open prize question.

---

## 1. Normalized definitions (ABF δ-radius convention)

Each definition cites its source note. `C = RS[F, L, k] ⊆ F^n`, `n = |L|`,
`ρ = k/n`, relative Hamming distance `Δ(·, ·)`.

**Proximity-gap error `ε_pg`.** (ABF Def; BCHKS Def 1.1, in `bchks.md`.) `C` has
a proximity gap up to radius `δ` if for every pair `f, g ∈ F^n`, either many
points of the affine line `{f + z·g : z ∈ F}` are δ-close to `C` *and* the pair
is δ-close to the 2-interleaved code, or only a few are. `ε_pg(C, δ)` is the
probability over `z` of the "bad" event — a line that is δ-close at a random `z`
without the pair being globally δ-close. In BCHKS's accounting this is `a/q`,
where `a = #{z : Δ(f+z·g, C) ≤ δ}` is the count of exceptional scalars and
`q = |F|`.

**Correlated-agreement error `ε_ca`.** (Crites-Stewart Def 2, in
`crites-stewart.md`; BCHKS Def 1.1.) `u₀,…,u_ℓ` have correlated agreement (CA)
with `C` of density `≥ 1−δ` if there is a common set `D′ ⊆ L`, `|D′|/n ≥ 1−δ`,
and codewords `v₀,…,v_ℓ ∈ C` with `u_i = v_i` on all of `D′` for every `i`.
`ε_ca(C, δ)` is the probability that a random line combination is δ-close while
no such common agreement set exists. `ε_ca = 1` means CA is *totally broken*
(every line point is close, yet no global structure) — the Crites-Stewart Cor 1
regime.

**Mutual-correlated-agreement error `ε_mca`.** (ABF Def 4.3; Bordage-Chiesa
Def 3.14, in `bordage-chiesa.md`.) The MCA strengthening: with combination
generator `G(x)·U`, MCA requires that whenever a δ-fraction agreement set `T`
makes `G(x)·U|_T ∈ C|_T` for a random `x`, then in fact *every* `u_j|_T ∈ C|_T`
(the agreement is "mutual" — all components, not just the combination, lie in
the code on `T`). `ε_mca(C, δ)` is the failure probability. MCA ⇒ CA
(Bordage-Chiesa Lem 3.22), so MCA is the stronger notion; the prize uses
`ε_mca`, which is **loss-free** (ABF Remark 4.4 — no proximity-loss slack).

**Proximity loss `ε*` (a.k.a. δint−δfld).** (BCHKS, in `bchks.md`; the symbol
clashes — see §2.) The additive gap between the per-line radius the hypothesis
sees and the conclusion radius the theorem delivers: `Δ([f,g], C²) ≤ δ + ε*`.
"Lossless" / `ε* = 0` means the conclusion holds at the *same* radius δ. ABF's
`ε_mca` carries **no** proximity loss by definition; BCHKS Thm 1.5 achieves
`ε* = 0` up to the Johnson radius.

**Line / correlated agreement combination.** (All notes.) The base case is the
**affine line** generator `x ↦ (1, x)`, i.e. the combination `f + z·g` for a
random scalar `z ∈ F`. More generally a **polynomial generator**
`G(x) = (P_j(x))_j` with `P_j` linearly independent (Bordage-Chiesa Def 3.1);
univariate powers `(1, x, …, x^{ℓ−1})` give the "curve" / power-combination of
degree ℓ−1 (ACFY / WHIR). FRI/STIR/WHIR folding reduces to the affine line, so a
line-level statement bounds the whole folding hierarchy.

**Interleaved / batched code `C^{≡m}`.** (Bordage-Chiesa §10, in
`bordage-chiesa.md`; ABF.) The `m`-fold interleaving: codewords are `m`-tuples
`(c₁,…,c_m)`, `c_i ∈ C`, viewed as a single code over the alphabet `F^m` on the
same domain `L`. Distance is the *coordinate-wise* (row) Hamming distance — a
position disagrees if any of the `m` components does. `C^{≡2}` is the pair code
appearing in the soundness expression. Bordage-Chiesa Lem 10.1:
`ε_mca(C^{≡m}, δ) ≤ m·ε_mca(C, δ)` (a benign linear union).

**Smooth domain.** (ABF Def 2.12; used in `kambire.md`, `bchks.md`.) `L ⊆ F` is
smooth if it is a multiplicative subgroup (or coset) of order a **power of two**,
`|L| = 2^r` — the FFT-friendly evaluation domain deployed in FRI/STIR/WHIR. (NB
the BCHKS / Kambiré negative constructions use multiplicative subgroups
`⟨ω⟩` of order `n = 2^t`; `𝔽_q^*` itself, order `q−1`, is *not* strictly smooth
— see the VERIFY note in §4.)

**List size `|Λ(C^{≡m}, δ)|`.** (ABF; Crites-Stewart Def 1, in
`crites-stewart.md`.) `Λ(C^{≡m}, δ)` is the list of interleaved codewords within
relative radius δ of a given (worst-case) word: for the underlying code,
`Λ(C, δ) = {v ∈ C : Δ(u, v) ≤ δ}` for the worst `u`. The grand list-decoding
challenge asks for the largest δ with `|Λ(C^{≡m}, δ)| ≤ ε*·|F| = 2⁻¹²⁸·|F|`.
Inside the Johnson radius the RS list size is poly(n)/poly(η) (Johnson bound).

---

## 2. Notation-clash table (ABF vs BCHKS vs Crites-Stewart vs Goyal-Guruswami)

Condensed from `literature/notes/_litA_synthesis.md`. **Carry the ABF convention
everywhere; the most dangerous clash is BCHKS's `δ` = minimum distance vs
everyone else's `δ` = radius, and BCHKS/GG's `γ` meaning opposite things.**

| Symbol | ABF (2026/680) | Crites-Stewart (2025/2046) | BCHKS (2025/2055) | Goyal-Guruswami (2025/2054) |
|---|---|---|---|---|
| proximity radius | **δ** | **δ = f/n** (`f` = abs. errors) | **γ** | **δ** (= 1−R−η) |
| code min-distance | δmin(C) (= "capacity") | 1−ρ (Singleton) | **δ** (= 1−k/n) | 1−R (Singleton) |
| soundness error | ε_pg / ε_ca / ε_mca | ε | **a/q** | **err** |
| proximity loss | (none in ε_mca) | — | **ε*** | **γ** (conclusion slack) |
| rate | ρ | ρ = k/n | ρ = k/n = 1−δ | R |
| line combination | `f + γ·f₂` | `u₀ + z·u₁` | `f + z·g` | power curve `Σ u_j α^j` |
| "capacity" | Singleton 1−ρ | **list-decoding** 1−H_q(ρ) | Singleton 1−δ | Singleton 1−R |

Two load-bearing consequences of the clash:
- When importing **BCHKS Thm 1.5** into the ABF convention: ABF radius `δ` =
  BCHKS `γ`; ABF `ρ` = BCHKS `ρ = 1−δ_BCHKS`; ABF error `n/(η⁵|F|)` = BCHKS
  `a/q` with `a = O(n/η⁵)`; ABF `η = 1−√ρ − δ` = BCHKS `η` (Johnson slack). The
  calculator's `bchks-johnson` bound performs exactly this translation.
- **Crites-Stewart's "capacity" is the *list-decoding* capacity `1−H_q(ρ)`**,
  strictly below the Singleton `1−ρ` that ABF/BCHKS/GG call capacity. The gap
  `H_q(ρ) − ρ ≤ 1/log₂q` (CS Claim 1) is itself a load-bearing quantity (≤ 1/31
  for 31-bit fields).

---

## 3. Target parameter family

### 3a. The three rate-only radii

UD, Johnson and Singleton capacity depend **only on ρ** (not the field):

| ρ | unique-dec (1−ρ)/2 | Johnson 1−√ρ | Singleton cap 1−ρ |
|---|---|---|---|
| 1/2 | 0.250000 | 0.292893 | 0.500000 |
| 1/4 | 0.375000 | 0.500000 | 0.750000 |
| 1/8 | 0.437500 | 0.646447 | 0.875000 |
| 1/16 | 0.468750 | 0.750000 | 0.937500 |

(Strict ordering `(1−ρ)/2 < 1−√ρ < 1−ρ` holds for all ρ ∈ (0,1); asserted in the
calculator tests.)

### 3b. List-decoding capacity radius `1−H_q(ρ)` (field-dependent)

The fourth "radius" — the generic q-ary list-decoding capacity, the
Crites-Stewart ceiling — *does* depend on the field via `q = |F|`. It sits
`H_q(ρ)−ρ ≤ 1/log₂q` below the Singleton capacity. Fields:
M31 = 2³¹−1 (≈31.0 bits), BabyBear = 15·2²⁷+1 (≈30.91), Goldilocks = 2⁶⁴−2³²+1
(≈64.0), a 128-bit prime, a 256-bit prime. The gap to Singleton is in
parentheses:

| ρ | Singleton | M31 / BabyBear (≈31b) | Goldilocks (64b) | 128-bit | 256-bit |
|---|---|---|---|---|---|
| 1/2 | 0.500000 | 0.4678 (gap 0.0322) | 0.4844 (0.0156) | 0.4922 (0.0078) | 0.4961 (0.0039) |
| 1/4 | 0.750000 | 0.7225 (0.0275) | 0.7370 (0.0130) | 0.7436 (0.0064) | 0.7468 (0.0032) |
| 1/8 | 0.875000 | 0.8558 (0.0192) | 0.8661 (0.0089) | 0.8707 (0.0043) | 0.8729 (0.0021) |
| 1/16 | 0.937500 | 0.9251 (0.0124) | 0.9319 (0.0056) | 0.9348 (0.0027) | 0.9362 (0.0013) |

The gap halves as `log₂q` doubles — the `O(1/log₂q)` rate (calculator
`list_decoding_capacity_radius`, tested). Domain sizes of interest: `n = 2^r`,
`r ∈ {16, 20, 24, 30}` (all `≤ 2⁴⁰`). Target `ε* = 2⁻¹²⁸`.

---

## 4. Exact known bounds as functions of (ρ, n, |F|, δ, η, m)

All source-tagged. The calculator (`calculator/bounds.py`) implements each;
`verified=True` = traceable to a note, `verified=False` = a flagged `VERIFY`
placeholder for a hidden constant.

### Positive (deployed smooth RS, all fields)

**[P1] Unique decoding** — classical BCIKS20 (via `crites-stewart.md`).
For `δ ≤ (1−ρ)/2`:
```
ε_mca(C, δ)  ≤  n / |F|.
```
*(verified.)* Leading constant 1, numerator `n`.

**[P2] Bordage-Chiesa Thm 9.2 — the headline positive MCA bound**
(`bordage-chiesa.md`). For `C = RS[F, D, k]` over **any** domain `D` (incl. smooth
multiplicative subgroups), the affine-line / univariate-powers generator,
message degree `d = k−1 = ⌊ρn⌋−1`, and any integer `m ≥ 3`:
```
ε_mca(C, δ)  ≤  (m + 1/2)⁷ · n² · d / ( 3 · ρ^{3/2} · |F| ),
                              valid for  δ ≤ 1 − (1 + 1/(2m))·√ρ.
```
*(verified.)* Loss-free (= ABF `ε_mca`). As `m → ∞` the validity window → the
Johnson radius `1−√ρ`; larger `m` pushes δ toward Johnson but inflates the error
via the `(m+1/2)⁷` factor — the **m-vs-δ trade** the solver sweeps. Matches
BCIKS20 correlated agreement and lifts it to MCA + all polynomial generators.

**[P3] BCHKS Thm 1.5 — better n-scaling, hidden ρ-constant** (`bchks.md`; = ABF
Thm 4.12 core). In the ABF convention, with `η = 1−√ρ − δ > 0`:
```
ε_mca(C, δ)  ≤  C_ρ · n / ( η⁵ · |F| ),        valid for  δ ≤ 1 − √ρ − η,   ε* = 0.
```
*(VERIFY: the hidden `O_ρ(·)` constant `C_ρ`; calculator placeholder `C_ρ = 1`.)*
The η-exponent 5 is firm (BCHKS Table 1, "this work", `a = Θ_δ(n/η⁵)`); the
leading constant is an unspecified `O_ρ(·)` (asymptotic leading term
`1/(48·ρ^{3/2})` only for `η ≪ √ρ`). **Linear in `n`** (vs P2's `n²·d ~ n³`), so
asymptotically the stronger bound — but kept unverified until `C_ρ` is pinned.

**[P-comp] Interleaving — Bordage-Chiesa Lem 10.1** (`bordage-chiesa.md`):
```
ε_mca(C^{≡m}, δ)  ≤  m · ε_mca(C, δ).
```
*(verified.)* A benign linear union; not the bottleneck below Johnson.

→ **δ_known_positive = Johnson `1 − √ρ`** for smooth RS over all fields (P2/P3).

### Negative (deployed smooth RS)

**[N1] Crites-Stewart Cor 1** — `ε_ca = 1` (`crites-stewart.md`). For `q ≥ 10`,
**`q ≥ n`**, in the k-window
`n(1−H_q(δ)) + 2 + √(n·H_q(δ)−f) ≤ k ≤ n−f−2`, there exist `u⁽⁰⁾, u⁽¹⁾` with
`Δ(u⁽¹⁾, C) > δ` yet `Δ(u⁽⁰⁾ + λu⁽¹⁾, C) ≤ δ` for **every** `λ ∈ F`. So CA fails
*outright* reaching down to the **list-decoding capacity `1−H_q(ρ)`**; `δ*_C`
is provably strictly below the Singleton `1−ρ`. **Claim 1:** Singleton →
list-dec capacity gap `≤ 1/log₂q`.

**[N2] Kambiré Thm 1** — the explicit smooth-domain prime-field counterexample
(`kambire.md`). For every `C > 0` and rate `ρ ∈ (0, 1/2)` there are infinitely
many `n`, `k` with **`C = RS[F_p, ⟨ω⟩, k]`, `D = ⟨ω⟩` a smooth multiplicative
subgroup, `n = 2^t`, prime `p ≡ 1 mod n`, `p < n^A` (`A = K·log 8`)**, and words
`f = X^{rm}`, `g = X^{(r−1)m}` such that at
```
δ = (1 − ρ) − 2/(K(ρ,C)·log n)          (i.e. η = 2/s = 2/(K log n) below capacity)
```
there are `≥ n^C` scalars `z` with `Δ(f+z·g, C) ≤ δ`, yet `Δ([f,g], C²) > δ` — a
direct proximity-gap **and** CA failure. **Prime fields only; ρ = 1/2 excluded**
(stated for the open interval). **Finite-field caveat (key number) — a *finite-field
inferred ceiling*, not a proven endpoint:** the construction needs `8^s ≤ |F|`, so
`s ≤ log₂|F|/3` and the gap is a **constant** `η ≥ 6/log₂|F|` (≈ 0.047 at 128-bit,
≈ 0.023 at 256-bit) — **NOT o(1)** at deployed sizes. The calculator encodes this as
the headline `δ_unsafe = (1−ρ) − 6/log₂|F|`, which is **inferred** arithmetic from the
`[4^s, 8^s]` Linnik window (VERIFY) — *not* a proven endpoint.

> **Assembled lemma (pointer; full statement in `technical-note.md` §3.2).**
> *Lemma (per-field instantiated negative ceiling — ASSEMBLED, R13).* For `F = F_p`
> with `b = log₂p` bits, smooth domain `n = 2^t | p−1`, quotient parameter `s = n/m`
> a **power of two** (R13 integrality): under (i) **distinctness** `p > φ(s)^{φ(s)}`
> (KK25; cited, not re-proved), (ii) **count** `N(s,ρ) > 2^{b−128}`, with `N(s,ρ)`
> the per-rate bad-scalar count (Kambiré `N0_sum` for ρ < 1/2; the proven antipodal
> `3^{s/2}` for ρ = 1/2, carrying the N1 conditionality flag), and (iii) the
> construction's standing distinctness/degree hypotheses — then
> `ε_mca(C, δ) > 2⁻¹²⁸` at `δ = (1−ρ) − 2/s`, hence
> `δ*_C ≤ (1−ρ) − 2/s_max(b)` with `s_max(b)` the largest power of two satisfying
> (i)+(ii): `s_max = 16/16/32` at `b = 31/64/128`; at `b = 256` **no valid `s`
> exists** (count `3³² = 2⁵¹ ≪ 2¹²⁸`) — the mechanism establishes **no threshold
> ceiling** there. *Status: **ASSEMBLED from cited components** (Kambiré Thm 1
> count + KK25 distinctness + exact threshold arithmetic in
> `experiments/small_rs_atlas/n2_prize_threshold.py` + R13 power-of-two
> integrality); soft spots: literal KK25 lemma cited, not re-proved; ρ = 1/2
> N1-conditional; deg-convention (`N0_fire` vs `N0_sum`) noted.* The continuum
> `δ_unsafe = (1−ρ) − 6/log₂|F|` remains an asymptotic reference only (a
> finite-field inferred ceiling, not part of the lemma); the proven negatives
> underneath it are the asymptotic Kambiré result and the unconditional M31
> `δ = 1/2` witness (N3).

**[N3] BCHKS Thm 1.13** — prime multiplicative-subgroup loss (`bchks.md`).
Conditional on additive-combinatorics Conj 1.12, prime-field
multiplicative-subgroup RS at `γ = δmin − Θ(1/log n)` forces proximity loss
`Θ(1/log n)`. **UNCONDITIONAL for M31, ρ ≈ 1/2 — over the NON-SMOOTH domain `𝔽_q^*`:**
`C = RS[F_{2³¹−1}, 𝔽_q^*,
(q−1)/2 + 2/log₂q]` has `f, g` that are 1/2-close for *all* `z` yet CA-distance
`≥ 1/2 + 1/62 ≈ 0.516` (loss ≈ 0.016) — a genuine failure at radius δ = 1/2 for
a deployed prime at rate ≈ 1/2. (`𝔽_q^*` has order `q−1 = 2(2³⁰−1)`, not a power
of two — not smooth in ABF's sense; whether the witness lands on a strictly
power-of-two subgroup is VERIFY, per `bchks.md`. Not a proven smooth-domain
prize witness as stated.)

**[N4] BCHKS Thm 1.9 — THE KEY BARRIER** (`bchks.md`; all RS). For
`γ = LDR_{F,D,q}(δ) + 2/n` (just past the list-decoding radius for list size q),
∃ `f, g` with `≥ q/(2n)` close line points but `Δ([f,g], C²) ≥ δ − 1/n`, so
`a/q ≥ 1/(2n)` independent of q. **Improving MCA / proximity gaps beyond the
Johnson radius for any RS code REQUIRES first improving its list-decoding radius
(for list size q) beyond Johnson** — a separately hard / long-open problem.

→ **δ_known_negative (prime smooth RS) = capacity − a small constant**: clean
asymptotic `(1−ρ) − 2/(K log n)` (N2, proven), finite-field
`(1−ρ) − 6/log₂|F|` (a **finite-field inferred ceiling** — inferred from the `8^s ≤ |F|`
window, VERIFY, asymptotic reference only; the per-field threshold ceiling
`(1−ρ) − 2/s_max(b)` is the assembled lemma above, R13), with the **proven**
unconditional M31 witness at δ = 1/2 (N3) and the list-dec-capacity ceiling
`1 − H_q(ρ)` (N1).

### Wrong-code-class positives (do NOT lift δ*_C for smooth RS)

Goyal-Guruswami 2025/2054 (capacity `1−R−η` for **folded / subspace-design /
random** RS), Gao-Cai 2025/870 (√-lossy list-decodability ⇒ proximity-gaps
reduction; strong RS instantiations need **random puncturing / folding**),
Yuan-Zhu (random linear codes only). These reach near capacity but for code
classes *other than* plain smooth RS, so they sharpen the contrast without
raising `δ*_C` for the deployed code.

---

## 5. Gap intervals `[δ_known_positive, δ_known_negative]` per rate family

Positive side `δ_known_positive = 1−√ρ` (P2/P3) is field-type-agnostic. The
negative side is **also field-agnostic** (corrected; see `n2-verdict.md`):

- **Prime fields:** `δ_known_negative = (1−ρ) − 6/log₂|F|` (N2, finite-field
  Kambiré) — a **finite-field inferred ceiling** (continuum value inferred from the
  `8^s ≤ |F|` window, VERIFY, asymptotic reference only; the per-field threshold
  version is the **assembled lemma in §4 (R13)**; **R13 correction:** `s` must be a
  power of two, so the honest per-field ceiling is `(1−ρ) − 2/s_max(b)` with
  `s_max(b)` the largest power-of-two `s` passing the KK25 distinctness calibration
  `p > φ(s)^{φ(s)}` — see the table below — and at 256-bit **no
  threshold-established ceiling below `R_cap` exists** — i.e., no Kambiré-type
  near-capacity, constant-η ceiling, from any known mechanism; the generic
  CS25/Elias ceiling at ≈ `R_cap` stands at every field size); the *proven* negative anchor is
  the unconditional M31 / BCHKS witness at ρ≈1/2, δ = 1/2. The gap is *closed at the
  top* by this inferred ceiling (proven at the M31 δ=1/2 point, inferred along the
  continuum).
- **Genuine odd-characteristic extension fields:** `δ_known_negative` is the
  **same** inferred ceiling `(1−ρ) − 6/log₂|F|` (same R13 power-of-two correction
  applies). The Kambiré-type smooth-domain counterexample is
  **field-agnostic** — its distinct-bad-scalar count is a characteristic-zero
  cyclotomic invariant (`e₁,e₂` are reductions mod the characteristic `p` of
  fixed elements of `ℤ[ξ_s]`, independent of the extension degree `e`), so genuine
  `GF(p^e)` realizes the same count as `GF(p)` (verified exactly, `n2-verdict.md`).
  The earlier "OPEN — Linnik is prime-specific" reading is superseded: the same
  near-capacity band is foreclosed for extensions too. (The construction is
  sub-threshold at 256-bit and does not refute CGHLL Conj 2, exactly as for
  primes — so it forecloses only the same constant-`η = 2/s` band, not more.)

### Prime *and* genuine odd-char extension fields (the gap is bracketed; negative endpoint is a finite-field inferred ceiling)

(δ_positive = **proven** Johnson floor; δ_negative = the **finite-field inferred
ceiling**, per the R13 correction `(1−ρ) − 2/s_max(b)` with `s_max(b)` the largest
power-of-two `s` passing the KK25 distinctness calibration `p > φ(s)^{φ(s)}` —
`s_max = 32` at b=128 — per the **assembled lemma in §4 (R13)**; soft spots: KK25
cited not re-proved, ρ=1/2 N1-conditional. At 256-bit the
mechanism is **infeasible at the prize threshold** (count `3³² = 2⁵¹ ≪ 2¹²⁸`), so
there is **no threshold-established 256-bit ceiling below `R_cap` (R13)** — i.e., no
Kambiré-type near-capacity, constant-η ceiling; the generic CS25/Elias ceiling at
≈ `R_cap` stands at every field size; the 256-bit
mechanism-asymptotic reference is `(1−ρ) − 2/64` in the R13 power-of-two
convention (`s = 64`, the largest distinctness-allowed power of two — primary,
matching the N1 figures), with the pre-R13 continuum `(1−ρ) − 6/256` retained
only as a secondary continuum reference. The
negative *mechanism* is field-agnostic.)

| ρ | δ_positive (Johnson, proven) | δ_negative @128-bit (assembled §4 lemma, R13: s_max=32) | δ_negative @256-bit (mechanism-asymptotic only, NOT threshold-established — R13 power-of-two `(1−ρ)−2/64`; pre-R13 continuum in brackets) | open gap @128-bit |
|---|---|---|---|---|
| 1/2 | 0.29289 | 0.4375 | (0.46875, reference only; continuum 0.47656) | [0.293, 0.438] |
| 1/4 | 0.50000 | 0.6875 | (0.71875, reference only; continuum 0.72656) | [0.500, 0.688] |
| 1/8 | 0.64645 | 0.8125 | (0.84375, reference only; continuum 0.85156) | [0.646, 0.813] |
| 1/16 | 0.75000 | 0.875 | (0.90625, reference only; continuum 0.91406) | [0.750, 0.875] |

(**R13 correction — the old M31 number IMPROVES.** At 31 bits the old
continuum value `(1−ρ) − 6/31` (ρ=1/2 → 0.30645) was **too pessimistic**: `s` must
be a power of two, and under the KK25 `φ(s)^{φ(s)}` calibration `s_max = 16` at
b=31, giving `δ_negative = (1−ρ) − 2/16`, e.g. ρ=1/2 → **0.375**.
**BabyBear/circle re-scoping caveat (R29/M4):** the direct 31-bit anchor is
**BabyBear** (`p = 15·2²⁷+1`, 2-adicity 27), where the assembled lemma applies
verbatim; literal M31 has `p−1 = 2·(2³⁰−1)` with `v₂(p−1) = 1`, hence **no smooth
domain of order ≥ 4**, and the assembled lemma is **vacuous over literal M31** —
the M31 figure refers to the deployed circle/`F_{p²}` smooth subgroup of order
`p+1 = 2³¹` via the field-agnostic N2 mechanism (distinctness calibrated at the
characteristic `p ≈ 2³¹`, giving the same `s_max = 16`), that row inheriting N2's
**ESTABLISHED-MODULO-FORMALIZATION** grade rather than the prime-field ASSEMBLED
grade. The prize's open window at M31, ρ=1/2 is therefore [0.293, 0.375] — wider
than previously stated, **not** "almost closed" (the re-scoping strengthens, not
weakens, the withdrawal: the direction of the correction is unchanged). Same
`s_max = 16` values at b=64 (Goldilocks: 0.40625 was too *optimistic*, corrected
to 0.375). See §6.)

### Extension fields (same bracket as prime — field-agnostic)

Genuine odd-characteristic extension fields use the **same** table as the prime
table above: `δ_negative = (1−ρ) − 2/s_max(b)` per R13 (continuum
`(1−ρ) − 6/log₂|F|` mechanism-asymptotic only), per the same **assembled §4 lemma
(R13)** (KK25 cited not re-proved, ρ=1/2 N1-conditional; the continuum value stays
VERIFY, asymptotic reference only; the Kambiré-type
counterexample *mechanism* is a characteristic-zero cyclotomic invariant —
`n2-verdict.md`). E.g. ρ=1/2: `[0.293, 0.438]` @128-bit (R13); at 256-bit no
threshold-established ceiling below `R_cap` — i.e., no Kambiré-type near-capacity,
constant-η ceiling (R13; the generic CS25/Elias ceiling at ≈ `R_cap` stands at every
field size; mechanism-asymptotic `[0.293, 0.477]` retained
for reference only) — identical to prime. There is **no**
separate "extension OPEN" window; extensions buy no advantage and are no more
vulnerable. The construction is sub-threshold at 256-bit and Conj-2-consistent for
both field types.

**Important caveats** (from `kambire.md`, `bchks.md`, `n2-verdict.md`):
- Kambiré states `ρ ∈ (0, 1/2)` *open*; exactly `ρ = 1/2` is not covered by N2
  as written (covered instead by the unconditional M31 BCHKS witness N3, also
  prime-only).
- BCHKS Thm 1.13's general prime statement is **conditional on Conj 1.12**; only
  M31 / (M31)⁴ are unconditional.
- N1 (Crites-Stewart) needs `q ≥ n`; its ε_ca=1 strip thins as `q ≫ n`.

---

## 6. δ*_C tables — the headline numeric deliverable (summary)

Full machine-generated tables: **`calculator/out/delta_star_tables.md`**
(`python3 calculator/cli.py delta-star`). For each (ρ, field, n):
`δ_safe` = largest δ where the best **verified** positive bound certifies
`ε_mca ≤ 2⁻¹²⁸` (capped at Johnson; **proven** floor); `t@safe` = query count of the
full composite at `δ_safe` — NOTE this is NOT "queries for 128 bits via the query
term": at `δ_safe` the certified `ε_mca` nearly saturates the `2⁻¹²⁸` budget, so the
query term must cover only the residual slack `2⁻¹²⁸ − ε_mca − |Λ|/|F| ≈ 2⁻¹³⁵·⁸`,
inflating `t@safe = 272` at prime256/`δ_safe = 0.29251`, where the query term alone
needs `t = 257` (`−log₂(1−δ) = 0.4992`, `128/0.4992 ≈ 256.4` — see the §7 table); `δ_unsafe` = the **finite-field inferred ceiling**, per
the R13 correction `(1−ρ) − 2/s_max(b)` with `s_max(b)` the largest power-of-two `s`
passing the KK25 distinctness calibration `p > φ(s)^{φ(s)}` (`s_max = 16` at b=31 —
anchored at **BabyBear** (`p = 15·2²⁷+1`, 2-adicity 27), where the lemma applies
verbatim; **literal M31 caveat (R29/M4):** `v₂(p−1) = 1` gives no smooth domain of
order ≥ 4, so the assembled lemma is vacuous over literal M31 — its 0.375 transfers
via the deployed circle/`F_{p²}` smooth subgroup of order `p+1 = 2³¹` at N2's
ESTABLISHED-MODULO-FORMALIZATION grade —
`16` at b=64, `32` at b=128; at b=256 the mechanism is **infeasible at the prize
threshold** — count `3³² = 2⁵¹ ≪ 2¹²⁸` — so there is **no threshold-established
256-bit ceiling below `R_cap`** — i.e., no Kambiré-type near-capacity, constant-η
ceiling; the generic CS25/Elias ceiling at ≈ `R_cap` stands at every field
size) — the smallest δ at which the finite-field Kambiré instantiation
gives `ε_mca > 2⁻¹²⁸` (per the **assembled §4 lemma, R13** — soft spots: KK25 cited
not re-proved, ρ=1/2 N1-conditional; the continuum `6/log₂|F|` reading stays VERIFY,
asymptotic reference only — the proven negative point is
the M31 δ=1/2 witness). Its *mechanism* is field-agnostic: the same value for prime
and genuine odd-char extension (`n2-verdict.md`). (The generated `δ_unsafe`
column in `calculator/out/delta_star_tables.md` is regenerated per R13
(2026-06-10): `bounds.kambire_unsafe_delta` now enforces power-of-two `s` and
returns the `(1−ρ)−2/s_max(b)` values — `0.375/0.375/0.4375` at 31/64/128-bit,
ρ=1/2 — with `None`/`[δ_safe, —)` at 256-bit, matching this table.)

**The honest punchline — `ε_mca ≤ 2⁻¹²⁸` from the single-code bound is a
FIELD-SIZE constraint.** The verified Bordage-Chiesa error has an `n²·d`
numerator costing ≈ `2·log₂n + log₂(ρn)` bits (≈ 60 bits at `n = 2²⁰`, ≈ 90 at
`2³⁰`). Certifying `ε_mca ≤ 2⁻¹²⁸` therefore needs roughly **`log₂|F| ≳ 256`**
for `n` up to `2³⁰`:

| ρ = 1/2, n = 2²⁰ | δ_safe (verified positive) | t@safe | bits | δ_unsafe (prime) | open gap |
|---|---|---|---|---|---|
| **Mersenne31** (31b, prime) | **infeasible** | — | — | **0.375** (R13; **improves** old 0.30645; via the deployed circle/`F_{p²}` subgroup of order `p+1 = 2³¹`, N2-grade ESTABLISHED-MODULO-FORMALIZATION — literal M31 has `v₂(p−1) = 1`, no smooth domain of order ≥ 4, lemma vacuous there; R29/M4) | [0.000, 0.375] |
| BabyBear (31b, prime) | infeasible | — | — | 0.375 (R13; improves old 0.30587; the direct 31-bit anchor — lemma applies verbatim, R29/M4) | [0.000, 0.375] |
| Goldilocks (64b, prime) | infeasible | — | — | 0.375 (R13; old 0.40625 was too optimistic) | [0.000, 0.375] |
| prime128 (128b, prime) | infeasible | — | — | 0.4375 (R13; old 0.45313 was too optimistic) | [0.000, 0.438] |
| **prime256** (256b, prime) | **0.29251** (≈ Johnson) | **272** (composite; query term alone: 257) | 128.0 | **no threshold-established Kambiré-type (near-capacity, constant-η) ceiling below `R_cap` (R13)** — the generic CS25/Elias ceiling at ≈ `R_cap` stands; mechanism-asymptotic 0.46875 (`(1−ρ)−2/64`, R13 power-of-two; pre-R13 continuum 0.47656) reference only | [0.293, —) |
| **ext256** (256b, genuine ext) | **0.29251** | **272** (composite; query term alone: 257) | 128.0 | same as prime (R13: no Kambiré-type below-`R_cap` ceiling threshold-established; 0.46875 reference only, continuum 0.47656) | [0.293, —) |

(**R13 headline:** the 31-bit `δ_unsafe` values **rose** from ≈0.306 to 0.375 — the
old continuum numbers were too pessimistic, so the small-field window is *less*
closed than previously stated — while the 64/128-bit values *fell* (old values too
optimistic) and 256-bit lost its established below-`R_cap` ceiling entirely — i.e.,
no Kambiré-type near-capacity, constant-η ceiling below `R_cap` remains established
at 256-bit, while the generic CS25/Elias ceiling at ≈ `R_cap` stands at every field
size; see the VERIFY register.)

Mersenne31 / BabyBear / Goldilocks **cannot** certify `ε_mca ≤ 2⁻¹²⁸` from the
single-code bound at all. Deployments over those fields reach 128-bit soundness
via the **`(1−δ)^t` query term** plus protocol repetition / proof-of-work and by
running the soundness argument in a **large extension** — not from a single-code
`2⁻¹²⁸` MCA certificate. (BCHKS Thm 1.5's linear-`n` bound would relax the
field-size requirement *if* its hidden `C_ρ` were pinned — currently `VERIFY`.)

### The ρ = 1/2 headline (the 64 → 128-bit jump)

At ρ = 1/2 the proven positive ceiling is Johnson `δ = 1−√(1/2) ≈ 0.29289`. The
single-query catch probability there is `1−δ = √(1/2)`, so the query-phase factor
is `(1−δ)^t = (1/√2)^t`, and

```
t = 128 queries  →  (1/√2)¹²⁸ = 2⁻⁶⁴      ← only 64 bits at the Johnson radius.
```

Reaching the Singleton capacity `δ → 1/2` would give `(1/2)^t`, i.e.
`(1/2)¹²⁸ = 2⁻¹²⁸` — the full 128 bits. **Closing Johnson → capacity at ρ = 1/2
is the headline 64 → 128-bit jump** that the prize is about. (The query-term
count is best-case; over small fields the `ε_mca + |Λ|/|F|` floor blocks it
entirely — see the per-ρ tables.)

| field | n | δ_safe | (1−δ_safe) | bits/query | t for 128 bits (query term only) |
|---|---|---|---|---|---|
| Mersenne31 | 2²⁰ | 0.29289 (Johnson; δ_safe infeasible) | 0.70711 | 0.5000 | 257 |
| prime256 | 2²⁰ | 0.29251 | 0.70749 | 0.4992 | 257 |
| prime256 | 2³⁰ | 0.28537 | 0.71463 | 0.4847 | 265 |

---

## 7. Which theorem moves SNARK parameters (the exit-criterion one-pager)

**The single thing that moves SNARK parameters: a positive MCA theorem for plain
smooth-domain RS *above* the Johnson radius at ρ = 1/2.** Everything below is
why, and what each candidate result would (or would not) buy.

1. **The lever is the query term, and the Johnson radius pins it at 64 bits.**
   The deployed soundness is `max(ε_mca + |Λ|/|F|, (1−δ)^t)`. With δ capped at
   Johnson (the proven positive ceiling), the query factor at ρ=1/2 is
   `(1/√2)^t`, so **128 queries buy only 2⁻⁶⁴**. To get 128 bits from t = 128
   you must operate at δ → capacity = 1/2, where the factor is `(1/2)^t`. Any
   theorem that does **not** raise the *provable-safe δ above Johnson at ρ=1/2*
   leaves the query count (hence proof size) where it is.

2. **What a positive theorem must beat — the BCHKS Thm 1.9 barrier [N4].**
   Improving MCA beyond Johnson for *any* RS code provably **requires** first
   improving the RS list-decoding radius (list size q) beyond Johnson — a
   separately hard, in-general-false-for-some-RS, long-open problem. So the
   "easy" route is blocked; this is why the prize poses BOTH the MCA and the
   list-decoding challenges, and why the project's core bet is **line-decoding**:
   show that many close codewords on a line for a *smooth* domain must align as a
   line of codewords (⇒ MCA with the right error). The smooth structure is the
   only handle not used by the general (barrier-bound) RS argument.

3. **What the negatives already foreclose — and what they don't.**
   - *Prime fields:* the top of the band is closed — Kambiré [N2] /
     Crites-Stewart [N1] / BCHKS [N3] prove `δ*_C < (1−ρ) − Θ(1/log n)`. The
     naive "up-to-capacity" target is **dead** for prime smooth RS. At small
     fields the window narrows (M31, ρ=1/2: open gap [0.293, 0.375] per the R13
     power-of-two correction; the older "almost closed" [0.293, 0.306] figure
     was too pessimistic).
   - *Genuine odd-characteristic extension fields:* the top is closed at the
     **same** band as prime — the Kambiré-type counterexample is **field-agnostic**
     (characteristic-zero cyclotomic invariant; `n2-verdict.md`), so the
     extension-field analogue of Kambiré is established (`δ*_C < (1−ρ) − Θ(1/log n)`),
     not open. Extensions buy no advantage; there is no separate extension opening.
     (The construction is sub-threshold at 256-bit and does not refute CGHLL Conj 2,
     for both field types.) The genuine open frontier is a *positive* MCA theorem
     above Johnson (sub-lemma P′). Note the asymmetry: the **negative endpoint is
     field-agnostic**, but the positive **P′** route should first be proved for
     **prime-field multiplicative subgroups**, with the extension-field analogue
     stated separately (open) — P′ is *not* itself field-agnostic (`p-prime-route.md`
     §9(e)).
   - *ρ = 1/2 exactly:* not covered by Kambiré as written; the unconditional
     obstruction there is the M31 BCHKS witness (loss ≈ 0.016). (The N1 ρ=1/2
     extension is field-agnostic, like the rest of the counterexample.)

4. **Wrong-code-class positives do not move plain-RS parameters.** Goyal-
   Guruswami / Gao-Cai / Yuan-Zhu reach near capacity for folded / subspace-
   design / random-punctured RS and random linear codes. If a deployment is
   willing to **change code family** to folded / subspace-design RS, those give
   capacity-radius MCA *today* — but at an alphabet-size / folding-overhead cost
   that Workstream F must check survives. For *plain* smooth RS they buy nothing.

5. **Decision-grade summary for designers.**
   - *Today, provably:* operate at δ ≤ Johnson; budget ~256-bit fields (or a
     large extension) for a single-code `2⁻¹²⁸` MCA certificate; expect ~64-bit
     query-term security per 128 queries at ρ=1/2 (→ need ~256 queries, or
     repetition / PoW, for 128-bit).
   - *The one theorem that changes this:* a proof that smooth-domain RS at
     ρ=1/2 has `ε_mca ≤ 2⁻¹²⁸` for some `δ > 0.293` (ideally → 0.5). It would cut
     the query count toward 128 and shrink FRI/STIR/WHIR proofs at the headline
     rate. Absent it, the bracket `[0.293, 0.4375]` at ρ=1/2, 128-bit (R13
     power-of-two ceiling) — proven floor `0.293`,
     per-field ceiling `0.4375` (per the assembled §4 lemma, R13; KK25 cited not
     re-proved, ρ=1/2 N1-conditional), with the negative *mechanism* field-agnostic, the same for prime
     and genuine extension (cyclotomic invariant; `n2-verdict.md`) — is exactly the
     prize.

---

## Provenance / VERIFY register

| Item | Source note | Status |
|---|---|---|
| UD error `n/|F|` | `crites-stewart.md` (BCIKS20) | verified |
| Bordage-Chiesa Thm 9.2 MCA `(m+½)⁷n²d/(3ρ^{3/2}|F|)`, window `1−(1+1/2m)√ρ` | `bordage-chiesa.md` | verified |
| Bordage-Chiesa Lem 10.1 interleaving `m·ε_mca` | `bordage-chiesa.md` | verified |
| Prime no-go (Kambiré Thm1 / CS Cor1 / BCHKS Thm1.13) | `kambire.md`, `crites-stewart.md`, `bchks.md` | verified |
| Kambiré finite-field gap `6/log₂|F|` (the headline `δ_unsafe` — a **finite-field inferred ceiling**) | `kambire.md` "INFERRED" (8^s ≤ |F|); assembled lemma (R13) stated §4 / `technical-note.md` §3.2 | **INFERRED, not proven — and now PER-FIELD RESOLVED at threshold (2026-06-10, from the existing `n2_prize_threshold` computation):** at `ε* = 2⁻¹²⁸` the prize condition is `count > 2^{b−128}`. (i) `b = 128`: any nontrivial witness passes (`2^{b−128} = 1`), the largest distinctness-feasible `s = b/3` gives exactly `η = 6/b` — the table value is **supported**, upgradeable to a theorem modulo the distinctness calibration (`8^s ≤ \|F\|` generic floor), `n \| p−1` smoothness, and the deg-convention (`N0_fire` strict vs `N0_sum` Kambiré-literal). (ii) `b = 256`: **unsupported** — the count needs `s ≳ 162` while the generic distinctness floor caps `s ≤ b/3 ≈ 85`; computed feasibility boundary `b* ≈ 140` bits (`results/n2_prize_threshold.json`). The 256-bit `δ_unsafe` column must be read as mechanism-asymptotic, NOT threshold-established; rescue path = per-prime distinctness certification (`p ∤` explicit subset-sum differences), which is now assessed as BLOCKED with current tools: enumeration at the required `s ≈ 162` means distinguishing `~C(81,40) ≈ 2^{77}` subset sums (computationally infeasible), and the height/norm route fails for a fixed prime (`\|Norm(Δ)\| ≤ (2s)^{φ(s)} ≈ 2^{676} > 2^{256}`, so individual divisibility cannot be excluded, and no union bound over the `~2^{154}` difference pairs is available for fixed `p`). The honest 256-bit status for a negative ceiling below `R_cap` (i.e., a Kambiré-type near-capacity, constant-η ceiling — the generic CS25/Elias ceiling at ≈ `R_cap` stands at every field size): OPEN, needs a new idea. Concrete lead (recorded 2026-06-10): full distinctness is not needed — only `count_p ≥ 2^{128} = p^{1/2}` — and the values are `r`-fold subset sums of the multiplicative subgroup `μ_s ⊂ F_p^*` (`s ~ 324 = p^{0.033}` would give char-0 count `2^{~257}` with margin; `r = ρs+2 ≈ 83` pinned by the rate). Whether `~83`-fold subset sums of a `p^{0.033}`-subgroup have image `≥ p^{1/2}` is a Bourgain–Glibichuk–Konyagin-type growth question — now SURVEYED (`literature/notes/subgroup-sumset-growth.md`, 2026-06-10): **OPEN and not close** — best known bounds give `~2^{15}` vs the `2^{128}` target (≈113-bit gap; covering known only at `~2^{63}` summands, Cipra–Cochrane; the `s ≍ log p` regime named open even for qualitative cancellation, Kowalski 2401.04756). Scoping note: this verdict is specific to the tiny-subgroup regime of the 256-bit ceiling; the P′ keystone's own regime uses much larger subgroups (`\|H\| = n = p^{0.12}-p^{0.65}` deployed), where GK iteration genuinely reaches `√p` with hundreds of summands — the negative verdict must not be over-propagated to P′. (iii) `ρ = 1/2` WIRED via the proven N1 antipodal count `3^{s/2}` (binomial theorem). (iv) **FORMAL LEMMA ASSEMBLED + s-integrality CORRECTED (R13):** `s` must be a power of two (`s = n/m`, `n = 2^t`), so the round-6b value `s=42` was invalid; honest per-field ceilings use the largest power-of-two `s` passing distinctness. Two calibrations: conservative Linnik window (`p ≥ 8^s`) and the literal KK25 distinctness hypothesis (`p > φ(s)^{φ(s)}`, `φ(s) = s/2`), the latter weaker at small `b`. Resulting per-field theorem values (`δ_unsafe = (1−ρ) − 2/s_max`): `b=31`: `s_max=16` (φ-calibration), `δ_unsafe(1/2) = 0.375` (table 0.3065 was too pessimistic); `b=64`: `s_max=16`, `δ_unsafe(1/2) = 0.375` (table 0.4062 too optimistic); `b=128`: `s_max=32`, `δ_unsafe(1/2) = 0.4375` (table 0.4531 and round-6b 0.4524 both too optimistic); `b=256`: distinctness allows `s=64` but count `3^{32} = 2^{51} ≪ 2^{128}` — INFEASIBLE at threshold (unchanged). Status: lemma assembled from cited components (Kambiré Thm 1 count + KK25 distinctness + exact threshold arithmetic); soft spots flagged: literal KK25 lemma cited not re-proved; ρ=1/2 carries the N1 conditionality; deg-convention (`N0_fire` vs `N0_sum`) tracked per row. The headline tables' continuum `6/b` column should be replaced by the power-of-two `2/s_max(b)` values — DONE (2026-06-10): `calculator/out/delta_star_tables.md` regenerated from the R13-corrected `bounds.kambire_unsafe_delta`. |
| BCHKS Thm 1.5 `C_ρ·n/(η⁵|F|)`, hidden `C_ρ` | `bchks.md` | **VERIFY** `C_ρ` (placeholder 1; η-exp 5 firm) |
| Extension-field capacity window (= same band as prime) | `n2-verdict.md` (exact cyclotomic-invariant computation) | **field-agnostic** — counterexample extends to genuine `GF(p^e)`; sub-threshold at 256-bit; not a Conj-2 refutation |
| List size `n·m/η` | Johnson bound; deep analysis separate | **VERIFY** constant + exponents |
| BCHKS Thm 1.13 general prime statement | `bchks.md` | **conditional on Conj 1.12** (M31 unconditional) |
| `ρ = 1/2` boundary in Kambiré | `kambire.md` | **VERIFY** (stated for open `(0,1/2)`) |
| Whether `𝔽_q^*` / construction lands on a strictly power-of-two smooth subgroup | `bchks.md`, `kambire.md` | **VERIFY** (matters for hitting the FFT domain) |
