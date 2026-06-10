# BCHKS — On Proximity Gaps for Reed–Solomon Codes

> Source via `curl` (WebFetch 403-blocked). Text via `pdftotext -layout` from
> `/tmp/eprint_2025_2055.pdf`. Formulas transcribed verbatim (ASCII rendering mine).

## Bibliographic info
- **Authors:** Eli Ben-Sasson, Dan Carmon, Ulrich Haböck (all StarkWare), Swastik Kopparty,
  Shubhangi Saraf (both Univ. of Toronto, Math + CS).
- **Title:** "On Proximity Gaps for Reed–Solomon Codes".
- **Id / date:** IACR eprint 2025/2055, https://eprint.iacr.org/2025/2055. **November 6, 2025.**
- **Cited in ABF survey as `[BCHKS25]`** (its Thm 1.3 = ABF Thm 4.9.2; Thm 4.6 = ABF Thm 4.12;
  Thm 1.9 = ABF Thm 5.2; Cor 1.7 = ABF Thm 4.18; the prime-field negative result feeds ABF Thm 4.16).
- Author list MATCHES the project brief (Haböck = "Habock"). Note ABF also cites Haböck's solo
  work separately as `[Hab25]`.

## ⚠ Notation convention (DIFFERENT from ABF / Crites–Stewart)
- **`δ` here = the code's MINIMUM relative distance** `= 1 - k/n` (= what ABF calls `δmin(C)`),
  NOT the proximity radius.
- **`γ` = the proximity (decoding) radius** (= what ABF/CS call `δ`).
- **`a/q` = the soundness error** (`a` = number of "exceptional"/bad `z` values; `q=|F|`).
- **`ε∗` = the proximity loss** (gap between the per-line radius `γ'` and the conclusion radius).
- `ρ := k/n = 1 - δ` (the paper notes this is "off-by-1/n from the rate"). Johnson radius
  `J(δ) = 1 - sqrt(1-δ)`; double-Johnson `J_2(δ)=1-(1-δ)^{1/4}`; 1.5-Johnson `J_{1.5}(δ)=1-(1-δ)^{1/3}`.

## Definition 1.1 (Proximity gaps).
`C ⊆ F_q^n` has proximity gaps up to radius `γ∈[0,1]` with parameters `a∈N`, `ε∗∈[0,1]` if for
every `f,g∈F_q^n` and every `γ' ∈ [0,γ]`:
```
|{ z ∈ F_q : ∆(f + z·g, C) ≤ γ' }| ≥ a   ⟹   ∆([f, g], C^2) ≤ γ' + ε∗.
```
(`C^2` = the 2-interleaved code; `∆([f,g], C^2)` = correlated-agreement distance.) "soundness
error" `= a/q`; "proximity loss" `= ε∗`. Applications want `γ` a constant in (0,1), `ε∗ → 0`,
sometimes `ε∗ = 0`.

---

## Main results — POSITIVE (exact statements)

### Theorem 1.3 (up to half-minimum-distance; = ABF Thm 4.9 item 2)
`C = RS[F_q, D, k]`, `n=|D|`, `δ = 1 - k/n`. For `γ ∈ [δ/3, δ/2 - 1/n]`. If
`S = {z : ∆(u_0 + z·u_1, C) ≤ γ}` has size `a ≥ (δ/γ - 1)·(1/(δ-2γ))`, then
`∆([u_0,u_1], C^2) ≤ (1 + 1/(a-1))·γ`. Equivalently, **for proximity loss `ε∗` it suffices to take**
```
a ≥ max( (δ/γ - 1)·1/(δ-2γ) ,  1 + γ/ε∗ ).
```
Coarser (Table 1): `a ≥ 1/(δ/2 - γ)`. **Headline: if `ε∗` and `δ/2-γ` are positive constants,
then `a = O(1)`** — saving a Θ(n) factor in soundness error over [BCIKS20] (which needed `a>n`).

### Corollary 1.4 (lossless, ε∗=0, up to half-distance).
For `δ ≥ 2/(3√n)` and `γ ∈ [δ/3, δ/2 - 3/(δn)]`: if `|S| = a > γ·n + 1` then
`∆([u_0,u_1], C^2) ≤ γ` (0 proximity loss). Optimal on its range (matches [AHIV17, RZ18]).

### Theorem 1.5 (up to the Johnson radius, ε∗=0; the KEY positive RS result; = ABF Thm 4.12 core)
`C = RS[F_q, D, k]`, `n=|D|`, `δ=1-k/n`, `ρ = k/n = 1-δ`. For `γ ∈ [0, 1 - sqrt(ρ))`, let
`η = 1 - sqrt(ρ) - γ`, `m = max(⌈sqrt(ρ)/(2η)⌉, 3)`. If `S = {z : ∆(u_0+z·u_1, C) ≤ γ}` has size
```
a >  ( 2(m+1/2)^5 + 3(m+1/2)γρ ) / (3 ρ^{3/2}) · n  +  (m+1/2)/sqrt(ρ)   =  O_ρ( n / η^5 ),   (eq.1)
```
then `∆([u_0,u_1], C^2) ≤ γ` (**ε∗ = 0**). For `η ≪ sqrt(1-δ)`, RHS asymptotics `O_δ(n/η^5)` with
leading constant `1 / (48(1-δ)^{3/2})`. **Improves [BCIKS20 Thm 5.1] by more than a factor n**
(was `O_δ(n^2/η^7)` — see Table 1). This is the bound that gives MCA error `O_ρ(n/(η^5|F|))` at
`γ = J(δ) - η` once converted to MCA form (via collinearity/Lemma 4.6).

### Table 1 (positive landscape, RS of distance δ; shaded = arbitrary linear codes)
| Ref | radius γ | a | ε∗ |
|---|---|---|---|
| [RVW13] | arbitrary | a≥2 | γ |
| [BKS18] | arbitrary | a≥2 | `(a/(a-1))(γ-γ^2)`; →`J^{-1}(γ)-γ` as a→∞ |
| trivial | <δ | a≥2^{Ω_δ(n)} | 0 |
| [AHIV17] | <δ/4 | a≥2 | `γ/(a-1)`; =0 for a>γn+1 |
| [RZ18,BKS18] | <δ/3 | a≥2 | `γ/(a-1)`; =0 for a>γn+1 |
| [BCIKS20] | <δ/2 | a>n | 0 |
| **This work** | <δ/2 | `a≥1/(δ/2-γ)` | `γ/(a-1)`; =0 for a>γn+1, `2δ-γ≥Ω(1/√n)` |
| [BKS18] | `<1-(1-δ)^{1/4}-η` (double-Johnson) | `Θ_δ(1/η^2)` | `O_δ(1/(η^2 a))·γ`; =0 for `a=Ω_δ(n/η^2)` |
| [BGKS20] | `<1-(1-δ)^{1/3}-η` (1.5-Johnson) | `Θ_δ(1/η)` | `O_δ(1/(ηa))·γ`; =0 for `a=Ω_δ(n/η)` |
| [BCIKS20] | `<1-(1-δ)^{1/2}-η` (Johnson) | `Θ_δ(n/η^7)` | 0 |
| **This work** | `<1-(1-δ)^{1/2}-η` (Johnson) | `Θ_δ(n/η^5)` | 0 |

---

## Main results — NEGATIVE (exact statements)

### Theorem 1.6 (refutes the n^τ-bounded proximity-gaps conjecture for ALL τ; CHAR 2)
Fix integer `τ>0`, `λ_τ = 2^{-(τ+2)}`, arbitrary const `ϵ>0`. Take `δ = 1-λ_τ`, `γ = 1-4λ_τ`.
Then for **all `F_q` of characteristic 2** there are RS codes `C = RS[F_q, D, (1-δ)n]`, domain
`D` with `n=|D|=O(q^{1/(τ(1+ϵ))})`, distance δ, and words `f,g` such that
```
|{z ∈ F_q : ∆(f+zg, C) ≤ γ}| ≥ (1-o(1))·q ≥ n^{τ(1-ϵ)},
```
yet `∆([f,g], C^2) ≥ 1 - 2λ_τ = (2/3)δ + (1/3)γ`. **So at radius `γ = δ - Ω_τ(1)`, soundness
error `a/q` can't beat `n^{-(τ-o(1))}` without large proximity loss** → refutes Conjecture 1.2
(below) for every constant τ. (Construction: D ≈ random F_2-subspace of F_q; f,g = monomials
`X^u, X^v`, u>v>k; subspace-polynomial coefficient structure + 2nd moment method.)

### Corollary 1.7 (tightness of Thm 1.5 at the Johnson radius; = ABF Thm 4.18; CHAR 2)
Arbitrary const `ϵ>0`, `δ = 15/16`, `γ = 1-sqrt(1-δ) = 3/4`. For all `F_q` of char 2 there are
RS codes `C=RS[F_q,D,(1-δ)n]`, `n=|D|=q^{(1+ϵ)/2}`, and `f,g` with
`|{z : ∆(f+zg,C) ≤ γ}| ≥ (1-o(1))q = n^{2(1-ϵ)}` yet `∆([f,g],C^2) ≥ 7/8 = γ + 1/8`. **So `a`
must jump from `O(n)` (Thm 1.5, just below J(δ)) to `Ω(n^{2-ϵ})` exactly at `γ = J(δ)`.**

### Conjecture 1.2 (the conjecture they refute) — n^τ-bounded proximity gaps.
For constant δ∈(0,1): ∃ const τ s.t. for every η>0, `C` has proximity gaps up to `γ = δ - η`
with `ε∗ = o_η(1)` and `a = O_η(n^τ)`. (Stronger form asks `ε∗ = 0`.) **Refuted for all τ=O(1)
by Theorem 1.6** (in char 2). The paper explicitly notes: "Versions of Conjecture 1.2 may still
be true ... even with τ=1, for fields of PRIME cardinality, or for well chosen evaluation
domains D over fields of characteristic 2." — i.e. the deployed prime-field smooth case is NOT
refuted by Thm 1.6.

### Theorem 1.9 (good proximity gaps ⇒ good list-decoding, RS only; = ABF Thm 5.2 source)
`C = RS[F_q, D, k]`, `|D|=n`, `k=(1-δ)n`. Let `γ = LDR_{F_q,D,q}(δ) + 2/n`. Then ∃ `f,g` with
`|{z : ∆(f+z·g, C) ≤ γ}| ≥ q/(2n)` but `∆([f,g], C^2) ≥ δ - 1/n`. **⇒ soundness error `a/q ≥
1/(2n)`, independent of q.** So pushing the proximity-gap radius beyond the list-decoding radius
(for list size q) with `a < q/(2n)` is impossible — **improved list-decoding bounds are a
PREREQUISITE for improved proximity gaps beyond Johnson.** (`LDR_{F_q,D,L}(δ)` = Def 1.8 = largest
γ with all radius-γ Hamming balls containing ≤ L codewords; `LDR_{F_q,n,L}(δ)` = min over D.)

### Theorem 1.13 (PRIME-FIELD / multiplicative-subgroup negative result — THE deployed-case one)
**Conditional on additive-combinatorics Conjecture 1.12.** Suppose `(q,a,b)` is *admissible*
(Def 1.11: ∃ mult. subgroup `G⊆F_q^*` with `|G|=b` and the `ℓ=b/2`-wise distinct-element sumset
`|G^{(+ℓ)}| ≥ a`), with `b` even. Let `H⊇G` be any mult. subgroup, `C=RS[F_q, D, k]` with `D=H`,
`n=|D|`, `k=(1/2 - 2/b)n`, relative distance `δ = 1/2 + 2/b`. Then ∃ `f,g` with
```
|{z ∈ F_q : ∆(f+z·g, C) ≤ δ - 2/b}| ≥ a,    yet  ∆([f,g], C^2) ≥ δ - 1/b.
```
- **Conjecture 1.12:** for infinitely many primes q, ∃ `b ≤ 10 log q` with `(q, q/10, b)` admissible.
  (Weaker than infinitude of Mersenne primes; for Mersenne `q=2^p-1`, `G=⟨-2⟩` makes `(q,q,2log_2(q+1))`
  admissible.)
- **Consequence (verbatim):** the conjecture implies that for RS over **prime fields**, for
  infinitely many n, **when `γ = δ - Θ(1/log n)` and `a ≥ q/10 ≥ n/10`, there must be a proximity
  loss of `Θ(1/log n)`.** (Table 2 row: `δ≈1/2, γ=δ-1/log n, a=n=q, ε∗ = 1/(2log_2 n)`, Mersenne q.)
- **Explicit instantiations (UNCONDITIONAL for these q):**
  - `q = M31 = 2^31-1`, `C = RS[F_q, F_q^*, (q-1)/2 + 2/log_2 q]`: ∃ `f,g` with `∆(f+z·g, C) ≤ 1/2`
    for ALL `z∈F_q`, yet `∆([f,g], C^2) ≥ 1/2 + 1/62 ≈ 0.516`.
  - `q = (M31)^4 ≈ 2^124`: `(q,q,b)` admissible with `b=8·31 ≈ 2log_2 q`; RS code of any length n
    (with `b | n`, `n | q-1`), distance `δ = 1/2 + 2/b ≈ 0.508`, proximity gaps at
    `γ = δ - 2/b = 1/2 ≈ δ - 0.008` have proximity loss `ε∗ = 1/b ≈ 0.004` with `a = q`.

### Theorems 1.15, 1.16 (phase transitions, δ = o(1) regime; prime fields)
- **Thm 1.15 (transition at δ/3):** ∃ infinitely many q, `C=RS[F_q,F_q,k]`, `n=q`, `k=n-c`
  (so `δ=c/n`), with `|{z : ∆(f+zg,C) ≤ δ/3}| ≥ (q-1)/c = O(1/δ)`, yet `∆([f,g],C) ≥ 2δ/3`.
  (Table 2: prime q, any const c.)
- **Thm 1.16 (transition at δ/2):** fix `0<ϵ<1/4`; ∃ infinitely many primes q,
  `C=RS[F_q,D,k]`, `n=|D|=O(q^{0.5+2ϵ})`, `k=n-Θ(n^ϵ)` (`δ=Θ(n^{-(1-ϵ)})`), with
  `|{z : ∆(f+zg,C) ≤ δ/2}| ≥ q-1 = Ω(n^{2-8ϵ})`, yet `∆([f,g],C^2) ≥ 3δ/4`. (Table 2: ε∗=δ/4.)
- Both are `δ=o(1)` ONLY (can't hold for `δ=Ω(1)`, since positive J/J_{1.5} results cover that).

### Theorem 1.17 (limit on STARK / DEEP-ALI soundness)
For the DEEP-ALI IOP of the (unsatisfiable) CYCLE-SUM CSP with `C=RS[F_q,D,k]`, `D` = union of t
cosets of an order-`a` subgroup `G`, `|D|=a·t=n`, `k=a=n/t`, `δ=1-1/t`: ∃ a prover strategy that
isn't rejected and produces `(h_1,…,h_c)` with `Pr[∆([h_1,…,h_c], C^c) ≤ (1+γ_q)/2] ≥ Ω(1/n)`,
where `γ_q = LDR_{F_q,D,q}(δ) + 1/n`. ⇒ the [BGKS20] STARK soundness-error bound `O(n/q)` cannot
be improved much; the cheating probability is governed by RS list-decodability at list size q.
(E.g. if `δ=0.84` and `γ_q=J(δ)`, then `(1+γ_q)/2 = 0.8` ≫ trivial.)

---

## Regime of validity
- **Code class:** `C = RS[F_q, D, k]`, distance `δ = 1-k/n`.
- **Positive (Thm 1.3, 1.5, Cor 1.4):** ALL RS codes (any domain D, including smooth multiplicative
  subgroups), `δ = Ω(1)` for the headline asymptotics. Valid up to the **Johnson radius**
  `γ < 1 - sqrt(1-δ) = 1 - sqrt(ρ)`, ε∗=0, `a = O_δ(n/η^5)`.
- **Negative Thm 1.6 & Cor 1.7:** **characteristic-2 fields ONLY**, special random-subspace domains.
  Do NOT apply to prime fields or to deployed smooth domains.
- **Negative Thm 1.13:** **prime fields, D = multiplicative subgroup (i.e. SMOOTH-type domain)** —
  the directly relevant deployed case — but **CONDITIONAL on Conjecture 1.12** (additive
  combinatorics); only the M31 and (M31)^4 instantiations are unconditional, and those have
  `δ ≈ 1/2` (rate ≈ 1/2) and only push `γ` to `≈ δ - O(1/log q)` with small (not zero) loss.
- **Negative Thm 1.9:** ALL RS codes; ties proximity-gap radius to the list-decoding radius.
- **Negative Thm 1.15/1.16:** prime fields, but `δ = o(1)` (irrelevant to our constant-rate targets).

---

## Relevance to our targets (ρ∈{1/2,1/4,1/8,1/16}, smooth domain, ε*=2^-128)
- **Positive: confirms `δ_known_positive` for smooth RS at the Johnson radius** with the BEST
  current error constant: at `γ = 1-sqrt(ρ)-η`, `a = O_δ(n/η^5)`, so MCA/CA soundness error
  `a/q = O_ρ(n/(η^5 q))` and `ε∗=0`. To hit `2^-128` you need `q ⪆ 2^128·n/η^5`. For ρ=1/2 the
  Johnson radius is `1-sqrt(1/2) ≈ 0.293`; ρ=1/4 → 0.5; ρ=1/8 → 0.646; ρ=1/16 → 0.75.
- **Negative (deployed smooth prime case):** Thm 1.13 (cond. on Conj 1.12) says you canNOT reach
  within `Θ(1/log n)` of `δ = δmin` with small proximity loss over prime-field multiplicative-
  subgroup domains. At ρ=1/2 (δ≈1/2), the M31 instantiation **unconditionally** exhibits `f,g`
  that are `1/2`-close for ALL z yet have CA-distance `≥ 1/2 + 1/62 ≈ 0.516` — i.e. a genuine
  proximity-gap failure (loss `≈ 0.016`) at radius `γ = 1/2` for a rate-≈1/2 smooth-domain RS code
  over a deployed prime. This is concrete evidence that for ρ=1/2 the achievable `δ*_C` sits
  well below `δmin`, near the Johnson radius, NOT near capacity.
- **Thm 1.9 reframes the prize:** improving MCA/CA beyond Johnson for smooth RS REQUIRES first
  improving the RS list-decoding radius (for list size q) beyond Johnson — which is itself a
  hard, long-open problem (and false for some RS codes [JH01, BSKR06]). This directly motivates
  the prize's "grand list-decoding challenge" as a prerequisite.
- The `ε∗ = 0` of Thm 1.5 matters because ABF's `ε_mca` is the loss-free notion (Remark 4.4).

---

## Placement of δ_known_positive / δ_known_negative (smooth multiplicative-subgroup RS)
- **δ_known_positive = `1 - sqrt(ρ) - η`** (Johnson radius − slack), `a=O_δ(n/η^5)`, ε∗=0 (Thm 1.5).
  Best known constant for the deployed code.
- **δ_known_negative:**
  - char-2 smooth-ish: `γ = J(δ)` already needs `a ≥ n^{2-ϵ}` (Cor 1.7), and `γ = δ-Ω_τ(1)` needs
    `a ≥ n^{τ-o(1)}` (Thm 1.6) — but char 2 only.
  - **prime + multiplicative-subgroup (deployed):** `γ = δ - Θ(1/log n)` forces `ε∗ = Θ(1/log n)`
    (Thm 1.13, cond. Conj 1.12); UNCONDITIONALLY for M31, ρ≈1/2: at `γ=1/2` there is a CA failure.
  - all RS: `γ > LDR_{F_q,D,q}(δ)` forces `a/q ≥ 1/(2n)` (Thm 1.9).
- **Open strip for prime smooth RS:** `[ 1-sqrt(ρ) (positive), δmin - Θ(1/log n) (negative, cond.) ]`.

---

## Open questions the paper states
1. **Does Conjecture 1.2 hold for PRIME fields, even with τ=1?** (Thm 1.6 only refutes char 2.)
   Explicitly flagged as "a very exciting direction for future research."
2. Do char-2 analogues of Thm 1.6 hold for **all constant-characteristic** fields? (not checked).
3. **Prove additive-combinatorics Conjecture 1.12** (or its admissibility statements for specific
   practical q) — would make the prime-field negative result (Thm 1.13) unconditional for inf. many n.
4. **Determine `LDR_{F_q,n,poly(n)}(δ)`** — the RS list-decoding radius for polynomial list size — a
   "very basic and well-studied" open problem that (via Thm 1.9) gates all proximity-gap progress
   beyond Johnson.
5. Extend Thm 1.9 (proximity-gaps ⇒ list-decoding) to **general linear codes** (only known for RS).

---

## INFERRED — VERIFY
- **Notation conflict is the big one:** BCHKS `δ` = min-distance (`δmin`), `γ` = radius, `a/q` =
  soundness error, `ε∗` = proximity loss. ABF/CS use `δ` = radius and `ε_mca/ε_ca` = soundness
  error (≈ `a/q` here), with NO proximity loss in ABF's `ε_mca`. When mapping BCHKS Thm 1.5 to
  ABF Thm 4.12: ABF's `δ` (radius) = BCHKS `γ`; ABF's `ρ` = BCHKS `ρ=1-δ`; ABF error `n/(η^5|F|)` =
  BCHKS `a/q` with `a=O(n/η^5)`. I have cross-checked this map and it is consistent, but the exact
  `(m+1/2)` constant placement should be eyeballed once against the source.
- The M31 numbers (`1/2 + 1/62`, `≈0.516`) and (M31)^4 numbers (`δ≈0.508`, `ε∗≈0.004`, `γ=1/2`)
  are quoted verbatim from §1.4.3; I have not re-derived them.
- Thm 1.13 is CONDITIONAL on Conj 1.12 for the "infinitely many n" statement; the M31 / (M31)^4
  instantiations are presented as unconditional (admissibility checkable). Treat the general
  prime-field negative result as conjectural until Conj 1.12 (or specific admissibility) is proven.
- Whether the M31 instantiation's domain `F_q^*` counts as "smooth" per ABF Def 2.12 (mult. coset
  of a power-of-two-order subgroup): `F_q^*` has order `q-1 = 2^31-2 = 2·(2^30-1)`, which is NOT a
  power of two, so `F_q^*` itself is NOT smooth in ABF's strict sense. Thm 1.13 allows `D=H` any
  mult. subgroup `⊇ G` with `b|n, n|(q-1)`; a power-of-two-order `H` (smooth) requires `b | 2^r`
  and `2^r | q-1`. **VERIFY whether the negative construction can be instantiated on a strictly
  power-of-two-order (smooth) subgroup** — this matters a lot for whether the deployed FFT domain
  is directly hit. (The Θ(1/log n) prime-field consequence is stated for multiplicative subgroups
  generally; the strict power-of-two smoothness constraint may or may not be compatible with the
  admissibility requirement.)
