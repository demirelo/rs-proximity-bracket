# Crites–Stewart — On Reed–Solomon Proximity Gaps Conjectures

> Source fetched via `curl` (WebFetch 403-blocked). Text via `pdftotext -layout` from
> `/tmp/eprint_2025_2046.pdf`. Formulas transcribed verbatim (ASCII rendering mine).

## Bibliographic info
- **Authors:** Elizabeth Crites, Alistair Stewart (Web3 Foundation).
- **Title:** "On Reed–Solomon Proximity Gaps Conjectures".
- **Id / date:** IACR eprint 2025/2046, https://eprint.iacr.org/2025/2046. **December 19, 2025.**
- **Cited in ABF survey as `[CS25]`** (its Cor 1 = ABF Thm 4.17; its Thm 2 = ABF Thm 5.3).
- **One-line:** A NEGATIVE-results paper. Disproves the up-to-capacity CA / MCA / list-decoding
  conjectures for RS, showing failure already at the **list-decoding capacity** `1 - H_q(ρ)`
  (strictly below the MDS capacity `1-ρ`). Second result: CA with small error ⇒ RS list
  decoding (a proximity-gaps ⇒ list-decodability reduction).

---

## Main results — exact theorem statements

### Notation
- `RS(F, D, k)` = evaluations on domain `D` (`|D|=n`) of polys of degree `≤ k-1`; rate `ρ=k/n`.
- `f` = an **absolute** error count (NOT a fraction); relative distance is `δ = f/n`.
- `H_q(x) = x·log_q(q-1) - x·log_q(x) - (1-x)·log_q(1-x)` (q-ary entropy).
- `(δ,L)`-list decodable (**Definition 1**): for every `u∈F_q^n`, `|{v∈C : ∆(u,v) ≤ δ·n}| ≤ L`.
- **Correlated agreement (Definition 2):** `u_0,…,u_ℓ` have CA with C of density `≥1-δ` if
  `∃ D'⊆D`, `v_0,…,v_ℓ∈C` with `|D'|/|D| ≥ 1-δ` and `u_i = v_i` on all of `D'` for all `i`.
  (This is the existence side; the proximity-gap "error" `ε` is the probability that a random
  line point is δ-close WITHOUT global CA — see Thm 1.4 below for the [BCIKS] error form.)

### List-decoding capacity (the bound everything is measured against)
**Theorem 7.4.1 (List-Decoding Capacity) [GRS14, Eli57].** `q≥2`, `0 ≤ δ < 1-1/q`, small `η>0`,
large `n`:
- (i) if `ρ ≤ 1 - H_q(δ) - η`, there exists a `(δ, O(1/η))`-list decodable code;
- (ii) if `ρ ≥ 1 - H_q(δ) + η`, every `(δ,L)`-list decodable code has `L ≥ q^{Ω(ηn)}`.

So **list-decoding capacity = `1 - H_q(δ)`**, which is strictly below the MDS/Singleton
capacity `1-δ` because `H_q(δ) > δ` (Claim 1).

**Claim 1 (the `1/log₂q` gap):** for `δ ≤ 1 - 1/(q-1)`,
```
0 ≤ H_2(δ)/log_2 q  -  δ/((ln 2)(q-1)(log_2 q))  ≤  H_q(δ) - δ  ≤  H_2(δ)/log_2 q  ≤  1/log_2 q.
```
i.e. `H_q(δ) - δ ≤ 1/log_2 q`. **This is the key practical constant: moving from "up to
capacity `1-ρ`" to "up to list-decoding capacity `1-H_q(ρ)`" costs at most `1/log₂q` in
distance.** Fields used in SNARKs have `log₂q ⪆ 31`, so the loss is `≤ 1/31`.

### Theorem 1 — proximity gap fails beyond list-decoding capacity (CORE NEGATIVE result)
> **Theorem 1.** Suppose `f < n-k` and `∆(u^{(1)}, RS(F_q, D, k)) > f`. Then there exists a
> `u^{(0)}` such that there are **at least** `  C(n,f) / ( q^{n-f-k} g(f(n-f)/q) C(n,f) ) `
> [printed as `C(n,f) / ( q^{n-f-k}·g(f(n-f)/q) )` over the binomial — see note] values of
> `λ ∈ F_q` such that `∆(u^{(0)} + λ u^{(1)}, RS(F_q, D, k)) ≤ f`, where
> ```
> g(x) = exp(x)                 when x ≤ 3/2,
>      = exp(2√x) / (√(2π)·⌊√x⌋) when x > 3/2.
> ```

(Underlying **Lemma 1**: for uniformly random `u∈F_q^n`, `Pr[∆(u, RS) ≤ f]` is sandwiched, with
upper bound `≈ q^{H_q(f/n)·n - n + k} = q^{(H_q(δ)-1+ρ)n}`. When `ρ > 1 - H_q(δ)`, i.e. beyond
list-decoding capacity, this probability is large, so a random line point is δ-close with high
probability — defeating the proximity-gap dichotomy.)

### Corollary 1 — `ε_ca = 1` regime (= ABF Theorem 4.17)
> **Corollary 1.** Given `n, q` with `q ≥ 10`, `q ≥ n`, if
> ```
> n(1 - H_q(f/n)) + 2 + sqrt( n·H_q(f/n) - f )  ≤  k  ≤  n - f - 2,
> ```
> then there exist `u^{(0)}, u^{(1)} ∈ F_q^n` with `∆(u^{(1)}, RS) > f` and **for EVERY `λ∈F_q`**,
> `∆(u^{(0)} + λ u^{(1)}, RS(F_q,D,k)) ≤ f`. (Indeed one `u^{(0)}` works for any far `u^{(1)}`.)

Meaning: in this k-window (k just above the list-decoding-capacity dimension `n(1-H_q(f/n))`,
up to near the Singleton dimension), the proximity-gap error is `ε = 1` — every point on the
line is close, yet there is no correlated agreement. Equivalent δ-form (from the abstract /
ABF Thm 4.17 reading): re-expressing via `1 - H_q(δ) + 2/n + sqrt((H_q(δ)-δ)/n) ≤ ρ ≤ 1-δ-2/n`
gives `ε_ca(C,δ)=1`; since `H_q(δ) ≈ δ - 1/log q` for large q, breakdown occurs at
`δ = 1 - ρ - η` with `η ≈ 1/sqrt(n·log q) + 2/n - 1/log q`.
Application note (their §): "when `k/n ≫ log₂q`, a similar result holds with
`q(1 - H_q(k/n)) ⪅ f ≤ n-k-1`" (typo for `n(1-H_q(k/n))`).

### Disproved conjectures (reproduced verbatim by the paper)
- **Conjecture 2.3 [BGKS20] — list-decodability up to capacity:** for every `ρ>0` there is `C_ρ`
  s.t. every RS code of rate ρ is list-decodable from `δ ≤ 1-ρ-η` with list size `(n/η)^{C_ρ}`.
  **FALSE** by Thm 7.4.1(ii) (since `H_q(δ)>δ`).
- **Conjecture 8.4 [BCI⁺23] — CA up to capacity:** ∃ `c_1,c_2` s.t. for all `η>0`, the
  [BCIKS] proximity-gap/CA theorems hold for `δ ≤ 1-ρ-η` with error
  `ε ≤ (1/(ηρ)^{c_1})·(n^{c_2}/q)` (and the curve version with `(ℓn)^{c_2}/q`). **FALSE** in the
  Corollary-1 regime, where `ε` cannot be `<1`; taking `q > (2n)^{c_1+c_2}` breaks it.
  - For reference, the **true** [BCIKS Thm 1.2/1.4] error (which Conj 8.4 tried to extend to
    capacity): unique-decoding `δ∈(0,(1-ρ)/2]` → `ε_U = n/q`; Johnson `δ∈((1-ρ)/2, 1-√ρ)`,
    `η:=1-√ρ-δ` → `ε_J = (k+1)^2 / (2·min{η, √ρ/20}^7 · q) = O( (1/(ηρ)^{O(1)})·(n^2/q) )`.
- **Conjecture 4.12 [ACFY24b / WHIR] — MCA for smooth RS** (Definition 7, proximity generator
  `Gen(ℓ;α)=(1,α,…,α^{ℓ-1})`, `C=RS[F,L,m]`, `ρ=2^m/|L|`):
  - (1) Up to Johnson `B*=√ρ`: `err = (ℓ-1)·2^{2m} / (|F|·(2·min{1-ρ-δ, √ρ/20})^7)`.
  - (2) Up to capacity `B*=ρ`: ∃ `c_1,c_2,c_3` s.t. for `0<δ<1-ρ-η`,
    `err* = (ℓ-1)^{c_2}·d^{c_2} / (η^{c_1}·ρ^{c_1+c_2}·|F|)`. **FALSE** (MCA⇒CA, and CA fails).

### Minimally modified conjectures (the paper's proposed "best standing" replacements)
- **Our Conjecture 1 (list-decoding, prime fields):** Conj 2.3 with `δ ≤ 1-ρ-η` replaced by
  `H_q(δ) ≤ 1-ρ-η`, for prime `q`.
- **Our Conjecture 2 (proximity gap / CA, prime fields):** Conj 8.4 with `δ ≤ 1-ρ-η` replaced by
  `δ ≤ 1 - H_q(δ) - 1/n - η`, for prime `q`. (Their construction's error ≈ `q^{(H_q(δ)-1+ρ)n}`,
  which under `ρ ≤ 1-H_q(δ)-1/n` is `≤ 1/q`.)
- **Our Conjecture 3 (MCA, prime fields):** Conj 4.12 with `0<δ<1-ρ-η` replaced by
  `0 < H_q(δ) < 1 - 1/n - ρ - η`, for prime `q`.

### Theorem 2 — CA implies RS list decoding (= ABF Theorem 5.3) (POSITIVE-direction reduction)
> **Theorem 2.** If `RS(F_q, D, k)` satisfies correlated agreement over lines with `f < n-k-1`
> errors with error parameter `ε < (q-n)/(kq)`, then `RS(F_q, D, k+1)` is `(f/n, L)`-list
> decodable, where
> ```
> L = ⌈ ε q (q-n) / (q - n - k ε q) ⌉.
> ```
> For `ε < (q-n)/(2kq)` this simplifies to `L ≤ 2εq`.

(Mechanism, Claim 3/4: divide pointwise by `(x-a)`; polynomial-remainder theorem maps a list of
`RS(k+1)` codewords near `u` into many line points near `RS(k)`; Schwartz–Zippel: distinct
degree-≤k polys collide at random `a` w.p. `≤ k/(q-n)`.)

**The paper's own caveat (verbatim, important):** "Correlated agreement results with `ε > 1/k`
are useful; indeed, the only known results for small fields and large k, e.g., `q≈2^32, k≈2^20`,
that are used in SNARKs in practice fall into this regime. For these parameters, **Theorem 2
does not demonstrate anything.**" Conversely, any conjecture with `ε = F(n,k,f)/q` implies a list-
decoding conjecture for `q > 2k·F(n,k,f)`. → Strong positive CA (`ε≪1/k`) would imply hard list-
decoding results, so the community should seek **large-error** (M)CA instead.

---

## Regime of validity
- **Code class:** RS over a general domain `D`; the constructions (Thm 1 / Cor 1) need a far word
  `u^{(1)}` (a deep hole, e.g. `x^k` or `1/(x-a)`). **Domain need NOT be smooth** — the results
  hold for arbitrary `D` including smooth multiplicative subgroups; the constructive deep holes
  used (`x^k`, `1/(x-a)`) exist over any domain.
- **Field regime for Corollary 1:** requires `q ≥ 10` and **`q ≥ n`** (field at least as large as
  the domain). This is the regime where the gap `H_q(δ)-δ ≈ 1/log q` is the binding loss. Note
  this `q ≥ n` is satisfied by SNARK params (`q≈2^31`, `n≤2^30`) but is a real constraint.
- **δ regime:** the negative results bite **near capacity**, specifically for `δ = f/n` with
  `k` in the window `[n(1-H_q(δ)) + 2 + sqrt(nH_q(δ)-f), n-f-2]` — i.e. `δ` between the list-
  decoding capacity `1-H_q^{-1}(1-ρ)` and the Singleton capacity `1-ρ`.
- **Theorem 2 regime:** needs `ε < (q-n)/(kq) ≈ 1/k` and `f < n-k-1`. Large-field / small-k
  only; does NOT cover deployed `q≈2^31, k≈2^20`.

---

## Relevance to our targets (ρ∈{1/2,1/4,1/8,1/16}, smooth domain, ε*=2^-128)
- **Negative side:** establishes that the "naive up-to-capacity" target `δ → 1-ρ` is
  **impossible** for plain RS — CA/MCA fail in a strip below `1-ρ`, already at the list-decoding
  capacity `1-H_q(ρ)`. So the prize's `δ*_C` is **strictly below `1-ρ`** for every rate. The
  reachable ceiling is at best `≈ 1 - H_q(ρ)` (list-decoding capacity), which is `~1/log₂q`
  below `1-ρ`. For a Mersenne-31/BabyBear field (`log₂q≈31`) this ceiling is `(1-ρ) - O(1/31)`.
- **At ρ=1/2:** Singleton capacity `1-ρ=0.5`; list-decoding capacity `1-H_q(0.5) = 0.5 - O(1/log q)`.
  CA is `= 1` (totally broken) once `δ` enters the Cor-1 window near `0.5`. So `δ*_C(ρ=1/2) < 0.5`,
  and is at most about `0.5 - 1/log₂q ≈ 0.5 - 0.032 ≈ 0.468` for `q≈2^31` (UPPER bound on the
  ceiling, not a positive achievability — positive achievability only known to the Johnson radius
  `0.293` per ABF Thm 4.12).
- **For the prize's threshold form (`ε* = 2^-128`):** Cor 1 gives `ε_ca=1 ≫ 2^-128` in its window,
  so that window is firmly outside `δ*_C`. The paper does NOT give a `2^-128`-calibrated positive
  δ; it only pushes the ceiling down to list-decoding capacity and refutes capacity claims.
- **Theorem 2 is mostly a barrier-flag for our SNARK params:** at `q≈2^31, k≈2^20` it says
  nothing, but it warns that any FUTURE small-error (`ε≪1/k`) positive CA result we might prove
  would automatically imply a (hard, decades-open) RS list-decoding result — i.e. the easy route
  is blocked; we should expect to need large-error (M)CA.

---

## Placement of δ_known_positive / δ_known_negative implied by this paper
- **δ_known_negative (this paper's contribution):** `δ` at the **list-decoding capacity**
  `1 - H_q(ρ)` (equivalently k in the Cor-1 window) already gives `ε_ca = 1`. So
  `δ_known_negative ≤ 1 - H_q(ρ)` in the sense that AT this radius CA is fully broken. Since
  `H_q(ρ) ≈ ρ + 1/log_2 q`... [careful: Cor 1 is phrased in terms of f and k; the cleanest
  statement is `ε_ca(C,δ)=1` whenever `1-H_q(δ)+2/n+sqrt((H_q(δ)-δ)/n) ≤ ρ ≤ 1-δ-2/n`]. The
  **upper ceiling for `δ*_C` is therefore `≈ 1 - H_q(ρ) = (1-ρ) - Θ(1/log q)`, NOT `1-ρ`.**
- **δ_known_positive:** this paper proves NO positive lower bound on `δ*_C`. (It provides the
  modified conjectures, but those are conjectures.) Positive knowledge stays at the Johnson
  radius from BCHKS / ABF Thm 4.12.
- **Net gap for smooth prime RS:** `[ 1-√ρ (positive, from elsewhere), 1-H_q(ρ) (this paper's
  ceiling) ]`, a strip of width `≈ (1-H_q(ρ)) - (1-√ρ) = √ρ - H_q(ρ) ≈ √ρ - ρ - 1/log q`.

---

## Open questions the paper states (§"Directions for future research")
1. Pursue (quantum) attacks on **bounded distance decoding** of RS (relevant since hash-based
   SNARKs are pitched as post-quantum; existing BDD-hardness results — [GV05] binary, [GGG18]
   exp-large fields, [CW04] full-domain — do NOT apply to SNARK-sized fields).
2. **Prove or disprove the three modified conjectures** (Our Conj 1–3) — up to list-decoding
   capacity for prime fields.
3. **Pursue large error probabilities** for (mutual) correlated agreement — because (Thm 2)
   small-error CA would imply list-decodability (a much harder, long-open problem). So the
   productive regime is large-error (M)CA.

---

## INFERRED — VERIFY
- **Notation conflict with ABF / others:** Crites–Stewart use `f` = ABSOLUTE error count and
  `δ = f/n`; ABF and BCHKS use `δ` directly as the relative radius. Also CS write list-decoding
  capacity `1 - H_q(δ)` whereas ABF's "capacity" = `δmin ≈ 1-ρ` (Singleton). When comparing
  "capacity," ABF means Singleton `1-ρ`; CS's central point is that the relevant ceiling is the
  SMALLER list-decoding capacity `1-H_q(ρ)`. Keep these straight in synthesis.
- Theorem 1's exact "number of λ" denominator rendered awkwardly in the PDF
  (`C(n,f) / ( q^{n-f-k} g(...) )` form). The structure (binomial over `q^{n-f-k}·g`) is clear
  and consistent with the proof, but the precise placement of the binomial should be re-checked
  against the source PDF before citing the bound numerically.
- My numeric `δ*_C(ρ=1/2) ⪅ 0.468` for `q≈2^31` is arithmetic (`0.5 - 1/31`) from Claim 1's
  `1/log₂q` bound; it is an UPPER bound on the achievable ceiling, not something the paper
  states as a number, and it is NOT a proven positive achievability.
- The `q ≥ n` hypothesis in Cor 1 means the cleanest negative results are stated for fields not
  too much larger than the domain. For our 256-bit-field targets with `n≤2^30`, `q ≥ n` holds
  easily, but the *strength* of the `ε_ca=1` window (its width in δ) shrinks with `1/log q`, so
  for very large fields the broken strip is thinner — VERIFY how Cor 1 degrades as `q ≫ n`.
- Whether Cor 1 / Thm 1 hold specifically for the deployed **smooth multiplicative-subgroup**
  domain: the constructions only need a deep hole, which exists for any domain, so I read this
  as YES (applies to smooth domains), but the paper does not single out smooth domains — VERIFY.
