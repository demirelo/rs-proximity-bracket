# Bordage, Chiesa, Guan, Manzur — All Polynomial Generators Preserve Distance with Mutual Correlated Agreement

> Content extracted verbatim from the source PDF (eprint 2025/2051, dated **May 19, 2026**;
> retrieved via curl 2026-06-02; 55pp). This is the leading POSITIVE result on MCA. It is
> the up-to-Johnson positive side; its near-capacity *negative* remark points to DG25/CS25
> (Diamond–Gruen / Carmon–Stewart) and is consistent with Kambiré.

## Bibliographic info

- **Title:** "All Polynomial Generators Preserve Distance with Mutual Correlated Agreement"
- **Authors:** Sarah Bordage, Alessandro Chiesa, Ziyi Guan, Ignacio Manzur (all EPFL)
- **Identifier:** Cryptology ePrint Archive, Paper 2025/2051. PDF title-page date: **May 19, 2026**
  (the "2025" is the ePrint slot; this is a later revision).
- **Keywords:** proximity testing, distance preservation, mutual correlated agreement.
- **Funding note (verbatim):** "Ethereum Foundation has allocated over $1M towards resolving
  certain conjectures about MCA [E25]." (This is the Proximity Prize.)

## Main results — exact statements + formulas (define symbols)

### Definitions (Section 3, verbatim)
- **Def 3.8 (RS code):** `RS[F, D, k] ⊆ F^{|D|}` = `[|D|, k, |D|−k+1]_F`-code of evaluations
  over `D ⊆ F` of polynomials of degree `≤ k−1`. Rate `ρ := k/|D|`. (n := |D|.)
- **Def 3.11 (zero-evading):** `G : S → F^ℓ` is zero-evading with error `ε_ZE` if
  `max_{v ∈ F^ℓ\{0}} Pr_{x∈S}[G(x)·vᵀ = 0] ≤ ε_ZE`.
- **Def 3.12 (MDS generator):** `M_G` = `(|S|×ℓ)` matrix with rows `{G(x)}`, `C_G` = code it
  generates; `G` is MDS if `C_G` is MDS.
- **Def 3.14 (MCA), verbatim:** `G : S → F^ℓ` has mutual correlated agreement for `C ⊆ Σ^n`
  with error `ε_MCA : [0,1]→[0,1]` if for all `u_1,…,u_ℓ ∈ Σ^n`, `γ ∈ [0,1]`, with
  `U := Mat(u_1,…,u_ℓ)`:
  `Pr_{x∈S}[ ∃ T ⊆ [n] : |T| ≥ n(1−γ) ∧ (G(x)·U)|_T ∈ C|_T ∧ ∃ j∈[ℓ], u_j|_T ∉ C|_T ] ≤ ε_MCA(γ)`.
- **Def 3.21 (CA):** `G` has CA with error `ε_CA(e,t)` if for all `1≤t<e≤n` and `U` with
  `Δ(U, C^ℓ) ≥ e`: `Pr_{x∈S}[Δ(G(x)·U, C) ≤ e−t] ≤ ε_CA(e,t)`.
- **Lemma 3.22 (MCA ⇒ CA):** if `G` has MCA with `ε_MCA` then it has CA with
  `ε_CA(e,t) := ε_MCA((e−1)/n)`.
- **Def 3.1 (Def 2, polynomial generator):** `G(x) = (P_j(x))_{j∈[ℓ]}`, `P_1,…,P_ℓ` linearly
  independent polynomials in `F[X_1,…,X_s]`.

### Theorem 1 (informal, all linear codes — Section 1.1, verbatim)
Let `G : S_1×…×S_s → F^ℓ`, `G(x)=(P_j(x))`, `d_i := max_j deg_{X_i}(P_j)`. For every linear
code `C ⊆ F^n` with relative distance `δ_C` and tradeoff `η ∈ (0,1)`, `G` has MCA with error
```
              max{n·γ, 1} · Σ_i (d_i/|S_i|)            if γ ≤ δ_C / (max_i d_i + 2)
ε_MCA(γ) :=
              O( (n/η) · Σ_i (d_i/|S_i|) )             if γ ≤ 1 − (1 − δ_C + η)^{1/(max_i d_i + 2)}
```
Both bounds = zero-evading error `ε_ZE ≤ Σ_i d_i/|S_i|` times a multiplicative overhead. The
unique-decoding (upper) bound is proven tight for univariate powers (App B).

### Theorem 6.1 (FORMAL, MDS generators for all linear codes — Section 6, verbatim)
Let `G : S → F^ℓ` be an MDS generator with `dim C_G = ℓ ≥ 2`. Let `C ⊆ Σ^n` be `F`-linear with
relative distance `δ_C`. Set `ρ_C := 1 − δ_C`. For every `η ∈ (0,1)`, `G` has MCA with error
```
              max{n·γ, 1} · (ℓ−1)/|S|                                              if γ < δ_C/(ℓ+1)
ε_MCA(γ) :=   (n·γ_ℓ / η)·(ℓ−1)/|S|
              + max{ 2(ℓ−1) / (η·((ρ_C+η)^{1/(ℓ+1)} − (ρ_C+η)^{1/ℓ})·|S|) ,
                     ℓ(ℓ+1)/(η·|S|) }                                              if γ ≤ 1 − (ρ_C+η)^{1/(ℓ+1)}
              1                                                                    otherwise
```
where `γ_ℓ := 1 − (ρ_C + η)^{1/ℓ}`. (Lemma 6.2 isolates the unique-decoding case:
`ε_MCA(γ) = max{n·γ,1}·(ℓ−1)/|S|` for `γ < δ_C/(ℓ+1)`.)

### Lemma 1 (affine space = affine line, Section 1.1, verbatim)
For any linear code `C ⊆ F^n`, the affine space generator `G(x_1,…,x_s)=(1,x_1,…,x_s)` has MCA
with error
```
              max{n·γ,1} · 1/|F|        if γ ≤ δ_C/3
ε_MCA(γ) :=
              O( (n/η) · 1/|F| )        if γ ≤ 1 − (1 − δ_C + η)^{1/3}
```
**Independent of the number of variables `s`** (matches `ε_ZE = 1/|F|`), and imposes NO lower
bound on `|F|` (unlike prior CA analyses, e.g. [GG25]).

### Theorem 2 / Definition 9.1 (FORMAL, RS codes up to Johnson — Sections 1.1 & 9, verbatim)
For `C = RS[F, D, k]`, `ρ := k/n`, polynomial generator `G(x)=(P_j(x))`, `d_i := max_j deg_{X_i}(P_j)`,
and any integer `m ≥ 3`:
```
                              (1/|F|) · ((m+1/2)^7 / (3ρ^{3/2})) · d · n^2     if γ ≤ 1 − (1 + 1/(2m))·√ρ
ε_MCA,RS,d(γ) :=
                              1                                               otherwise
```
**Theorem 9.2:** `G` has MCA for `RS[F,D,k]` with error `Σ_{i∈[s]} ε_MCA,RS,d_i`.
Informal Theorem 2 (single-`d` version): `ε_MCA(γ) = O( (m^7 n^2 Σ d_i)/(ρ^{3/2} |F|) )` for
`γ ≤ 1 − (1 + 1/(2m))·√ρ`. This **matches the BCIKS 2020 correlated-agreement error** and
extends it to MCA and all polynomial generators. Lemma 9.3 (univariate powers
`x ↦ (1,x,…,x^d)`) is the base case, modifying Haböck [Hab25]'s affine-line proof and the
BCIKS Guruswami–Sudan curve argument (`D_X < (m+1/2)√ρ·n`, `D_Y < (m+1/2)/√ρ`,
`D_{YZ} < (m+1/2)^3/(6√ρ)·nd`).

### Lemma 10.1 (interleaving — Section 10, verbatim)
If `G` has MCA for `C ⊆ Σ^n` with error `ε_MCA`, then for every `k ∈ N`, `G` has MCA for the
`k`-interleaving `C^k` with error **`k · ε_MCA`**. (Direct sum: Lemma 10.3; concatenation:
Lemma 10.5.) — Directly relevant to `C^{equiv m}` in our soundness expression.

### Lower bound (Appendix B, verbatim)
**Lemma B.1.** Let `d ≥ 1`, `F` with `d | |F|−1`. For a linear code `C ⊆ F^n` with distance
`∆_C`, `0 < e ≤ (∆_C − 1)/2`, `|F| ≥ e·d + 1`. Then there exist `u_0,…,u_d ∈ F^n` with
`max_j Δ(u_j, C) ≥ e` and `Pr_{x∈F}[Δ(Σ_j x^j u_j, C) < e] ≥ e · d/|F|`.
Construction: `ω ∈ F^×` of order `d`, `G = ⟨ω⟩`, `e` cosets `G_i = a_i·G`;
`u_0 = (a_1^d,…,a_e^d, 0^{n−e})`, `u_j = 0` for `1 ≤ j ≤ d−1`, `u_d = (−1,…,−1,0^{n−e})`.
**This lower bound is a TIGHTNESS result at the unique-decoding radius** (matches the upper
`ε_MCA` within `δ_C/(d+2)`), NOT a near-capacity failure.

### Near-capacity NEGATIVE remarks (Section 1.1 "Open questions", verbatim)
- "[ACFY25] conjectures MCA for the univariate powers generator for `γ` close to the code's
  capacity `1−ρ` ... Yet **[DG25; CS25] shows that this is too much to hope for: there is a
  family of Reed–Solomon codes for which distance preservation fails for the affine line
  generator `x ↦ (1,x)` when `γ` is close to the code's capacity `1−ρ`.**"
- "[BGKS20] constructs a linear code `C` and words `u_1, u_2` far from `C` such that their
  linear combination via the affine line generator drops with probability almost one." — gives
  intuition that any all-linear-codes result holds only up to some `γ` and errors → 1 "when
  the blocklength `n` of `C` approaches `|F|`."

## Regime of validity — δ range; code class; fields

- **δ range (positive):** up to the **Johnson radius** `γ ≤ 1 − (1+1/(2m))√ρ → 1 − √ρ` as
  `m→∞` for RS (Thm 2/9.2); up to `1 − (1−δ_C+η)^{1/(ℓ+1)}` for general MDS generators
  (list-decoding regime, Thm 6.1). The all-linear-codes Theorem 1 is the same shape with
  `δ_C` in place of RS structure.
- **Code class:** (i) ALL `F`-linear codes (Thm 1, Thm 6.1 for MDS generators); (ii) Reed–
  Solomon `RS[F, D, k]` for **arbitrary** evaluation domain `D ⊆ F` (Thm 2 / 9.2) — this
  **includes smooth multiplicative-subgroup domains** since `D` is unrestricted. So the RS
  Johnson-bound positive result DOES apply to our smooth-domain case.
- **Fields:** general finite fields `F`; **no lower bound on `|F|`** for the affine-space MCA
  (Lemma 1). Works for prime and extension fields alike.
- **Generators covered:** all polynomial generators `(P_j)` linearly independent — includes
  identity, univariate powers `(1,x,…,x^d)`, affine space `(1,x_1,…,x_s)`, multilinear, and the
  low-degree-poly generators of BCL22/BCGL22/DP24a/MZ25/BMMS25b. Also one explicit
  NON-polynomial zero-evading generator [AGHP92] is shown to have MCA (Section 2.6) — first
  result beyond polynomial generators.

## Relevance to our targets — ρ∈{1/2,1/4,1/8,1/16}, smooth domain, ε*=2^-128

- This pins `δ_known_positive` for smooth-domain RS at the **Johnson radius** `1 − √ρ` with a
  concrete, MCA-grade error: `ε_MCA(γ) ≤ (m^7 n^2 d)/(O(ρ^{3/2})|F|)` for `γ ≤ 1−(1+1/(2m))√ρ`,
  for the line generator (`ℓ=2`, `d=1`). For `ε* = 2^-128` we need
  `(m^7 n^2)/(ρ^{3/2}|F|) ≤ 2^-128` and choose `m` to push `γ` toward `1−√ρ` (larger `m` ⇒
  closer to Johnson but worse error, so there is a concrete `m`-vs-`δ` trade — Workstream B).
- Concrete Johnson radii: ρ=1/2 → `1−√0.5 ≈ 0.293`; ρ=1/4 → `0.5`; ρ=1/8 → `1−0.354 ≈ 0.646`;
  ρ=1/16 → `0.75`. These are the largest δ where we currently have a PROVEN small `ε_mca` on
  smooth domains.
- Interleaving (Lemma 10.1) gives the `C^{equiv m}` term: `ε_MCA(C^m, γ) ≤ m·ε_MCA(C, γ)` — a
  benign linear factor, so the interleaving bound is NOT the bottleneck below the Johnson
  radius.

## Placement of δ_known_positive / δ_known_negative implied

- `delta_known_positive` (smooth RS, all fields): **`1 − √ρ`** (Johnson), with MCA error
  `~ m^7 n^2 / (ρ^{3/2}|F|)`. This is the strongest positive smooth-domain anchor in the lit.
- `delta_known_negative` (this paper, via DG25/CS25 + Kambiré): distance preservation fails
  "close to capacity `1−ρ`" for the affine line generator. This paper does not give the
  explicit `δ` (it cites DG25/CS25); Kambiré makes it explicit at `(1−ρ) − 2/(K log n)` on
  smooth prime-field domains.
- Between `1−√ρ` and `(1−ρ)−o(1)`: OPEN.

## Open questions stated (Section 1.1, verbatim/paraphrased)

- Full characterization of which (zero-evading) generators have MCA — proven beyond polynomial
  generators, so the characterization remains open.
- MCA behavior under tensor products and distance-amplification operations on codes.
- For RS: "Can one bound the MCA error for distances `γ` beyond `1 − √ρ`?" Whether MCA holds
  in the band `(1−√ρ, 1−ρ)` between Johnson and capacity is open (and DG25/CS25/Kambiré show
  it cannot hold all the way to capacity).

## INFERRED — VERIFY

- The RS error `(m+1/2)^7/(3ρ^{3/2})·d·n^2/|F|` is for the univariate-powers/affine-line case
  (`ℓ=2`, `d=1` for the line). For our line generator I read `d=1`, giving
  `ε_MCA ≈ (m+1/2)^7 n^2 / (3ρ^{3/2}|F|)`. *Verify the exact `d` and `m` plugging for the
  specific generator we deploy (FRI uses the affine line / univariate powers).* The
  `m`-vs-`γ` trade `γ ≤ 1−(1+1/(2m))√ρ` is the lever for Workstream B's query-count calc.
- DG25 and CS25 are cited as the source of the near-capacity affine-line failure; this note
  records only Bordage–Chiesa's one-sentence summary. *The actual `δ` and code class in
  DG25/CS25 must be read from those papers* (CS25 = Carmon–Stewart 2025/2046, in the project
  lit list; DG25 = Diamond–Gruen). Bordage–Chiesa's phrasing "family of Reed–Solomon codes ...
  affine line generator ... close to capacity" is consistent with Kambiré but does not by
  itself say whether DG25/CS25 use smooth domains or prime fields.
- Lemma B.1 lower bound is at unique decoding `e ≤ (∆_C−1)/2`; it is a *tightness* witness for
  the constant in the unique-decoding `ε_MCA`, NOT evidence about the Johnson-to-capacity band.
