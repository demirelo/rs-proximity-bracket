# Kambiré — Proximity Gaps Conjecture Fails Near Capacity over Prime Fields

> CRITICAL NEGATIVE RESULT. This is the single most important paper for `delta_known_negative`
> on the deployed smooth-domain case. All content below is extracted verbatim from the
> source PDF (arXiv:2604.09724v1, retrieved via curl 2026-06-02). The paper is a 6-page
> note; I read it in full.

## Bibliographic info

- **Title:** "Proximity Gaps Conjecture Fails Near Capacity over Prime Fields"
- **Author:** Antonio Kambiré (single author)
- **Identifier:** arXiv:2604.09724v1 [cs.IT] (cross-listed cs.CR, math.IT), CC-BY 4.0
- **Dated:** "April 1st 2026" on title; arXiv stamp "9 Apr 2026"
- **Provenance:** "flesh out a sketch by Krachun and Kazanin" — ref [4] = Dmitry Krachun
  and Stepan Kazanin, "Failure of the proximity gap conjecture for Reed-Solomon code close
  to the capacity regime", *Personal communications*, 2026.
- **Key cited works:**
  - [1] Ben-Sasson, Carmon, Haböck, Kopparty, Saraf, "On proximity gaps for Reed–Solomon
    codes," ECCC Report No. 169 (2025) — the construction "follow[s] the same
    multiplicative-subgroup and sumset template as **Theorem 7.1 of [1]**."
  - [2] Ben-Sasson, Carmon, Ishai, Kopparty, Saraf, "Proximity gaps for Reed-Solomon
    codes," ePrint 2020/654 (the original proximity-gaps conjecture).
  - [3] Carmon, Goldberg, Haböck, Lerer, Lesokhin, Papini, Samocha, "S-two whitepaper,"
    ePrint 2026/532, 2026 — "Appendix A.5 of [3] used this result to formalize a
    conjecture on list- and curve-decodability properties of Reed-Solomon codes over prime
    fields, up to the information-theoretic limit."

## Main results — exact statements + formulas (define symbols)

### Conjecture being refuted (Introduction, verbatim)
> "The proximity gaps conjecture, introduced in [2], asserts that if many points on an
> affine line `f + z·g` are each close to a Reed-Solomon code, then the line itself must be
> explained by a nearby codeword pair, that is, the pair `[f, g]` is close to the
> corresponding interleaved Reed-Solomon code, a condition called correlated agreement."

### Headline claim (Introduction, verbatim)
> "We construct block lengths `n` and dimensions `k`, a prime field `F_p`, and words
> `f, g ∈ F_p^D` such that for `δ = (1 − k/n) − Ω(1/log n)` there are at least `n^C` distinct
> scalars `z ∈ F_p` with `Δ(f + z·g, C) ≤ δ`, while simultaneously `Δ([f,g], C^2) > δ`, where
> `C = RS[F_p, D, k]`."

### Theorem 1 (verbatim, Section 2)
> **Theorem 1.** For every constant `C > 0` and rate `ρ ∈ (0, 1/2)`, there exist infinitely
> many block lengths `n`, dimensions `k`, such that with `δ = (1 − k/n) − Ω(1/log n)`, the
> following hold:
> - There exists a prime `p < n^A` with `p ≡ 1 (mod n)`, for some constant `A = A(ρ, C)`.
> - Let `ω` be a primitive `n`-th root of unity in `F_p^×`, set `D = ⟨ω⟩`, and `C = RS[F_p, D, k]`.
>   Then there exist `f, g ∈ F_p^D` such that
>   `|{ z ∈ F_p : Δ(f + z·g, C) ≤ δ }| ≥ n^C`,    `Δ([f,g], C^2) > δ`.

Symbols: `Δ(·, C)` = relative Hamming distance to the code; `C^2` = the 2-wise interleaved
RS code (so `Δ([f,g], C^2) > δ` means **no** δ-correlated-agreement / no common agreement set
of density `≥ 1−δ` for the pair `(f,g)`); `D = ⟨ω⟩` is the multiplicative subgroup of order
`n`; rate `ρ = k/n`.

This is a violation of BOTH the (single-line) proximity-gap statement AND correlated
agreement: a line with `≥ n^C` near-codeword points (`n^C >> 1`, super-polynomially many of
the `q = p < n^A` scalars) yet `[f,g]` is δ-far from the interleaved code.

### The exact gap below capacity (Setting-parameters block, verbatim)
The "Ω(1/log n)" is made fully explicit. Choose integers:
- `C > 0` (controls how many distinct sums).
- `ρ := u / 2^v ∈ (0, 1/2)` with `u, v ∈ Z_{≥0}`, `u < 2^{v-1}` (the code rate).
- `L(ρ,C) = max{ C / (ρ·log(1/(2ρ))) , (9/2)·log 8 }`, and pick `K = K(ρ,C)` a power of 2 with
  `L(ρ,C) ≤ K ≤ 2·L(ρ,C)`.
- `s := 2^α` for integer `α` large enough that `α ≥ v`, `α ≥ log_2 K`, and `K ≤ 2^α` with
  `2^α / K ≥ α` (so `K | 2^α`). "This controls how close we get to capacity."
- `r := ρs + 2 = u·2^{α−v} + 2` (a positive integer; "controlling the relative distance").
- `m := 2^{ 2^α / K − α }` (a power of 2).
- **`n := s·m`,  `k := (r−2)·m`.**

Verified identities (verbatim):
- (1) `ρ = (r−2)/s = (r−2)m / (sm) = k/n`.
- (2) `K·log_2 n = K·log_2(sm) = K·log_2(2^α · 2^{2^α/K − α}) = K·(2^α/K) = 2^α = s`.

Therefore (verbatim):
> "Then we set `δ := 1 − r/s` which is `η := (1 − ρ) − δ = 2/s ∈ Ω(1/log(n))` away from
> capacity."

**`η = (1−ρ) − δ = 2/s = 2/(K·log n)`.** Since `K = K(ρ,C)` is a constant depending only on
`(ρ,C)`, the gap below capacity is `η = Θ(1/log n)` — it vanishes as `n → ∞`. Equivalently
`s = K·log n` (natural log; the doc sets `log ≡ ln` for the number-theory part).

### Construction of the bad line (Proof, "Constructing the Counterexample", verbatim)
- `ξ` = a primitive `s`-th root of unity in `F_p^×`; `H := ⟨ξ⟩ ⊂ D` (subgroup of order `s`).
- `H^{(+r)} := { Σ_{i=1}^r e_i | e_1,…,e_r ∈ H distinct }` (the `r`-fold distinct-element sumset).
- **`f := X^{rm}`,  `g := X^{(r−1)m}`**, line `L := { f + λ·g | λ ∈ F_p } ⊂ F_p^D`.
- Claim: for each `λ = ξ_1 + ξ_2 + … + ξ_r ∈ H^{(+r)}`, `Δ(X^{rm} + λ·X^{(r−1)m}, C) ≤ δ`.
- Number of cosets of `H` in `D` is `|D|/|H| = n/s = m`. Pick `r` cosets:
  `H_j := { a ∈ D | a^m = ξ_j }` for `j = 1,…,r`.
- Polynomial identity (verbatim):
  `∏_{a ∈ H_1 ∪…∪ H_r} (X − a) = ∏_{j=1}^r (X^m − ξ_j) = X^{rm} − (ξ_1+…+ξ_r) X^{(r−1)m} + R(X)`
  `= X^{rm} − λ X^{(r−1)m} + R(X)`, where `deg R ≤ (r−2)m`.
- Hence `X^{rm} − λ X^{(r−1)m}` agrees with `R(X)` (a degree-`< k` polynomial, `k=(r−2)m`) on
  the set `H_1 ∪…∪ H_r`, which has size `rm = (1−δ)n`. So `Δ(X^{rm} + λX^{(r−1)m}, C) ≤ δ`.
  (Note `1 − r/s = 1 − rm/n = 1 − (1−δ)`... i.e. agreement on `rm` points gives distance `≤ δ`
  with `δ = 1 − r/s`.)
- **No correlated agreement (verbatim):** "assume it did. Then we would have some `D' ⊂ D` of
  size `(1−δ)n` for which every point on `L` agrees with some polynomial of degree at most `k`.
  In particular, `X^{(r−1)m}` agrees with a polynomial `q(X)` of degree at most `k` on `D'`, but
  this implies `|D'| ≤ k = (r−2)m` as `q(X)` can have at most `k` roots in `F_p`. This
  contradicts `|D'| = rm`." (Because `X^{(r−1)m}` itself is degree `(r−1)m > k`, so it cannot
  agree with any degree-`<k` poly on more than `k` points.)

### Counting many distinct scalars (Proof, "Counting the Number of Sums")
- Need `|H^{(+r)}|` distinct values of `λ` to remain distinct in `F_p` (no collisions mod `p`).
- Uses a **quantitative Linnik theorem** to find a prime `p ≡ 1 (mod n)` in `[4^s, 8^s]` with
  `p < 8^s = 8^{K log n} = n^{K log 8}`, so `A = K·log 8`. Counts `T ≥ 8^s/(n^{3/2} log(8^s))`
  such primes (uses `n = 2^t ⇒ φ(n)=n/2`).
- "Bad" primes are those dividing `Res(Φ_s, Q)` where `Q(x)=x^{i_1}+…+x^{i_r}−(x^{j_1}+…+x^{j_r})`
  (i.e. primes causing two distinct `r`-tuples to give the same sum). Bound: each tuple-pair
  yields `B ≤ log_4(s)` bad primes; total bad triples `≤ B·binom(s,r)^2 = log_4(s)·(2^s)^2 << T`
  for large `s` (using `K > 9/(2 log 8)`). So a good prime exists.
- **List/multiplicity size (verbatim):**
  `a := |H^{(+r)}| = binom(s/2, r) ≥ (s/(2r))^r`, and rewriting `r = ρs + 2`:
  `a ≥ (1/(2ρ))^{ρs+2} ≈ (1/(2ρ))^{ρ K log n + 2} = n^{ρ K log(1/(2ρ))} · (1/(2ρ))^2`.
  "Since `K > C/(ρ log(1/(2ρ)))` and `1/(2ρ) > 1`, we get `a > n^C` as needed."

So the number of scalars `z` with `Δ(f+zg, C) ≤ δ` is `≥ n^C` for an arbitrarily large
constant `C` (choose `K` large). The fraction of bad scalars is `≥ n^C / p ≥ n^C / n^A`,
which is a nontrivial (super-`1/poly`) fraction since `C` can be pushed up.

## Regime of validity — δ range; code class; fields

- **δ range / radius of failure:** `δ = (1 − ρ) − η` with `η = 2/s = 2/(K log n)`. The failure
  is at radius **`Θ(1/log n)` below capacity `(1−ρ)`**, i.e. asymptotically AT capacity. It
  says NOTHING about failure at or below the Johnson radius `1 − √ρ` (where positive results
  hold). The construction needs `r = ρs+2` and `rm = (1−δ)n` agreement points, so `δ` is
  pinned to `1 − r/s` exactly; it does not give a family of failures at smaller `δ`.
- **Code class:** Reed–Solomon `C = RS[F_p, D, k]` with **`D = ⟨ω⟩` a multiplicative subgroup
  of order `n`** (a smooth, FFT-friendly domain — exactly the FRI/STIR/WHIR deployed case).
  `n = sm` is a product of powers of 2 ⇒ `n = 2^t` smooth (the proof explicitly uses
  `n = 2^t ⇒ φ(n) = n/2`). It is NOT a random/punctured domain and NOT the full field.
- **Fields:** **PRIME fields only** — `F_p` with `p ≡ 1 (mod n)`, `p < n^A`, `A = K log 8`.
  Nothing about extension fields; the prime structure (Linnik, primes `≡ 1 mod n`) is essential
  to make the subset sums distinct. The field is `≈ n^A`-sized (polynomial in `n`), i.e. SMALL
  relative to `n` (NOT a 128-bit / 256-bit field for the asymptotic family — see below).
- **Rate:** `ρ ∈ (0, 1/2)` (covers ρ ∈ {1/4, 1/8, 1/16} and the open boundary toward 1/2; the
  stated range is the open interval `(0, 1/2)`, so `ρ = 1/2` is excluded as stated, but
  `ρ → 1/2^-` is allowed).

## Relevance to our targets — ρ∈{1/2,1/4,1/8,1/16}, smooth domain, ε*=2^-128

- This is the FIRST explicit, self-contained counterexample to the up-to-capacity proximity-
  gap / correlated-agreement statement **on a smooth multiplicative-subgroup domain `D=⟨ω⟩`
  with `n=2^t`** — i.e. directly the deployed FRI/STIR/WHIR setting. It applies for our rates
  ρ ∈ {1/4, 1/8, 1/16} (all in `(0,1/2)`) and ρ → 1/2.
- **What it rules out:** any theorem of the form "`ε_mca(C, δ) ≤ 2^-128` for `δ` up to
  `(1−ρ) − o(1)` (gap shrinking with `n`) over prime fields with smooth domains" is FALSE.
  At `δ = (1−ρ) − 2/(K log n)` the proximity gap *fails outright* (a constant/non-negligible
  fraction of scalars are bad), so `ε_mca` and the bad-combining probability are `≥ n^C/p`,
  nowhere near `2^-128`. Equivalently the interleaved list `Λ(C^{equiv 2}, δ)` ... actually the
  failure is that `[f,g]` is δ-FAR yet the line has `n^C` near points — a direct proximity-gap
  failure, which forces any honest `ε_mca` bound to be large at this δ.
- **CRUCIAL CAVEAT for our `ε* = 2^-128` targets (field-size mismatch):** the counterexample
  family uses `p < n^A = n^{K log 8}` (a SMALL field, polynomial in `n`). For deployed
  parameters we use `|F| ≈ 2^128`–`2^256` with `n = 2^r` up to ~`2^30`. In that regime
  `log_2 n ≤ ~30`, so `K log n` is a CONSTANT-sized `s`, and `η = 2/(K log n)` is a CONSTANT,
  not `o(1)`. The asymptotic statement "fails at `o(1)` below capacity" is an `n → ∞`
  statement; for fixed finite `(n, |F|)` it provides a counterexample at a *constant* additive
  gap `η = 2/s` below capacity, with `s` constrained by the requirement that a good prime
  `p ≡ 1 (mod n)`, `p ∈ [4^s, 8^s]`, exists and is `≤ |F|`. For `|F| ≈ 2^128`: need
  `8^s ≤ 2^128 ⇒ s ≤ 128/3 ≈ 42`, giving `η = 2/s ≥ ~0.047` below capacity at best. So over a
  128-bit prime field the construction (as written) only forces failure roughly `η ≳ 0.05`
  below capacity, NOT arbitrarily close to it. **This is the key number to pin down precisely
  in Workstream A/B.** (See INFERRED below.)
- The construction is for the 2-wise line (affine line generator `x ↦ (1,x)`), which is the
  weakest generator and the base case for all batch/folding reductions — so it bounds the
  whole hierarchy (a line counterexample propagates to FRI folding).

## Placement of δ_known_positive / δ_known_negative implied

- `delta_known_negative` (smooth-domain RS, prime field): proximity gaps / CA / MCA
  **FAIL** at `δ = (1−ρ) − 2/(K(ρ,C)·log n)` for the line generator. As `n→∞` this →
  capacity `1−ρ`. For finite `|F| = p`, the achievable gap is `η = 2/s` with `s` bounded by
  `8^s ≤ p`, i.e. `s ≤ log_8 p = (log_2 p)/3`, so the proven failure point is at most
  `η_min ≈ 2 / ((log_2 |F|)/3) = 6 / log_2 |F|` below capacity (≈ 0.047 for 128-bit). I.e.
  `δ_known_negative ≤ (1−ρ) − 6/log_2|F|` for prime-field smooth RS (concretely; see INFERRED).
- `delta_known_positive` (smooth-domain RS, any field): the matching POSITIVE bound is the
  Johnson radius `δ < 1 − √ρ` (Bordage–Chiesa Theorem 9.2; BCIKS 2020), where `ε_mca` is
  `O(n^2/|F|)`-small. Nothing positive is known between `1−√ρ` and the failure point.
- So the gap interval `[δ_known_positive, δ_known_negative]` for smooth-domain RS at, e.g.,
  ρ=1/4 is roughly `[1−√(1/4), (1−1/4) − η] = [0.5, 0.75 − η]`, with `η` constant for finite
  fields. The Johnson radius and capacity differ by `√ρ − ρ`; for ρ=1/4 that's `0.5 − 0.25 = 0.25`,
  and the negative result only closes the top `η ≈ 0.05` of it (for 128-bit), leaving a wide
  open band.

## Open questions stated

The note itself is terse and states no explicit open-problem list, but implies:
- It addresses behavior "above the Johnson bound," which is "still less clear and the subject
  of active research." The whole region `(1−√ρ, (1−ρ)−η)` between the Johnson radius and the
  proven failure point is left open.
- Via [3] (S-two whitepaper, App A.5): this result is used to "formalize a conjecture on
  list- and curve-decodability properties of Reed-Solomon codes over prime fields, up to the
  information-theoretic limit" — i.e. the precise conjectural threshold is open.

## INFERRED — VERIFY

- **(High value, must verify) Concrete finite-field failure gap.** The "Ω(1/log n)" gap is
  asymptotic. For deployed `|F| ≈ 2^b` (b∈{128,256}), the construction needs a prime
  `p ≡ 1 (mod n)` in `[4^s, 8^s]` with `p ≤ |F|`, forcing `s ≤ b/3` (from `8^s ≤ 2^b`). Then
  `η = 2/s ≥ 6/b`. So I INFER the construction (as literally written) yields a smooth-domain
  prime-field counterexample at gap **`η ≈ 6/log_2|F|`** below capacity — about `0.047`
  (4.7 percentage points) for 128-bit, `0.023` for 256-bit. NOT vanishingly close to
  capacity at deployed sizes. *Verify by working through whether `s` can be decoupled from
  `p` (e.g. larger `n`, or whether `m` can absorb the slack) and recomputing the minimal
  achievable `η(|F|, n)`.* This number directly bounds how much of the "Johnson-to-capacity"
  band the negative result actually forecloses for us.
- **(Verify) `ρ = 1/2` boundary.** Theorem 1 states `ρ ∈ (0, 1/2)` (open). Whether the
  construction extends to exactly `ρ = 1/2` (our top target rate) is not stated. The
  parametrization `ρ = u/2^v`, `u < 2^{v-1}` enforces `ρ < 1/2` strictly. Need to check if a
  separate argument covers `ρ = 1/2`, or if `ρ=1/2` smooth-domain is genuinely untouched by
  this paper.
- **(Verify) Whether the failure is of MCA specifically or only CA/proximity-gap.** The
  theorem violates the proximity-gap statement and `Δ([f,g],C^2) > δ` (no CA). CA failure ⇒
  MCA failure a fortiori (MCA is strictly stronger). So `ε_mca` is large here too — but
  confirm the MCA definition (Bordage–Chiesa Def 3.14 / ACFY Def 4.9) is violated, not just CA.
- **(Verify) Extension-field reach.** Paper is prime-only. Whether an analogous construction
  exists over `F_{p^e}` with smooth subgroup domains (relevant for M31/BabyBear extensions) is
  not addressed. The number-theoretic engine (primes `≡1 mod n`) is prime-specific.
- **(Verify) Relation to BCHKS Thm 7.1 [1] and DG25/CS25.** Kambiré says the line construction
  follows "the same multiplicative-subgroup and sumset template as Theorem 7.1 of [1]"
  (Ben-Sasson–Carmon–Haböck–Kopparty–Saraf, ECCC 169/2025). Cross-check what Thm 7.1 of [1]
  states (likely the same/closely related negative result, possibly with a tighter or
  different parameter regime) and how DG25 (Diamond–Gruen, MDS near-capacity failure) and
  CS25 (Carmon–Stewart) relate. These are separate papers in the project's lit list.
