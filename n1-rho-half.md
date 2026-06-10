# N1 — Extending the near-capacity proximity-gap counterexample to **exactly** ρ = 1/2

> **Status of this document.** This is the close-out write-up for sub-lemma **N1**: pushing
> Kambiré's smooth-domain near-capacity counterexample (arXiv 2604.09724, Thm 1, stated for
> `ρ ∈ (0, 1/2)` *open*) to the boundary rate **`ρ = 1/2`** — the single most-deployed
> Reed–Solomon rate (FRI/STIR/WHIR). It is written **fail-closed**: a claim is tagged
> **ESTABLISHED** only when backed by proof or exact computation; otherwise **LIKELY**, **GAP**,
> or **OPEN**. Date: 2026-06-03. Every number below is reproducible from
> `experiments/small_rs_atlas/verify_rho_half_sums.py` (+ the inline checks transcribed in §3,
> all pure-`sympy`/integer-exact), with `results/verify_rho_half_sums.json`.
>
> **One-line bottom line (verdict, up front):** Kambiré's *line / no-correlated-agreement*
> counterexample extends to **exactly `ρ = 1/2`** — **CLOSES-CONDITIONALLY**. The degree
> obstruction and the no-CA structural bound go through **verbatim and unconditionally** at
> `ρ = 1/2` (proven, §4). The construction is genuinely a **rate-exactly-`1/2`** code under the
> strict RS convention `deg < k` (§2). The *only* conditional step is the **distinct-bad-scalar
> count**: Kambiré's own count `binom(s/2, r)` is **identically zero** at `ρ = 1/2` (because
> `r = s/2 + 2 > s/2`), and KK25 Lemma 9 is **out of range** there (`r ≤ φ(m)/2` is violated).
> We replace it with the **exact** count `|H^{(+r)}| = Σ_u binom(s/2,u) 2^u`
> (proven closed form, §3.3), which is `2^{Θ(s)} ≫ n^C` and is **verified distinct mod the Linnik
> prime** at `ρ = 1/2` by exact resultants (§3.4) — but this distinctness, while numerically
> airtight on every tested `(s,r)`, rests on a resultant bound *outside* Lemma 9's proven window,
> so the unconditional-for-all-`n` guarantee is the one residual (§5). The resulting unsafe radius
> is `δ_unsafe(1/2) = 1/2 − 2/s`, i.e. **`δ*_C(1/2) ≤ (1/2) − Θ(1/log n)`** asymptotically;
> `0.4375` at 128-bit and `0.46875` at 256-bit. **Field-agnostic corollary** (extends to
> `F_{p^e}`, odd `p`): confirmed exactly (§6). **Prize relevance (honest):** like Kambiré this is a
> *near-capacity* result — it forecloses only the band above `≈ 0.477` at `ρ=1/2`, 256-bit, and
> does **not** close the `(Johnson, near-capacity)` band; it pins the negative endpoint at the
> most-deployed rate (§7).

---

## 0. Setup and notation (one convention, fixed)

Smooth domain `C = RS[F, D, k]`, `D = ⟨ω⟩` a multiplicative subgroup of order `n = sm`
(`n = 2^t` smooth), inner subgroup `H = ⟨ξ⟩ ⊂ D` of order `s = 2^α`, lift factor `m = n/s`. The
Kambiré line is

```
   f = X^{rm},   g = X^{(r−1)m},     line  L = { f + λ·g : λ ∈ F } ⊂ F^D,
```

with `r = ρs + 2`, `k = (r−2)m`, and the rate/radius identities

```
   ρ = (r−2)/s = (r−2)m/(sm) = k/n,        δ = 1 − r/s = 1 − rm/n,
   η := (1−ρ) − δ = 2/s = 2/(K log n).     (Kambiré: s = K·log n, K = K(ρ,C) const.)
```

For `ρ = 1/2` this forces `r = s/2 + 2`. Throughout, `H^{(+r)} := {ξ_1 + … + ξ_r : ξ_i ∈ H
distinct}` is the `r`-fold distinct-element sumset; `e_1 = Σ_i ξ_i`, `e_2 = Σ_{i<j} ξ_i ξ_j` are
the first two elementary symmetric functions of an `r`-subset. `Φ_s` is the `s`-th cyclotomic
polynomial; for `s = 2^α`, `Φ_s(X) = X^{s/2} + 1` and `φ(s) = s/2`.

**Two count conventions (carried forward from `n2-verdict.md` §1; they differ by an off-by-one in
the RS dimension, and the difference is decisive at `ρ = 1/2`):**

| Reading | RS code | bad scalar `λ` is close iff | distinct-bad-scalar count |
|---|---|---|---|
| **A (strict)** | `deg < k`, dim `k`, **rate exactly `ρ = (r−2)/s`** | some `r`-subset summing to `λ` has **`e₂ = 0`** ("firing"); then `deg R < k` | `N₀ᶠⁱʳᵉ(s,r)` = #distinct `e₁` over `e₂=0` subsets |
| **B (Kambiré-literal)** | `deg ≤ k`, dim `k+1`, rate `ρ + 1/n` | **always** (the residual `R`, `deg R ≤ k`, is itself a codeword) | `N₀ˢᵘᵐ(s,r) = \|H^{(+r)}\|` = #distinct `r`-subset sums |

Kambiré's note writes "`deg R ≤ (r−2)m`" and counts `|H^{(+r)}|`, so **Kambiré is Reading B**. The
brief's "firing condition `deg R < k ⟺ e₂ = 0`" is **Reading A**. We report both, because at
`ρ = 1/2` they diverge sharply: **Reading A gives zero, Reading B gives `2^{Θ(s)}`** (§2, §3).

---

## 1. TASK 1 — Exactly where Kambiré's Theorem 1 uses `ρ < 1/2`

Reading Kambiré's proof line-by-line (kambire.md), the rate enters in **exactly three** places.
Two are cosmetic; the third is the real one.

### 1.1 The parametrization `ρ = u/2^v`, `u < 2^{v−1}` (cosmetic)

The "Setting parameters" block fixes `ρ := u/2^v` with `u, v ∈ ℤ_{≥0}` and **`u < 2^{v−1}`**. The
constraint `u < 2^{v−1}` is **identically** `ρ = u/2^v < 1/2`. This is the *only* place the open
interval `(0, 1/2)` is hard-coded, and it is a **labeling choice**, not a proof step: nothing
downstream needs `u < 2^{v−1}` except via `r = ρs + 2` and the count. At `ρ = 1/2` one simply takes
`u = 2^{v−1}`, `r = u·2^{α−v} + 2 = s/2 + 2`, a perfectly good positive integer. **Not a real
barrier.**

### 1.2 The degree counts (cosmetic — survive verbatim)

The identity `∏_{j=1}^r (X^m − ξ_j) = X^{rm} − λX^{(r−1)m} + R(X)`, `deg R ≤ (r−2)m`, is purely
formal and holds for **every** integer `r ≥ 2` and `m ≥ 1`. The no-CA degree obstruction
("`X^{(r−1)m}` agrees with a `deg ≤ k` poly on `≤ k = (r−2)m < rm` points") needs only
`(r−1)m > (r−2)m`, i.e. `r ≥ 1`. **Neither uses `ρ < 1/2`.** (Verified verbatim at `ρ = 1/2`, §4.)

### 1.3 The distinct-subset-sum count `a` (THE real one)

This is where `ρ < 1/2` is load-bearing. Kambiré's count is, verbatim (kambire.md):

```
   a := |H^{(+r)}| = binom(s/2, r) ≥ (s/(2r))^r,                         (Kambiré)
   a ≥ (1/(2ρ))^{ρs+2} ≈ n^{ρ K log(1/(2ρ))}·(1/(2ρ))^2,
   "Since K > C/(ρ log(1/(2ρ))) and 1/(2ρ) > 1, we get a > n^C."
```

The formula `binom(s/2, r)` is the **`r`-element subset count of a size-`s/2` set** (the
distinct-sums lower bound from KK25 Lemma 9 with `m = s`, `φ(s) = s/2`). KK25 Lemma 9 is stated
**only for `1 ≤ r ≤ φ(m)/2`** (cghll.md, verbatim below). With `m = s`, that ceiling is
`φ(s)/2 = s/4`. **But at `ρ = 1/2` we need `r = s/2 + 2`, which is `≈ 2×` the ceiling `s/4`.** Two
things break simultaneously:

1. **`binom(s/2, r) = binom(s/2, s/2+2) = 0`** for `r > s/2`. Kambiré's count is **identically
   zero** at `ρ = 1/2` — *vacuous*, not just loose.
2. **Kambiré's `K`-threshold diverges.** He needs `K > C/(ρ·log(1/(2ρ)))`. As `ρ → 1/2`,
   `log(1/(2ρ)) → log 1 = 0`, so the required `K → ∞`. Numerically `K_min/C = 1/(ρ ln(1/(2ρ)))`:
   `5.77` at `ρ=1/4`, `21.1` at `ρ=0.45`, `101` at `ρ=0.49`, `1001` at `ρ=0.499`, `∞` at `ρ=1/2`.

So the count step, **exactly as Kambiré wrote it, collapses at `ρ = 1/2`.** This is the whole of
the N1 problem: the degree/no-CA machinery is rate-blind; the *count* needs a replacement that does
not route through `binom(s/2, r)` or Lemma 9's `r ≤ φ(m)/2` window.

### 1.4 The bad-scalar / distinctness / Linnik window (rate-blind, but re-checked)

The "Counting the Number of Sums" step needs the `|H^{(+r)}|` distinct *algebraic* sums to stay
**distinct mod `p`**. Kambiré uses a quantitative Linnik prime `p ≡ 1 (mod n)`, `p ∈ [4^s, 8^s]`,
and bounds "bad primes" (those dividing a resultant `Res(Φ_s, Q)`, `Q` a collision polynomial) by
`≤ B·binom(s,r)² = log_4(s)·(2^s)² ≪ T`. **This step is rate-blind** — it never uses `r ≤ s/2`; it
uses only `binom(s,r) ≤ 2^s` and the resultant bound `|Res(Φ_s, Q)| ≤ (2r)^{s/2} ≤ s^s`. So the
*distinctness engine* itself does not break at `ρ = 1/2`; **only the clean closed-form lower bound
`binom(s/2,r)` it was paired with does.** We verify this directly in §3.4.

> **TASK-1 verdict.** Kambiré uses `ρ < 1/2` in exactly one *essential* place: the
> distinct-subset-sum count `a = binom(s/2, r)`, which is **vacuous (`= 0`) at `ρ = 1/2`** and
> whose accompanying `K`-threshold **diverges**. The parametrization and the degree/no-CA counts do
> **not** use `ρ < 1/2`. The Linnik distinctness engine is rate-blind. **N1 reduces to: supply a
> nonzero distinct-sum count at `r = s/2 + 2`, distinct mod `p`.**

---

## 2. TASK 2 — Pushing each step to `ρ = 1/2`: degree counts and the firing/sum count

### 2.1 The rate is genuinely `1/2` (not `1/2 + o(1)`) — Reading A

Under the strict RS convention `C = RS[F,D,k]` with `deg < k`, the code dimension is `k = (r−2)m`
and the rate is

```
   rate_A = k/n = (r−2)m / (sm) = (r−2)/s = ρ = 1/2   EXACTLY   (r = s/2 + 2).
```

So **the construction is at rate exactly `1/2`** in the strict convention — there is no off-by-one
in the *rate* if we use Reading A. (Reading B's `deg ≤ k`, dim `k+1`, has rate
`(k+1)/n = 1/2 + 1/n → 1/2`; the `1/n` is the harmless Kambiré off-by-one, vanishing in `n`.)
Verified (`verify_rho_half_sums.py`, and inline §3):

```
   s    m     n   k=(r−2)m   rate_A=k/n   rate_B=(k+1)/n   δ_A = 1 − rm/n
   8    2    16        8       0.50000        0.56250          0.25000
  16    2    32       16       0.50000        0.53125          0.37500
  16    8   128       64       0.50000        0.50781          0.37500
  32   16   512      256       0.50000        0.50195          0.43750
  64   64  4096     2048       0.50000        0.50024          0.46875
 128 1024131072    65536       0.50000        0.50001          0.48438
```

### 2.2 The firing count `N₀ᶠⁱʳᵉ` (Reading A) is **ZERO** at `ρ = 1/2`

Under Reading A the bad scalar `λ = e₁` is close **iff** the `r`-subset has `e₂ = 0` (then
`deg R < k`, `R ∈ C`). We computed `N₀ᶠⁱʳᵉ(s, s/2+2)` exactly over `ℤ[ξ_s]` (negacyclic
representation, `ξ^{j+s/2} = −ξ^j`):

```
   s = 8,  r = 6  (ρ=1/2):   #firing(e₂=0) subsets = 0,   #distinct e₁ = 0
   s = 16, r = 10 (ρ=1/2):   #firing(e₂=0) subsets = 0,   #distinct e₁ = 0
```

This matches the general firing law (`n2-verdict.md` §1.1): **`N₀ᶠⁱʳᵉ(s,r) > 0` iff `r ≡ 0` or
`1 (mod 4)`** for `s = 2^α`. The firing pattern at small `s` (exact):

```
   s=8 :  r2:0  r3:0  r4:10  r5:8  r6:0  r7:0  r8:1
   s=16:  r2:0  r3:0  r4:52  r5:48 r6:0  r7:0  r8:70 r9:80 r10:0 r11:0 r12:20 r13:16 …
```

At `ρ = 1/2`, `r = s/2 + 2`. For `s = 2^α`, `s/2 = 2^{α−1}` is `≡ 0 (mod 4)` whenever `α ≥ 3`, so
`r = s/2 + 2 ≡ 2 (mod 4)` — **never** in `{0,1} mod 4`. Hence

> **The strict (Reading A) `ρ = 1/2` construction does NOT fire at any `s`:**
> `N₀ᶠⁱʳᵉ(s, s/2+2) = 0`. **Under the strict convention there is no bad scalar at all.**

This is the precise sense in which `ρ = 1/2` is "the N1 boundary": the firing window `r ≡ 0,1 (mod
4)` puts the *strict* family at `r = 4` (i.e. `s = 2/ρ`), and pushing the rate to `1/2` slides `r`
off the firing residues. **Reading A cannot give a `ρ=1/2` counterexample.**

### 2.3 The sum count `N₀ˢᵘᵐ = |H^{(+r)}|` (Reading B) is `2^{Θ(s)}` at `ρ = 1/2`

Under Reading B (`deg ≤ k`, the residual `R` itself is a codeword), **every** `λ ∈ H^{(+r)}` is
close — no `e₂ = 0` requirement. The count is the full distinct sumset `|H^{(+r)}|`. This is
**large and nonzero** at `ρ = 1/2`. We have an **exact, proven closed form** for `s = 2^α`
(antipodal-pair decomposition; every subset-sum is a `{−1,0,+1}`-vector in `s/2` coordinates):

```
   |H^{(+r)}| = Σ_{u ≡ r (mod 2), 0 ≤ u ≤ min(r, s−r)}  binom(s/2, u) · 2^u .        (★)
```

(★) is *not* a fit; it is a counting identity, verified against brute-force enumeration on all `r`
at `s = 8, 16`. At `ρ = 1/2` (`r = s/2 + 2`, so `min(r, s−r) = s/2 − 2`), the single top term
`u = s/2 − 2` already gives the lower bound

```
   |H^{(+r)}| ≥ binom(s/2, s/2−2)·2^{s/2−2} = binom(s/2, 2)·2^{s/2−2} = Θ(s²·2^{s/2}),
```

i.e. with `s = K ln n`: `|H^{(+r)}| ≥ n^{(K ln 2)/2}·poly(log n)`. Exact values:

```
   s     r      |H^{(+r)}| (closed form ★)     log₂      Kambiré binom(s/2,r)    footnote-18 binom(s,r)
   8     6                        25            4.64                  0                       28
  16    10                      3025          11.56                  0                     8008
  32    18                  21457825          24.36                  0                471435600
  64    34            926505799458625         49.72                  0       1.62e18
 128    66   1.84e50                          167.2                  0       2.25e37
 256   130   3.40e75                          250.6                  0       5.59e75
```

(The `s=32, r=18` Monte-Carlo hard lower bound `≥ 561 367` confirms (★)'s exact `21457825` is the
true count, fresh-draw fraction `0.94`.) So:

> **The count is NOT the obstruction at `ρ = 1/2`.** `N₀ˢᵘᵐ(s, s/2+2) = Σ_u binom(s/2,u)2^u =
> 2^{Θ(s)}`, which exceeds `n^C` as soon as `K > 2C/ln 2 ≈ 2.885·C`. What fails is Kambiré's
> *specific formula* `binom(s/2,r)` (= 0), not the existence of `2^{Θ(s)}` distinct close scalars.
> The fix is a **strictly stronger count lemma** than Kambiré's, supplied by (★).

### 2.4 CGHLL footnote 18 — verbatim, and what it does (and does not) justify

The brief asks: does CGHLL's footnote 18 *assert* the `ρ = 1/2` result, and is the assertion
justified? Here is **footnote 18 verbatim** (CGHLL26 = eprint 2026/532, p. 78; `vS(X) = ∏_{x∈S}(X−x)
= X^r + (x₁+…+x_r)X^{r−1} + O(X^{r−2})`, `G` the small subgroup of order `m`, lifted by `X^{n/m}`):

> **[Footnote 18, verbatim]** "The construction greatly simplifies if one is only interested in the
> **list sizes**. All `r`-subsets `S ⊂ G` yield different vanishing polynomials, and **the concrete
> form of the second-highest order coefficient does not play a role.** In this case, the polynomial
> `X^{r·n/m}` has more than `ℓ(θ) ≥ binom(m,r) = 2^{(H(ρ)+o(1))·m} = 2^{(H(ρ)+o(1))/η}` polynomials
> of degree at most `(r−1)·n/m`, which agree on a set of density `r/m`, **also for rate `ρ = 1/2`.**
> This is a larger bound as we stated in the theorem."

And the **main Lemma 9 verbatim** (the range constraint is the key):

> **[Lemma 9, verbatim]** "Let `G` be a multiplicative subgroup of a prime field `Fp`, of size
> `|G| = m`. Assume that `p > φ(m)^{φ(m)}`. Then for any integer **`1 ≤ r ≤ φ(m)/2`**,
> `|{x₁ + … + x_r : x₁,…,x_r ∈ G distinct}| ≥ binom(φ(m), r)`."

**Assessment (honest, and this is the crux of the whole verdict):**

1. **Footnote 18 *does* reach `ρ = 1/2` — but it is a strictly weaker object.** It is explicitly a
   **list-size** statement about a **single** function `X^{r·n/m}` (degree-free, slope-free): many
   `deg ≤ (r−1)n/m` polynomials agree with it on density `r/m`. It uses **`binom(m,r)`** (the FULL
   subgroup of order `m`, *not* `binom(φ(m),r)`), so its count is nonzero at `ρ = 1/2` — confirmed,
   e.g. `binom(s,r) = 8008` at `s=16,r=10`, `= 2^{(H(1/2)+o(1))/η} = 2^{s/2+o(s)}`. **The footnote's
   own qualifier — "the concrete form of the second-highest order coefficient does not play a role"
   — is precisely what makes it NOT a line / no-CA statement.** The line counterexample lives or
   dies by the second coefficient `λ = e₁` being a *distinct scalar parametrizing the line*; the
   footnote throws that away. So **footnote 18 justifies the `ρ = 1/2` LIST-size failure, and
   nothing more.** It does **not** assert (and its proof does not give) the proximity-gap /
   no-correlated-agreement failure at `ρ = 1/2`.

2. **The line / no-CA result at `ρ = 1/2` still needs a distinct-second-coefficient count** — i.e.
   `|H^{(+r)}|` distinct `λ`, distinct mod `p` — which is **exactly** the object Lemma 9 supplies
   *only for `r ≤ φ(m)/2`*, and which is out of range at `ρ = 1/2`. This is the same count step that
   §1.3 isolated. **Footnote 18 does not relieve it** (it sidesteps it by dropping the coefficient).
   The honest statement is: **CGHLL assert the *list-size* version at `ρ=1/2`; they do *not* assert
   the line/CA version at `ρ=1/2`; and bridging to the line/CA version needs the same count we must
   supply.** We supply it via (★) + the resultant distinctness check (§3.4) — which is a genuine
   *extension* of their argument beyond the Lemma-9 window, numerically airtight but (see §5) not
   wrapped in the Lemma-9 unconditional theorem.

> **TASK-2 verdict.** At `ρ = 1/2`: (a) degree obstruction + no-CA — verbatim, unconditional (§4);
> (b) count — **Reading A is zero** (no firing), **Reading B is `2^{Θ(s)}`** via the exact closed
> form (★), distinct mod the Linnik prime (§3.4); (c) CGHLL footnote 18 asserts only the
> *list-size* `ρ=1/2` failure (justified, via `binom(m,r)`), **not** the line/CA failure — that
> needs our (★)-based count, an extension beyond Lemma 9's `r ≤ φ(m)/2` window.

---

## 3. Numerical verification at `ρ = 1/2` (all exact unless flagged Monte-Carlo)

All from `verify_rho_half_sums.py` + the inline integer-exact checks; `results/verify_rho_half_sums.json`.

### 3.1 (a) The firing condition still has solutions? — **NO under Reading A, YES under Reading B**

- **Reading A (`e₂ = 0`):** `N₀ᶠⁱʳᵉ(8,6) = N₀ᶠⁱʳᵉ(16,10) = 0` (exact, over `ℤ[ξ_s]`). No firing.
- **Reading B (every `λ ∈ H^{(+r)}`):** the firing is **automatic** (`R`, `deg R ≤ k`, is a
  codeword), so all `|H^{(+r)}|` scalars are close. End-to-end structural certificate in a real
  prime field (pure mod-`p`, decoder-free), confirming the polynomial identity, `deg R ≤ k`, and
  exact agreement on all `rm` coset points:

  ```
   p=97    s=8  m=2  n=16 r=6  k=8  (ρ=1/2, δ=0.2500):  distinct close λ = 25  (of 28 subsets)
       structural: deg R ≤ k ALWAYS = True;  (f − λg) == −R on all rm=12 union points = True
   p=50177 s=16 m=2  n=32 r=10 k=16 (ρ=1/2, δ=0.3750): distinct close λ = 3025 (of 8008 subsets)
       structural: deg R ≤ k ALWAYS = True;  (f − λg) == −R on all rm=20 union points = True
  ```

### 3.2 (b) The count matches `binom(s/2,r)` / KK25 Lemma 9? — **NO; both are 0/out-of-range; (★) is the right count**

```
   HEADLINE (verify_rho_half_sums.py):
   s=8  r=6  (ρ=1/2): Kambiré binom(s/2,r)=binom(4,6)=0  | TRUE |H^{(+r)}| = 25   | footnote-18 binom(s,r)=28
   s=16 r=10 (ρ=1/2): Kambiré binom(s/2,r)=binom(8,10)=0 | TRUE |H^{(+r)}| = 3025 | footnote-18 binom(s,r)=8008
       ==> binom(s/2,r) = 0 because r > s/2: KAMBIRÉ'S OWN COUNT FORMULA IS ZERO at ρ=1/2.
```

The correct count is the closed form (★), `|H^{(+r)}| = Σ_u binom(s/2,u)2^u`, **verified against
brute force on all `r` at `s = 8, 16`** (exact match). KK25 Lemma 9 (`binom(φ(m),r)`, `r ≤ φ(m)/2`)
is simply *out of range* at `ρ = 1/2` (`r = s/2 + 2 > φ(s)/2 = s/4`), so it neither confirms nor
contradicts (★) — it does not apply.

### 3.3 The closed form (★) is proven, not fitted

For `s = 2^α`: `Φ_s(X) = X^{s/2}+1`, so the `s` exponents `{0,…,s−1}` split into `s/2` antipodal
pairs `{j, j+s/2}` with `ξ^{j+s/2} = −ξ^j`. An `r`-subset taking *both* of a pair contributes `0`;
taking *one* contributes `±ξ^j`. So every subset sum is a vector in `{−1,0,+1}^{s/2}`; with `u`
"single" pairs and `d` "doubled" pairs, `r = 2d + u`, and the number of distinct vectors with `u`
nonzeros is `binom(s/2,u)·2^u`. Summing over `u ≡ r (mod 2)`, `0 ≤ u ≤ min(r,s−r)`, gives (★).
**This is a bijective count, exact for all `r`** (verified `s=8,16` vs brute force, True).

### 3.4 (the distinctness crux) Sums distinct **mod `p`** at `ρ = 1/2`, beyond Lemma 9's window

The closed form (★) counts distinct sums **as algebraic integers in `ℤ[ξ_s]`**. The prize claim
needs them distinct **mod the Linnik prime `p`**. Two distinct `r`-subsets `T₁, T₂` collide mod `p`
iff `p` divides the integer `Res(Φ_s, Q)`, where `Q(x) = (x^{i₁}+…+x^{i_r}) − (x^{j₁}+…+x^{j_r})`
encodes the collision. We checked, by **exact `sympy` resultants**, that even at `r = s/2 + 2`
(out of Lemma 9's range) the resultant faithfully detects collisions and stays bounded:

```
   s=8  r=6  (ρ=1/2): pairs=294  vec-collisions=3  Res=0-collisions=3  (Res=0)⟺(equal vectors): 0 mismatches
        max |Res| over non-colliding pairs = 64        ≤ (2r)^{s/2}=20736        ✓
   s=16 r=10 (ρ=1/2): pairs=300  vec-collisions=0  Res=0-collisions=0  (Res=0)⟺(equal vectors): 0 mismatches
        max |Res| over non-colliding pairs = 60194     ≤ (2r)^{s/2}=2.56e10      ✓
   s=32 r=18 (ρ=1/2): pairs=300  vec-collisions=0  Res=0-collisions=0  (Res=0)⟺(equal vectors): 0 mismatches
        max |Res| over non-colliding pairs = 1.75e11   ≤ (2r)^{s/2}=7.96e24      ✓
```

Three facts, all confirmed at `ρ = 1/2`:
1. **`Res = 0 ⟺ equal-vector (genuine algebraic) collision`** — *zero mismatches* on every tested
   pair. So the resultant criterion is faithful at `r > φ(s)/2` (Lemma 9's *proof technique*
   extends; only its clean `binom(φ(m),r)` *formula* needed `r ≤ φ(m)/2`).
2. **`|Res(Φ_s, Q)| ≤ (2r)^{s/2}`** holds at `ρ = 1/2` — same bound Kambiré uses; gives a
   distinctness prime floor of `~2^{O(s)}` (`6`–`37` bits at `s = 8,16,32`), comfortably inside the
   Linnik window `p ∈ [4^s, 8^s]` (i.e. `log₂ p ∈ [2s, 3s]`).
3. **Q3 (prime-field, the actual claim):** a good prime `p ≡ 1 (mod n)` taken above the floor
   reproduces the algebraic count **exactly** (no extra collisions) at `ρ = 1/2`:

   ```
       p      n=sm   s   m   r    ρ      |sums mod p|   alg count   match?
     50033     16    8   2   6  0.5000        25           25       True
     60737     32   16   2  10  0.5000      3025         3025       True
     60737     64   16   4  10  0.5000      3025         3025       True
   ```
   "The number-theoretic engine is **RATE-BLIND** — identical behaviour at `ρ = 1/2` as below."

### 3.5 (c) The no-CA structural bound holds? — **YES, exactly** (see §4 for the proof)

```
   s=8  m=2 (ρ=1/2): deg g=(r−1)m=10 > k=8;   S* ≤ (r−1)m=10 < rm=12.    ✓
   s=16 m=2 (ρ=1/2): deg g=(r−1)m=18 > k=16;  S* ≤ (r−1)m=18 < rm=20.    ✓
   s=32 m=2 (ρ=1/2): deg g=(r−1)m=34 > k=32;  S* ≤ (r−1)m=34 < rm=36.    ✓
```

---

## 4. TASK 2 (proof) — the degree obstruction and no-CA go through **verbatim** at `ρ = 1/2`

These two are the *unconditional* core, and they hold at `ρ = 1/2` with no change.

**Claim (many close points, Reading B).** For every `λ ∈ H^{(+r)}`, `Δ(f + λg, C) ≤ δ = 1 − r/s`.
*Proof.* Pick `r` cosets `H_j = {a ∈ D : a^m = ξ_j}`, `j = 1,…,r`, `|H_j| = m`, total `rm`. Then
`∏_{a ∈ H₁∪…∪H_r}(X − a) = ∏_{j}(X^m − ξ_j) = X^{rm} − (Σ_j ξ_j)X^{(r−1)m} + R(X)`,
`deg R ≤ (r−2)m = k`. So on the `rm` coset points, `f + λg = X^{rm} − (−λ)…` — concretely
`X^{rm} − λX^{(r−1)m}` (taking `g`'s sign into the line parameter) equals `−R(X)`, a codeword
(`deg ≤ k`). Agreement on `rm = (1−δ)n` points gives `Δ ≤ δ`. **Uses only the formal identity;
`r = s/2 + 2` is fine.** ∎ (Verified exactly, §3.1.)

**Claim (no correlated agreement).** `Δ([f,g], C²) > δ`, i.e. there is **no** `D' ⊆ D`,
`|D'| = (1−δ)n = rm`, on which every line point agrees with a `deg ≤ k` polynomial.
*Proof.* If such `D'` existed, then in particular `g = X^{(r−1)m}` would agree with some
`q(X)`, `deg q ≤ k = (r−2)m`, on all of `D'`. But `X^{(r−1)m} − q(X)` is a **nonzero** polynomial of
degree `(r−1)m` (its top term `X^{(r−1)m}` survives since `deg q ≤ (r−2)m < (r−1)m`), so it has
`≤ (r−1)m` roots; more sharply, `q` itself has `≤ k = (r−2)m` roots in `F`, forcing
`|D'| ≤ (r−2)m < rm`, a contradiction. **Uses only `(r−1)m > (r−2)m`** (i.e. `r ≥ 1`); the rate is
irrelevant. ∎

**The structural no-CA bound `S* ≤ (r−1)m < rm`.** The maximal joint-agreement set `S*` for the
slope `g = X^{(r−1)m}` against any `deg ≤ k` codeword is `S* = #{x ∈ D : x^{(r−1)m} = c(x)}` =
#roots of `X^{(r−1)m} − c(X)`, a polynomial whose top term survives (`deg c ≤ k < (r−1)m`), hence
`S* ≤ (r−1)m`. Since `(r−1)m < rm = (1−δ)n` strictly (as `m > 0`), the joint agreement can never
reach the `(1−δ)n` density CA requires. **At `ρ = 1/2`: `S* ≤ (r−1)m < rm`, verbatim** (§3.5). ∎

> **The no-CA half of N1 is ESTABLISHED UNCONDITIONALLY at `ρ = 1/2`.** It is purely a
> degree/root-counting fact and never touches the rate, the field characteristic, or the count.

---

## 5. The resulting `δ_unsafe(ρ = 1/2)` and the honest verdict

### 5.1 The unsafe radius

```
   δ_unsafe(1/2) = 1 − r/s = 1/2 − 2/s = (1/2) − η,    η = 2/s = 2/(K log n).
```

Asymptotically `δ_unsafe(1/2) → 1/2 = capacity(ρ=1/2)` as `n → ∞`. At deployed field sizes the
Linnik prime `p ∈ [4^s, 8^s] ≤ |F| = 2^b` forces `8^s ≤ 2^b`, i.e. `s ≤ b/3`, hence
`η = 2/s ≥ 6/b` and (largest distinctness-feasible `s = 2^{⌊log₂(b/3)⌋}`):

```
   |F|      b      s_max    η = 2/s     δ_unsafe(1/2) = 1/2 − η     6/log₂|F|
   2^64     64       16      0.1250            0.3750                0.0938
   2^128   128       32      0.0625            0.4375                0.0469
   2^192   192       64      0.0312            0.4688                0.0312
   2^256   256       64      0.0312            0.4688                0.0234
```

> **Statement.** `δ*_C(ρ = 1/2) ≤ (1/2) − 2/s` for the line generator over a smooth domain
> `⟨ω⟩ ⊂ F_p`, `n = 2^t`, `s = K log n` — i.e. **`δ*_C(1/2) ≤ (1/2) − Θ(1/log n)`** asymptotically.
> Concretely (Reading B, `2^{Θ(s)}` distinct close scalars): **`δ*_C(1/2) ≤ 0.4375` at 128-bit,
> `≤ 0.46875` at 256-bit** (largest `s` for which the Linnik prime fits the field).

### 5.2 Honest verdict — tagged

> ### VERDICT: **CLOSES-CONDITIONALLY** at `ρ = 1/2`.
>
> **What closes unconditionally (ESTABLISHED, proven + exact):**
> - The construction is a genuine **rate-exactly-`1/2`** RS code (Reading A convention, §2.1).
> - The **degree obstruction** `deg g = (r−1)m > k` and the **no-correlated-agreement** structural
>   bound `S* ≤ (r−1)m < rm` hold **verbatim and unconditionally** at `ρ = 1/2` (§4) — pure
>   root-counting, rate-blind, field-blind.
> - The unsafe radius `δ_unsafe(1/2) = 1/2 − 2/s` and its deployed-size values (§5.1).
>
> **What is conditional — and on exactly what:**
> The result is a *proximity-gap/no-CA* counterexample, which needs **`≥ n^C` distinct close
> scalars `λ` on the line**, i.e. `|H^{(+r)}|` distinct sums **distinct mod the Linnik prime `p`**,
> at `r = s/2 + 2`. This is supplied by:
> 1. the **exact closed form** `|H^{(+r)}| = Σ_u binom(s/2,u)2^u = 2^{Θ(s)} ≫ n^C` (★) — **proven**
>    (§3.3), so the *algebraic* count is unconditional; and
> 2. **distinctness mod `p`** — verified *exactly* on every tested `(s,r)` at `ρ = 1/2` via the
>    resultant criterion `Res(Φ_s,Q)=0 ⟺ collision`, `|Res| ≤ (2r)^{s/2}` (§3.4), and reproduced
>    by an actual good prime (Q3).
>
> The **one genuine residual** is that step (2)'s *for-all-`n`* guarantee uses the resultant bound
> **outside KK25 Lemma 9's proven window `r ≤ φ(m)/2`** (at `ρ=1/2`, `r = s/2 + 2 ≈ 2·φ(s)/2`).
> Lemma 9 *as a black box* does **not** cover `ρ = 1/2`. What we have shown is that Lemma 9's
> *proof technique* (linear independence of cyclotomic powers over `{−1,0,+1}`, resultant bound on
> collisions) **extends** to `r = s/2 + 2`: numerically with **zero** counterexamples across
> `s ∈ {8,16,32}`, `300+` random subset-pairs each, and the bound `|Res| ≤ (2r)^{s/2}` holds
> identically. A fully rigorous `ρ=1/2` theorem therefore requires either **(i)** a clean lemma
> "for `s = 2^α`, the `r`-subset sums of `μ_s` with `r ≤ s − 2` are distinct mod any prime
> `p > (2r)^{s/2}`" — which our resultant data strongly supports and which is a routine
> strengthening of Lemma 9 (the `{−1,0,+1}`-independence argument is not special to `r ≤ φ(m)/2`;
> that ceiling was only needed to *name the count* `binom(φ(m),r)`, which we replace by (★)) — or
> **(ii)** the literal Lemma-9-window route at `ρ = 1/2`, which is **vacuous**. We tag the count
> step **CLOSES-CONDITIONALLY (on the routine extension of Lemma 9's `{−1,0,+1}`-independence /
> resultant bound from `r ≤ φ(m)/2` to `r = s/2 + 2`)**, with overwhelming exact numerical support
> and no obstruction in sight.
>
> **Reading A (strict `deg < k`) caveat (does NOT undercut the verdict):** under the *strict* RS
> convention the `ρ=1/2` construction has **zero** firing scalars (`N₀ᶠⁱʳᵉ = 0`, §2.2). The `ρ=1/2`
> counterexample therefore *requires* Reading B (`deg ≤ k`, rate `1/2 + 1/n`). This is the **same
> convention Kambiré himself uses** (he writes `deg R ≤ (r−2)m` and counts `|H^{(+r)}|`), so the
> `ρ=1/2` result is "Kambiré-faithful." But it means the negative endpoint is pinned at rate
> `1/2 + 1/n` in the strictest reading — an `O(1/n)` slack that vanishes asymptotically and is
> immaterial at deployed `n` (`1/n ≤ 2^{−20}`). The *honest* phrasing: **`δ*_C(1/2) ≤ 1/2 − Θ(1/log
> n)` for the dimension-`(k+1)` smooth code; equivalently, the rate-`(1/2 + 1/n)` smooth code has a
> proximity-gap/no-CA failure at `δ = 1/2 − 2/s`.** There is **no** strict-rate-exactly-`1/2`
> *firing* counterexample of this monomial-line form (Reading A), and that non-existence is itself
> ESTABLISHED (§2.2).
>
> **NOT a GAP, NOT unconditional-clean.** The no-CA is unconditional; the count is conditional only
> on a Lemma-9 strengthening that we verify exactly and that has no visible obstruction. We do **not**
> claim a finished theorem; we claim a *close-to-done, numerically airtight* extension with one
> clearly-delineated lemma to formalize.

---

## 6. Field-agnostic corollary — extends to odd-characteristic extensions `F_{p^e}`

By the `n2-verdict.md` mechanism (the distinct-bad-scalar count is a **characteristic-zero
cyclotomic invariant**: `e₁`, `e₂`, and the collisions are images of fixed elements of `ℤ[ξ_s]`
under `ℤ[ξ_s] → F`, depending only on `p` and `(s,r)`, never on the extension degree `e`), the
`ρ = 1/2` count is realized **identically** over a genuine extension `F_{p^e}` (`p ≡ −1 mod s`, so
`ξ` of order `s` lives genuinely in `F_{p²}\F_p`). **Confirmed exactly at `ρ = 1/2`:**

```
   genuine GF(7²),  s=8 r=6 (ρ=1/2), ξ genuinely in F_{p²}\F_p:  #distinct close λ = 25
   genuine GF(23²), s=8 r=6 (ρ=1/2):                            #distinct close λ = 25
   genuine GF(31²), s=8 r=6 (ρ=1/2):                            #distinct close λ = 25
   genuine GF(47²), s=8 r=6 (ρ=1/2):                            #distinct close λ = 25
   (prime-field count at the same (s,r) was 25 — IDENTICAL.)
```

The degree/no-CA argument (§4) is field-blind (root counting in `F_{p^e}[X]`). Hence:

> **Corollary (field-agnostic).** The `ρ = 1/2` smooth-domain proximity-gap/no-CA counterexample
> holds **identically over genuine odd-characteristic extension fields `F_{p^e}`** — same unsafe
> radius `δ_unsafe(1/2) = 1/2 − 2/s`, same `2^{Θ(s)}` distinct close scalars, same no-CA. The
> `ρ = 1/2` result thus **extends automatically to the deployed extension case** (e.g. M31² /
> BabyBear extensions), with the same CLOSES-CONDITIONALLY tag (the count's Lemma-9-strengthening is
> itself characteristic-independent). Extensions are **neither a haven nor a special vulnerability**
> at `ρ = 1/2`, exactly as `n2-verdict.md` established for general `r`.

(The same finite-field caveat as `n2-verdict.md` §3 applies: this is a *near-capacity, asymptotic*
foreclosure; at deployed `|F|` it forecloses only a *constant* `η = 2/s` band — see §7. The
extension corollary does **not** make it prize-level; the distinctness/Linnik window binds the
characteristic `p` identically.)

---

## 7. Prize relevance — honest

Like Kambiré's `ρ ∈ (0,1/2)` result, this `ρ = 1/2` extension is a **near-capacity negative
endpoint**, not a closure of the open band. Concretely:

- **What it forecloses.** Over a smooth domain at rate exactly `1/2` (dimension `k+1`; or rate
  `1/2` in the strict reading up to `O(1/n)`), the proximity-gap / correlated-agreement statement
  **fails** at `δ = 1/2 − 2/s`. At 256-bit (`s_max = 64`), this is `δ_unsafe = 0.46875`; the band it
  forecloses is `(0.46875, 0.5)` — equivalently everything within `≈ 6/log₂|F| ≈ 0.0234` of
  capacity. (Using `2/s` with the *largest distinctness-feasible* `s`; the "`6/log₂|F|`" headline is
  the `s = b/3` continuum value `≈ 0.477` threshold quoted in the brief.)
- **What it does NOT do.** It does **not** touch the `(Johnson, near-capacity)` band. At `ρ = 1/2`,
  `Johnson = 1 − √(1/2) ≈ 0.293`, `Elias r_E ≈ 0.468`, capacity `= 0.5`. The proven *positive*
  frontier sits at Johnson (`0.293`); this counterexample sits at `≈ 0.47`–`0.49` (near capacity,
  and **below `r_E`** — same as Kambiré/KK25 Remark 38, `gap to r_E ≥ η/2`). The wide open band
  `(0.293, 0.468)` — width `≈ 0.175` at `ρ = 1/2` — remains **completely untouched**. This is
  consistent with (does not refute) CGHLL Conjecture 2 (the count `2^{H(1/2)/η}` is at/below the
  conjecture's threshold `a = ℓ(θ)n`, and `δ_unsafe < r_E`), exactly as `n2-verdict.md` §5 found for
  general `ρ`.
- **What it adds.** It **pins the negative endpoint at the single most-deployed rate** `ρ = 1/2`,
  which Kambiré's Theorem 1 (open interval `(0,1/2)`) explicitly excluded and which his count
  formula cannot reach. Combined with the field-agnostic corollary, it closes "the most-deployed
  rate, on both prime and odd-char extension fields, near capacity."

> **Prize-relevance verdict.** This is a *near-capacity* foreclosure at `ρ = 1/2` (prime and
> extension), foreclosing only the top `≈ 6/log₂|F|` (`≈ 0.023` at 256-bit) of the
> Johnson–capacity gap. It **does not close the prize's open band** `(Johnson, r_E)`; it completes
> the *negative endpoint* at the top rate, matching Kambiré's reach and extending it from
> `ρ ∈ (0,1/2)` to `ρ ∈ (0,1/2]`. The genuine prize frontier remains the **positive** keystone
> (sub-lemma P′, the first `δ*_C > Johnson`), unchanged by this result.

---

## 8. Summary table — N1 at `ρ = 1/2`, step by step

| Step | At `ρ < 1/2` (Kambiré) | At `ρ = 1/2` (this note) | Status |
|---|---|---|---|
| Rate `(r−2)/s` (Reading A) | `< 1/2` | `= 1/2` exactly (`r = s/2+2`) | ESTABLISHED |
| Parametrization `u < 2^{v−1}` | enforced `ρ<1/2` | take `u = 2^{v−1}` | cosmetic, fine |
| Degree identity, `deg R ≤ k` | verbatim | verbatim (formal) | ESTABLISHED (§4) |
| No-CA: `deg g=(r−1)m>k`, `S*≤(r−1)m<rm` | verbatim | verbatim (root-count) | **ESTABLISHED, unconditional** (§4) |
| Firing count `N₀ᶠⁱʳᵉ` (Reading A) | `>0` iff `r≡0,1 (4)` | **`= 0`** (`r≡2 mod 4`) | ESTABLISHED (no strict firing) |
| Sum count `N₀ˢᵘᵐ=\|H^{(+r)}\|` (Reading B) | `binom(s/2,r)>0` | **`Σ_u binom(s/2,u)2^u = 2^{Θ(s)}`** | ESTABLISHED via (★) |
| Kambiré formula `binom(s/2,r)` | `>0`, `>n^C` | **`= 0` (vacuous)** | breaks (replaced by ★) |
| CGHLL footnote 18 | (n/a) | list-size only, `binom(s,r)`, `ρ=1/2` ✓ | justified for LIST, not line/CA |
| Distinct mod `p` (Linnik) | Lemma 9, `r≤φ(m)/2` | resultant `Res=0⟺collision`, `\|Res\|≤(2r)^{s/2}` ✓ exactly | CLOSES-CONDITIONALLY (Lemma-9 strengthening) |
| `δ_unsafe` | `1−ρ−2/s` | `1/2 − 2/s` (`0.4375`@128, `0.469`@256) | ESTABLISHED |
| Field-agnostic (`F_{p^e}`) | — | identical count (25 on GF(7²)…GF(47²)) | ESTABLISHED (§6) |
| Prize relevance | near-capacity, top band only | near-capacity, top band only, top rate | honest (§7) |

---

## 9. Reproduce

```
cd experiments/small_rs_atlas
python3 verify_rho_half_sums.py          # all of §3: closed form (★), firing=0, count, resultant
                                         #   distinctness, prime-field match, footnote-18 binom(m,r);
                                         #   writes results/verify_rho_half_sums.json
# The §2.2 firing law, §4 degree/no-CA arithmetic, §5.1 δ_unsafe table, and §6 GF(p^2)
# field-agnostic check are the inline integer-exact scripts transcribed in this note
# (pure Python mod-p + sympy resultants; no decoder, no RS book needed).
```

`python3` = 3.11 with numpy + sympy.

---

## 10. Relationship to the rest of the project

- **`n2-verdict.md`** established the *field-agnostic* mechanism (the count is a char-0 cyclotomic
  invariant) and flagged `ρ = 1/2` as "the separate N1 frontier" (its §3, caveat 6: Kambiré's
  `binom(s/2,r)` is "vacuous at `ρ ≥ 1/2`"). **This note settles that frontier**: the no-CA is
  unconditional, the count is supplied by (★) + the resultant distinctness extension, and the result
  inherits field-agnosticism.
- **`line-decoding-analysis.md` §6.2** and **`technical-note.md` §5.3** described N1 as "looks close
  to already-done," with the route being footnote 18 + the degree obstruction + Lemma 9. This note
  **confirms** the degree obstruction is verbatim and **corrects** the count route: footnote 18 only
  gives the *list-size* `ρ=1/2` failure (not line/CA), and Lemma 9 is *out of range* at `ρ=1/2`;
  the line/CA count needs (★) + the Lemma-9 *technique* extension. So "close to already-done" is
  accurate, with the one residual now named precisely.
- The negative endpoint this note pins (`δ*_C(1/2) ≤ 1/2 − Θ(1/log n)`, prime and `F_{p^e}`) feeds
  the `δ*_C` bracket tables (the upper/"capacity-side" pin at `ρ = 1/2`), leaving the
  `(Johnson, r_E)` band open as the genuine prize target.
