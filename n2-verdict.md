# N2 — Definitive verdict: does the near-capacity proximity-gap counterexample extend to odd-characteristic extension fields?

> **Status of this document.** This is the settling verdict for sub-lemma **N2** after the
> flip-and-reflip history (Wave-4: "does NOT extend" → Wave-5 + orchestrator probe: that
> was an `r=3`-only artifact). It is written **fail-closed**: a claim is tagged
> **ESTABLISHED** only when backed by proof or exact computation; the central
> field-agnosticity claim is tagged **ESTABLISHED-MODULO-FORMALIZATION** — captured by the
> Theorem of §2.1 (cyclotomic-invariant argument + exact computation as verification), with
> the one open formalization step being the explicit finiteness bound for the bad-prime set
> `B(s,r)`; weaker claims are **LIKELY** or **OPEN**. Date: 2026-06-03. All numbers below are reproducible from
> `experiments/small_rs_atlas/{n2_char0_count.py, n2_count_laws.py, n2_prize_threshold.py,
> n2_crosscheck.py}` (+ the existing `n2_hardening.py`, `counterexample_extension.py`),
> with `results/{n2_prize_threshold.json, n2_crosscheck.json}`.

---

## 0. One-paragraph bottom line

The counterexample **mechanism is field-agnostic and extends structurally to genuine
odd-characteristic extension fields** (ESTABLISHED-MODULO-FORMALIZATION — the
field-agnosticity **Theorem of §2.1**, with the exact computation as supplementary
verification; the one formalization gap is the explicit finiteness bound for the bad-prime
set `B(s,r)`). The number of
distinct bad scalars is a **characteristic-zero cyclotomic quantity** — identical for a
prime field and a genuine extension `GF(p^2)` up to finite-field saturation — so "do
extensions get fewer?" is answered **no** (they match the prime; sometimes the extension
is *less* saturated, never systematically fewer). However, the extension counterexample is
**NOT prize-level at 256-bit** (no `(s,r)` with `s ≤ 2^40`, `ρ ∈ {1/2,1/4,1/8,1/16}` makes
the distinct-bad-scalar count exceed `2^-128·|F|` over a `≥ 141`-bit field — the field
large enough for the subset sums to stay distinct is too small to hold `2^128` of them).
It **is** "prize-level" at **128-bit** in the trivial sense that the threshold there is
`2^0 = 1`. And — the decisive correction to the prior framing — the extension finding
**does NOT refute CGHLL Conjecture 2's extension clause**: the failure radius
`δ = capacity − 2/s` sits *below* the Elias radius `r_E`, and the bad-scalar count
(`≈ 2^{H(2ρ)/(2η)}`) is *below* Conjecture 2's exponentially-large line-decoding threshold
`a = ℓ(θ)·n = 2^{H(ρ)/η}·n`, so it is **consistent with** (indeed calibrated to) the
conjecture, on extensions exactly as on primes. What it refutes is the *naive
up-to-MDS-capacity* proximity-gap statement and the *hope that odd-char extensions are
safer than primes* — both now closed on the negative side over extensions too.

---

## 1. The firing condition and the two count conventions (exact, decisive)

Smooth domain `C = RS[F, L, k]`, `D = ⟨ω⟩` of order `n = sm`, inner subgroup `H = ⟨ξ⟩` of
order `s`, `k = (r−2)m`, `ρ = (r−2)/s`, line `f = X^{rm}`, `g = X^{(r−1)m}`,
`δ = 1 − r/s = capacity − 2/s`. For an `r`-subset `{ξ₁,…,ξ_r} ⊆ H`, the vanishing-polynomial
identity gives

```
∏_j (X^m − ξ_j) = X^{rm} − e₁·X^{(r−1)m} + e₂·X^{(r−2)m} − e₃·X^{(r−3)m} + …
R(X) := ∏ − (X^{rm} − λ·X^{(r−1)m}),   λ = e₁ = Σξ_i,   deg R ≤ (r−2)m = k,
leading term of R = e₂·X^{(r−2)m}  (degree exactly k unless e₂ = 0).
```

A bad scalar `λ` makes `X^{rm} − λX^{(r−1)m}` agree with a **codeword** on the `rm`-point
coset union, i.e. `Δ ≤ δn`. Whether `R` is a codeword depends on the RS dimension
convention — and the two readings give two genuinely different counts (**both exactly
verified**, `n2_crosscheck.py`/`n2_count_laws.py`):

| Reading | RS code | bad-scalar `λ` is close iff | distinct-bad-scalar count |
|---|---|---|---|
| **A (strict)** | `deg < k` (dimension `k`, rate `ρ=(r−2)/s`) | some `r`-subset summing to `λ` has **`e₂ = 0`** (firing) | `N₀ᶠⁱʳᵉ(s,r)` = #distinct `e₁` over `e₂=0` subsets |
| **B (Kambiré-literal)** | `deg ≤ k` (dimension `k+1`, rate `(r−2)/s + 1/n`) | **always** (`R` itself is a codeword) | `N₀ˢᵘᵐ(s,r) = \|H^{(+r)}\|` = #distinct `r`-fold subset sums |

**Exactly verified** (prime `GF(17)` and genuine extension `GF(7²)`, full
`rs.dist_to_code` enumeration): under Reading A, `{distinct close λ} = {distinct firing
λ} = N₀ᶠⁱʳᵉ`, with **zero** close `λ` outside `H^{(+r)}`; under Reading B (code `deg ≤ k`),
**all** `|H^{(+r)}|` of them are close. So Kambiré's note carries a harmless off-by-one
(`δ`, `ρ` differ by `O(1/n)` between the two), but the *count* differs by a large factor.
The brief's "firing condition `deg R < k ⟺ e₂ = 0`" is **Reading A**; Kambiré's claimed
`binom(s/2,r)` bad scalars is **Reading B**. We report both.

### 1.1 Why the prior `r=3` flip was an artifact

`N₀ᶠⁱʳᵉ(s,r) > 0` **iff `r ≡ 0` or `1 (mod 4)`** for `s = 2^t` (exactly verified, `s=16,32`;
fires at `r = 4,5,8,9,12,13,25,28,29`, zero at `r = 3,6,7,10,11,…`). In particular
**`N₀ᶠⁱʳᵉ(s,3) = 0` in characteristic zero — for primes too.** Wave-4 tested only `r=3`
and saw firing **only in tiny fields** (saturation "bad-prime" artifacts), concluding
"extension does not extend". The truth: at `r=3` the strict construction does not fire
*generically in any field*. The genuine signal lives at `r ≡ 0,1 (mod 4)` (e.g. `r=4`),
and there it is field-agnostic — which is exactly what Wave-5 + the orchestrator probe
found at `r=4`.

---

## 2. General-`r` count — extension vs prime (TASK 1 crux)

### 2.1 The count is a characteristic-zero cyclotomic invariant (the mechanism)

Both `e₁` and `e₂` are images of **fixed elements of the cyclotomic ring `ℤ[ξ_s]`** under
the reduction `ℤ[ξ_s] → F`. Two firing subsets collide (same `λ`) in `F` iff the
characteristic `p` divides an integer resultant-type quantity (`Res(Φ_s, ·)`) measuring
their difference in `ℤ[ξ_s]` — a quantity that depends **only on `p` and `(s,r)`, never on
the extension degree `e`**. A primitive `s`-th root of unity exists in `F_p` iff
`p ≡ 1 (mod s)` and **genuinely in `F_{p²}\F_p` iff `p ≡ −1 (mod s)`** (its order `s` does
not divide `p−1`). In both cases the reduction is injective on the relevant differences
once `p` is large, so **both realize the same characteristic-zero count
`N₀(s,r) := #distinct e₁ in ℤ[ξ_s]`.** (`s = 2^t` ⇒ `Φ_s = X^{s/2}+1`, so each `ξ^j` is a
signed unit vector and `N₀` is computed exactly with integer arithmetic; see
`n2_char0_count.py`.)

This is the mechanism; we now state it as a theorem (with the exact computation below
becoming supplementary verification of it).

> **Theorem (Field-agnosticity of the Kambiré subset-sum mechanism).** Fix `s = 2^t` and
> `r` in the admissible range. There is an explicit finite set of "bad primes" `B(s,r)` —
> the primes dividing the relevant cyclotomic resultants `Res(Φ_s, ·)` governing `e₁/e₂`
> collisions — such that for every prime `p ∉ B(s,r)`, the number of distinct bad scalars
> (firing scalars with no correlated agreement) over `F_p` equals that over every genuine
> odd-characteristic extension `F_{p^e}`.

**PROOF SKETCH (the cyclotomic-invariant argument).** Fix `s = 2^t` and an admissible `r`.
The firing data are the elementary symmetric functions `e₁ = Σ_{i} ξ_{j_i}` and
`e₂ = Σ_{i<i'} ξ_{j_i}ξ_{j_{i'}}` of an `r`-subset of the `s`-th roots of unity, i.e.
`e₁, e₂` are the images of **fixed elements of `ℤ[ξ_s]`** (independent of the field) under
the unique ring map `ℤ[ξ_s] → F` sending `ξ_s` to a chosen primitive `s`-th root of unity
in `F`. A degeneracy or collision — two distinct firing subsets giving the same `λ = e₁`,
or a subset becoming non-firing (`e₂ = 0`) — is the **vanishing in `F` of a fixed difference
element `D ∈ ℤ[ξ_s]`** (a difference `e₁−e₁'`, or an `e₂`). Since `ξ_s` is a root of the
fixed cyclotomic polynomial `Φ_s` (`= X^{s/2}+1` for `s = 2^t`), `D` vanishes under the
reduction iff the characteristic `p` divides the integer **resultant** `Res(Φ_s, D̃)` (where
`D̃` is the fixed polynomial representing `D`) — equivalently, iff `p` divides the norm
`N_{ℚ(ξ_s)/ℚ}(D)`, a **fixed nonzero integer** depending only on `s` and `(s,r)`. Crucially,
this divisibility condition is a statement about `p` **alone**: it does not see the extension
degree `e`, because the entire computation takes place inside the single ring `ℤ[ξ_s]` and
its reduction is determined by the image of `ξ_s` (which exists in `F_p` when `p ≡ 1 mod s`
and genuinely in `F_{p^e}` — e.g. `F_{p²}`, `p ≡ −1 mod s` — when `s ∤ p−1`), not by which
field `F_{p^e}` we land in. Collect into `B(s,r)` every prime dividing any of the finitely
many such fixed nonzero norms (over the finitely many subset pairs / `e₂` elements). For
`p ∉ B(s,r)` **no** fixed difference vanishes, so the firing/collision pattern over `F_p` and
over every `F_{p^e}` is literally the *same combinatorial pattern lifted from `ℤ[ξ_s]`* —
hence the distinct-bad-scalar count is a **characteristic-only invariant**, identical across
all genuine extensions of `F_p`, for every `p` off the finite bad-prime set. ∎(sketch)

**What remains to fully formalize (HONEST).** The argument above reduces field-agnosticity to
one residual: an **explicit description and finiteness bound for `B(s,r)`** — i.e. a clean
statement of the resultant/norm condition `p ∤ Res(Φ_s, ·)` together with a bound on the bad
primes (e.g. `|B(s,r)|` and `max B(s,r)` in terms of `s, r`, matching the empirically observed
saturation "bad-prime" artifacts of §1.1 and the field-size saturation of §2.2). The
*structure* of the argument (vanishing of a fixed `ℤ[ξ_s]`-element ⇔ `p | Res`, `e`-independent)
is complete; the missing piece is the quantitative resultant bound delimiting `B(s,r)`. We
therefore tag the field-agnosticity **ESTABLISHED-MODULO-FORMALIZATION** (not bare PROVEN):
the mechanism and its `e`-independence are established, the bad-prime set is described but its
explicit finiteness bound is the formalization gap.

**The exact computation is VERIFICATION supporting the Theorem.** The tables of §2.2 — produced
by `n2_char0_count.py` / `n2_general_r_count.py` (the Newton/elementary-symmetric DP), each
**validated against brute-force enumeration**, and showing `char0 == prime(ξ∈F_p) ==
genuine-ext(ξ∈F_{p²})` exactly in every unsaturated case — are *supplementary verification* of
this Theorem at concrete `(s,r,p)`, not the primary evidence: they confirm the count is the
predicted characteristic-zero invariant off the bad-prime set (the only deviations being the
field-SIZE saturation of §2.2, which hits prime and extension comparably, exactly as the
finite-`B(s,r)` picture predicts).

### 2.2 The exact tables (extension MATCHES prime)

Smart counting throughout: a **Newton/elementary-symmetric DP** on `(e₁,e₂)` field-element
states (adjoining `x` sends `(e₁,e₂) → (e₁+x, e₂+x·e₁)`), width `≤ min(binom(s,k), q²)`,
validated against brute force — **never** a `binom(s,r)` enumeration. Char-0 done with the
negacyclic cyclotomic DP/brute.

**`N₀ᶠⁱʳᵉ` (Reading A) — EXACT, all match `char0 == prime(ξ∈F_p) == genuine-ext(ξ∈F_{p²})`:**

| `s` | `r` | `N₀ᶠⁱʳᵉ` char-0 | prime (`p≡1`) | genuine ext `GF(p²)` (`p≡−1`) | match | note |
|---|---|---|---|---|---|---|
| 8 | 4 | 9 | 9 | 9 | ✓ | `=(s/2−1)²` |
| 8 | 5 | 8 | 8 | 8 | ✓ | |
| 16 | 4 | 49 | 49 | 49 | ✓ | `=(s/2−1)²` |
| 16 | 5 | 16 | 16 | 16 | ✓ | |
| 16 | 8 | 49 | 49 | 49 | ✓ | |
| 16 | 9 | 48 | 48 | 48 | ✓ | |
| 16 | 12 | 17 | 17 | 17 | ✓ | |
| 16 | 13 | 16 | 16 | 16 | ✓ | |
| 32 | 4 | 225 | 225 | 225 | ✓ | `=(s/2−1)²` |
| 32 | 5 | 32 | 32 | 32 | ✓ | |
| 32 | 28 | 33 | 33 | 33 | ✓ | |
| 32 | 29 | 32 | 32 | 32 | ✓ | |

`N₀ᶠⁱʳᵉ` **matches exactly in every case** (it is small enough never to saturate any of the
fields used). Closed form at `r=4`: `N₀ᶠⁱʳᵉ(s,4) = (s/2−1)²` (verified `s=8,16,32`).

**`N₀ˢᵘᵐ = |H^{(+r)}|` (Reading B) — matches when the field is large enough; the only
deviations are field-SIZE saturation, hitting prime and extension comparably:**

| `s` | `r` | `N₀ˢᵘᵐ` char-0 | prime (large `p≡1`) | genuine ext `GF(p²)` (large `p≡−1`) | match |
|---|---|---|---|---|---|
| 16 | 8 | 3281 | 3281 | 3281 | ✓ |
| 16 | 9 | 3280 | 3280 | 3280 | ✓ |
| 16 | 12 | 1233 | 1233 | 1233 | ✓ |
| 16 | 13 | 464 | 464 | 464 | ✓ |
| 32 | 4 | 29601 | 28545 (sat.) | **28097** (sat.) | both short |
| 32 | 5 | 144288 | 124576 (sat.) | **129088** (sat.) | both short |

At `s=32` the chosen fields (`~2·10⁵` prime, `GF(479²)/GF(863²)`) are slightly too small to
hold the full `≈ 3·10⁴ – 1.4·10⁵` distinct sums, so **both** prime and extension under-count
by saturation — and in the second row the **genuine extension is closer to char-0 than the
prime** (`129088 > 124576`). There is **no** regime where the extension systematically gets
fewer.

### 2.3 Does it reach Kambiré's `binom(s/2,r)`?

Yes, and **exceeds it**: `N₀ˢᵘᵐ(s,r) ≥ binom(s/2,r)` (Kambiré/KK25 Lemma 9, rigorous) and
the measured value is strictly larger (e.g. `s=32,r=4`: `N₀ˢᵘᵐ = 29601` vs
`binom(16,4)=1820`; `s=32,r=10`: `1.2·10⁷` vs `8008`). Empirically `log₂N₀ˢᵘᵐ ≈ 0.75·s` at
`ρ=1/4` (i.e. `2^{Θ(s)}`, exponential in `s`, faster than the `binom` lower bound). `binom`
is therefore a safe under-estimate; `N₀ˢᵘᵐ` is `2^{Θ(s)}`.

> **General-`r` count verdict.** For **both** count conventions the distinct-bad-scalar
> count is a characteristic-zero cyclotomic invariant — it **MATCHES between prime and
> genuine odd-char extension** for general `r` (as it did for `r=4`), with deviations only
> from field-size saturation that hit both sides comparably. **Extensions do not get
> fewer.** (The haven does NOT re-open.)

---

## 3. Prize-level threshold (TASK 2)

Condition (brief): a genuine `ε_mca` violation needs `#distinct bad scalars > 2^-128·|F|`,
at `δ = capacity − 2/s`, with `ρ = (r−2)/s ∈ {1/2,1/4,1/8,1/16}` and `s ≤ 2^40`. Because the
count is field-agnostic, the same `(s,r)` works for prime **and** genuine extension —
**subject to the field being large enough that the subset sums stay distinct** (Kambiré's
Linnik window `p ∈ [4^s, 8^s]`, i.e. `log₂|F| ≳ 3s`; for a genuine extension `GF(p²)` the
same window applies to the characteristic `p`, with `|F| = p²` only buying headroom — so a
"no witness" conclusion under the `3s` floor is robust for extensions). A witness needs
**both**: (i) count `> 2^{b−128}`, and (ii) distinctness feasible (`b ≥ 3s`, i.e.
`s ≤ b/3`). Using the rigorous `N₀ˢᵘᵐ ≥ binom(s/2,r)`:

| field | `ρ` | witness? | detail |
|---|---|---|---|
| **128-bit** (`b−128 = 0`, threshold count `> 2^0 = 1`) | 1/16 | **YES** | `s=2^4=16, r=3, δ = 0.9375 − 0.125 = 0.8125`, count `≥ 2^{5.8} > 2^0`, distinctness floor `2^{48} ≤ 2^{128}` |
| | 1/8 | **YES** | `s=2^3=8, r=3, δ = 0.875 − 0.25 = 0.625`, count `≥ 2^2 > 2^0`, floor `2^{24}` |
| | 1/4 | **YES** | `s=2^4=16, r=6, δ = 0.75 − 0.125 = 0.625`, count `≥ 2^{4.8} > 2^0`, floor `2^{48}` |
| | 1/2 | **NO** | `binom(s/2,r)=0` for `ρ≥1/2` (`r>s/2`); Kambiré's count bound is **vacuous** at `1/2` (needs the footnote-18 list-size variant, a different object) |
| **256-bit** (threshold count `> 2^128`) | 1/16, 1/8, 1/4 | **NO** | INCOMPATIBLE window: count `> 2^128` needs `s ≥ 2^9 = 512`, but distinctness needs `s ≤ 2^6 = 64` (`8^s ≤ 2^256`). **The field large enough for the sums to be distinct is too small to hold `2^128` of them.** |
| | 1/2 | **NO** | vacuous `binom` (as above) |

**Field-bit boundary** (`ρ=1/4`, rigorous `binom` LB, `8^s` distinctness floor):
**a witness exists iff `b ≤ b* = 140` bits.** So **128-bit is (trivially) prize-level,
256-bit (and any `b ≥ 141`) is sub-threshold.** (Analytically: empty window iff
`2(b−128) > b/3` ⇔ `b > 154`; the discrete `s=2^t` + exact-`binom` boundary lands at 140.)

> **Convention caveat on the 128-bit witnesses (important).** All three Reading-B 128-bit
> witnesses above have **`N₀ᶠⁱʳᵉ = 0`** under the **brief's strict Reading A** (firing
> `e₂=0`, code `deg < k`): `(16,3)`,`(8,3)` are `r=3` (never fires) and `(16,6)` is
> `r ≡ 2 (mod 4)` (verified `N₀ᶠⁱʳᵉ(16,6)=0`) — the *strict* construction gives **zero** bad
> scalars there. The strict design construction `r = ρs+2` fires (under Reading A) **only
> when `r ≡ 0,1 (mod 4)`**, which at the target rates forces essentially **`r = 4`**, i.e.
> `s = 2/ρ`: the strict 128-bit witnesses are exactly
> `ρ=1/4 → (s,r)=(8,4)` (`N₀ᶠⁱʳᵉ=9`), `ρ=1/8 → (16,4)` (`49`), `ρ=1/16 → (32,4)` (`225`),
> all with distinctness floor `≤ 2^{96} ≤ 2^{128}` (verified). For a *given* rate, the strict
> construction fires only at this single small `s` (larger `s` pushes `r = ρs+2` off
> `{0,1} mod 4`); so under Reading A the strict bad-scalar supply at a target rate is the
> *single* family `(s,r)=(2/ρ, 4)`, count `(1/ρ−1)²`, which **never reaches `2^128`** (it is
> `O(1/ρ²)`) — i.e. strict Reading A is sub-threshold at 256-bit by an even wider margin.
> **`ρ=1/2` has NO strict firing witness at any `s`** (`r=s/2+2 ≡ 2 (mod 4)` ⇒ `N₀ᶠⁱʳᵉ=0`,
> verified `s=8,16`) — the strict `ρ=1/2` construction does not fire at all (the N1
> boundary). The 256-bit "no witness" is convention-robust: it fails for the *larger*
> Reading-B count, a fortiori for Reading A.

**Strict firing count** `N₀ᶠⁱʳᵉ(s,4) = (s/2−1)²` (rate `ρ=2/s`, not a target rate):
128-bit witness `s=8`; **256-bit never** — it grows only `~s²`, so at the largest
distinctness-feasible `s≈2^85` the count is `~2^{11} ≪ 2^{128}`.

> **Prize-level verdict.** A witness over a **256-bit** field does **NOT exist** in the
> prize regime (`s ≤ 2^40`, target rates), for genuine extensions exactly as for primes —
> the same finite-field obstruction Kambiré's asymptotics hide (`δ` only a *constant*
> `η = 2/s ≳ 6/log₂|F|` below capacity at deployed sizes). A witness **does** exist at
> **128-bit**, but only because the threshold there degenerates to `2^0 = 1`; the realized
> bad-scalar fraction is then a tiny `≥ 2^{−128}`, i.e. *exactly at* the boundary, not
> deep into a violation. **The extension counterexample is therefore sub-threshold at the
> headline 256-bit target and only marginal at 128-bit — it is NOT a prize-level
> near-capacity counterexample at cryptographic field size.** This matches the prime-field
> conclusion: the construction is asymptotic, and at fixed `(n,|F|)` it only forecloses a
> *constant* `η = 2/s` band below capacity.

---

## 4. Exact decoder-free cross-check (TASK 3)

No genuine-extension *firing* case with `r ≥ 5` fits `q^k ≤ 3·10⁶` (`k=(r−2)m ≥ 6` forces
`q^k` huge), so per the brief we verify `r=4` at genuine extensions, decoder-free (we never
decode — closeness is exact Hamming distance to the **full** codeword set, or an exact
in-field polynomial-agreement count; `S*` is exact branch-and-bound or exact-by-degree).
Results in `results/n2_crosscheck.json`:

- **(a) Full `rs.dist_to_code` + exact branch-and-bound `S*`** at genuine `GF(3²)`, `s=4`,
  `r=4` (the only genuine case with `q^k = 6561 ≤ 3·10⁶`): witness `{1,6,2,3}=μ₄` is
  genuinely-extension, firing `λ=0` is close (`dist=0 ≤ δn`), exact `S*=5 < rm=8` (no CA).
  **Caveat: `δ = 1 − r/s = 0` here (degenerate)** — it validates the machinery and
  genuineness but is not a near-capacity radius.
- **(b) Codeword-free EXACT structural certificate** at NON-degenerate genuine extensions
  (`δ > 0`): build `R(X)` in-field, confirm `deg R < k`, evaluate `f − λg` and the codeword
  `−R` on the **full domain** and confirm they agree on **exactly `rm = (1−δ)n` points**
  (so `dist ≤ δn`, exact, no decoder); bound the joint CA **exactly by polynomial degree**
  (`g = X^{(r−1)m}` agrees with any `deg<k` codeword on `≤ (r−1)m` points, a degree-`(r−1)m`
  polynomial having `≤ (r−1)m` roots), giving `S* ≤ (r−1)m < rm`:
  - `GF(31²)`, `s=8`, `δ=0.5`: witness genuine, `deg R = 2 < 4`, agree `8/16`, **9 distinct
    close scalars (10/10 firing subsets genuine)**, `S* ≤ 6 < 8` → **CERTIFIED BAD**.
  - `GF(127²)`, `s=8`, `δ=0.5`: same, **9 distinct close scalars** → **CERTIFIED BAD**.
  - `GF(31²)` & `GF(127²)`, `s=16`, `δ=0.75`: **49 distinct close scalars** (`=(s/2−1)²`),
    `S* ≤ 6 < 8` → **CERTIFIED BAD**.
- **(c) Exact prime control** `GF(17)`, `s=8`, `δ=0.5` (full enumeration + exact `S*`):
  `close_count = 9 = predicted`, **exact `S* = 5 < 8`**, `is_bad_line = True` — anchors that
  the field-agnostic firing line is genuinely bad.

The one residual: an *exact branch-and-bound* `S*` at a *non-degenerate* genuine extension
is infeasible (`q^k` too big), but the no-CA bound there is **exact-by-polynomial-degree**
(rigorous), and (c) gives the exact `S*` on the size-matched prime. **A real, decoder-free-
verified, near-capacity bad line on a genuine odd-char extension is confirmed (`δ = 1/2`
and `3/4`, multiple genuine extensions).**

---

## 5. CGHLL Conjecture 2 reconciliation (TASK 4) — the decisive correction

### 5.1 Conjecture 2, stated precisely (verbatim, eprint 2026/532 App. A.5)

> **Conjecture 2.** *Any Reed–Solomon code with domain of definition `D` in a prime field
> `Fp`, is line-decodable up to the Elias radius `r_E = r_E(ρ)`, with threshold
> `a = ℓ(θ)·n + o(n)`, for concurrency number `n = |D|`, and where `ℓ(θ)` is as in
> Conjecture 1. **The same is true if the alphabet of the code is from an extension field
> of `Fp`.***

with `ℓ(θ) ≤ c₁·2^{c₂·H(ρ)/η}` (Conj 1, calibrated `c₁=c₂=1` to KK25), at `θ = 1−ρ−η`, and
the **Elias radius** `r_E` defined by `ρ = 1 − H_p(r_E)`, satisfying
`1−ρ−1/log₂p ≤ r_E < 1−ρ`. Line-decodability (Def. 24, `M=1`): if **more than `a`** of the
proximate codewords on a line are `θ`-close, then `> b` of them lie on a single
codeword-line. Via the A.2 machinery this yields `ε_mca ≤ ℓ(θ)·n/|F| + o(n/|F|)`.

### 5.2 Radius placement (exact)

Our failure radius is `δ = capacity − 2/s = (1−ρ) − η` with `η = 2/s`. The construction's
field is `p ≈ 8^s` (Linnik), so `log₂p ≈ 3s` and `r_E ≈ (1−ρ) − 1/(3s)`. Hence

```
gap to capacity, ours:  η      = 2/s
gap to capacity, r_E:   1/log₂p ≈ 1/(3s)
```

Since `2/s > 1/(3s)`, **`δ` lies strictly BELOW `r_E`** (computed for `ρ=1/4,1/8,1/16`,
`s=16…1024`; `r_E − δ > 0` throughout, `≈ +0.10` at `s=16` down to `+0.0016` at `s=1024`).
This **reproduces CGHLL Remark 38 exactly**: the KK25/Kambiré failure is "clearly below the
Elias radius", gap `≥ η/2`.

### 5.3 Does our finding refute Conjecture 2 (extension clause)? — **NO**

This is the crux, and the answer corrects the prior framing. Two things must both hold for
a refutation, and **neither does**:

1. **The count must exceed the conjecture's threshold `a`.** Conjecture 2 does *not* claim
   "the proximity gap is good (small `ε`) at every `θ ≤ r_E`". It claims line-decodability
   with the **exponentially-large** threshold `a = ℓ(θ)·n = 2^{H(ρ)/η}·n` — i.e. you need
   **more than `a`** close combinations before collinearity (CA) is forced. The KK25/Kambiré
   bad-scalar count (which we confirmed is field-agnostic) is `≈ 2^{H(2ρ)/(2η)}`, and by
   `H(2ρ) ≤ 2H(ρ)` (CGHLL Rmk 38) this is `≤ 2^{H(ρ)/η}` — i.e. **at or below `a`, not
   above it.** Theorem 37 itself states the failure needs "more than `a = 2^{H(2ρ)/(2η)}`"
   combinations — the threshold is **tight, not violated**. Conjecture 2 is *calibrated*
   (`c₁=c₂=1`) precisely so KK25 is its worst case.
2. **Even the radius is consistent.** `δ < r_E` puts us inside the region Conjecture 2
   speaks about, but with the large threshold `a`; having `≈ a` close non-collinear scalars
   is exactly what "need `> a` for collinearity" permits.

Our contribution to this picture is that **the entire argument carries over verbatim to
genuine odd-char extensions** — the count `2^{H(2ρ)/(2η)}` is realized identically (§2), so
the **extension clause of Conjecture 2 is consistent-and-tight in exactly the same way as
the prime clause.** We give the extension clause its first concrete supporting evidence on
the *negative* side: the worst known example (KK25) extends, and it sits right at the
conjecture's threshold, on extensions as on primes.

### 5.4 What our finding **does** refute

- The **naive up-to-MDS-capacity** proximity-gap / correlated-agreement statement (the
  original [BCI⁺20] conjecture) — at `δ = capacity − 2/s` it fails **over genuine odd-char
  extension fields too**, not only over primes. (ESTABLISHED structurally; sub-threshold at
  crypto field size — see §3.)
- The **hope that odd-characteristic extension fields are *safer* than prime fields**
  (that the absence of the char-2 subspace-polynomial obstruction makes `F_{p^e}` immune):
  **false.** The obstruction here is **multiplicative** (cyclotomic subset sums), not
  additive, and it lives in `F_{p^e}` identically to `F_p`. CGHLL's "absence of additive
  subspaces in prime fields" hope does not buy extra safety for the *multiplicative*
  counterexample, and `F_{p^e}` does contain `F_p`-subspaces besides.

> **CGHLL Conj 2 verdict.** **Radius gap, not a refutation.** Our `δ = capacity − 2/s` is
> below `r_E`, and the bad-scalar count is at/below Conjecture 2's exponential threshold
> `a`, so the finding is **consistent with Conjecture 2's extension clause** (it confirms
> the clause is tight, identically to the prime clause). It refutes only the *naive
> MDS-capacity* gap and the *"extensions are safer"* hope.

---

## 6. The three tagged verdicts (fail-closed)

### (a) The counterexample mechanism extends structurally to odd-char extensions — **ESTABLISHED-MODULO-FORMALIZATION** (field-agnosticity Theorem §2.1; computation is verification)

The field-agnosticity is captured by the **Theorem of §2.1** (field-agnosticity of the
Kambiré subset-sum mechanism): off an explicit finite bad-prime set `B(s,r)`, the
distinct-bad-scalar count over `F_p` equals that over every genuine extension `F_{p^e}`.
The firing condition `e₂=0`, the bad scalar `λ=e₁`, and the count are **characteristic-zero
cyclotomic invariants** — `e₁,e₂` are reductions mod `p` of fixed elements of `ℤ[ξ_s]`, and a
collision is the vanishing of a fixed resultant mod `p` (independent of `e`). The exact
computation (§2.2: `char0 == prime == genuine-ext` exactly for `N₀ᶠⁱʳᵉ` in all tested `(s,r)`;
for `N₀ˢᵘᵐ` when the field is large enough, with saturation hitting both sides comparably) is
**supplementary verification** of the Theorem, and a real, decoder-free-verified, non-degenerate
(`δ=1/2, 3/4`) bad line exists on genuine `GF(31²)`, `GF(127²)` (§4). **The Wave-4 "does not
extend" was an `r=3` saturation artifact; the mechanism extends.** Tagged
**ESTABLISHED-MODULO-FORMALIZATION** (not bare PROVEN): the mechanism and its `e`-independence
are established and exactly verified; the one residual is the explicit finiteness bound for the
bad-prime set `B(s,r)` (the resultant condition) — see §2.1 "What remains to fully formalize".

### (b) It is prize-level (`ε_mca > 2^-128` at 256-bit within `s ≤ 2^40`) — **OPEN → effectively NO (ESTABLISHED sub-threshold) at 256-bit; trivially YES at 128-bit**

By exact threshold computation: **no `(s,r)`** in the prize regime makes the
distinct-bad-scalar count exceed `2^-128·|F|` over a field of `≥ 141` bits — the window is
empty (distinctness needs `s ≤ b/3`, the count needs `s ≳ 2(b−128)`). So **NOT prize-level
at 256-bit** (ESTABLISHED: no witness). At **128-bit** the threshold is `2^0=1` and
witnesses exist (`ρ=1/16,1/8,1/4`), so it is "prize-level" only in that degenerate boundary
sense (realized fraction exactly `≈ 2^{−128}`). Tagged **OPEN** *only* in the sense that a
*cleverer* extension-specific construction beating KK25's count is not ruled out — but the
**KK25/Kambiré construction itself is sub-threshold at 256-bit on extensions exactly as on
primes** (ESTABLISHED). This is field-agnostic: extensions buy *no* advantage here.

### (c) It refutes CGHLL Conjecture 2's extension clause — **ESTABLISHED FALSE (it does NOT refute it)**

`δ = capacity − 2/s < r_E`, and the bad-scalar count `≈ 2^{H(2ρ)/(2η)} ≤ a = ℓ(θ)n`, so the
finding is **consistent with** Conjecture 2 (prime *and* extension clause) — it confirms the
conjecture's threshold is tight, it does not break it. **Radius/threshold gap, not a
refutation.** (ESTABLISHED by the exact radius computation + the verbatim threshold form.)

---

## 7. Downstream claims now CONTRADICTED (LIST ONLY — do not edit; the human will fix)

The following assertions across the repo predate this verdict and are now **contradicted**
by the established findings above (mechanism extends ⇒ "extension is the *opening*/*OPEN
negative*/*recommended attack*" is wrong; and the extension counterexample does **not**
refute CGHLL Conj 2 and is sub-threshold at 256-bit). **I did not edit these.**

### `technical-note.md`
- **L10** (abstract, "Third"): "*isolating the odd-characteristic **extension-field** case
  as the genuine, under-explored opening*". — Now: the extension case is *settled* (mechanism
  extends, field-agnostic, sub-threshold at 256-bit), not an open opening.
- **L89**: "*for **extension fields** [`δ_known-negative`] is **OPEN** (every unconditional
  negative is prime-specific)*". — Now: the negative extends to genuine odd-char extensions
  (ESTABLISHED structurally), at the same constant-`η` band as primes.
- **L188 (§5.2 heading)**: "**The real opening: odd-characteristic extension fields**". —
  Contradicted (not "the real opening").
- **L190**: "*touched by **no** unconditional counterexample … **This is the sharpest gap in
  the literature and the recommended primary attack.***". — Contradicted (the multiplicative
  counterexample *does* touch it; field-agnostic).
- **L203 (Sub-lemma N2 box)**: "*extend Kambiré to an odd-characteristic extension field —
  **highest value***". — Contradicted as "highest value / open"; the extension extension is
  done and is sub-threshold at 256-bit.
- **L220 (Caveat 3)** & **L236 (table row "`δ_known-negative` for extension fields | OPEN
  (no no-go)")**: both assert the extension negative is open/untested. — Contradicted.
- **L253** ("What remains … (2) MCA, negative"): "*N1 … and **N2 (extension fields)** … N2 is
  the decisive, under-explored, recommended primary attack*". — Contradicted (N2 settled).

### `submission-checklist.md`
- **L48 (row A6)**: "*Matching impossibility for **extension** fields … the genuine opening
  (sub-lemma N2) | **OPEN***". — Contradicted.
- **L68 (row D3)**: "*Extend the no-go to odd-characteristic extension fields … **the
  recommended primary attack** | **OPEN (precise target)***". — Contradicted (the no-go
  extends; sub-threshold at 256-bit).
- **L81**: "*N2 (the odd-characteristic extension-field negative — **the recommended primary
  attack**)*". — Contradicted.
- **L83 ("What we explicitly do NOT claim")**: "*The **odd-characteristic extension-field
  negative is UNTESTED** (sub-lemma N2 — the genuine opening …)*". — Contradicted (now
  tested + established).
- **L127**: "*isolates the odd-characteristic extension-field frontier (sub-lemma N2) as the
  concrete next attack*". — Contradicted.

### `README.md`
- **L3**: "*(ii) we **map the extension-field frontier**, isolating the odd-characteristic
  `F_{p^e}` negative (sub-lemma N2) as the genuine, untouched opening … the extension-field
  negative is untested*". — Contradicted (settled, field-agnostic, sub-threshold at
  256-bit). Also the per-rate bracket rows tag `δ_known-negative` as prime-only — the
  negative now extends to genuine odd-char extensions (same band).

### `line-decoding-analysis.md`
- **L339 / L406 / L466**: the discussion framing the *smooth subgroup of an odd-char
  extension field* as "**the most genuinely-open and highest-leverage**" scenario and the N2
  target box. — Contradicted (resolved; not the highest-leverage opening). Note any line
  asserting CGHLL Conj 2's extension clause is "unsupported / the open opening" should be
  updated to: the extension clause is *consistent with* the KK25 extension (our finding),
  **and the KK25 counterexample does not refute it** (radius/threshold gap).

### `n2-extension-experiment.md`
- **L19–21 / L35 / L321 / L410** (Wave-4/Wave-5 TL;DR): "*the … multiplicative counterexample
  **DOES extend** … **the prize opening is REAL** … N2 is reopened on the negative side*". —
  The "DOES extend" half is **correct** (and now strengthened to field-agnostic + the
  general-`r` law), but the framing "**prize opening is REAL**" is **contradicted**: the
  extension counterexample is **sub-threshold at 256-bit** and only marginal at 128-bit, and
  it **does not refute CGHLL Conj 2**. The "opening" language overclaims prize-relevance.
- **L560–562**: "*… CGHLL26 Conj. 2's **extension-safety clause** … `GF(p^e)` extension
  alphabets **are not safe** from the smooth-[subgroup counterexample]*". — **Imprecise/now
  corrected**: the extension alphabets are "not safe" only in the same asymptotic,
  sub-threshold, below-`r_E` sense as prime fields; this does **not** contradict Conj 2's
  extension clause (which allows the failure below `r_E` at count `≤ a`). The clause is
  *consistent-and-tight*, not refuted.

> Net correction for the human: replace every "extension fields are the genuine
> opening / N2 is the recommended positive attack / extension negatives are OPEN/UNTESTED"
> with: **"the smooth-subgroup counterexample is field-agnostic and extends to genuine
> odd-char extensions (ESTABLISHED), but it is sub-threshold at 256-bit (no witness with
> `s ≤ 2^40`) and does NOT refute CGHLL Conjecture 2's extension clause (radius `δ < r_E`,
> count `≤` the conjecture's threshold `a`)."** N2 is **closed**, on the side of "extends but
> not prize-level and not a Conj-2 refutation".

---

## 8. Honest caveats / limits of this verdict

1. **Extension degree.** All genuine-extension exact work is `e = 2` (`GF(p²)`, the only
   proper subfield being `GF(p)`). The cyclotomic field-agnosticism argument (§2.1) is
   `e`-independent (it depends only on the characteristic `p`), so `e ≥ 3` is **LIKELY** the
   same, but is not separately enumerated here. M31² (`e=2`) is the deployed case, so this
   is the relevant one.
2. **`N₀ˢᵘᵐ` growth law `2^{Θ(s)}`** is an empirical fit (slope `≈0.75` at `ρ=1/4`,
   `s ≤ 32`) plus the rigorous lower bound `binom(s/2,r)`. The prize-threshold conclusion
   uses **only the rigorous `binom` lower bound** for the "no 256-bit witness" claim, so it
   does not depend on the fitted slope. (A larger true count only *raises* the count, which
   does **not** create a 256-bit witness because the binding constraint is *distinctness*
   `s ≤ b/3`, not the count.)
3. **Distinctness floor `3s` (`8^s`).** This is Kambiré's Linnik window for prime `p`. For a
   genuine extension `GF(p²)` the binding requirement is on the characteristic `p` (so
   `|F| = p²` is *easier*, more headroom). Using the `|F| ≥ 8^s` floor on `|F|` is therefore
   **conservative for extensions** — the "no 256-bit witness" conclusion is robust (an
   extension can only do better on distinctness, and it still fails the count/distinctness
   incompatibility). A precise extension-specific resultant bound could only *enlarge* the
   feasible `s` slightly; the `b* ≈ 140` boundary would shift modestly but `256 ≫ 140`
   leaves the 256-bit verdict unchanged.
4. **Exact `S*` at non-degenerate genuine extensions** is infeasible (`q^k` too large); the
   no-CA there is exact-by-polynomial-degree (rigorous) and the exact branch-and-bound `S*`
   is shown on the size-matched prime and on the degenerate (`δ=0`) genuine `GF(3²)`.
5. **Off-by-one convention.** Kambiré's note uses "`deg ≤ (r−2)m`" while `RS[F,L,k]` is
   "`deg < k`"; the two differ by `1/n` in rate and a large factor in count (§1). We compute
   **both** and the prize-threshold conclusion holds for both (Reading B/`N₀ˢᵘᵐ` is the more
   favorable one and still fails at 256-bit).
6. **`ρ = 1/2`.** Kambiré's `binom(s/2,r)` count is **vacuous** at `ρ ≥ 1/2` (this is the N1
   question). The `ρ=1/2` proximity-gap counterexample needs the footnote-18 *list-size*
   variant (`X^{rn/m}` with `≥ binom(m,r)` agreeing polynomials), a different object not
   analyzed in this N2 count; `ρ=1/2` is therefore **not** covered by the witnesses here and
   remains the separate N1 frontier.

---

## 9. Reproduce

```
cd experiments/small_rs_atlas
python3 n2_char0_count.py        # char-0 cyclotomic firing count (DP vs brute, validated)
python3 n2_count_laws.py         # the TWO count laws (N0_fire, N0_sum) vs binom(s/2,r)
python3 n2_prize_threshold.py    # prize threshold; writes results/n2_prize_threshold.json
python3 n2_crosscheck.py         # decoder-free cross-check; writes results/n2_crosscheck.json
python3 probe_count_scaling.py   # the orchestrator r=4 probe (baseline, GF_p2 vs prime)
# field-agnosticism (char0 == prime(p==1 mod s) == genuine ext(p==-1 mod s)) is in the
# inline runs cited in §2.2; n2_general_r_count.py holds the DP + brute validators.
```
(`python3` = 3.11 with numpy+sympy; the older `python` 3.9 runs the numpy-only modules.)
