# Sub-problem 2: the interleaved list-size challenge

Proximity Prize, ABF survey (eprint 2026/680, §1, boxed p.5). Analysis + a
self-contained resolver (`calculator/listsize_resolution.py`).

Date: 2026-06-02. Notation throughout is the ABF convention: `delta` = relative
radius, `rho = k/n`, Johnson radius `J = 1 - sqrt(rho)`, Singleton/MDS capacity
`= 1 - rho`, unique-decoding radius `UD = (1-rho)/2`, `q = |F|`, `eps* = 2^-128`.

> **One-paragraph thesis.** Sub-problem 2 is governed almost entirely by a single
> quantity — the *budget* `B := eps*·|F| = 2^-128·|F|` — and not by the hard
> coding theory that gates the MCA challenge. The list-size challenge **decouples**
> from the BCHKS Thm 1.9 barrier: that barrier ties *MCA beyond Johnson* to RS
> list-decoding **at list size `q`**, whereas sub-problem 2 only needs the list
> `≤ 2^-128·q`. The consequence is a clean trichotomy in the field size, with the
> **crossover at exactly `|F| = 2^128`**. For `|F| > 2^128` we establish a **proven
> two-sided bracket** `delta*_C^(2) ∈ [J - o(1), R_cap]`, where Johnson + interleaving
> (smooth domains) give a proven lower floor `J - eta_min`,
> `eta_min = 1/(2 rho B^{1/m})` — Johnson proves `J - eta_min` at finite budget,
> **not** `J` itself; `J - o(1)` is our explicit lower-endpoint convention — and
> `R_cap := H_q^{-1}(1-rho) ≈ (1-rho) - 1/log2 q` (the exact inverse-entropy
> crossing, the proven object — R41 naming, previously written `r_E`) is a proven
> upper ceiling (Elias/CS25; the nearby formula `1 - H_q(rho)` is its labeled
> **large-q approximation**, the convention in which the 5dp ceiling values here
> are quoted, the two radii differing beyond the fifth decimal at some deployed
> rates, §2.3 — every code, incl. this one via deep holes, has list
> `≥ q^Ω(n) ≫ B`). The value lies in `[J - o(1), R_cap]` — **not** at `R_cap`. Reaching
> `≈ R_cap` from below is **conjectural**: it needs a worst-case large-list
> (list `≤ B^{1/m}`) smooth-domain RS list-decoding-beyond-Johnson theorem, which
> is unproven (ABF §7.9 / sub-lemma P′/(D2)), weaker than and not gated by Thm 1.9,
> but still open at cryptographic field size. For `|F| = 2^128` the problem
> collapses to unique decoding (`delta* = (1-rho)/2`, proven). For `|F| < 2^128` it
> is *vacuous* — no `delta*` exists. The single open input — a proven worst-case
> small-list bound for smooth-domain RS in `(J, R_cap)` — governs both the upper
> reach to `≈ R_cap` for large fields and the strong form of the `|F| = 2^128` case.

---

## 1. Definitions: the interleaved list and the columnwise metric

### 1.1 The objects

* `C = RS[F, L, k]`: Reed–Solomon code, evaluation domain `L ⊆ F`, `|L| = n`,
  message degree `< k`, rate `rho = k/n`. Minimum relative distance
  `delta_min(C) = 1 - rho + 1/n` (MDS / Singleton, exact for RS).
* **Interleaved code (ABF Definition 2.9):**
  `C^{≡m} := { (u_1, …, u_m) : each u_j ∈ C } ⊆ (F^m)^n`. A codeword of
  `C^{≡m}` is an `m`-tuple of RS codewords; we view it as an `n`-column word whose
  `i`-th column is the vector `(u_1[i], …, u_m[i]) ∈ F^m`.
* **Columnwise / mixed metric.** For two `m`-row words `W = (w_1,…,w_m)` and
  `c = (c_1,…,c_m)`,
  ```
  d(W, c) := #{ i ∈ [n] : ∃ j ∈ [m],  w_j[i] ≠ c_j[i] }.
  ```
  A column is an "error" iff **any** of the `m` rows disagrees there. Relative
  distance is `d(W,c)/n`. (This is the metric in which `C^{≡m}` has minimum
  distance `delta_min(C)` — the same as the single code, since the rows are
  independent codewords of the *same* `C` and a single row at minimum distance
  already forces `delta_min(C)` differing columns.)

### 1.2 The list and the prize quantity

* **Single-code list (ABF §1):** for a target `f ∈ F^n`,
  `Λ(C, delta, f) := { c ∈ C : Δ(c, f) ≤ delta }`, and the *list size* is
  `Λ(C, delta) := max_{f} |Λ(C, delta, f)|`.
* **Interleaved list:** for a target `F^m`-word `T = (t_1,…,t_m)`,
  `Λ(C^{≡m}, delta, T) := { (u_1,…,u_m) ∈ C^{≡m} : d(T, (u_1,…,u_m))/n ≤ delta }`,
  and `Λ(C^{≡m}, delta) := max_T |Λ(C^{≡m}, delta, T)|`.
* **The prize success criterion (ABF §1, boxed):** for given `eps*` and constant
  `m`, determine the **largest `delta*_C ∈ [0,1]` with
  `|Λ(C^{≡m}, delta*_C)| ≤ eps*·|F|`**, together with a proof that for every
  `delta > delta*_C` the bound fails. No efficient decoder is required — only the
  value. We write this prize value `delta*_C^{(2)}(rho, |F|, m)`.

### 1.3 Why this term appears, and why `|F|` is *inside* the constraint

The toy-protocol soundness (ABF Lemma 6.6) is
`max( eps_mca(C,delta) + |Λ(C^{≡2},delta)|/|F| , (1-delta)^t )`. The list-size term
enters as `|Λ(C^{≡2},delta)| / |F|`. Demanding this be `≤ eps* = 2^-128` is exactly
`|Λ(C^{≡2},delta)| ≤ 2^-128·|F|`. **Critically, the right-hand side scales with `|F|`.**
This is the lever (Section 3).

### 1.4 Relation to single-RS list-decoding

The interleaved list is squeezed between the single-code list and its `m`-th power
(ABF Definition 2.9):
```
|Λ(C, delta)|  ≤  |Λ(C^{≡m}, delta)|  ≤  |Λ(C, delta)|^m.                  (★)
```
*Lower:* any single-code list of size `L` for a target `f` yields `L` interleaved
codewords near `T = (f, c_2, …, c_m)` (rows `2..m` exact). *Upper:* the columnwise
metric means each interleaved codeword near `T` projects, row by row, into a
single-code list near `t_j`; a tuple in the interleaved list forces *every* row to
be `delta`-close (a column is good only if all rows agree), so the tuple lies in
the product of the `m` per-row lists, of size `≤ ∏_j |Λ(C, delta, t_j)| ≤ Λ(C,delta)^m`.

There is also the **`m`-independent** GGR11 bound (ABF Lemma 2.10): for
`delta ∈ [0, delta_min(C))`, with `η := delta_min(C) - delta`, `b := ⌈delta/η⌉`,
`r := ⌈log_2(delta_min(C)/η)⌉`,
```
|Λ(C^{≡m}, delta)|  ≤  ⌈(b+r)/r⌉ · |Λ(C, delta)|^r       (any m ≥ 1).      (★★)
```
The exponent `r` is independent of `m` but blows up as `delta → delta_min(C)`. For
our regimes the simple product bound `(★)` is never worse by more than the last
significant digit, and is cleaner; we use `(★)` as the working interleaving relation
and note `(★★)` as the `m`-uniform refinement that the experiment's "m=3 matched m=2"
finding reflects (both give the same `delta*` to the resolved precision).

**Takeaway:** sub-problem 2 reduces to *bounding the single-RS list `Λ(C, delta)`*,
then raising to the `m`-th power (a benign constant-factor-in-the-exponent step,
since `m` is constant). All the action is in `Λ(C, delta)` and in the budget `B`.

---

## 2. The single-RS list as a function of `delta`: three regimes, three sourced bounds

The resolver uses exactly these three bounds; each is stated with its source and
validity window. We never extrapolate a bound past its proven range.

### 2.1 Unique-decoding regime — list `= 1` (rigorous, all fields)

For `delta < delta_min(C)/2 = (1 - rho + 1/n)/2`, every Hamming ball of radius
`delta` contains **at most one** codeword (triangle inequality + minimum distance;
MDS for RS). Hence
```
|Λ(C, delta)| = 1   and, by (★),   |Λ(C^{≡m}, delta)| = 1,   for delta < UD.
```
*Source:* Singleton/minimum-distance bound. Rigorous, all `F`, all `m`, all domains.

### 2.2 Johnson regime — list a *constant* (rigorous, all RS incl. smooth)

**ABF Corollary 3.3** (MDS specialization of Johnson 1962, ABF Theorem 3.2): for
`delta = 1 - sqrt(rho) - η` with `η > 0`,
```
|Λ(C, delta)|  ≤  1 / (2·η·rho).
```
This is a **constant independent of `n`**. By `(★)`,
`|Λ(C^{≡m}, delta)| ≤ (1/(2·η·rho))^m`. *Source:* ABF Cor 3.3. Rigorous for **all**
RS codes — in particular smooth multiplicative-subgroup domains — and all fields.
This is the rigorous positive anchor for sub-problem 2: **up to the Johnson radius
the interleaved list is a constant**, hence trivially below any budget `B ≥ 1`.

### 2.3 Capacity regime — Elias worst-case explosion (the negative ceiling)

Between Johnson and capacity the worst-case single-RS list is governed by the
list-decoding-capacity volume bound. **Crites–Stewart 2025/2046, Theorem 1 +
Lemma 1** (built on Elias 1957 / GRS list-decoding capacity, CS25 Thm 7.4.1): for a
worst-case (deep-hole) target, the single-RS list at radius `delta` satisfies
```
|Λ(C, delta)|  ≈  q^{ ( H_q(delta) - (1 - rho) ) · n }          for delta > 1 - H_q(rho),
```
and is `≤` a constant (the Johnson constant) for `delta ≤ 1 - H_q(rho)`. Here
`H_q` is the `q`-ary entropy and `1 - H_q(rho)` is the **list-decoding-capacity
radius**. The model the resolver uses:
```
log_q |Λ(C, delta)| = max( 0,  (H_q(delta) - (1 - rho)) · n ),  interleaved × m.
```
The crucial features:

* Below `1 - H_q(rho)` the exponent is `0`: the worst-case list is a constant
  (consistent with §2.2). Above it the list goes **super-polynomial in `n`**.
* `1 - H_q(rho) = (1 - rho) - (H_q(rho) - rho)`, and by CS25 Claim 1
  `0 ≤ H_q(rho) - rho ≤ 1/log_2 q`. So the explosion onset sits **`~1/log2 q`
  below Singleton capacity** — exactly what the resolver reports
  (`cap − delta*[Cap] ≈ 1/log2 q`). For `q = 2^31` the gap is `~0.032`; for
  `q = 2^256` it is `~0.004`.
  *(Notation — the `R_cap` convention, load-bearing: two nearby radii appear here —
  the **exact inverse-entropy crossing** `R_cap := H_q^{-1}(1 - rho)`, i.e. the
  `delta` solving `H_q(delta) = 1 - rho`, which is what the resolver's volume
  model computes and is the **proven ceiling object** (R41 naming; previously
  written `r_E`), and its **large-q approximation**, the formula value
  `1 - H_q(rho)`. These are **not equal at any deployed
  rate** (verified by direct computation): the crossing `R_cap` sits **above** the
  approximation
  at `rho = 1/2` (by `+9.6e-5` at 31-bit, `+1.7e-7` at 256-bit) and **below** it
  at `rho ≤ 1/4`; the solve-vs-formula gap is `≤ 0.0017` across deployed rates and
  field sizes (max `≈ 0.00166` at 31-bit, `rho = 1/8`). We keep the R42 quoting
  convention — every 5dp ceiling value quoted here is a value of the large-q
  approximation `1 - H_q(rho)` (the resolver's printed label is still `r_E`),
  not of the crossing; no claim here is sensitive to the
  `≤ 0.0017` difference, and wherever a proven *ceiling* is asserted, the
  mathematically exact proven object is the crossing `R_cap = H_q^{-1}(1 - rho)`.)*

*Source:* CS25 Thm 1/Lemma 1, Thm 7.4.1, Claim 1 (= ABF Thm 4.17 / Thm 3.6 inputs).

**This is a CEILING ONLY — a worst-case / negative-side statement, not a positive
achievability.** Read carefully what CS25 provides: Lemma 1 sandwiches the
close-count for a *uniformly random* center (an **average**, a first-moment
quantity), and Thm 1 / Thm 7.4.1(ii) are **lower** bounds (they exhibit *a* word /
*a* code with a *large* list). An average count and a lower bound **cannot
upper-bound the worst case**. So the volume bound rigorously tells us where the
worst-case list *must* blow up — giving the **upper ceiling `delta*_C^(2) ≤ R_cap`**
— but it gives **no positive (lower) reach** for `delta*` beyond Johnson. The only
proven worst-case *upper* bound on this code's list is the Johnson bound (§2.2),
valid up to `J`. It applies to RS over **any** domain (the deep holes `x^k`,
`1/(x-a)` exist over every domain), so the `R_cap` ceiling is genuine for
smooth-domain RS too. **The ceiling `R_cap` (= the crossing
`H_q^{-1}(1 - rho)`, per the §2.3 notation) is therefore NOT an achieved
`delta*`; treating it as one (reading the volume-bound crossing as the positive
answer) is exactly the step that is unproven for smooth domains** — see §4.

The matching **information-theoretic lower bound** (CS25 Thm 7.4.1(ii) / Elias):
for `rho ≥ 1 - H_q(delta) + η`, *every* `(delta, L)`-list-decodable code has
`L ≥ q^{Ω(ηn)}`. So past the list-decoding-capacity radius the worst-case list is
super-polynomial **for any code** — smooth RS cannot escape this by structure.
That makes `1 - H_q(rho)` a hard upper ceiling for any `delta*` once the budget `B`
is sub-`q^{Ω(n)}` (which it always is: `B = 2^{logF-128}` is at most `q`). It does
**not**, however, certify that `delta*` *reaches* `R_cap`.

---

## 3. The field-size lever, quantified

Everything turns on the budget `B := eps*·|F| = 2^-128·|F|`, i.e.
`log2 B = log2|F| - 128`. The resolver's `print_lever_summary()` tabulates it:

```
       field  logF     budget = 2^-128|F|                                   forces
         M31    31                  2^-97 2^-97 < 1  => delta*_C DOES NOT EXIST (degenerate)
       M31^2    62                  2^-66 2^-66 < 1  => delta*_C DOES NOT EXIST (degenerate)
  Goldilocks    64                  2^-64 2^-64 < 1  => delta*_C DOES NOT EXIST (degenerate)
       M31^4   124                   2^-4 2^-4 < 1  => delta*_C DOES NOT EXIST (degenerate)
     128-bit   128                   2^+0 = 1  => |Lambda|<=1 (unique decoding of interleaved obj)
     192-bit   192                  2^+64        = 2^64  => |Lambda|<=2^64 (loose)
     256-bit   256                 2^+128      = 2^128  => |Lambda|<=2^128 (loose)
```

This yields a **trichotomy with crossover at `|F| = 2^128`**:

### Regime A — `|F| < 2^128`: DEGENERATE (no `delta*` exists)

`B = 2^{logF-128} < 1`. Any target that *is* an interleaved codeword has
`|Λ(C^{≡m}, delta)| ≥ 1` for all `delta ≥ 0` (the codeword itself is at distance 0).
So `|Λ| ≤ B < 1` is violated for **every** `delta`. Hence **`delta*_C^{(2)}` does
not exist** — precisely ABF's "assuming `|F|` sufficiently large so that such a
`delta*_C` exists." For M31, M31², Goldilocks, and even M31⁴ (`≈2^124`), the
list-size term `|Λ|/|F|` *cannot* be driven to `2^-128` on its own; soundness at
these fields must lean on the per-query term `(1-delta)^t` or accept a larger `eps*`.
**This is not "open"; it is vacuous.** The honest statement is: sub-problem 2 as
posed has no solution for `|F| < 2^128`.

### Regime B — `|F| = 2^128`: BINDING (`delta* = (1-rho)/2`)

`B = 1`. The constraint is `|Λ(C^{≡m}, delta)| ≤ 1`: **unique decoding of the
interleaved object**. The largest `delta` with the interleaved list provably `= 1`
for *every* target is the unique-decoding radius:
```
delta*_C^{(2)} = delta_min(C)/2 = (1 - rho + 1/n)/2  →  (1-rho)/2.
```
*Why not higher?* Above `UD` a worst-case target has two single-RS codewords in its
ball (two messages whose evaluations agree on `> (1-delta_min)n` points exist once
`delta > delta_min/2`), hence two interleaved codewords (extend with identical other
rows), so `|Λ| ≥ 2 > B`. *Why the Johnson constant doesn't help here:* at `B = 1`
the Johnson bound `(1/(2ηρ))^m ≤ 1` forces `1/(2ηρ) ≤ 1`, i.e. `η ≥ 1/(2ρ) ≥ 1 > J`
for `rho ≤ 1/2` — infeasible. So the Johnson model yields nothing beyond `UD` at the
knife-edge field; the rigorous answer is `delta* = (1-rho)/2`. This is the **most
demanding deployed case** (128-bit prime / Goldilocks² / M31⁴-rounded-up).

### Regime C — `|F| > 2^128`: PROVEN BRACKET `[J - o(1), r_E]` (upper reach conjectural)

`B = 2^{logF-128} ≥ 2^64` is a large budget.

* **PROVEN lower floor (Johnson model).** `(1/(2ηρ))^m ≤ B` is satisfiable
  with `η → 0`, so `delta*_C^{(2)} ≥ J - o(1) = 1 - sqrt(rho) - o(1)` rigorously,
  with a tiny slack `η_min = 1/(2ρ·B^{1/m})` that is astronomically small
  (`< 2^-31` already at 192-bit). I.e. for any `|F| > 2^128` the interleaved list
  is provably a small constant `≤ B` *all the way up to the Johnson radius*, with
  room to spare. This floor is field-robust and holds for smooth domains.
* **PROVEN upper ceiling (capacity/Elias model).** The worst-case list provably
  exceeds `B` once `delta` passes the list-decoding-capacity radius `r_E = 1 - H_q(rho)`,
  because above it the list is `q^{Ω(n)} ≫ B` (CS25 lower bound, §2.3, for *every*
  code via deep holes). So `delta*_C^{(2)} ≤ r_E = (1-rho) - Θ(1/log2 q)`. The
  resolver reports the formula `r_E = 1 - H_q(rho)` at the field's true `q` as the
  bracket ceiling (R29/M2 correction); its `delta*[Cap]` column is the finite-`n`
  volume-bound *crossing* (`≈ H_q^{-1}(1 - rho)`), whose offset from the formula
  is **rate-dependent in sign** (§2.3 notation): **above** the formula at
  `rho = 1/2` (`+1.7e-7` at 256-bit; `+9.6e-5` at the 31-bit anchor), and
  `≈ 2·10^-5` **below** it at `rho ≤ 1/4`, `n = 2^20` (the `n` and `m` dependence
  only moves the last digit). The mathematically exact proven ceiling object is
  the crossing; the bracket column quotes the `r_E` formula convention, and the
  crossing must not be read as an *achieved* `delta*`.

**The value therefore lies in the proven bracket `[J - o(1), r_E]`** — e.g. `[0.293, 0.496]`
at `rho = 1/2`, 256-bit — **not at the single value `r_E`.**

**Reaching `≈ r_E` from below is CONJECTURAL, not proven.** The Elias/volume bound
(§2.3) is an *average + lower* bound, so it certifies the ceiling `delta* ≤ r_E` but
gives **no positive reach** beyond Johnson: it does not upper-bound the worst-case
list of *this specific smooth code* in the gap `(J, r_E)`. To prove
`delta*_C^{(2)} ≈ r_E` one needs a worst-case large-list (list `≤ B^{1/m}`)
RS list-decoding-beyond-Johnson theorem for smooth domains — the §4 conjecture
(ABF open problem §7.9 / sub-lemma P′/(D2)). That is unproven; it is weaker than and
not gated by Thm 1.9 (§4.1), but still open at cryptographic field size. So for
large fields **sub-problem 2 is proven to `J - o(1)` and bracketed below `r_E`, with the
upper reach to `r_E` conjectural** — `r_E` is within `1/log2 q` of Singleton
capacity, but it is the *ceiling*, not a certified answer.

---

## 4. The central claim to investigate (stated precisely)

The rigorous facts above leave exactly one gap — the **worst-case list of
*smooth*-domain RS in the band `(J, r_E)`** — and it is the *single open input* that
governs **both** (a) the upper reach to `≈ r_E` for `|F| > 2^128` (Regime C) and
(b) the strong form of the knife-edge `|F| = 2^128` case (Regime B). State it as a
conjecture.

> **Conjecture (smooth-domain interleaved list stays small up to capacity).**
> Let `C = RS[F, L, k]` with `L` a smooth (power-of-two multiplicative-subgroup)
> domain, `|L| = n`, `rho = k/n`. There is a function `P(rho, η)` — ideally
> `P = poly(1/η)` or even `O(1/η)`, independent of `n` — such that for every
> constant `m` and every `delta ≤ (1 - rho) - η`,
> ```
> |Λ(C^{≡m}, delta)|  ≤  P(rho, η)^m.
> ```
> In particular the interleaved list is `poly(n)`-bounded (indeed `n`-independent)
> for all `delta` up to `1 - rho - O(1/n)`.

* This conjecture is the open input the bracket needs in **two** ways. (i) For
  `|F| > 2^128`: it is exactly what would certify the *positive reach* from the
  Johnson floor `J` up to the Elias ceiling `r_E`, turning the proven bracket
  `[J - o(1), r_E]` into a single value `≈ r_E`. (ii) If true with `P = poly`, then for
  `|F| ≥ 2^128 · P(rho,η)^m` the list-size constraint holds up to
  `delta = (1-rho) - η`, lifting the binding `|F| = 2^128` case from `UD` up toward
  capacity (modulo the `poly` factor in the field requirement). Absent it, the
  proven object is the bracket and the binding-case value is `UD`.
* **What the experiment says.** Wave-1 exact enumeration (small fields `q ≤ 64`)
  found the *max interleaved list size stays exactly 1* from below Johnson, through
  Johnson, and well beyond, departing from 1 only at `≈ capacity − 0.06` and
  exploding *at* capacity; `m=3` matched `m=2`; smooth domains are comparable to
  random up to a small additive excess.
  This is strong **experimental support** for the conjecture in its strongest form
  (`P` a small constant, list essentially 1 until the rate-driven explosion at
  `⌈(1-delta)n⌉ → k`, i.e. `delta → 1 - rho`). **Caveat:** tiny fields and per-target
  sampling; it does *not* prove anything for cryptographic `|F|`, where the deep-hole
  constructions that power the Elias/CS25 worst case need a large field to exist.
* **What proving it requires (the OPEN input, not currently available for smooth
  domains).** A smooth-domain list-decoding theorem: that multiplicative-subgroup RS
  has worst-case list size `poly(1/η)` at radius `(1-rho)-η`. This is exactly the
  *derandomization* of the random-RS result (ABF Thm 3.6 / [AGL24]:
  `|Λ(C, 1-rho-η)| ≤ 2(1-rho-η)/η` w.h.p. over a **random** domain, needing
  `|F| ≥ n + k·2^{20(1-rho-η)/η}`) to the **structured** smooth domain — ABF open
  problem §7.9. **This is the same bound that `line-decoding-analysis.md` §3
  (INPUT 2 / sub-lemma (D2)) flags as "a theorem only for random domains," listed
  there as the OPEN sub-lemma P′/(D2) for smooth domains. It is NOT currently
  available for smooth domains** — it is the open input whose proof would close the
  `[J - o(1), r_E]` bracket to `≈ r_E`. No such theorem is known; CS25 + Kambiré show the
  analogous *MCA/CA* statement *fails* near capacity over prime fields, but those are
  CA failures, not list-size lower bounds (see §5), so they neither prove nor refute
  the list-size conjecture.

### 4.1 Does sub-problem 2 dodge the BCHKS Thm 1.9 barrier? — **Yes.**

BCHKS Theorem 1.9 (= ABF Thm 5.2): if RS has good proximity gaps at radius
`γ = LDR_{F,D,q}(δ) + 2/n` (the list-decoding radius **for list size `q`**) then a
counterexample forces soundness error `≥ 1/(2n)`. Its consequence (BCHKS §1, ABF
§7): **improving MCA beyond Johnson for smooth RS requires first beating the RS
list-decoding radius at list size `q = |F|`** — a hard, long-open problem.

Sub-problem 2 is **strictly easier**, for a precise reason: it asks for the list to
be `≤ eps*·|F| = 2^-128·|F|`, *not* `≤ |F|`.

* The BCHKS barrier is about `LDR` at list size `q` (the radius up to which balls
  hold `≤ q` codewords). Sub-problem 2's constraint is the *list size itself*, capped
  at `2^-128·q`, which for `q = 2^256` is `2^128 ≪ q`, and for `q = 2^128` is `1`.
* More importantly, the **list-size question never invokes proximity gaps / MCA at
  all.** It is a pure list-decoding (combinatorial) question about `Λ(C^{≡m}, δ)`. The
  Johnson bound resolves it (list = small constant) up to `J` with *zero* coding-theory
  difficulty; the Elias bound caps it at `1 - H_q(rho)`. Neither step needs an
  MCA/CA theorem, and neither is obstructed by Thm 1.9 — Thm 1.9 only bites when one
  tries to convert a *proximity-gap* statement into list-decoding, the opposite
  direction.
* Concretely: for `|F| = 2^256`, sub-problem 2 is **proven to `J - o(1)` and bracketed
  below `r_E = 1 - H_q(rho)`** by the rigorous Johnson + Elias bounds, *with no
  progress on the MCA barrier whatsoever*. The MCA challenge (sub-problem 1) at the
  same field has its positive frontier stuck at the Johnson radius. **This is the
  cleanest demonstration that the two grand challenges are not equivalent, and that
  sub-problem 2 is strictly easier (proven to `J - o(1)`, bracketed below `r_E`, not gated
  by Thm 1.9) — as Wave-1 flagged.** (Whether the bracket closes to `≈ r_E` is the
  separate, conjectural large-list smooth-domain question of §4.)

The one place the conjecture in §4 *would* brush the barrier: proving the *strong*
form (smooth-domain list `poly(1/η)` up to capacity-`η`) is a smooth-domain
list-decoding theorem, and CS25 Thm 2 (= ABF Thm 5.3) shows that *small-error CA* for
RS implies RS list-decodability — so a list-decoding theorem and a CA theorem are
linked in one direction. But sub-problem 2 does **not** require small-error CA: it
requires only the *list size bound*, which is the hypothesis CS25 Thm 2 *concludes*,
not one it needs. So even the strong-form conjecture is a list-decoding statement,
not a proximity-gap statement, and is not gated by Thm 1.9. It is gated instead by
the (also hard, but different) open problem of RS list-decoding on structured
domains.

---

## 5. Separating the list-size question from the MCA negative results

It is tempting to read Kambiré / BCHKS Thm 1.13 / CS25 Cor 1 ("CA fails near
capacity over smooth prime fields") as *also* killing sub-problem 2 near capacity.
**It does not, and the distinction is the heart of why sub-problem 2 is more
tractable.** Those are statements about `eps_ca`/`eps_mca` (the *probability that a
random line point is close without correlated agreement*), not about `|Λ|`:

* **Kambiré Thm 1 / BCHKS Thm 1.13:** construct `f,g` with many `z` such that
  `f + z·g` is `delta`-close to `C` yet `[f,g]` has large *correlated-agreement*
  distance. That is a *proximity-gap* failure (the dichotomy breaks); it bounds
  `eps_ca` from below. It says nothing directly about how many codewords sit in a
  *single* ball.
* **CS25 Cor 1 (ABF Thm 4.17):** `eps_ca = 1` in a strip below capacity. Again a CA
  statement.

The list-size quantity `|Λ(C^{≡m}, delta)|` is *upstream* of all of these. The only
list-side lower bound that bites is the **Elias volume bound** (§2.3): past
`1 - H_q(rho)` the worst-case list is `q^{Ω(n)}` for *every* code, smooth or not.
That is what caps `delta*_C^{(2)}` at `1 - H_q(rho)` for large fields — and it is a
*list* statement, cleanly matching the *list* question. So:

* The MCA negatives (Kambiré, BCHKS 1.13, CS25 Cor 1) **do not lower** `delta*_C^{(2)}`
  below `1 - H_q(rho)`. They lower the *MCA* threshold (sub-problem 1), not the
  *list-size* threshold (sub-problem 2).
* The list-size threshold is pinned, top and bottom, by *list* facts: Johnson (lower,
  rigorous, smooth) and Elias (upper, rigorous, all codes). The gap between them
  (`J` to `1 - H_q(rho)`) is where the §4 conjecture lives, and the experiment
  suggests smooth RS fills it (list stays `O(1)`), but it is **unproven for
  cryptographic fields**.

---

## 6. The verdict

**Sharp statement.** *Sub-problem 2 (interleaved list size, `eps* = 2^-128`,
`rho ∈ {1/2,1/4,1/8,1/16}`, smooth domain, constant `m`) is:*

| `|F|` regime | status | `delta*_C^{(2)}` | basis |
|---|---|---|---|
| `|F| < 2^128` (M31, Goldilocks, M31², M31⁴) | **DEGENERATE** | does not exist | `2^-128|F| < 1`; list `≥ 1` always |
| `|F| = 2^128` (128-bit prime, Goldilocks², M31⁴↑) | **PROVEN (binding)** | `(1-rho)/2` (UD radius), exact | budget `= 1` ⇒ unique decoding; MDS |
| `|F| > 2^128` (192-bit, 256-bit) | **PROVEN BRACKET** `[J - o(1), r_E]`; upper reach CONJECTURAL | proven `J - o(1) ≤ delta*_C^{(2)} ≤ r_E = (1-rho) - Θ(1/log2 q)`; `≈ r_E` conjectural | Johnson floor `J - eta_min` (proven, smooth) + Elias/CS25 ceiling (proven); closing to `r_E` needs the §4 open input |

**Crossover field size.** `|F| = 2^128` exactly. Below it the problem is vacuous; at
the knife edge `delta*_C^{(2)} = (1-rho)/2` (proven); above it the proven object is
the two-sided bracket `[J - o(1), r_E]` (proven floor `J - eta_min`, proven ceiling `r_E`), with the
positive reach to `≈ r_E` conjectural (§4). The deployed-relevant open pieces are
**both** (a) whether the large-field bracket closes to `≈ r_E`, and (b) the
knife-edge `|F| = 2^128` strong form — governed by the same open input.

**For `rho = 1/2` specifically (the brief's headline ask), `n = 2^20`, `m = 2`:**

* **256-bit field:** the **proven bracket is `delta*_C^{(2)} ∈ [0.29289, 0.49609]`**
  (lower endpoint in the `J - o(1)` convention: Johnson proves `J - eta_min` with
  `eta_min = 1/(2 rho B^{1/m}) ≤ 2^-29` here, not `J` itself — Johnson gives list
  `=` small constant `≤ 2^128 = B` up to `J - eta_min`) and
  `r_E = 1 - H_{2^256}(1/2) ≈ 0.49609` is a proven ceiling (Elias gives
  list `= q^{Ω(n)} ≫ B` above `r_E`). The value lies in this bracket, **not** at
  `0.496`: reaching
  `≈ 0.496` from below is **conjectural**, requiring a worst-case smooth-RS small-list
  bound throughout `(J, r_E)` (§4, the open input). So **sub-problem 2 at 256-bit is
  proven to `J - o(1)` and bracketed below `r_E` — strictly easier than the MCA challenge,
  but not "solved at capacity."**
* **M31 (`2^31`):** `delta*_C^{(2)}` **does not exist** — `2^-128·2^31 = 2^-97 < 1`,
  so no radius satisfies the constraint. The list-size term cannot be made `2^-128`
  at M31; this is a degeneracy of the problem statement at small fields, not an open
  research question. (If one instead asked for `|Λ| ≤ 1`, the rigorous answer would be
  the UD radius `0.25`; but that is `eps* = 2^-31`, not `2^-128`.)
* **Crossover:** `|F| = 2^128`. At exactly `2^128`, `delta*_C^{(2)} = (1-rho)/2 = 0.25`
  (rigorous, unique decoding). Just above (e.g. 192-bit) the proven ceiling `r_E` is
  `≈ 0.4948`; at 256-bit `≈ 0.4961` — with the floor `J = 0.293` in both cases.

**What is genuinely OPEN.** A single open input — *can one prove the worst-case
interleaved list stays `≤` a small constant (or `poly(n)`, i.e. `≤ B^{1/m}`) for
smooth-domain RS at radii between the Johnson radius `0.293` and the ceiling
`r_E = 1 - H_q(rho) ≈ 0.492`, for `rho = 1/2`?* — governs **both** the large-field
upper reach (would close the `[J - o(1), r_E]` bracket to `≈ r_E`) and the knife-edge
`|F| = 2^128` case (the §4 conjecture would push it from `0.25` toward `0.49`). The
experiment supports this at tiny fields; no proof exists for cryptographic fields.
**It dodges the BCHKS Thm 1.9 barrier**
(§4.1) because it is a list-decoding statement, not a proximity-gap statement —
making it strictly more tractable than the MCA grand challenge.

---

## 7. Honest provenance: proven / conjectured / experimental

**Proven (rigorous, citable):**

* `(★)` `|Λ(C,delta)| ≤ |Λ(C^{≡m},delta)| ≤ |Λ(C,delta)|^m` — ABF Def 2.9.
* List `= 1` for `delta < delta_min/2` — Singleton/MDS.
* Johnson list `≤ 1/(2ηρ)` at `delta = J - η`, all RS incl. smooth — ABF Cor 3.3.
* Elias/CS25 **lower bound + average** (NOT a worst-case upper bound): worst-case
  list `≥ q^{Ω(ηn)}` for every code past list-decoding capacity, and Lemma 1's
  random-center count `≈ q^{(H_q(delta)-(1-rho))n}` — CS25 Thm 1/Lemma 1/Thm 7.4.1,
  Claim 1. This proves the **ceiling** `delta*_C^{(2)} ≤ r_E = 1 - H_q(rho)`; it does
  **not** give any positive reach beyond Johnson (§2.3).
* The trichotomy in `|F|` and the crossover at `2^128` — direct arithmetic on
  `B = 2^-128·|F|` (this document, §3; resolver).
* `delta*_C^{(2)} = (1-rho)/2` at `|F| = 2^128` (exact, binding/MDS).
* **For `|F| > 2^128`: the proven two-sided bracket `J - o(1) ≤ delta*_C^{(2)} ≤ r_E`**
  — proven floor (Johnson + interleaving, smooth) and proven ceiling (Elias/CS25),
  combine the above (this document, §3; resolver). The value lies in `[J - o(1), r_E]`.
* Sub-problem 2 decouples from BCHKS Thm 1.9 and is **strictly easier than
  sub-problem 1** — §4.1 (the barrier is a proximity-gap→list-decoding statement;
  sub-problem 2 is a pure list statement, whose bracket endpoints are pinned without
  any MCA input for large fields).

**Conjectured (plausible, unproven for cryptographic fields):**

* **Closing the `[J - o(1), r_E]` bracket to `≈ r_E` for `|F| > 2^128`** is CONJECTURAL: it
  requires a worst-case large-list (list `≤ B^{1/m}`) smooth-domain RS
  list-decoding-beyond-Johnson theorem in the band `(J, r_E)`. The volume-bound
  crossing the resolver reports as `delta*[Cap] ≈ r_E` is this **conjectural upper
  reach**, not a proven `delta*`.
* The §4 smooth-domain conjecture: worst-case interleaved list `≤ P(rho,η)^m`
  (`P` poly or constant) up to `delta = (1-rho)-η`. This is the single open input
  governing both the large-field upper reach and the `|F| = 2^128` strong form.
  Requires a structured-domain RS list-decoding theorem (derandomize ABF Thm 3.6 /
  [AGL24]); ABF open problem §7.9. **This is the same bound `line-decoding-analysis.md`
  §3 flags as a theorem only for random domains (INPUT 2 / sub-lemma (D2)) — NOT
  currently available for smooth domains.**

**Experimentally suggested (Wave-1, tiny fields `q ≤ 64`, per-target):**

* Max interleaved list `= 1` from below Johnson through and well past it, departing
  from 1 only at `≈ capacity − 0.06`, exploding at capacity; `m=3 ≡ m=2`; smooth is
  comparable to random up to a small additive excess. Supports the *strong* form of the §4 conjecture, but is **not
  evidence for cryptographic `|F|`** (the deep-hole worst-case constructions need a
  large field to instantiate, so small-field enumeration cannot see them).

**Caveats / things to VERIFY before citing numerically:**

* CS25 Thm 1's exact deep-hole list multiplicity (the binomial/`g(·)` denominator)
  came through PDF extraction imperfectly; the *scaling* `q^{(H_q(delta)-(1-rho))n}`
  (= Lemma 1 / Elias) is solid and is what the resolver uses, but the precise leading
  constant should be eyeballed against the source if a tighter `delta*` near the
  explosion is needed (it would only move `delta*[Cap]` by `O(1/n)`).
* The resolver caps the alphabet fed to `H_q` at `2^320` — strictly above every
  field size resolved, so no deployed field is ever capped. (R29/M2 correction:
  the earlier cap of `2^200` was **not** harmless — `r_E` depends on `q` at the
  `1/log2 q` scale, so the capped resolver reported the 200-bit ceiling `0.49500`
  instead of the true 256-bit `0.49609` at `rho = 1/2`; the old "stable to >6 dp
  for any cap ≥ 2^160" claim was false for the capacity column.)
* "256-bit field" in the table is a stand-in for any `|F| ≈ 2^256` (large prime or
  extension); "128-bit" likewise covers M31⁴ (`≈2^124`, rounded), a 128-bit prime, or
  Goldilocks². The verdict depends only on `log2|F|` relative to `128`.

---

## 8. Resolver and sample output

`calculator/listsize_resolution.py` is self-contained (imports **only**
`proximity_parameters`; does not touch `bounds.py`). It implements the three models of
§2, computes `delta*_C^{(2)}` by bisection / closed form, and prints the table over
`rho ∈ {1/2,1/4,1/8,1/16} × {M31, Goldilocks, 128-bit, 256-bit} × n ∈ {2^16,2^20,2^24}`
for `m = 2` and `m = 3`. Run: `python3 calculator/listsize_resolution.py`.

### 8.1 The field-size lever

```
THE FIELD-SIZE LEVER:  eps* * |F| = 2^-128 * |F|  is the entire list-size budget.
       field  logF     budget = 2^-128|F|                                   forces
         M31    31                  2^-97 2^-97 < 1  => delta*_C DOES NOT EXIST (degenerate)
       M31^2    62                  2^-66 2^-66 < 1  => delta*_C DOES NOT EXIST (degenerate)
  Goldilocks    64                  2^-64 2^-64 < 1  => delta*_C DOES NOT EXIST (degenerate)
       M31^4   124                   2^-4 2^-4 < 1  => delta*_C DOES NOT EXIST (degenerate)
     128-bit   128                   2^+0 = 1  => |Lambda|<=1 (unique decoding of interleaved obj)
     192-bit   192                  2^+64        = 2^64  => |Lambda|<=2^64 (loose)
     256-bit   256                 2^+128      = 2^128  => |Lambda|<=2^128 (loose)
Crossover: |F| = 2^128 is the knife-edge.  Below it sub-problem 2 is
vacuous (no delta*); at it, delta* = (1-rho)/2 (exact); above it, delta*
lies in the PROVEN BRACKET [J - o(1), r_E] with r_E = 1 - H_q(rho) ~
(1-rho) - 1/log2|F| (Johnson proves J - eta_min at finite budget, not J;
reaching ~r_E from the Johnson floor is CONJECTURAL, not proven).
```

### 8.2 The prize table (m = 2; m = 3 is identical to the resolved precision)

```
  rho       field  logF      n |     UD      J    cap |  d*[UD]  d*[Joh] d*[Cap*] |    PROVEN d* (bracket / value)    regime
  1/2         M31    31 2^20  | 0.2500 0.2929 0.5000 |     -        -        -   |                 does not exist DEGENERATE
  1/2  Goldilocks    64 2^20  | 0.2500 0.2929 0.5000 |     -        -        -   |                 does not exist DEGENERATE
  1/2     128-bit   128 2^20  | 0.2500 0.2929 0.5000 | 0.25000  0.25000  0.49219 |            0.25000 (exact, UD)   BINDING
  1/2     256-bit   256 2^20  | 0.2500 0.2929 0.5000 | 0.25000  0.29289  0.49609 |             [0.29289, 0.49609]   BRACKET
  1/4         M31    31 2^20  | 0.3750 0.5000 0.7500 |     -        -        -   |                 does not exist DEGENERATE
  1/4     128-bit   128 2^20  | 0.3750 0.5000 0.7500 | 0.37500  0.37500  0.74358 |            0.37500 (exact, UD)   BINDING
  1/4     256-bit   256 2^20  | 0.3750 0.5000 0.7500 | 0.37500  0.50000  0.74681 |             [0.50000, 0.74683]   BRACKET
  1/8     128-bit   128 2^20  | 0.4375 0.6464 0.8750 | 0.43750  0.43750  0.87066 |            0.43750 (exact, UD)   BINDING
  1/8     256-bit   256 2^20  | 0.4375 0.6464 0.8750 | 0.43750  0.64645  0.87285 |             [0.64645, 0.87288]   BRACKET
 1/16     128-bit   128 2^20  | 0.4688 0.7500 0.9375 | 0.46875  0.46875  0.93478 |            0.46875 (exact, UD)   BINDING
 1/16     256-bit   256 2^20  | 0.4688 0.7500 0.9375 | 0.46875  0.75000  0.93616 |             [0.75000, 0.93618]   BRACKET
```

(Columns: `d*[UD]` = unique-decoding model (rigorous, list ≤ 1); `d*[Joh]` = Johnson
model (PROVEN floor, smooth RS); `d*[Cap*]` = Elias/CS25 volume-bound *crossing* — the
CONJECTURAL upper *reach* if read as achieved, since the volume bound is an average +
lower bound, not a proven worst-case bound. The PROVEN upper *ceiling*
(above it the worst-case list is `q^Ω(n) ≫ B`) is reported in the bracket column
as the formula `r_E = 1 - H_q(rho)` at the field's true `q` (R29/M2); the exact
proven object is the inverse-entropy crossing `H_q^{-1}(1 - rho)`, whose offset
from the formula is rate-dependent in sign — **above** at `rho = 1/2` (`+1.7e-7`
at 256-bit; `+9.6e-5` at the 31-bit anchor), `≈ 2·10^-5` **below** at `rho ≤ 1/4`
(§2.3 notation). The "PROVEN d*" column gives the honest answer: `does not exist` (DEGEN),
the exact UD value (BIND), or the proven bracket `[J - o(1), r_E]` (BRACKET). Full
`n ∈ {2^16,2^20,2^24}` rows in the script output — `n` moves only the last digit of
`d*[Cap*]`.)

### 8.3 The sharp `rho = 1/2` verdict

```
SHARP VERDICT FOCUS:  rho = 1/2, n = 2^20, m = 2
  Johnson J = 0.29289, capacity = 0.50000, UD = 0.25000
        M31 (2^ 31): regime=DEGENERATE  delta*_C^(2) = DOES NOT EXIST
  Goldilocks (2^ 64): regime=DEGENERATE  delta*_C^(2) = DOES NOT EXIST
    128-bit (2^128): regime=   BINDING  delta*_C^(2) = 0.25000  (exact, UD radius)
    192-bit (2^192): regime=   BRACKET  delta*_C^(2) in PROVEN BRACKET [0.29289, 0.49479]  (upper reach ~0.49479 CONJECTURAL)
    256-bit (2^256): regime=   BRACKET  delta*_C^(2) in PROVEN BRACKET [0.29289, 0.49609]  (upper reach ~0.49609 CONJECTURAL)
```
