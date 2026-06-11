# KKH 2026/782 verification against the per-field negative ceiling

**Verdict: PARTIALLY DISCHARGES the KK25 citation debt.** The extracted
2026/782 text states and proves a cyclotomic subset-sum distinctness lemma
(`Lemma 1`) and proves prime-field near-capacity list/proximity-gap
counterexamples. It does **not** fully discharge the repo debt as currently
written: the printed lemma is not numbered `Lemma 9`, its stated prime bound is
`p > s^(s/2)` rather than `p > phi(s)^phi(s)`, and its stated range still does
not cover the N1 `rho = 1/2` out-of-range case `r = s/2 + 2`.

Source read: `/tmp/eprint_2026_782.txt` in full (965 extracted-text lines).
The extracted text itself identifies the paper as **"Failure of proximity gaps
close to capacity"** by **Dmitry Krachun, Stepan Kazanin, and Ulrich Haboeck**;
this differs from the user-facing shorthand "A Note on Mutual Correlated
Agreement" / "Han". I treat the local extracted text as the source under
review.

Evidence tags:

- **[T1]** Direct statement/proof in `/tmp/eprint_2026_782.txt`, quoted or
  line-referenced.
- **[T2]** Arithmetic or inference combining [T1] with the repo statements in
  `assembled-vs-proven.md`, `negative-endpoint-ledger.md`, or
  `n1-rho-half.md`.

## 0. Repo-side obligations checked first

[T1, repo] `assembled-vs-proven.md` identifies the single load-bearing import:

> "KK25 Lemma 9 distinctness (the `{-1,0,+1}` cyclotomic-independence /
> resultant-distinctness bound) -- CITED, NOT RE-PROVED. This is the
> load-bearing import. Where we use it inside its stated range (`r <= phi(s)/2`,
> prime calibration `p > phi(s)^{phi(s)}`), the only gap is that we have not
> re-derived its proof. Where the rho = 1/2 close-out (sub-lemma N1) needs it
> past its stated range, the claim is explicitly CONDITIONAL..."

[T1, repo] `negative-endpoint-ledger.md` Section 2.2 states the assembled lemma
needs:

> "(i) distinctness: `p > phi(s)^{phi(s)}` (the KK25 distinctness hypothesis;
> cited, not re-proved here)"

and uses tabled `s_max = 16` at `b = 31`, `s_max = 16` at `b = 64`,
`s_max = 32` at `b = 128`, with the `rho = 1/2` row carrying the N1
conditionality flag.

[T1, repo] `negative-endpoint-ledger.md` Section 3.3 says N1 needs a named
extension:

> "A routine, numerically-airtight-but-not-yet-formalized strengthening of KK25
> Lemma 9's `{-1,0,+1}`-cyclotomic-independence / resultant-distinctness bound
> past its stated `r <= phi(m)/2` range to `r=s/2+2`..."

## 1. Distinctness lemma in 2026/782

### 1.1 Does it state one?

[T1] Yes. The distinctness/count lemma is **Lemma 1**, not "Lemma 9":

> **Lemma 1.** "Let `G` be a multiplicative subgroup of a prime field `Fp`, of
> size `s`, a power of two. If `p > ss/2`, then for any integer `1 <= r <= s/2`,
> `|{x1 + ... + xr : distinct x1, ... , xr in G}| >= 2^r * binom(s/2, r)`."

The extraction renders the exponent as `ss/2`; the proof immediately below uses
`(2r)^phi(s) <= s^(s/2)`, so I read the hypothesis as `p > s^(s/2)`.

### 1.2 Does it prove one?

[T1] Yes. The proof is the cyclotomic/resultant proof we need, for the lemma's
stated family and range:

> "The key idea is that the first `phi(s) = s/2` powers of an `s`-th primitive
> root of unity are linear independent over `Q`, and this independence is
> preserved 'modulo' large primes `p`, that is when substituting the root with a
> generator `g` of `G`."

[T1] It restricts to antipodal-free subset sums and encodes them by
sum-polynomials:

> "We consider specific sums of `r` elements from `G`: Those which stem from
> subsets `S subset G` of size `|S| = r` with `S cap (-S) = empty`. Since
> `-g^i = g^{s/2+i}`, each such sum corresponds to a sum-polynomial..."

[T1] The proof's resultant bound is:

> "`|Res(R(X), Phi_s(X))| = product_{Phi_s(z)=0} |R(z)| <= (2r)^phi(s) <=
> s^(s/2)`"

and it concludes:

> "Since `Phi_s(X)` is irreducible over the rationals, and `deg P(X) < phi(s)`,
> we therefore must have `R(X) = 0`, in other words, `P(X) = Q(X)`. This
> completes the proof."

[T1] The paper also proves an asymptotic prime-selection result, **Lemma 2**:

> **Lemma 2.** "For `beta > max{2tau + 1, 12/5}`. If `s` is taken sufficiently
> large in our construction from Section 2.1, then there exists a prime
> `p in [n^beta, 2n^beta]` so that `Res(P(X) - Q(X), Phi_s(X)) not equiv 0 mod p`,
> for every pair of distinct sum-polynomials `P(X), Q(X)` as defined in the
> proof of Lemma 1."

### 1.3 Exact hypothesis range

[T1] The printed `Lemma 1` range is:

- field: prime field `Fp`;
- group: multiplicative subgroup `G <= Fp^*`;
- group size: `s`, a power of two;
- prime lower bound: `p > s^(s/2)`;
- subset size: `1 <= r <= s/2`;
- conclusion: at least `2^r * binom(s/2, r)` distinct `r`-element sums.

[T2] This is not the same statement as the repo's imported "Lemma 9" text in
`n1-rho-half.md`, which quotes:

> "Assume that `p > phi(m)^{phi(m)}`. Then for any integer `1 <= r <= phi(m)/2`,
> `|{x1 + ... + x_r : x1,...,x_r in G distinct}| >= binom(phi(m), r)`."

The 2026/782 `Lemma 1` has a stronger count and wider `r` range than that quote,
but a stronger printed prime-size hypothesis.

### 1.4 Coverage of assembled lemma hypothesis (i)

[T2] Fail-closed verdict by exact printed statement:

| field bits `b` | repo `s_max` | repo calibration `phi(s)^phi(s)` | 2026/782 Lemma 1 calibration `s^(s/2)` | exact-statement coverage |
|---:|---:|---:|---:|---|
| 31 | 16 | `8^8 = 2^24` | `16^8 = 2^32` | **No** for a 31-bit prime/BabyBear-sized characteristic |
| 64 | 16 | `8^8 = 2^24` | `16^8 = 2^32` | **Yes** on the prime-size inequality |
| 128 | 32 | `16^16 = 2^64` | `32^16 = 2^80` | **Yes** on the prime-size inequality |

[T2] Range coverage for the deployed non-N1 rates is better:
for `r = rho*s + 2`, the printed range `r <= s/2` covers `rho = 1/4, 1/8, 1/16`
at `s = 16` and `s = 32`. It does **not** cover `rho = 1/2`, where
`r = s/2 + 2`.

[T2] There is a proof-level sharpening visible in the printed proof:
the proof only needs `p > (2r)^phi(s)` for a fixed `r`, because that is the
displayed resultant bound before it is relaxed to `s^(s/2)`. Under this
proof-level inference, a BabyBear-sized 31-bit prime can pass the deployed
non-N1 `s=16` rows:

| `s` | `rho` | `r = rho*s + 2` | `(2r)^phi(s)` | log2 bound |
|---:|---:|---:|---:|---:|
| 16 | 1/4 | 6 | `12^8` | 28.680 |
| 16 | 1/8 | 4 | `8^8` | 24.000 |
| 16 | 1/16 | 3 | `6^8` | 20.680 |

This is **not** a named theorem statement in the paper, so I tag it [T2], not a
full T1 discharge of the repo's exact `p > phi(s)^phi(s)` citation.

### 1.5 Coverage of N1 (`rho = 1/2`)

[T1] No. `Lemma 1` is stated only for `r <= s/2`; N1 needs `r = s/2 + 2`.

[T2] The paper does contain a final Theorem 1 for arbitrary target
`delta in (0,1)` via the appendix quotient construction, so it gives a
near-`rho = 1/2` asymptotic prime-field proximity-gap line. But that is not the
repo's N1 statement: it does not prove the monomial-line `r=s/2+2` distinctness
extension, and it does not give the exact per-field `rho=1/2` table row.

[T2] Therefore N1 remains **conditional** as stated in the repo. The paper
partially supports the method but does not discharge the named out-of-range
Lemma-9 extension.

## 2. Counterexample parameters and table cross-check

### 2.1 Main theorem parameters

[T1] **Theorem 1** states:

> "For `tau >= 1`, `delta in (0,1)`, and `beta > max{tau + 1, 12/5}`. There
> exist infinitely many `n = 2^b` and primes `p equiv 1 mod n`, `p = Theta(n^beta)`,
> together with Reed-Solomon codes `C = RS[Fp, H, k]` with distance
> `delta* = delta + o(1/log n)` and evaluation domain `H subset Fp^*`, a
> multiplicative subgroup of size `n`..."

with

> "`eta = c/(tau log n) * (1 + o(1))`, where `c in (0,1)` is a constant which
> depends only on `delta*`."

[T1] Its list-size and bad-scalar conclusions are:

> "`ell(delta* - eta) = |{c in C : Delta(c,u) <= delta* - eta}| >= 2^{c/eta}
> = n^{tau-o(1)}`"

and

> "`a(delta* - 2eta) = |{lambda in Fp : Delta(u0 + lambda u1, C) <= delta* - 2eta}|
> >= 2^{c/eta} = n^{tau-o(1)}`, and yet `Delta(u1,C) >= delta* - eta`."

### 2.2 Monomial-line construction parameters

[T1] Section 2.1 is the monomial-line construction matching the repo's Kambire
line:

> "`u0(x) = x^{rm}`, `u1(x) = x^{(r-1)m}`, for `x in H`, and the Reed-Solomon
> code `C` is generated by the polynomials `p(X) in Fp[X]_{<= (r-2)m}`."

[T1] The radii are:

> "`Delta(u0 + c1 u1, C) <= delta* - 2/s` but `Delta(u1, C) >= delta* - 1/s`,
> where `delta* = 1 - (r - 2)/s` is the distance of the code `C`."

[T1] The Section 2 bad-scalar count from Lemma 1 is:

> "`a = 2^r * binom(s/2, r)`"

and asymptotically:

> "`a sim 2^{s(rho + 1/2 H_2(2rho) - o(1))} = 2^{s*(c-o(1))} = 2^{(c-o(1))/eta}`."

[T1] **Proposition 1** covers the direct monomial-line statement for
`delta in [1/2,1)` and `beta > 2tau + 1`:

> "There exist functions `u0, u1 in Fp^H` with
> `|{lambda in Fp : Delta(u0 + lambda u1, C) <= delta* - 2eta}| >=
> 2^{(c-o(1))/eta} = n^{tau-o(1)}`, and yet `Delta(u1,C) >= delta* - eta`."

[T1] **Remark 3** says the polynomial-field-size condition can be improved:

> "Lemma 2 can be improved to `beta > max{tau + 1, 12/5}` by allowing that few
> resultants may be zero modulo `p`."

### 2.3 Appendix quotient-line parameters

[T1] Appendix A proves the final theorem without Lemma 1:

> "In this section we show how to prove Theorem 1 without Lemma 1."

[T1] It starts from a bad list center, then uses quotienting. The line is:

> "`u0(x) = u(x)/(x^m - z^m)`, `u1(x) = 1/(x^m - z^m)`, `x in H`"

and the code is:

> "`C = RS[Fp, H, k]` ... generated by the polynomials of degree at most
> `k = (r - 2)m`."

[T1] The appendix radii are again:

> "`Delta(u1,C) >= 1 - (r-1)/s = delta* - 1/s`"

and

> "`Delta(u0 + lambda u1, C) <= 1 - r/s = delta* - 2/s`."

[T1] **Proposition 3** gives the final all-`delta` count:

> "For `tau >= 1`, `delta in (0,1)`, and `beta > max{tau + 1, 12/5}`... There
> exist functions `u0, u1 in Fp^H` with
> `|{lambda in Fp : Delta(u0 + lambda u1, C) <= delta* - 2eta}| gtrsim
> 2^{(c-o(1))/eta} = n^{tau-o(1)}`, and yet `Delta(u1,C) >= delta* - eta`."

[T1] **Proposition 4** gives an extreme variant:

> "For any `delta in (0,1)` and `beta > 12/5`... 
> `|{lambda in Fp : Delta(u0 + lambda u1, C) <= delta* - 2eta}| >= p/(2n)`,
> and yet `Delta(u1,C) >= delta* - eta`."

### 2.4 Cross-check against `delta_unsafe = (1-rho) - 2/s_max`

[T2] The monomial and quotient constructions both have
`delta* = 1 - (r-2)/s`. With actual rate `rho = (r-2)/s`, the close radius is

`delta* - 2/s = 1 - r/s = (1-rho) - 2/s`.

So the paper's radius formula agrees with the repo table formula.

[T2] The table arithmetic itself is consistent:

| `s_max` | `rho=1/2` | `rho=1/4` | `rho=1/8` | `rho=1/16` |
|---:|---:|---:|---:|---:|
| 16 | 0.3750 | 0.6250 | 0.7500 | 0.8125 |
| 32 | 0.4375 | 0.6875 | 0.8125 | 0.8750 |
| 64 reference | 0.46875 | 0.71875 | 0.84375 | 0.90625 |

[T2] Inconsistencies or caveats to carry:

1. The exact printed distinctness lemma does **not** use the repo's
   `p > phi(s)^phi(s)` calibration. It states `p > s^(s/2)`. This matters at
   the 31-bit `s=16` row.
2. The exact printed distinctness lemma does **not** cover `rho=1/2` because
   `r=s/2+2 > s/2`.
3. The paper's theorem is asymptotic over infinitely many primes
   `p equiv 1 mod n`, `p = Theta(n^beta)`. It does not itself state the finite
   deployed `s_max = 16/16/32` rows.
4. The appendix proves a broader quotient-line counterexample without Lemma 1.
   That is useful evidence, but it is not the repo's monomial-line N1
   distinctness extension.
5. Boundary-text caveat: Theorem 1 phrases `c in (0,1)`, while Appendix A sets
   `c = H_2(rho)` in Proposition 3. At the exact boundary `rho = 1/2`, this
   entropy value is `1`, not strictly below `1`. This has no effect on the
   table identity `delta_unsafe = (1-rho)-2/s`, but it is another reason to keep
   exact-`rho=1/2` claims fail-closed.

## 3. Interleaved list-size and field-agnosticity

### 3.1 Ordinary list-size evidence

[T1] The paper has a strong ordinary RS list-size lower bound. The abstract
says:

> "The same construction gives a slightly stronger list-decoding lower bound."

[T1] Theorem 1 Item 1 gives:

> "`ell(delta* - eta) = |{c in C : Delta(c,u) <= delta* - eta}| >= 2^{c/eta}
> = n^{tau-o(1)}`."

[T1] Section 3, **Proposition 2**, gives the direct list version for the Section
2 construction:

> "`ell(delta* - eta) = |{c in C : Delta(c,u1) <= delta* - eta}| sim
> 2^{(c-o(1))/eta} = n^{tau-o(1)}`."

[T1] **Remark 5** says this list lower bound needs no additive combinatorics:

> "We do not need any additive combinatorics to see that different subsets `S'`
> yield different polynomials `p_S'(X)`. This follows from that different
> subsets have different vanishing polynomials!"

[T2] Bearing on our interleaved LIST-SIZE sub-problem: this is ordinary
single-code RS list-size evidence near capacity. The extracted text contains no
interleaved-code theorem and does not settle our interleaved list-size bracket.
It is consistent with the repo's statement that near-capacity ordinary lists can
be large, while the interleaved sub-problem remains governed by the separate
smooth-domain interleaved-list conjecture.

### 3.2 Field-agnosticity / extension fields

[T1] The paper's scope is prime fields. The abstract says:

> "for Reed-Solomon codes over multiplicative subgroups of prime fields"

[T1] The main purpose sentence says:

> "for Reed-Solomon codes over smooth multiplicative subgroups of prime fields
> `F = Fp`, both Conjecture 1 and Conjecture 2 fail near capacity."

[T1] Theorem 1 uses:

> "primes `p equiv 1 mod n`, `p = Theta(n^beta)`"

and an evaluation domain:

> "`H subset Fp^*`, a multiplicative subgroup of size `n`."

[T1] Lemma 1 is also explicitly prime-field:

> "Let `G` be a multiplicative subgroup of a prime field `Fp`..."

[T2] Therefore 2026/782 does **not** state or prove the repo's odd-characteristic
extension-field transfer. The cyclotomic-integer/norm proof is compatible with
the repo's field-agnostic mechanism, but that transfer is not a theorem in this
paper. Keep the extension-field result at the repo's existing EMF/internal
status, not as a T1 result from 2026/782.

## 4. Final debt ledger

[T2] What 2026/782 discharges:

- It gives a T1 proof of a cyclotomic/resultant distinctness lemma for prime
  fields, `s` a power of two, `1 <= r <= s/2`, and `p > s^(s/2)`.
- It gives T1 prime-field near-capacity proximity-gap and list-size
  counterexamples with radius `delta* - 2/s`, asymptotically written
  `delta* - 2eta`.
- It confirms the algebraic shape behind the repo's `delta_unsafe =
  (1-rho) - 2/s` formula.

[T2] What remains not discharged:

- The exact repo calibration `p > phi(s)^phi(s)` is not the printed 2026/782
  lemma statement.
- The 31-bit `s=16` row is not covered by the printed `p > s^(s/2)` statement,
  though the proof-level fixed-`r` bound covers deployed non-N1 rates.
- N1 at `rho = 1/2` remains outside the printed distinctness range.
- Extension-field field-agnosticity is not stated in the paper.
- Interleaved list-size is not addressed directly.
