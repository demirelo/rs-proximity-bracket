# Sumset Growth of Small Multiplicative Subgroups — What Known Theorems Give at (s, r, p) = (324, 83, ~2^256)

> Literature reconnaissance, 2026-06-10. Question: for H = μ_s ⊂ F_p^*, p ≈ 2^256 prime,
> s = |H| = 324 ≈ p^0.0326, do KNOWN theorems give |image of r-fold sums| ≥ p^{1/2} ≈ 2^128
> for r ≈ 83? Both the multiset sumset rH and the distinct-element subset sumset
> H^{(+r)} = {h_1+...+h_r : h_i ∈ H pairwise distinct} are considered.
>
> **HEADLINE VERDICT: OPEN.** The best unconditional lower bound any published theorem gives
> for either image at r = 83 is ≈ 3.9·10^4 ≈ 2^15.3 (multiset) / 2.0·10^4 ≈ 2^14.3
> (distinct), versus the target 2^128. The shortfall is a factor ≈ 2^113 (needed growth
> exponent 15.35 in |H|; known exponent ≈ 2). No surveyed result, including the
> sums-of-products (BGK) line, the Stepanov/Shkredov line, the Waring-mod-p line, and the
> subset-sums-over-subgroups line, comes within 110 bits of the target at r = 83.

## 0. Setup, conventions, and the numerical targets

- p ≈ 2^256 prime, H = μ_s the multiplicative subgroup of order s = 324 (s | p−1).
  log2 s = 8.34, so s = p^δ with **δ = 0.0326**.
- Target: image ≥ p^{1/2} = 2^128 = |H|^{15.35}. So a growth exponent ≥ 15.35 is needed.
- **−1 ∈ H** (s even, H cyclic ⇒ contains the unique order-2 element). Hence ±-signed
  sums of n elements of H are contained in the plain n-fold sumset nH. This matters: several
  theorems below are stated for signed sums / differences and transfer to nH for free.
- Notation: nA = {a_1+...+a_n : a_i ∈ A} (multiset/with repetition); A^{(+r)} = sums of r
  pairwise-distinct elements; N·A^n = N-fold sums of n-fold products. These are three
  DIFFERENT objects; claims in this area routinely blur them. For a subgroup, products are
  free (H·H = H), so the sums-of-products literature applies verbatim to plain sumsets of H —
  but the distinct-element restriction is NOT free (see §3).
- Counting ceilings (computed here): |rH| ≤ C(s+r−1, r), log2 C(406,83) = 292.3;
  |H^{(+r)}| ≤ C(s,r), log2 C(324,83) = 261.7. Both ≥ 128 bits, so the target is not
  obstructed by counting. But note C(324,83) ≈ 2^261.7 > p ≈ 2^256: **full injectivity of
  83-subset sums is impossible by pigeonhole at these exact parameters** (relevant to §5).

---

## 1. Sums-of-products / covering line (BGK, Glibichuk–Konyagin and successors)

### 1a. Bourgain–Glibichuk–Konyagin exponential sums

**Statement** (BGK, J. London Math. Soc. 73 (2006) 380–398; expositions: Kurlberg
arXiv:0705.4573, Kowalski arXiv:2401.04756 Thm 1.1, quoted verbatim): *"Let γ > 0. There
exists ν > 0, depending only on γ, such that for any prime p and any subgroup H ⊂ F_p^×
with |H| ≥ p^γ, we have Σ_{x∈H} e(ax/p) ≪ |H| p^{−ν} for any a ∈ F_p^×."*

**Quantitative dependence ν(γ):**
- Kurlberg's writeup of the BGK proof yields ν ≫ exp(−exp(C/γ)) (doubly exponentially
  small; arXiv:0705.4573, end of §4).
- The sharpest explicit version cited by Kowalski (Remark 1.2(3)) is Shkredov, *Some
  remarks on the asymmetric sum-product phenomenon*, Moscow J. Comb. Number Th. 8 (2019)
  15–41, arXiv:1705.09703, **Corollary 16**: for |Γ| ≥ p^δ,
  max_{ξ≠0} |Γ̂(ξ)| ≪ |Γ| · p^{−δ/2^{7+2/δ}} — single-exponentially small in 1/δ.

**Standard consequence for sumsets** (completion argument, standard): if
Φ_H := max_{a≠0}|Σ_{x∈H} e_p(ax)| ≤ |H| p^{−ν}, then rH = F_p once (r−2)ν ≥ 1−δ, i.e.
r ≈ 2 + (1−δ)/ν.

**At our parameters** (computed here): δ = 0.0326 gives, via Shkredov's Cor. 16,
ν = δ·2^{−(7+2/δ)} = 2^{−73.3}, hence covering at **r ≈ 2^73**. Kurlberg's ν is
astronomically smaller still.

**Range-of-validity caveat (important for the asymptotic family):** in Kambiré-type
constructions s ≍ K·log n ≍ log p, i.e. H is POLYLOG-sized along the family. BGK requires
|H| ≥ p^γ for a *fixed* γ and so never applies asymptotically there. Kowalski, Remark
1.2(2), verbatim: *"One can wonder about even smaller subgroups, but some restriction is
certainly needed... It would be interesting to see if one could say something interesting
for subgroups H of size ≍ (log p)^C for some constant C > 0."* — i.e., for polylog-size
subgroups even nontrivial exponential-sum cancellation is explicitly OPEN. Pre-BGK
thresholds were even larger: nontrivial bounds with explicit ν were known for
|H| ≫ p^{3/7+ε} (Shparlinski 1991), p^{3/8+ε} (Konyagin–Shparlinski, CUP 1999 book),
p^{1/3+ε} (Heath-Brown–Konyagin 2000), p^{1/4+ε} (Konyagin). All ≫ our 324.

**VERDICT: OPEN at r = 83.** The line does prove FULL covering (much stronger than the
p^{1/2} target) — but with r ≈ 2^63–2^73 summands at δ = 0.0326 (see also §1b), a factor
≈ 2^56 more summands than available. At r = 83 it yields nothing beyond §2. Note this line
needs no products for subgroups (H^n = H), so "sums of products" is not a loophole here.

### 1b. Explicit covering counts (Glibichuk–Konyagin, Cipra–Cochrane)

- Glibichuk–Konyagin, *Additive properties of product sets in fields of prime order*
  (arXiv:math/0702729, in "Additive Combinatorics", AMS 2007): for any A with |A| > p^δ,
  N·A^n = F_p with n ≤ δ^{−C}, N ≤ exp(δ^{−C}) — qualitative shape, exponential in 1/δ.
- **Cipra–Cochrane**, *Sum-product estimates applied to Waring's problem over finite
  fields*, INTEGERS 11 (2011) #A68, Theorem 2 (verbatim): *"If A is a multiplicative
  subgroup of F_q^* for which γ*(A,q) is defined and |A| > 1, then with k = (q−1)/|A|,
  γ*(A,q) ≤ 633·(2k)^{log4/log|A|}."* Here γ*(A,q) = min r with rA = F_q. They note: if
  |A| = p^ℓ this is γ* ≪ 4^{1/ℓ}, vs Glibichuk's similar result with 6^{1/ℓ} (Glibichuk,
  Izv. RAN 2011, [14] therein). Theorem 3 (signed sums): δ*(A,q) ≤ (40/3)·k^{log4/log|A|};
  since −1 ∈ H, plain sums inherit this.
- These hold UNCONDITIONALLY for any subgroup size (down to |A| = 2-ish), including ours —
  this is the only line that applies verbatim at s = 324, p = 2^256.

**At our parameters** (computed here): k = (p−1)/324, log4/log 324 = 0.2398, so
γ* ≤ 633·(2k)^{0.2398} ≈ 2^{68.9} and δ* ≤ 2^{63.1}. So: **every element of F_p is a sum
of ≈ 2^63 elements of μ_324 — a theorem — but at r = 83 this gives nothing.**

**VERDICT: OPEN at r = 83** (KNOWN-SUFFICES only if r were allowed to be ≈ 10^19).
For r = 83 to suffice via this route one would need |H| ≥ p^{1/log_4 83} ≈ p^{0.31} ≈ 2^80.

---

## 2. Fixed-number-of-summands sumsets of subgroups (Stepanov/Shkredov line)

These hold for ALL subgroups with |Γ| below the stated p-power — no lower bound on |Γ|,
absolute constants (Stepanov/Mit'kin polynomial method) — so they DO apply at s = 324.

| Result | Statement | Validity | Source |
|---|---|---|---|
| Heath-Brown–Konyagin 2000 | E^+(Γ) ≪ \|Γ\|^{5/2} ⇒ \|2Γ\| ≫ \|Γ\|^{3/2} | \|Γ\| ≪ p^{2/3} | Q. J. Math. 51 (2000) 221–235 |
| Shkredov–Vyugin | \|Γ±Γ\| ≫ \|Γ\|^{5/3} log^{−1/2}\|Γ\| | \|Γ\| < p^{1/2} | arXiv:1102.1172 (abstract, verbatim range) |
| Hart | \|2Γ\| ≫ \|Γ\|^{8/5−ε} (extends range, weaker exp.) | \|Γ\| ≪ p^{5/9} | arXiv:1303.2729 |
| **Shkredov tripling** | \|3Γ\| ≫ \|Γ\|^2 / log\|Γ\| | \|Γ\| < p^{1/2} | arXiv:1504.04522 Thm 1 (verbatim) |
| Glibichuk | 8A ⊇ Z_p^* | \|A\| ≥ 2p^{1/2} | Mat. Zametki 79 (2006) |
| Schoen–Shkredov | 7Γ ⊇ Z_p^* | \|Γ\| > p^{1/2} | arXiv:1008.0723 |
| Shkredov / Hart | 6Γ ⊇ Z_p^* | \|Γ\| > p^{55/112+ε} / p^{11/23+ε} | arXiv:1303.2729 |
| Shkredov | 5Γ ⊇ F_p^* (−1∈Γ) | \|Γ\| ≫ p^{1/2} log^{1/3} p | arXiv:1311.5726 |

The |Γ|^2/log|Γ| tripling bound is, as far as this survey found, **the largest published
growth exponent for any bounded number of summands of a subgroup below p^{1/2}**; nothing
published pushes a fixed-r sumset beyond exponent 2 + o(1) in this regime. (All r ≤ 8
covering results require |Γ| ≳ p^{1/2}.)

### Does growth iterate? Best known |rΓ| as a function of r

Two unconditional mechanisms exist:

1. **Cauchy–Davenport/Kneser chain** (Cipra–Cochrane Lemma 10, verbatim): *"If A is a
   multiplicative subgroup of F_q^* containing f linearly independent points over F_p and
   |A| > 2, then for any positive integer n, |nA| ≥ min{q, n|A|}."* — linear in n only.
2. **Glibichuk–Konyagin iteration** (GK Lemma 5.2 for F_p; Cipra–Cochrane Lemma 14
   generalization): define A_0 = A, A_l = 2A_{l−1} − 2A_{l−1} + A − A. Then
   |A_l| ≥ (3/8)·min{|A|^l, (p−1)/2} (GK's F_p form, as quoted by Cipra–Cochrane), and
   A_l sits inside the t_l-fold ±sumset with t_l = 4t_{l−1}+2, i.e. **t_l = (5/3)4^l − 2/3**
   (t_1=6, t_2=26, t_3=106, t_4=426, ...). Since −1 ∈ H, A_l ⊆ t_l·H. So the exponent grows
   like **log_4 r**: |rH| ≳ |H|^{⌊log_4(3r/5)⌋}.

**At our parameters** (computed here): r = 83 admits l = 2 (t_2 = 26 ≤ 83 < t_3 = 106), so
|83H| ≥ (3/8)|H|^2 = 39366 ≈ 2^{15.3}. The CD/Kneser chain gives 83·324 = 26892 ≈ 2^{14.7};
tripling gives ≈ |H|^2/log = 12587·c. **Best known: |83H| ≥ 39366 ≈ 2^{15.3}.**
To reach 2^128 = |H|^{15.35} via the GK iteration needs l = 16, i.e. t_16 ≈ 7.2·10^9
summands. Conversely, r = 83 buys exponent 2 of the needed 15.35.

**VERDICT: OPEN.** Gap ≈ 2^{112.7}. No published iteration does better than exponent
~log_4 r for tiny subgroups; whether |rΓ| ≥ |Γ|^{cr} (or even |Γ|^{r/2}) holds for
|Γ| ≪ p^ε is not addressed anywhere we found.

---

## 3. Distinct-element subset sums H^{(+r)}

This is what the Kambiré/BCHKS template actually needs (sums of r distinct s-th roots of
unity), and it is where claims most easily overreach: **none of §1–§2 transfers
automatically**, because those bounds are for multiset sums. (A naive fix-r−2-elements
transfer costs O(r·s) ≈ 2^{14.7} excluded values — same order as the §2 bounds themselves,
so it destroys them; derived here, flagged as our observation.)

What IS known:

1. **Dias da Silva–Hamidoune** (Erdős–Heilbronn conjecture; Bull. London Math. Soc. 26
   (1994) 140–146; also Alon–Nathanson–Ruzsa polynomial method): for ANY A ⊆ F_p,
   |A^{(+r)}| ≥ min{p, r(|A|−r) + 1}.
   At (324, 83): **|H^{(+83)}| ≥ 83·241+1 = 20004 ≈ 2^{14.3}.** This is the best
   unconditional distinct-sum bound we found — and it uses no multiplicative structure
   at all. No published result exploits the subgroup structure of H for distinct subset
   sums in our size regime.
2. **Zhu–Wan**, *Counting subset sums over multiplicative subgroups of finite fields*
   (arXiv:1101.0289), Theorem 1.1: |M_H(k,b) − (1/q)·C((q−1)/m, k)| ≤ (2/√q)·C((q+k+m√q)/m, k)
   where m = index of H; Corollary 1.2 (verbatim): *"There is an effectively computable
   absolute constant 0 < c < 1 such that if m < c√q and 6 ln q < k ≤ (q−1)/2m, then
   M_H(k,b) > 0 for all b ∈ F_q."* I.e., every element is a k-subset-sum — but only for
   **index m < c√q, i.e. |H| ≫ √q**. Our index is ≈ 2^{247.7}; also k = 83 < 6 ln q ≈ 1065.
   Doubly inapplicable. (Same for Wang–Wan-type results on quadratic residues, index 2,
   arXiv:1702.03028.)

**VERDICT: OPEN.** Best known |H^{(+83)}| ≥ 2^{14.3} vs target 2^128; gap ≈ 2^{113.7}.
The subgroup-aware subset-sum literature (Zhu–Wan line) starts at |H| ≳ √q.

---

## 4. Waring-mod-p / value sets of x^k + ... (Heilbronn, Konyagin–Shparlinski, Cochrane–Pinner)

The image of r-fold sums of H is exactly the value set of the diagonal form
x_1^k + ... + x_r^k with k = (p−1)/s — so this literature is the same question.

- Cauchy (1813!): γ(k,p) ≤ k. Cipra–Cochrane Thm 1: γ*(A,q) ≤ k+1 for |A| > 2.
- Heilbronn's conjecture I, proven by Cipra–Cochrane–Pinner (2007), explicit form
  **γ(k,p) ≤ 83·√k** by Cochrane–Pinner (the constant 83 here is a coincidence, unrelated
  to our r): useless for us since k ≈ 2^{247.7}.
- Heilbronn's conjecture II, proven by Konyagin (cited as [18] in Cipra–Cochrane):
  γ(k,p) ≪_ε k^ε once |A| > c(ε) — ineffective c(ε), same exponential-in-1/δ content as §1b.
- Bovey-type lattice bounds for |A| = O(1): γ(A,q) ≤ c_1(t)·k^{1/φ(t)} (Cipra–Cochrane
  Thm 4). At t = 324, φ(324) = 108: k^{1/108} ≈ 2^{2.3} — but c_1(t) depends on cyclotomic
  coefficient sizes and is not computed for t = 324; the theorem is for FIXED t as
  p → ∞ with unspecified constants. Worth flagging: this is the only result we found whose
  SHAPE (k^{1/φ(t)} ≈ 5) is small at our parameters, but it bounds the COVERING number
  γ, with an unevaluated constant c_1(324); it gives no lower bound on the image at r = 83.
  [INFERRED — if c_1(324) were explicit and ≤ ~16, this would cover F_p with
  c_1·k^{1/108} ≤ 83 summands. Verify against Bovey, Acta Arith. 1977, and
  Cipra–Cochrane–Pinner 2007 before relying on it; we did NOT verify c_1(324).]
- Konyagin–Shparlinski (*Character sums with exponential functions*, CUP 1999): the
  Gauss-sum reduction gives |Σ_{x∈H} e_p(ax)| ≤ √p always (Weil), nontrivial only for
  |H| > √p. For |H| = 324 ≪ 2^128 it is vacuous.

**VERDICT: OPEN at r = 83** (with the c_1(t)·k^{1/φ(t)} covering bound flagged as the one
unverified potential exception for COVERING — note even it says nothing for image-size at
r below its threshold, and its constant is unevaluated).

---

## 5. Sums of s-th roots of unity for small s: distinctness / equidistribution (CFKLLS–Kambiré line)

- **Pigeonhole impossibility** (computed here): C(324,83) ≈ 2^{261.7} > p ≈ 2^256, so the
  83-subset sums of μ_324 CANNOT all be distinct mod any 256-bit prime. Any distinctness-
  based argument must restrict the family (Kambiré restricts to C(s/2, r) tuples;
  log2 C(162,83) = 157.9 < 256, so distinctness of THAT family is combinatorially possible).
- **Resultant + Linnik method** (BCHKS ECCC TR-169 (2025) Thm 7.1 template; Kambiré
  arXiv:2604.09724, see literature/notes/kambire.md): distinctness of the C(s/2,r) subset
  sums is proven for SOME prime p ≡ 1 (mod n) in the window [4^s, 8^s], by counting bad
  primes dividing resultants Res(Φ_s, Q). For s = 324 the window is [2^648, 2^972] —
  **p ≈ 2^256 is far outside it**, and the union bound fails at 2^256: bad-prime count
  ≈ C(s,r)^2·log_4 s ≈ 2^{524+} vastly exceeds the ≈ 2^{240} available primes ≡ 1 mod n
  below 2^256 (computed here). The method proves existence of good LARGE primes; it says
  nothing about any fixed deployed 256-bit prime, and nothing about image ≥ 2^128 when
  full distinctness fails.
- **CFKLLS** (Canetti–Friedlander–Konyagin–Larsen–Lieman–Shparlinski, *On the statistical
  properties of Diffie–Hellman distributions*, Israel J. Math. 120 (2000) 23–46): the
  exponential-sum machinery there (and the whole pre-BGK ladder, §1a) needs |H| ≥ p^{1/4+ε}
  at best. Not applicable at s = 324.
- **Untrau** (arXiv:2112.05441, *Equidistribution of exponential sums indexed by a subgroup
  of fixed cardinality*): for FIXED order d, as q → ∞ over q ≡ 1 (mod d), the sums
  Σ_{x∈μ_d} e_q(ax) equidistribute (random-like behavior on average over the modulus).
  Confirms the heuristic that for "most" p our sums behave randomly — and contains nothing
  for an individual fixed p.
- Char-0 structure of collisions: vanishing sums of roots of unity are classified
  (Lam–Leung, J. Algebra 224 (2000) 91–109: lengths lie in N·2 + N·3 for s = 2^2·3^4);
  this controls collisions for p outside the (astronomically large) bad set, which is
  exactly what cannot be certified for a fixed p.

**VERDICT: OPEN for fixed p.** Distinctness/equidistribution is KNOWN on average over p
(and full distinctness is pigeonhole-impossible at our exact (s,r,p)); for a fixed 256-bit
prime nothing is known beyond §2–§3.

---

## 6. Summary table and the honest gap

| Tool | Object | Gives at (324, 83, 2^256) | Needs | Verdict |
|---|---|---|---|---|
| BGK + completion (§1a) | rH covering | nothing for r < ~2^73 | r ≈ 2^73 | OPEN at r=83 |
| Cipra–Cochrane Thm 2/3 (§1b) | rH covering | rH = F_p at r ≈ 2^63 | r ≈ 2^63 | OPEN at r=83 |
| GK iteration / CC Lemma 14 (§2) | \|rH\| | ≥ (3/8)\|H\|^2 = 39366 ≈ 2^15.3 | \|H\|^15.35 | OPEN (gap 2^112.7) |
| CD/Kneser chain (§2) | \|rH\| | ≥ 83·324 = 26892 ≈ 2^14.7 | — | OPEN |
| Shkredov tripling (§2) | \|3H\| ≤ \|rH\| | ≳ \|H\|^2/log ≈ 2^13.6 | — | OPEN |
| Dias da Silva–Hamidoune (§3) | \|H^{(+r)}\| | ≥ 20004 ≈ 2^14.3 | 2^128 | OPEN (gap 2^113.7) |
| Zhu–Wan (§3) | H^{(+r)} covering | inapplicable (needs \|H\| ≳ √q, k > 6 ln q) | — | N/A |
| Bovey/CC Thm 4 (§4) | covering | shape k^{1/φ(324)} ≈ 5, constant c_1(324) unevaluated | verify c_1 | UNVERIFIED |
| Resultant+Linnik (§5) | distinctness | only for some p ∈ [2^648, 2^972] | fixed p ~ 2^256 | OPEN for fixed p |

**No KNOWN-SUFFICES anywhere; no CLOSE-GAP either — every applicable unconditional bound
is ≥ 110 bits short.** The question "|83-fold (subset) sums of μ_324| ≥ 2^128 for a fixed
256-bit prime" is, per this survey, genuinely open, and moreover sits squarely in the
regime (polylog-size subgroups, fixed modulus) that Kowalski's 2024 survey explicitly
names as open even for the much weaker question of any exponential-sum cancellation.

### What WOULD suffice (conditional targets, computed here — not in the literature)

Let Φ_H = max_{a≠0} |Σ_{x∈H} e_p(ax)| for the specific deployed p. Standard completion
arguments give:
- **Image ≥ √p with r = 83** already follows from Φ_H ≤ |H|^{1−c} with modest c: with
  square-root cancellation Φ_H ≈ √s, image ≥ √p needs only r ≥ 16, and FULL covering
  rH = F_p needs only r ≥ 62 (< 83). Even Φ_H ≤ s^{0.55} gives covering at r ≈ 71.
- So the entire question reduces to: does μ_324 mod the deployed p have random-like Gauss
  sums? True for almost all p (§5, heuristically and in Untrau's averaged sense), provable
  for none of the deployed ones by current technology. A single 324-point FFT-style
  computation CHECKS Φ_H numerically for a given p (Φ_H is computable: max over a of a sum
  of 324 roots of unity — p−1 values of a is infeasible, but Φ_H relates to the largest
  Fourier coefficient, which can be bounded via E^+(H), itself computable in ~s^3 ops;
  E^+(H) ≤ s^{5/2} is guaranteed by Heath-Brown–Konyagin, and a numerically certified
  E_r(H) for moderate r gives unconditional image bounds via |rH| ≥ |H|^{2r}/E_r(H)).
  This computational route — certify higher additive energies of the concrete μ_324 —
  appears to be the only path to closing the gap without new theorems.

## Reference list (Tier 1 = read in full text here; Tier 2 = abstract/secondary)

- T1 Cipra–Cochrane, INTEGERS 11 (2011) #A68 (full text read; Thms 1–4, Lemmas 10–15).
- T1 Shkredov, arXiv:1504.04522, tripling |3Γ| ≫ |Γ|²/log|Γ|, |Γ| < √p (pp. 1–4 read).
- T1 Hart, arXiv:1303.2729 (full text read; survey of 6A/7A/8A covering ladder).
- T1 Kowalski, arXiv:2401.04756 (pp. 1–6 read; BGK statement, polylog-open remark, Shkredov Cor. 16 pointer).
- T1 Kurlberg, arXiv:0705.4573 (text extracted; Thm 1.1, ν ≫ exp(−exp(C/η))).
- T1 Shkredov, arXiv:1705.09703 = Moscow J. Comb. NT 8 (2019) 15–41 (Cor. 16 extracted).
- T1 Zhu–Wan, arXiv:1101.0289 (Thm 1.1, Cor. 1.2 extracted).
- T2 Shkredov–Vyugin, arXiv:1102.1172 (abstract verbatim: |Γ±Γ| ≫ |Γ|^{5/3}log^{−1/2}, |Γ| < p^{1/2}).
- T2 Bourgain–Glibichuk–Konyagin, JLMS 73 (2006) 380–398 (via Kowalski/Kurlberg).
- T2 Glibichuk–Konyagin, arXiv:math/0702729 (abstract + secondary).
- T2 Heath-Brown–Konyagin, Q. J. Math. 51 (2000) 221–235 (via Hart Thm 2).
- T2 Schoen–Shkredov, arXiv:1008.0723 (7Γ covering, via Hart).
- T2 Shkredov, arXiv:1311.5726 (medium-size subgroups; 5-basis at √p·log^{1/3}p).
- T2 Cochrane–Pinner, INTEGERS 8 (2008) #A46 (γ ≤ 83√k, via Cipra–Cochrane intro).
- T2 Dias da Silva–Hamidoune, Bull. LMS 26 (1994) 140–146 (standard; statement well-known).
- T2 Untrau, arXiv:2112.05441 (abstract).
- T2 CFKLLS, Israel J. Math. 120 (2000) 23–46 (abstract + secondary).
- T2 Lam–Leung, J. Algebra 224 (2000) 91–109 (background, char-0 vanishing sums).
- T1 Kambiré, arXiv:2604.09724 (via repo note literature/notes/kambire.md, read in full).
- T2 BCHKS, ECCC TR-169 (2025) Thm 7.1 (via repo notes).

Flags: everything labeled "computed here" or [INFERRED] is this survey's arithmetic on top
of quoted theorems, not a literature claim. The single unverified potential exception to
the OPEN verdict is the Bovey-type constant c_1(324) in §4 — it concerns covering, not the
r = 83 image bound, and should be checked against Bovey 1977 / CCP 2007 if it ever matters.
