# Two Grand Challenges Are Not One: A Proven Bracket for the Interleaved List-Size Challenge at Large Fields

**Ömer Demirel** — Snowfall Finance — omer@progfi.xyz

**Proximity Prize — short paper.**
**Date:** 2026-06-10. **Convention:** ABF $\delta$-radius throughout ($\delta$ = relative proximity *radius*; $\rho = k/n$; Johnson radius $J = 1-\sqrt{\rho}$; Singleton/MDS capacity $1-\rho$; list-decoding capacity radius $R_{\mathrm{cap}}(q,\rho) := H_q^{-1}(1-\rho)$, the inverse-entropy crossing — see §1 for the convention).

*Companion artifact:* the per-field negative-ceiling assessments, the M31/BabyBear circle-domain interpretation, the 256-bit rescue-path assessments, the N1/N2 close-outs, and the experimental data are recorded in the companion research ledger, `negative-endpoint-ledger.md`, in the public repository at <https://github.com/demirelo/rs-proximity-bracket>. This paper is self-contained for everything it asserts.

---

## Abstract

The Ethereum Foundation Proximity Prize (ABF, eprint 2026/680) poses two grand challenges for the deployed smooth-domain Reed–Solomon code: the MCA challenge (the threshold radius for mutual correlated agreement at error $2^{-128}$) and the interleaved list-size challenge ($|\Lambda(C^{\equiv m},\delta)|\le 2^{-128}|F|$). We make five points. (1) The two challenges are not equivalent. (2) The list-size challenge is governed by the budget $B=2^{-128}|F|$: for $|F|<2^{128}$ it is vacuous; at $|F|=2^{128}$ the threshold is the unique-decoding supremum $(1-\rho)/2$, not attained (discreteness). (3) For $|F|>2^{128}$ we establish a proven bracket $\delta^\ast\in[J-o(1),\,R_{\mathrm{cap}}]$: floor $J-\eta_{\min}$ with $\eta_{\min}=1/(2\rho B^{1/m})$ from the Johnson bound; ceiling $R_{\mathrm{cap}}(q,\rho):=H_q^{-1}(1-\rho)$ from the Elias volume bound. (4) Because the bracket uses only standard list bounds, it dodges the BCHKS Theorem 1.9 barrier, which gates the MCA challenge but not the list-size challenge. (5) The near-capacity MCA negative mechanism is field-agnostic: a characteristic-zero cyclotomic invariant, identical over prime and odd-characteristic extension fields. The positive beyond-Johnson program (keystone P$'$, prime-field multiplicative subgroups) is not field-agnostic and remains open, as does the beyond-Johnson smooth-domain list question.

---

## 1. The prize problem and the two grand challenges

The deployed object is the Reed–Solomon code $C=\mathrm{RS}[F,L,k]\subseteq F^n$ with $L\subseteq F$ a **smooth** evaluation domain (a multiplicative subgroup or coset of order a power of two, $n=|L|=2^r$; ABF Def. 2.12), constant rate $\rho=k/n\in\{\tfrac12,\tfrac14,\tfrac18,\tfrac1{16}\}$, in the regime $k\le 2^{40}$, $|F|<2^{256}$, with target soundness $\varepsilon^\ast=2^{-128}$ (ABF, eprint 2026/680, §1). The two challenges are stated verbatim by ABF (boxed, pp. 4–5):

> **Grand MCA challenge.** Determine the largest $\delta^\ast_C\in[0,1]$ such that $\varepsilon_{\mathrm{mca}}(C,\delta^\ast_C)\le\varepsilon^\ast$, *with a proof that for all $\delta>\delta^\ast_C$, $\varepsilon_{\mathrm{mca}}(C,\delta)>\varepsilon^\ast$.*

> **Grand list-decoding challenge.** For a constant interleaving $m$, determine the largest $\delta^\ast_C\in[0,1]$ such that $|\Lambda(C^{\equiv m},\delta^\ast_C)|\le\varepsilon^\ast\cdot|F|$, *with the analogous proof that for all $\delta>\delta^\ast_C$ it fails.* No efficient list-decoder is required — only the value.

(We refer to the MCA challenge as **sub-problem 1** and the interleaved list-size challenge as **sub-problem 2** throughout.) Both feed the round-by-round knowledge-soundness error of the deployed FRI/STIR/WHIR protocols, which has the exact composite form (ABF Lemma 6.6, for $\delta\in(0,\delta_{\min}(C))$):
$$
\mathsf{soundness}(C,\delta,t)\;=\;\max\!\Big(\underbrace{\varepsilon_{\mathrm{mca}}(C,\delta)+\tfrac{|\Lambda(C^{\equiv 2},\delta)|}{|F|}}_{\text{combine term}},\;\underbrace{(1-\delta)^t}_{\text{query term}}\Big). \tag{1}
$$
The MCA challenge controls the first summand of the combine term; the list-size challenge controls the second. Three reference radii depend only on $\rho$: unique decoding $\mathrm{UD}=(1-\rho)/2$, Johnson $J=1-\sqrt\rho$, Singleton capacity $1-\rho$, with the strict ordering $\mathrm{UD}<J<1-\rho$ for all $\rho\in(0,1)$. A fourth radius is field-dependent: the **list-decoding capacity radius**
$$
R_{\mathrm{cap}}(q,\rho)\;:=\;H_q^{-1}(1-\rho),
$$
the exact inverse-entropy crossing — the radius $\delta$ solving $H_q(\delta)=1-\rho$. Every proven upper-ceiling statement in this paper is a statement about $R_{\mathrm{cap}}$.

**The large-$q$ approximation (labeled shorthand only).** The literature often writes the nearby formula value $1-H_q(\rho)$ for this radius; we call it the **large-$q$ approximation** and never use it as a proven ceiling. The two are *not equal at any deployed rate* (verified by direct computation): they differ by at most $0.0017$ across the deployed rates and field sizes (maximum $\approx0.00166$ at 31-bit, $\rho=\tfrac18$), with rate-dependent sign — at $\rho=\tfrac12$ the approximation sits *below* the true crossing $R_{\mathrm{cap}}$ (by $+9.6\cdot10^{-5}$ at 31-bit, $+1.7\cdot10^{-7}$ at 256-bit), while at $\rho\le\tfrac14$ it sits *above* it, so quoting the approximation as a ceiling could overstate the proven object by up to $0.0017$. The approximation sits a distance $H_q(\rho)-\rho\le 1/\log_2 q$ below Singleton capacity (CS25 Claim 1), and $R_{\mathrm{cap}}$ sits within $0.0017$ of it; no claim in this paper is sensitive to the difference, but every 5dp ceiling value quoted here is a value of $R_{\mathrm{cap}}$ itself.

The headline stakes at $\rho=\tfrac12$ (ABF §6.3.1): at the Johnson radius the per-query catch probability is $1-\delta=\sqrt{1/2}$, so the query term in (1) is $(1/\sqrt2)^t$ and $t=128$ queries yield only $(1/\sqrt2)^{128}=2^{-64}$ — **64 bits**. Operating at capacity $\delta\to\tfrac12$ would give $(1/2)^{128}=2^{-128}$ — the full **128 bits**. Closing $J\to$ capacity at $\rho=\tfrac12$ is the **$64\to128$-bit jump** the prize is about.

---

## 2. Definitions and conventions

**The notation clash (carry the ABF convention).** ABF/CS25/GG25 use $\delta$ for the proximity *radius*; **BCHKS and CGHLL use $\delta$ for the code's minimum distance** ($=1-\rho$, MDS) and $\gamma$ (resp. $\theta$) for the radius. When importing BCHKS into ABF: radius $\delta\equiv\gamma$, rate $\rho\equiv 1-\delta_{\text{BCHKS}}$, error $a/q$ with $a$ a count. **CS25's "capacity" is the *list-decoding* capacity**, strictly below the Singleton $1-\rho$ that ABF/BCHKS/GG call capacity. All statements below are in the ABF convention.

**Correlated agreement (CA), ABF Def. 4.1.** $u_0,\dots,u_\ell$ have CA with $C$ of density $\ge 1-\delta$ if there is a common $D'\subseteq L$, $|D'|/n\ge1-\delta$, and codewords $v_i\in C$ with $u_i=v_i$ on all of $D'$. The error $\varepsilon_{\mathrm{ca}}(C,\delta)$ is the probability over $\gamma$ that a random line point $f_1+\gamma f_2$ is $\delta$-close to $C$ while no such common $D'$ exists.

**Mutual correlated agreement (MCA), ABF Def. 4.3 (loss-free).**
$$
\varepsilon_{\mathrm{mca}}(C,\delta):=\max_{f_1,f_2\in F^n}\Pr_{\gamma\leftarrow F}\!\Big[\exists\,S\subseteq[n],\,|S|\ge(1-\delta)n,\;\Delta_S(f_1+\gamma f_2,C)=0\;\wedge\;\Delta_S\big((f_1,f_2),C^{\equiv 2}\big)>0\Big].
$$
MCA is loss-free by definition (ABF Rem. 4.4) and is the strongest of the three: $\varepsilon_{\mathrm{pg}}\le\varepsilon_{\mathrm{ca}}\le\varepsilon_{\mathrm{mca}}$ (ABF Fact 4.5).

**Interleaved code $C^{\equiv m}$, ABF Def. 2.9.** $C^{\equiv m}:=\{(u_1,\dots,u_m):u_j\in C\}\subseteq(F^m)^n$, with the *columnwise* metric (a column is an error iff any row disagrees). Its minimum distance equals $\delta_{\min}(C)$.

**Lists.** $\Lambda(C,\delta,f):=\{c\in C:\Delta(c,f)\le\delta\}$, $\Lambda(C,\delta):=\max_f|\Lambda(C,\delta,f)|$; analogously $\Lambda(C^{\equiv m},\delta)$ in the columnwise metric. The interleaving relation (ABF Def. 2.9 / Lemma 2.10 inputs):
$$
|\Lambda(C,\delta)|\;\le\;|\Lambda(C^{\equiv m},\delta)|\;\le\;|\Lambda(C,\delta)|^m. \tag{$\star$}
$$

---

## 3. The budget trichotomy and the proven bracket $[J-o(1),\,R_{\mathrm{cap}}]$

The list-size term in (1) is $|\Lambda(C^{\equiv 2},\delta)|/|F|$; demanding it be $\le\varepsilon^\ast=2^{-128}$ is exactly $|\Lambda(C^{\equiv m},\delta)|\le 2^{-128}\cdot|F|$. **The right-hand side scales with $|F|$.** Define the **budget**
$$
B\;:=\;\varepsilon^\ast\cdot|F|\;=\;2^{-128}\cdot|F|,\qquad \log_2 B=\log_2|F|-128.
$$
The entire challenge is decided by three citable ingredients.

**Step (i) — Johnson floor [PROVEN, all RS incl. smooth].** ABF Cor. 3.3 (MDS Johnson, Joh62 / ABF Thm. 3.2): for $\delta=1-\sqrt\rho-\eta$ with $\eta>0$,
$$
|\Lambda(C,\delta)|\;\le\;\frac{1}{2\eta\rho},
$$
a *constant independent of $n$*, valid for **all** RS codes including smooth multiplicative-subgroup domains. By $(\star)$, $|\Lambda(C^{\equiv m},\delta)|\le(1/(2\eta\rho))^m$. This is the rigorous positive floor: at fixed slack below Johnson, the interleaved list is bounded by a constant independent of $n$. It is below the target budget whenever $(1/(2\eta\rho))^m\le B$, equivalently $\eta\ge\eta_{\min}$.

**Step (ii) — volume/list-capacity ceiling [PROVEN lower bound, all codes].** CS25 Thm. 1 + Lemma 1 (built on Elias 1957 / GRS list-decoding capacity, CS25 Thm. 7.4.1): the average/volume argument gives
$$
|\Lambda(C,\delta)|\;\ge\;\frac{\operatorname{Vol}_q(n,\lfloor\delta n\rfloor)}{q^{n-k}},
\qquad\text{so}\qquad
\log_q|\Lambda(C,\delta)|\;\ge\;(H_q(\delta)-(1-\rho))\,n-o(n)
$$
whenever $H_q(\delta)>1-\rho$. Above the inverse-entropy crossing $R_{\mathrm{cap}}=H_q^{-1}(1-\rho)$ — the radius where this exponent turns positive by fixed slack — some received word has list size $q^{\Omega(n)}$; this is a **lower bound on the worst-case list, not an exact list formula**. Below the inverse-entropy crossing, the volume/average lower bound no longer forces exponential worst-case lists. It does not provide a positive worst-case upper bound; the only proven uniform upper bound used here remains Johnson (Step (i)). The list-size threshold therefore cannot exceed $R_{\mathrm{cap}}$ in any asymptotic family where this lower bound exceeds $B$. The matching information-theoretic lower bound (CS25 Thm. 7.4.1(ii)): for $\rho\ge1-H_q(\delta)+\eta$, *every* $(\delta,L)$-list-decodable code has $L\ge q^{\Omega(\eta n)}$ — so past the crossing the worst-case list is super-polynomial **for any code**, smooth RS included. This is a *list* statement matching a *list* question.

**Step (iii) — interleaving [PROVEN].** $(\star)$: $|\Lambda(C^{\equiv m},\delta)|\le|\Lambda(C,\delta)|^m$; raising to the $m$-th power is benign since $m$ is constant. (The $m$-independent GGR11 refinement, ABF Lemma 2.10, gives the same $\delta^\ast$ to the displayed precision.)

### 3.1 The trichotomy and the main Proposition

Composing (i)–(iii) yields a clean trichotomy in $\log_2 B=\log_2|F|-128$, with crossover at $|F|=2^{128}$.

> **Proposition (the interleaved list-size challenge: a proven bracket).** *Let $C=\mathrm{RS}[F,L,k]$ with $L$ smooth, $\rho=k/n$, constant interleaving $m$, $\varepsilon^\ast=2^{-128}$, $B=2^{-128}|F|$. Then:*
>
> *(A) **Degenerate**, $|F|<2^{128}$ ($B<1$): no $\delta^\ast_C{}^{(2)}$ exists. (Any target that is itself an interleaved codeword has $|\Lambda(C^{\equiv m},\delta)|\ge1$ for all $\delta\ge0$, so $|\Lambda|\le B<1$ fails everywhere.) This is **vacuous, not open**.*
>
> *(B) **Binding**, $|F|=2^{128}$ ($B=1$): $\delta^\ast_C{}^{(2)}=\delta_{\min}(C)/2\to(1-\rho)/2$ (the unique-decoding radius), rigorous via Singleton/MDS — **read as a supremum**: $|\Lambda(C^{\equiv m},\delta)|\le1$ for all relative radii $\delta<d_{\min}/(2n)$, while any $\delta\ge\lceil d_{\min}/2\rceil/n$ already admits a target within distance $\lceil d_{\min}/2\rceil$ of two interleaved codewords (take $c_1\ne c_2\in C$ at the MDS distance $d_{\min}=n-k+1$ and a midpoint word agreeing with each on half the differing columns), so list $\ge2>B$. With closed Hamming balls the half-distance endpoint itself is therefore **not attained** (discreteness); the largest attained grid radius is $(\lceil d_{\min}/2\rceil-1)/n$, and $\delta^\ast_C{}^{(2)}=\sup\{\delta:|\Lambda(C^{\equiv m},\delta)|\le1\}=d_{\min}/(2n)\to(1-\rho)/2$. (At $B=1$ the Johnson model gives nothing beyond UD: $(1/(2\eta\rho))^m\le1$ forces $\eta\ge1/(2\rho)\ge1>J$ for $\rho\le\tfrac12$, infeasible.)*
>
> *(C) **Proven bracket $[J-o(1),\,R_{\mathrm{cap}}]$**, $|F|>2^{128}$ (for arbitrary $B>1$ the proven positive side is $J-\eta_{\min}$; in the large-budget deployed regime $B\ge2^{64}$ — the deployed $\ge192$-bit sizes — this is $J-o(1)$ at all displayed precision, and we use $B\ge2^{64}$ below only as that deployed-regime convenience, not as a theorem hypothesis; **endpoint convention:** the proven floor is exactly $J-\eta_{\min}$ with $\eta_{\min}=1/(2\rho B^{1/m})$, and we write $J-o(1)$ for it — Johnson proves $J-\eta_{\min}$ at finite budget, not $J$ itself): it is **PROVEN** that $\delta^\ast_C{}^{(2)}\ge J-\eta_{\min}$ (Johnson floor (i) + interleaving: the list is $\le(1/(2\eta\rho))^m\le B$ all the way up to $J-\eta_{\min}$ — numerically $\eta_{\min}=2^{-32}$ at $\rho=\tfrac12$ and $\le2^{-29}$ for all deployed $\rho$, with equality at $\rho=\tfrac1{16}$, already at 192-bit; field-robust, smooth domains), and **PROVEN** that $\delta^\ast_C{}^{(2)}\le R_{\mathrm{cap}}(q,\rho)=H_q^{-1}(1-\rho)$ (Elias/CS25 ceiling (ii): above the crossing the average-list lower bound gives a received word with list $q^{\Omega(n)}\gg B$, for every code including this one). The value therefore lies in the **proven bracket $[J-o(1),\,R_{\mathrm{cap}}]$** — e.g. $[0.293,\,0.496]$ at $\rho=\tfrac12$, $256$-bit ($\eta_{\min}$ invisible at this precision) — **not** the single value $R_{\mathrm{cap}}$. Reaching $\approx R_{\mathrm{cap}}$ from below is **CONJECTURAL**: it requires a worst-case large-list (list $\le B^{1/m}$) RS list-decoding-beyond-Johnson theorem for smooth domains — unproven, ABF §7.9 / the keystone sub-lemma P$'$ of §6 — which is weaker than and **not** gated by BCHKS Thm. 1.9, but still open at cryptographic field size. (The worst-case small list throughout $(J,R_{\mathrm{cap}})$ for this specific smooth code is the conjectural input — see §4.)*

**Hypotheses and provenance, made explicit.** Steps (i) and (ii) are *standard, citable list bounds* (ABF Cor. 3.3; CS25 Thm. 1 / Thm. 7.4.1) used strictly inside their proven validity windows; we never extrapolate a bound past its range. Crucially, the Elias/CS25 volume bound (ii) is an **average + lower** bound (a random-center first-moment count / an existence of a large-list word), so it rigorously supplies the **upper ceiling** $\delta^\ast_C{}^{(2)}\le R_{\mathrm{cap}}$ but gives **no positive (lower) reach** beyond Johnson; the only proven worst-case *upper* bound on this code's list is Johnson, valid to $J$. Step (iii) and the trichotomy arithmetic on $B$ are **our composition** (verified in the accompanying repository's exact-arithmetic calculator). The Proposition asserts **no new list-decoding theorem** — that is the point of §4.

**The $\rho=\tfrac12$ verdict ($n=2^{20}$, $m=2$); ceiling values are $R_{\mathrm{cap}}$ at the field's true $q$, 5dp:**

| field — regime | proven $\delta^\ast_C{}^{(2)}$ | basis |
|:-------------------------|:-------------------------------------------|:----------------|
| M31 ($2^{31}$) — degenerate, Prop. (A) | does not exist | $B=2^{-97}<1$ |
| Goldilocks ($2^{64}$) — degenerate, Prop. (A) | does not exist | $B=2^{-64}<1$ |
| 128-bit ($2^{128}$) — binding, Prop. (B) | $0.25000$ (UD supremum; endpoint not attained — discreteness) | $B=1$, MDS |
| 192-bit ($2^{192}$) — bracket, Prop. (C) | $[0.29289,\,0.49479]$; conjectural upper reach $\approx0.49479$ ($R_{\mathrm{cap}}$) | Johnson + Elias |
| 256-bit ($2^{256}$) — bracket, Prop. (C) | $[0.29289,\,0.49609]$; conjectural upper reach $\approx0.49609$ ($R_{\mathrm{cap}}$) | Johnson + Elias |

(At $\rho=\tfrac12$ the crossing $R_{\mathrm{cap}}$ and the large-$q$ approximation agree at the displayed precision — they differ by $+4.1\cdot10^{-7}$ at 192-bit and $+1.7\cdot10^{-7}$ at 256-bit.) The same trichotomy holds for $\rho\in\{\tfrac14,\tfrac18,\tfrac1{16}\}$: the 256-bit brackets are $[J-o(1),\,R_{\mathrm{cap}}]$ with $R_{\mathrm{cap}}$ at $0.74681,\,0.87285,\,0.93616$ (exact inverse-entropy crossings at $q=2^{256}$; the large-$q$ approximation would give $0.74683/0.87288/0.93618$, differing at the fifth decimal). For $\rho=\tfrac12$, 256-bit, the proven bracket is $[0.293,\,0.496]$; the upper end $\approx0.496$ is a **conjectural** reach, not a proven $\delta^\ast$, since the worst-case smooth-RS list throughout $(J,R_{\mathrm{cap}})$ is unproven: **sub-problem 2 is proven to $J-o(1)$ and bracketed below $R_{\mathrm{cap}}$ — decoupled from Thm. 1.9, with a standard list-bound bracket (§4), but not "solved at capacity."**

---

## 4. Why the two challenges are not equivalent (the Thm 1.9 separation)

**[PROVEN] (BCHKS Thm. 1.9 = ABF Thm. 5.2; all RS).** For $\gamma=\mathrm{LDR}_{F,D,q}(\delta_{\min})+2/n$ (just past the list-decoding radius *for list size $q$*), there exist $f,g$ with $\ge q/(2n)$ close line points but $\Delta([f,g],C^2)\ge\delta_{\min}-1/n$, so $a/q\ge1/(2n)$ independent of $q$. **Consequence:** improving MCA / proximity gaps beyond the Johnson radius for *any* RS code **requires** first improving its list-decoding radius (at list size $q$) beyond Johnson — a separately hard, in-general-false-for-some-RS, long-open problem. The barrier is *conditional, not absolute*: it does not forbid $\delta^\ast_C>J$; it proves that any *small-error* ($\varepsilon_{\mathrm{mca}}\le\mathrm{poly}(n)/q$) certificate of $\delta^\ast_C>J$ for smooth RS is itself a list-decoding-beyond-Johnson result. The only crack is **large-error (M)CA** (error $a/q$ not $\ll1/n$ but still $\le2^{-128}$), which both CS25 and CGHLL flag as the realistic target. This is exactly why the prize poses *both* challenges.

**The list-size bracket dodges this barrier — and the reason is exactly that it uses no MCA / proximity-gap input.**

1. **Different quantity.** Thm 1.9 concerns $\mathrm{LDR}_{F,D,q}(\delta)$ — the radius up to which balls hold $\le q$ codewords (list size $q=|F|$). Sub-problem 2 caps the *list size itself* at $2^{-128}q$, which for $q=2^{256}$ is $2^{128}\ll q$ and for $q=2^{128}$ is $1$. These are different objects.

2. **No proximity-gap content.** The list-size question never invokes MCA/CA. It is a *pure list-decoding (combinatorial)* question about $\Lambda(C^{\equiv m},\delta)$. Steps (i)–(ii) bracket it with zero coding-theory machinery beyond the Johnson and Elias bounds, neither of which is an MCA theorem, and neither of which Thm 1.9 obstructs. Thm 1.9 bites only when one tries to convert a *proximity-gap* statement *into* list-decoding — the opposite direction.

3. **The MCA negatives do not lower $\delta^\ast_C{}^{(2)}$.** Kambiré / BCHKS Thm 1.13 / CS25 Cor. 1 bound $\varepsilon_{\mathrm{ca}}$ (the probability a random line point is close without correlated agreement) from *below*; they say nothing about how many codewords sit in a single ball. The only list-side lower bound that bites is the Elias volume bound (Step (ii)), which caps $\delta^\ast_C{}^{(2)}$ at $R_{\mathrm{cap}}$ — a list fact matching the list question. So the MCA negatives lower the *MCA* threshold (sub-problem 1), not the *list-size* threshold (sub-problem 2).

**The cleanest demonstration the two challenges are not equivalent:** for $|F|=2^{256}$, sub-problem 2 is *proven* to $J-o(1)$ and *bracketed* below $R_{\mathrm{cap}}$ by rigorous Johnson + Elias bounds, *with no progress on the MCA barrier whatsoever* — while the MCA challenge (sub-problem 1) at the same field is itself only bracketed, with its positive frontier stuck at the Johnson radius $0.293$ (§5). The list-size challenge is the more tractable of the two: it is **decoupled from Thm 1.9, with a standard list-bound bracket; for the deployed regime $|F|>2^{128}$ it is proven to $J-o(1)$ and bracketed below $R_{\mathrm{cap}}$ — and pinned (proven) at the knife-edge $|F|=2^{128}$ (the UD supremum $(1-\rho)/2$; endpoint not attained, discreteness — Prop. (B)) — with the upper reach to $\approx R_{\mathrm{cap}}$ a separate (large-list, smooth-domain) list-decoding question that is open at cryptographic scale.**

---

## 5. The MCA bracket in brief

For sub-problem 1 only the *bracket* $[\delta_{\text{known-positive}},\,\delta_{\text{known-negative}}]$ is known; the true $\delta^\ast_C$ lies inside it, and pinning it is open. **This is a bracket, not a theorem-complete determination.** Its known *endpoints* are field-type insensitive — the positive endpoint is proven over all fields, and the negative *mechanism* is field-agnostic (below) — but this must not be read as a field-agnostic solution strategy: the positive beyond-Johnson program is *not* field-agnostic (§6).

**Positive side — up to the Johnson radius, all fields [PROVEN, radius-level].** Bordage–Chiesa–Guan–Manzur (eprint 2025/2051, Thm. 9.2) give, for any domain including smooth subgroups and integer $m\ge3$,
$$
\varepsilon_{\mathrm{mca}}(C,\delta)\;\le\;\frac{(m+\tfrac12)^7\,n^2\,d}{3\,\rho^{3/2}\,|F|},\qquad\text{valid for }\delta\le1-\big(1+\tfrac1{2m}\big)\sqrt\rho,
$$
loss-free, with validity window $\to J$ as $m\to\infty$; BCHKS (eprint 2025/2055, Thm. 1.5) give a linear-in-$n$ bound to $J-\eta$ whose leading constant $C_\rho$ is an unspecified $O_\rho(\cdot)$ (we treat it as unverified until pinned). So $\delta_{\text{known-positive}}=J$ **as a radius-level theorem**. Whether the error meets $\varepsilon^\ast=2^{-128}$ is a separate *field-size* question: the verified $n^2d$ numerator costs $\approx2\log_2 n+\log_2(\rho n)$ bits, so only $\ge256$-bit fields certify $2^{-128}$ from the single-code bound; 31/64/128-bit deployments reach 128-bit soundness via the query term plus repetition/PoW and large-extension soundness arguments, not from a single-code MCA certificate.

**Negative side — near capacity [PROVEN constructions; assessment per field].** Crites–Stewart (eprint 2025/2046, Cor. 1) prove $\varepsilon_{\mathrm{ca}}=1$ in a window reaching down to the list-decoding capacity radius $R_{\mathrm{cap}}$, so $\delta^\ast_C$ is provably strictly below Singleton capacity at every field size. Kambiré (arXiv:2604.09724, Thm. 1; fleshing out the Krachun–Kazanin sketch) constructs, over smooth subgroups of prime fields $p\equiv1\bmod n$, explicit monomial lines $f=X^{rm}$, $g=X^{(r-1)m}$ with $\ge n^C$ bad scalars at $\delta=(1-\rho)-2/s$, $s=\Theta(\log n)$ — a direct proximity-gap and CA failure near capacity. At deployed field sizes the gap $2/s$ is a *constant*, not $o(1)$; the per-field instantiation of this ceiling, with its exact threshold arithmetic and status tags, is the subject of the companion ledger.

**Field-agnosticity of the negative mechanism [ESTABLISHED-MODULO-FORMALIZATION].** The Kambiré/KK25 counterexample mechanism is **field-agnostic**: the firing condition, the bad scalar, and the distinct-bad-scalar count are *characteristic-zero cyclotomic invariants* — reductions mod the characteristic $p$ of fixed elements of $\mathbb{Z}[\xi_s]$, independent of the extension degree $e$ — so genuine odd-characteristic extension fields $F_{p^e}$ realize the *same* count as the prime field $F_p$, for all primes outside the finite bad set excluded by the distinctness calibration (the KK25 hypothesis $p>\varphi(s)^{\varphi(s)}$). This is verified exactly ($10/10$ unsaturated cases; decoder-free certified bad lines on genuine $GF(31^2)$, $GF(127^2)$). The obstruction is **multiplicative** (cyclotomic subset sums), not additive, so the "no additive subspace" hope buys extension fields nothing: extensions are neither a haven nor a special vulnerability. The statement and its proof, in compact form:

> **Theorem (cyclotomic field-agnosticity of the negative mechanism) [ESTABLISHED-MODULO-FORMALIZATION].** *Fix the Kambiré configuration with quotient parameter $s$. There is a finite set of primes $B(s,r)$ — containing the primes dividing the norms of all nonzero firing-test values $e_2(T)$ and of all nonzero scalar differences $e_1(T)-e_1(T')$, over the finitely many index subsets $T,T'$ — such that for every prime $p\notin B(s,r)$ and every $e\ge1$, the firing condition, the bad scalar, and the distinct-bad-scalar count of the construction over $F_{p^e}$ coincide with their characteristic-zero values, independently of $e$.*
>
> *Proof.* The firing data are elementary symmetric functions $e_1,e_2$ of $s$-th roots of unity: fixed elements of $\mathbb{Z}[\xi_s]$, evaluated in $F_{p^e}$ via reduction modulo a prime $\tilde{\mathfrak{p}}$ above $p$. A fixed nonzero element $\sigma\in\mathbb{Z}[\xi_s]$ vanishes under reduction only if $\tilde{\mathfrak{p}}\mid\sigma$, hence only if $p\mid N_{\mathbb{Q}(\xi_s)/\mathbb{Q}}(\sigma)$ — a fixed nonzero rational integer. Applying this to each nonzero firing-test value $e_2(T)$ preserves the firing condition (non-firing subsets cannot become firing), and applying it to each nonzero difference $e_1(T)-e_1(T')$ preserves distinctness of the bad scalars (distinct firing scalars cannot collide). There are finitely many subsets and pairs, hence finitely many norms; $B(s,r)$ is the set of their prime divisors. Outside $B(s,r)$, every characteristic-zero distinction survives reduction; and since $\xi_s$ already lives in $F_p(\xi_s)\subseteq F_{p^e}$, the configuration — hence each count — is independent of the extension degree $e$. $\square$

The formalization residual (why EMF rather than PROVEN): a clean, citable uniform effectivity bound on $B(s,r)$ at every deployed $s$. The exact computational verification ($10/10$ unsaturated cases; certified bad lines on genuine $GF(31^2)$, $GF(127^2)$), the per-field ceiling values (31/64/128/256-bit), the M31/BabyBear circle-domain interpretation, and the 256-bit rescue-path assessments are recorded in the companion ledger (§2.2, §3.2).

**Scoping rule (load-bearing).** Three different "field" statements must not be conflated: (a) the known bracket *endpoints* are field-type insensitive; (b) the *negative mechanism* is field-agnostic in the precise sense above; (c) the *positive* beyond-Johnson program is **not** field-agnostic — it must be attacked first for prime-field multiplicative subgroups, with the extension-field analogue a separate open problem (§6).

---

## 6. Open problems

**The MCA grand challenge.** We do not claim a solution: $\delta^\ast_C$ remains open inside $[J,\,R_{\mathrm{cap}})$ for all fields. The genuine positive frontier is the following keystone, which simultaneously advances both grand challenges.

> **Sub-lemma P$'$ (smooth-domain small list just beyond Johnson — the positive keystone).** *For $L=\langle\omega\rangle\subset F_p$, $|L|=n$, rate $\rho$, there is $\varepsilon_0=\varepsilon_0(\rho)>0$ with*
> $$|\Lambda(\mathrm{RS}[F_p,L,\rho n],\,J+\varepsilon_0)|\;\le\;O_\rho(1/\varepsilon_0)\quad(\text{constant in }n).$$

Via the GG25 line-decoding bootstrap (GG25 Thm. 3.5: line-decodability implies $\varepsilon_{\mathrm{mca}}\le a/|F|$), P$'$ would supply the first provable smooth-domain MCA point strictly beyond Johnson — weaker than the full CGHLL Conjecture 2 (line-decodability to $R_{\mathrm{cap}}$), hence dodging the hardest part of the barrier. **Field scope: P$'$ is *not* field-agnostic.** As stated it is a **prime-field multiplicative-subgroup** statement: the char-2 analogue is *false at Johnson* (BCHKS Cor. 1.7, subspace polynomials), and the genuine odd-characteristic extension-field analogue is **open** and must be stated separately (extensions contain $F_p$-subspaces, so the "no additive subspaces" lever the prime route relies on is unavailable). Only the negative endpoint is field-agnostic (§5); the positive program is prime-fields-first.

**The beyond-Johnson smooth-domain list question.** A single open input governs the list-size bracket's upper reach to $\approx R_{\mathrm{cap}}$ for $|F|>2^{128}$: can the interleaved list be proven $\le B^{1/m}$ in the *worst case* for *smooth-domain* RS at radii between $J$ and $R_{\mathrm{cap}}$? (It does **not** bear on the knife-edge $|F|=2^{128}$ *value*, which is already pinned at the UD supremum, Prop. (B): at $B=1$ any list of size $2$ breaks the budget, and size-$2$ lists exist just past UD by the MDS midpoint argument.) This is the smooth-domain interleaved-list conjecture — equivalent to derandomizing the random-RS list bound (ABF Thm. 3.6) to the structured smooth domain — ABF open problem §7.9. It is itself a *list-decoding* statement, so it too dodges Thm 1.9; it is gated instead by the different (also hard) open problem of RS list-decoding on structured domains. Exact small-field enumeration (recorded in the companion ledger) supports its strong constant-list form at tiny fields; no proof exists for cryptographic fields, and the tiny-field experiments cannot see the cryptographic-scale deep-hole regime.

---

## Acknowledgments

This work was carried out with substantial assistance from AI research tools — primarily Anthropic's Claude, with independent adversarial review passes by OpenAI's GPT models. All load-bearing claims are backed by exact-arithmetic verification scripts in the accompanying repository (<https://github.com/demirelo/rs-proximity-bracket>); the status tags reflect audited verification states, not model assertions. Repository artifacts: `problem-ledger.md` (formal bracket ledger); `listsize-analysis.md` + `calculator/listsize_resolution.py` (sub-problem-2 bracket analysis); `negative-endpoint-ledger.md` (the companion research ledger to this paper); `experiments/small_rs_atlas/` (exact enumeration).

---

## References

- **[ABF]** G. Arnon, D. Boneh, G. Fenzi, *Open Problems in List Decoding and Correlated Agreement*, IACR eprint **2026/680** (Apr 8 2026). The survey framing the Proximity Prize; Defs. 2.9, 2.12, 4.1, 4.3; Cor. 3.3; Lemma 6.6; Thms. 4.12, 4.17, 4.21, 5.2, 5.3; open problem §7.9.
- **[BCGM / Bordage–Chiesa]** S. Bordage, A. Chiesa, Z. Guan, I. Manzur, *All Polynomial Generators Preserve Distance with Mutual Correlated Agreement*, IACR eprint **2025/2051** (rev. May 19 2026). Thm. 9.2 (positive MCA to Johnson, all polynomial generators); Lemma 10.1 (interleaving).
- **[BCHKS]** E. Ben-Sasson, D. Carmon, U. Haböck, S. Kopparty, S. Saraf, *On Proximity Gaps for Reed–Solomon Codes*, IACR eprint **2025/2055** (Nov 6 2025). Thm. 1.5 (positive, zero-loss to Johnson); Thm. 1.9 (the barrier); Cor. 1.7 (char-2 failure at $J$); Thm. 1.13 + Conj. 1.12 (prime mult.-subgroup; M31 unconditional).
- **[CS25]** E. Crites, A. Stewart, *On Reed–Solomon Proximity Gaps Conjectures*, IACR eprint **2025/2046** (Dec 19 2025). Thm. 1 / Lemma 1 / Cor. 1 (Elias volume bound; $\varepsilon_{\mathrm{ca}}=1$ strip); Thm. 7.4.1, Claim 1 (list-decoding capacity, $1/\log_2 q$ gap); Thm. 2 (CA $\Rightarrow$ list-decoding).
- **[GG25]** R. Goyal, V. Guruswami, *Optimal Proximity Gaps for Subspace-Design Codes and (Random) Reed–Solomon Codes*, IACR eprint **2025/2054** (rev. Mar 2026). Def. 3.1 (curve-decodability); Thm. 3.5 (line-decode $\Rightarrow$ MCA); Cors. 4.9/4.10, Thm. 5.16 (capacity for folded/subspace-design/random RS); Lemma 8 (the list-size bootstrap).
- **[Kambiré]** A. Kambiré, *Proximity Gaps Conjecture Fails Near Capacity over Prime Fields*, **arXiv:2604.09724** (Apr 2026). Thm. 1 (explicit smooth-prime counterexample at $(1-\rho)-\Theta(1/\log n)$); the $f=X^{rm},g=X^{(r-1)m}$ line; flesh-out of the Krachun–Kazanin (KK25) sketch.
- **[CGHLL26]** D. Carmon, L. Goldberg, U. Haböck, L. Lerer, I. Lesokhin (+ S. Papini, S. Samocha), *S-two Whitepaper*, IACR eprint **2026/532** (2026). App. A.5: the list-decoding capacity radius; Conjecture 1 (list-decodability to it); **Conjecture 2** (line-decodability to it, incl. extension-field clause); Thm. 37 (= KK25); Lemma 9 (distinct-sums count).
- **[Companion ledger]** Ö. Demirel, *The Proximity Prize Negative-Endpoint Ledger*, `negative-endpoint-ledger.md`, <https://github.com/demirelo/rs-proximity-bracket> (2026-06-10). Per-field negative-ceiling assessments, M31/BabyBear circle-domain interpretation, 256-bit rescue paths, N1/N2 close-outs, experimental data.
