# Literature Set A — Cross-Paper Synthesis (Proximity Prize)

Four papers, extracted verbatim into sibling files: `abf-openproblems.md`, `crites-stewart.md`,
`bchks.md`, `goyal-guruswami.md`. All PDFs fetched via `curl` (WebFetch was 403-blocked by
eprint) and parsed with `pdftotext -layout`.

---

## The prize's exact success criterion (per the ABF survey, eprint 2026/680, Apr 8 2026)

The Ethereum Foundation **Proximity Prize** (https://proximityprize.org/) is framed by ABF as
**two "grand challenges"** for the deployed code `C = RS[F, L, k]` with **L a smooth domain**
(multiplicative subgroup / coset of order a power of two), constant rate
`ρ = k/|L| ∈ {1/2, 1/4, 1/8, 1/16}`, and a target `ε* = 2^-128`, in the regime
`L ⊆ F` smooth, `k ≤ 2^40`, `|F| < 2^256`:

1. **Grand MCA challenge:** determine the **largest `δ*_C` such that `ε_mca(C, δ*_C) ≤ ε*`**,
   *with a proof that for all `δ > δ*_C`, `ε_mca(C, δ) > ε*`*. (`ε_mca` = mutual-correlated-
   agreement error w.r.t. lines, ABF Def 4.3, loss-free.)
2. **Grand list-decoding challenge:** for a constant interleaving `m`, determine the **largest
   `δ*_C` such that `|Λ(C^{≡m}, δ*_C)| ≤ ε*·|F|`**, with the analogous "for all `δ>δ*_C` it
   fails" proof. (No efficient decoder required — only the value `δ*_C`.)

Both matter because the deployed protocols' (FRI/STIR/WHIR) round-by-round knowledge soundness
error has the exact form (ABF Lemma 6.6) `max( ε_mca(C,δ) + |Λ(C^{≡2},δ)|/|F| , (1-δ)^t )`.
The known answer sits **at the Johnson radius** `1-√ρ` (where positive RS knowledge runs out,
ABF Thm 4.12 / BCHKS Thm 1.5, error `O_ρ(n/(η^5|F|))`); the prize asks whether `δ*_C` can be
pushed toward the Singleton capacity `1-ρ`. ABF's §6.3.1 quantifies the stakes: at ρ=1/2, plain
interleaved RS stuck at the Johnson radius gives a per-query factor `(1-δ)^t ≈ (1/√2)^128 = 2^-64`
— so **provably only ≈64 bits** of security — whereas reaching capacity `δ→1-ρ=1/2` would give
`(1/2)^128 = 2^-128`. Closing that 64→128-bit gap for the deployed smooth-domain code is the prize.

---

## Cross-paper comparison table

| Paper | Best POSITIVE result (δ + error form + code class) | Best NEGATIVE / lower bound (δ + code class) | Applies to smooth multiplicative-subgroup RS? | Key constants |
|---|---|---|---|---|
| **ABF** (2026/680, survey) | Aggregates: smooth RS up to Johnson `1-√ρ-η`, `ε_mca ≤ O_ρ(n/(η^5\|F\|))` (Thm 4.12); folded/subspace-design & random RS up to capacity `1-ρ-η` (Thm 4.13–4.15) | Aggregates: prime smooth RS within `Θ(1/log n)` of capacity has `ε_ca ≥ n^c/\|F\|` (Thm 4.16); char-2 RS just past Johnson `ε_ca ≥ n^{2(1-ε)}/\|F\|` (Thm 4.18); `ε_ca=1` strip (Thm 4.17) | **Y** — this is the survey defining the smooth-RS target | `ε*=2^-128`; ρ∈{1/2,1/4,1/8,1/16}; k≤2^40; \|F\|<2^256; Johnson `1-√ρ`; capacity `1-ρ` |
| **Crites–Stewart** (2025/2046, NEG) | None (a negative paper). Proposes modified conjectures up to **list-decoding capacity** `1-H_q(δ)` (loss `≤1/log₂q`) | **CA fully fails (`ε_ca=1`)** in the window `1-H_q(δ)+2/n+√((H_q(δ)-δ)/n) ≤ ρ ≤ 1-δ-2/n` (Cor 1); CA/MCA/list-decoding up-to-capacity (`1-ρ`) conjectures all FALSE; ceiling is `1-H_q(ρ)` not `1-ρ` | **partial** — uses general domain D + deep holes (exist on any domain incl. smooth); needs `q ≥ n`; smoothness not singled out — VERIFY | `1/log₂q` gap (≤1/31 for 31-bit fields); `H_q(δ)≈δ-1/log q`; Thm 2 needs `ε<(q-n)/kq≈1/k` |
| **BCHKS** (2025/2055) | smooth RS up to Johnson `1-√ρ-η`, `ε_mca=0`-loss with `a=O_δ(n/η^5)` ⇒ error `O_ρ(n/(η^5\|F\|))` (Thm 1.5) — best constant for deployed RS | char-2 RS: `n^τ`-conj FALSE, `γ=δ-Ω_τ(1)` needs `a≥n^{τ-o(1)}` (Thm 1.6), `γ=J(δ)` needs `a≥n^{2-ε}` (Cor 1.7); **prime mult.-subgroup RS: `γ=δ-Θ(1/log n)` forces loss `Θ(1/log n)`** (Thm 1.13, cond. on Conj 1.12; M31 unconditional) | **Y (both directions)** — Thm 1.5 positive for any domain; Thm 1.13 negative for prime + multiplicative subgroup (deployed-type) | M31: `f,g` 1/2-close ∀z yet CA-dist `≥1/2+1/62≈0.516`; (M31)^4: `δ≈0.508, γ=1/2, ε∗≈0.004`; Thm 1.9: `γ>LDR ⇒ a/q≥1/2n` |
| **Goyal–Guruswami** (2025/2054, POS) | **Capacity** `1-R-η`, MCA: folded RS `err·\|F\|≳n/η+1/η^3`, `s≳1/η^2` (Thm 1.2/Cor 4.10); subspace-design (Cor 4.9); random RS `err·q≳n/η+1/η^5` (Thm 1.3/5.16) — field LINEAR in n (subspace-design) | None (a positive paper) | **N** — covers folded RS, univ. multiplicity, subspace-design, RANDOM RS (random points), RLC/LDPC/AEL; NOT plain smooth RS | line-decode⇒MCA (Thm 3.5); folded RS alphabet `s≳1/η^2`; random RS needs `q>n^{C2}` (formal) |

---

## Where δ*_C sits for the DEPLOYED smooth-domain RS (the prize's actual question)

Combining the four papers, for **plain RS over a smooth multiplicative-subgroup domain over a
prime field** (the deployed case):

- **`δ_known_positive = 1 - √ρ - η`** (the Johnson radius minus a slack), with
  `ε_mca ≤ O_ρ(n/(η^5·|F|))` and **zero proximity loss** — BCHKS Theorem 1.5 (= ABF Thm 4.12).
  This is the BEST proven positive threshold for the deployed code, and it is *below* the
  Singleton capacity. To certify `ε_mca ≤ 2^-128` you need `|F| ⪆ 2^128 · n/η^5`; with 256-bit
  fields and `n ~ 2^{20..30}` there is room, so a small η suffices and `δ*_C ≥ 1-√ρ-η`.
- **`δ_known_negative`** for the deployed code:
  - **Ceiling at the list-decoding capacity** `1 - H_q(ρ) ≈ (1-ρ) - 1/log₂q`: Crites–Stewart
    Corollary 1 gives `ε_ca = 1` in a window reaching down to `1-H_q(ρ)`. So `δ*_C` is provably
    **strictly below `1-ρ`** — the naive "up to capacity" target is dead.
  - **Within `Θ(1/log n)` of capacity, prime smooth RS has `ε_ca ≥ n^c/|F|` ≫ 2^-128**: ABF
    Theorem 4.16 (= BCHKS + KK25) — captures prime fields AND smooth domains explicitly.
  - **Conditional/explicit BCHKS Theorem 1.13:** prime multiplicative-subgroup RS at
    `γ = δmin - Θ(1/log n)` must carry proximity loss `Θ(1/log n)`; the **Mersenne-31, ρ≈1/2
    instantiation is UNCONDITIONAL** and exhibits a genuine proximity-gap failure exactly at
    radius `1/2` (`f,g` that are 1/2-close for every `z`, yet CA-distance `≥ 0.516`).

### At ρ = 1/2 specifically (the headline case)
- Johnson radius `1-√(1/2) ≈ 0.293`; Singleton capacity `1-ρ = 0.5`; list-decoding capacity
  `1-H_q(1/2) ≈ 0.5 - 1/log₂q ≈ 0.468` for q≈2^31.
- **`δ_known_positive ≈ 0.293`** (Johnson, BCHKS Thm 1.5, ε_mca = O(n/η^5/|F|), loss 0).
- **`δ_known_negative`:** an *unconditional* CA failure at `δ = 0.5` (BCHKS Thm 1.13, M31), and a
  general ceiling at `≈ 0.468` (list-decoding capacity, Crites–Stewart). So for ρ=1/2 the open
  gap `[δ_known_positive, δ_known_negative]` is roughly **`[0.293, 0.468]`** (positive side firm;
  the upper side has the M31 failure pinned at 0.5 and the asymptotic Θ(1/log n) ceiling from
  ABF Thm 4.16).
- **Crucially, NO paper proves a positive `δ*_C > 0.293` for the deployed code.** Goyal–Guruswami
  reaches `≈0.5-η` but ONLY for FOLDED RS / random RS, not plain smooth RS. So whether plain
  smooth RS's `δ*_C` is at the Johnson radius (≈0.293), somewhere in the strip, or up near the
  list-decoding capacity (≈0.468) is **exactly the open prize question.**

---

## The decisive tension (what the prize must resolve)

- **Goyal–Guruswami** proves RS-*like* codes (folded, random-evaluation-point) reach the Singleton
  capacity `1-R` for MCA — there is NO information-theoretic barrier for RS-flavored codes in
  general, and **line-decodability ⇒ MCA** is the clean route (Thm 3.5).
- **Crites–Stewart + BCHKS** prove the *specific* structured/deployed RS cannot reach `1-ρ`:
  CA fails at the list-decoding capacity `1-H_q(ρ)` (CS), and prime smooth-subgroup RS has
  unavoidable loss within `Θ(1/log n)` of capacity (BCHKS Thm 1.13, M31 unconditional).
- **BCHKS Theorem 1.9** ties the knot: improving MCA/proximity-gaps beyond the Johnson radius for
  *any* RS code REQUIRES improving its list-decoding radius (for list size q) beyond Johnson — a
  long-open, in-general-false-for-some-RS problem. This is why ABF poses BOTH grand challenges and
  asks whether good list-decoding implies good MCA (open; only the sqrt-lossy ABF Thm 5.1 known).
- **The prize's core bet** (per `RESEARCH_PLAN.md`) — characterize when close RS codewords on a
  line must align as a *line of codewords* for smooth domains — is precisely asking whether plain
  smooth RS is **line-decodable** near capacity. Goyal–Guruswami is the template for what a
  positive answer buys (MCA up to capacity); BCHKS Thm 1.13 / CS Cor 1 are the candidate
  counterexamples that would instead pin `δ*_C` down. See also CGHLL26 (Carmon–Goldberg–Haböck–
  Lerer–Lesokhin, cited by ABF for line-decoding conjectures) — NOT in this set, worth fetching.

---

## Notation conflicts across the four papers (carry into all downstream work)

| Symbol | ABF (2026/680) | Crites–Stewart (2025/2046) | BCHKS (2025/2055) | Goyal–Guruswami (2025/2054) |
|---|---|---|---|---|
| proximity radius | **`δ`** | **`δ = f/n`** (`f` = abs. errors) | **`γ`** | **`δ`** (`= 1-R-η`) |
| code min-distance | `δmin(C)` (= "capacity") | `1-ρ` (Singleton) | **`δ`** (= `1-k/n`) | `1-R` (Singleton) |
| soundness error | `ε_pg/ε_ca/ε_mca` | `ε` (proximity-gap error) | **`a/q`** | **`err`** |
| proximity loss | (none in `ε_mca`; CA has `δint-δfld`) | — | **`ε∗`** | **`γ`** (conclusion slack) |
| rate | `ρ` | `ρ = k/n` | `ρ = k/n = 1-δ` | `R` |
| combination | line `f1+γf2` (`F_lines`) | line `u0+z·u1` | line `f+zg` | power curve `Σ u_j α^j`, degree `ℓ` |
| "capacity" | Singleton `1-ρ` | **list-decoding** `1-H_q(ρ)` (the key point) | Singleton `1-δ` | Singleton `1-R` |

**The single most dangerous conflict:** BCHKS's `δ` is the *minimum distance*, while ABF/CS/GG's
`δ` is the *radius*; and BCHKS/GG's `γ` mean opposite things (radius vs loss). Crites–Stewart's
"capacity" is the *list-decoding* capacity `1-H_q(ρ)` — strictly below the Singleton `1-ρ` that
the other three call capacity; this `1/log₂q ≤ 1/31` gap is itself a load-bearing quantity.

## Caveats / human-verification flags (see per-paper "INFERRED — VERIFY" sections)
- BCHKS Thm 1.13 general prime-field result is **conditional on additive-combinatorics
  Conjecture 1.12**; only M31 / (M31)^4 are unconditional, and `F_q^*` (order `q-1`, not a power
  of 2) is not *strictly* smooth per ABF Def 2.12 — VERIFY whether the negative construction
  survives on a strictly power-of-two-order subgroup (matters for hitting the FFT domain exactly).
- Goyal–Guruswami random-RS Thm 5.16 formally needs `q > n^{C2}` (poly larger than n), which the
  deployed Mersenne-31/BabyBear FFT field does NOT satisfy — the "field linear in n" claim is for
  the folded/subspace-design families, not random RS as formally stated.
- Per-rate numeric δ values (Johnson `1-√ρ`, capacity `1-ρ`, list-dec. capacity `1-H_q(ρ)`) are my
  arithmetic; none of the papers tabulate them per rate.
