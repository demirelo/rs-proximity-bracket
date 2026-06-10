# Why "ASSEMBLED" Is Not "PROVEN": the Per-Field Negative Ceiling's Dependency Graph

**One-page referee companion** to `negative-endpoint-ledger.md` §2.4 (the per-field
negative-ceiling assessment). Purpose: make the gap between the ASSEMBLED tag and a
PROVEN tag exact, so a red-team review can target it directly.

## The assembled statement

> For deployed field size `b = log2|F| ∈ {31, 64, 128}` and rate `ρ`, the MCA
> threshold satisfies `δ*_C ≤ (1−ρ) − 2/s_max(b)` with `s_max(b)` a power of two
> (`16/16/32`), realized by the Kambiré monomial-line construction on a smooth
> subgroup of order `s·m`. At `b = 256` no threshold-established ceiling of this
> type exists.

## The four inputs and their individual grades

1. **Kambiré's construction and bad-scalar count** (arXiv:2604.09724, Thm. 1) —
   **PROVEN (cited)**. Used verbatim inside its stated hypotheses. Not a risk
   surface beyond ordinary citation trust; the construction is also reproduced
   computationally at small scale (`experiments/small_rs_atlas/counterexample_kambire.py`).
2. **KK25 Lemma 9 distinctness** (the `{−1,0,+1}` cyclotomic-independence /
   resultant-distinctness bound) — **CITED, NOT RE-PROVED**. This is the load-bearing
   import. Where we use it *inside* its stated range (`r ≤ φ(s)/2`, prime calibration
   `p > φ(s)^{φ(s)}`), the only gap is that we have not re-derived its proof.
   Where the ρ = 1/2 close-out (sub-lemma N1) needs it *past* its stated range,
   the claim is explicitly **CONDITIONAL** — that extension is named, numerically
   airtight in every tested case (zero counterexamples through `s = 32`), and
   unformalized. **This is the single place the assembly could break.**
3. **Threshold arithmetic** (does the bad-scalar count clear `2^{b−128}`?) —
   **OURS, PROVEN** (exact integer arithmetic, regenerated tables,
   `calculator/bounds.py::kambire_smax`, 40/40 tests). The R13 s-integrality
   constraint (`s` a power of two on smooth domains) is elementary and proven.
4. **The p-decoupling VERIFY** — whether `s` can be chosen fully independent of `p`
   at every deployed size is flagged VERIFY (`kambire.md`); the per-field rows are
   stated at calibrations where the dependence is satisfied, so this affects
   *generality*, not the tabled rows.

## Exactly what would upgrade ASSEMBLED → PROVEN

- Re-derive (or have a referee verify) KK25 Lemma 9 within its stated range: the
  31/64/128-bit rows then become PROVEN as stated.
- Formalize the named out-of-range extension: sub-lemma N1's conditional closure at
  ρ = 1/2 then becomes unconditional.
- Neither upgrade changes any tabled number; the values are already computed at the
  calibrations the lemma's stated range supports (the conditional extension affects
  only the ρ = 1/2 close-out's *status*, not its value).

## What cannot break

The Johnson floor (positive side), the budget trichotomy, the list-size bracket, and
the Thm-1.9 separation are independent of all four inputs above — they cite only
ABF Cor. 3.3 and CS25 Thm. 1/7.4.1 inside stated ranges. A total failure of the
assembled ceiling would widen the negative-side open gap; it would not touch the
short paper's centerpiece claims.
