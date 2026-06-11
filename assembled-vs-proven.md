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
   resultant-distinctness bound) — **CITED, NOT RE-PROVED → PARTIALLY DISCHARGED
   (2026-06-11 update)**. The published successor **KKH eprint 2026/782 (Krachun–
   Kazanin–Haböck), Lemma 1** now *proves* the cyclotomic/resultant distinctness
   statement in print: prime field, `G ≤ F_p^*` of power-of-two order `s`,
   `1 ≤ r ≤ s/2`, printed prime bound `p > s^{s/2}`, count `2^r·C(s/2, r)`
   (stronger count than the KK25 sketch's `C(φ(s), r)`). Verified against our
   usage row-by-row (`literature/notes/kkh-2026-782-verification.md`):
   the **64-bit and 128-bit rows at ρ ∈ {1/4, 1/8, 1/16} are now backed by the
   printed statement** (`2^64 > 2^32 = 16^8`; `2^128 > 2^80 = 32^16`); the
   **31-bit row** (`p ≈ 2^31 < 2^32`) is covered by the *proof's* sharper
   intermediate bound `p > (2r)^{φ(s)}` (`≤ 2^28.7` at the deployed non-N1
   rates) but **not by the printed statement** — a proof-level inference, one
   grade below a verbatim citation. Where the ρ = 1/2 close-out (sub-lemma N1)
   needs `r = s/2 + 2`, *past* both the KK25 and the KKH-782 printed ranges,
   the claim remains explicitly **CONDITIONAL** — named, numerically airtight
   (zero counterexamples through `s = 32`), unformalized; KKH-782's appendix
   proves an all-δ quotient-line counterexample *without* Lemma 1
   (asymptotic), which supports the method but does not supply the per-field
   `ρ = 1/2` row. **The residual risk surface is therefore narrowed to: the
   31-bit printed-statement gap (proof-level covered) and the N1 range
   extension.**
3. **Threshold arithmetic** (does the bad-scalar count clear `2^{b−128}`?) —
   **OURS, PROVEN** (exact integer arithmetic, regenerated tables,
   `calculator/bounds.py::kambire_smax`, 40/40 tests). The R13 s-integrality
   constraint (`s` a power of two on smooth domains) is elementary and proven.
4. **The p-decoupling VERIFY** — whether `s` can be chosen fully independent of `p`
   at every deployed size is flagged VERIFY (`kambire.md`); the per-field rows are
   stated at calibrations where the dependence is satisfied, so this affects
   *generality*, not the tabled rows.

## Exactly what would upgrade ASSEMBLED → PROVEN

- ~~Re-derive (or have a referee verify) KK25 Lemma 9 within its stated range~~
  **DONE in print for the 64/128-bit non-N1 rows** (KKH 2026/782 Lemma 1,
  2026-06-11 verification read). Remaining for full upgrade: a stated (not
  proof-mined) distinctness theorem covering 31-bit primes at `s = 16` — i.e.,
  the fixed-`r` calibration `p > (2r)^{φ(s)}` promoted from KKH-782's proof to
  a citable statement (or a one-page re-derivation, which the proof makes
  routine).
- Formalize the named out-of-range extension (`r = s/2 + 2`): sub-lemma N1's
  conditional closure at ρ = 1/2 then becomes unconditional.
- Neither upgrade changes any tabled number; the values are already computed at the
  calibrations the lemma's stated range supports (the conditional extension affects
  only the ρ = 1/2 close-out's *status*, not its value).

## What cannot break

The Johnson floor (positive side), the budget trichotomy, the list-size bracket, and
the Thm-1.9 separation are independent of all four inputs above — they cite only
ABF Cor. 3.3 and CS25 Thm. 1/7.4.1 inside stated ranges. A total failure of the
assembled ceiling would widen the negative-side open gap; it would not touch the
short paper's centerpiece claims.
