# rs-proximity-bracket

Verification artifacts for the IACR ePrint paper

> **Two Grand Challenges Are Not One: A Proven Bracket for the Interleaved
> List-Size Challenge at Large Fields** — Ömer Demirel.

This repository contains the paper, its companion research ledger, the
supporting analysis documents, the parameter calculator, and the exact
small-field experiments that the paper and ledgers cite by name. It is a
curated artifact set: everything cited by the short paper is here; the
authors' ongoing research notes are not.

## Contents

| Path | What it is |
|------|------------|
| `technical-note-short.md` / `.pdf` | The short paper (ePrint submission). |
| `negative-endpoint-ledger.md` / `.pdf` | Companion research ledger: per-field negative-ceiling assessments, 31/64/128/256-bit tables, M31/BabyBear circle-domain story, 256-bit rescue-path assessments, N1/N2 close-outs, experimental data pointers. |
| `problem-ledger.md` | The formal bracket ledger (claims with status tags and sources). |
| `listsize-analysis.md` | List-size challenge analysis (budget trichotomy, bracket derivation). |
| `assembled-vs-proven.md` | What is assembled from cited components vs proven from scratch. |
| `n1-rho-half.md`, `n2-verdict.md`, `snark-impact.md` | N1/N2 close-out and SNARK-impact analyses. |
| `literature/notes/` | Reading notes for the cited literature (Kambiré, BCHKS, Bordage–Chiesa, Crites–Stewart, subgroup-sumset growth survey, synthesis). |
| `calculator/` | The proximity-parameter calculator (bounds, δ\* tables, list-size resolution) with its test suite and committed outputs (`calculator/out/`). |
| `experiments/small_rs_atlas/` | Exact small-field RS experiments cited by the paper/ledgers, with committed `results/` JSON/CSV. See its README. |

**Status-tag glossary** (PROVEN / ASSEMBLED / EMF / CONDITIONALLY ESTABLISHED /
INFERRED / EXPERIMENTAL / CONJECTURED–OPEN): see `negative-endpoint-ledger.md`,
"Status-tag glossary" near the top of the document.

References to `p-prime/...` in the ledgers point to the authors' research
notes, available on request.

## Running the verification artifacts

Requirements: Python 3 with `numpy`, `pytest`, and `sympy` (only for
`verify_rho_half_sums.py`).

Calculator test suite (40 tests):

```sh
python3 -m pytest calculator/
```

Key experiments (from `experiments/small_rs_atlas/`):

```sh
python3 n2_prize_threshold.py          # prize-threshold arithmetic, b* ≈ 140 boundary
python3 counterexample_kambire.py      # exact Kambiré coset-union bad list
python3 falsification_open_band.py --selftest
python3 high_moments.py --selftest
```

Each script writes its outputs under `experiments/small_rs_atlas/results/`;
the committed JSON/CSV there are the runs cited in the paper and ledgers.

## Author / contact

Ömer Demirel — Snowfall Finance — omer@progfi.xyz

## License

- **Documents** (`*.md`, `*.pdf`): Creative Commons Attribution 4.0
  International (CC BY 4.0), <https://creativecommons.org/licenses/by/4.0/>.
- **Code** (`calculator/`, `experiments/`): MIT License — see `LICENSE`.
