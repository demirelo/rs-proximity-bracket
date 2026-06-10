# Small-Field Reed–Solomon Proximity Atlas (public subset)

Exact experimental probes of whether *smooth* multiplicative-subgroup
Reed–Solomon codes admit "bad lines" (lines with many close points but no
correlated agreement) or large interleaved lists, and how their behavior
compares to **random-subset** and **full** evaluation domains over the same
field at the same rate and radius.

Everything is computed with **exact full-codeword enumeration** (no decoder,
no list-decoder), so the measured distances and list sizes are ground truth.
Hand-rolled finite fields, `numpy` + `multiprocessing`; no other dependencies
(except `sympy` for `verify_rho_half_sums.py`).

This directory is the public subset of the project's experimental atlas: it
contains exactly the scripts (and their `results/` outputs) cited by the
short paper (`technical-note-short.md`) and the ledgers
(`negative-endpoint-ledger.md`, `problem-ledger.md`, `listsize-analysis.md`).

## Library modules

```
ff.py                 finite fields: GF(p) and GF(2^m); generator + subgroup tools
ff_ext.py             odd-characteristic prime-power fields GF(p^m)
rs.py                 RS code: domains, encode, exact dist_to_code, MDS check
search_bad_lines.py   bad-line search (close-count vs exact common-agreement S*)
search_interleaved.py interleaved C^{≡m} list-size search (m=2,3)
n2_hardening.py       N2 structural certificate + exact prime r=4 control
```

## Experiment scripts (each writes JSON/CSV under `results/`)

```
run_atlas.py                  orchestrates the smooth/random/full domain battery (parallel)
scaleup_smooth_vs_random.py   larger-field smooth-vs-random domain comparison
counterexample_kambire.py     exact reproduction of the Kambiré coset-union bad list
counterexample_extension.py   the same mechanism over odd-char extension fields GF(p^m)
interleaved_onset.py          onset of interleaved (m=2,3) list growth vs radius
singlelist_past_johnson.py    exact single-code list sizes past the Johnson radius
n2_char0_count.py             characteristic-zero distinct-value count (N2), exact DP
n2_count_laws.py              count laws for the N2 mechanism (uses n2_char0_count)
n2_general_r_count.py         general-r distinct-count over prime and extension fields
n2_prize_threshold.py         prize-threshold arithmetic: count vs budget 2^{b-128}, b* boundary
n2_crosscheck.py              independent N2 cross-check (structural certificate + exact control)
verify_rho_half_sums.py       cyclotomic/resultant verification of the rho=1/2 antipodal sums (sympy)
cluster_certificate.py        non-coset "cluster" word certificates in the open band   [--selftest]
exact_center.py               exact center search around cluster certificates          [--selftest]
falsification_open_band.py    exact open-band list-size falsification battery          [--selftest]
high_moments.py               exact high-moment / Markov probe of the open-band list   [--selftest]
```

Scripts marked `[--selftest]` have a fast built-in self-test
(`python3 <script> --selftest`). The others run their full (small) battery
when invoked directly, e.g. `python3 n2_prize_threshold.py`.

`results/` holds the committed outputs the paper and ledgers cite by name
(e.g. `results/n2_prize_threshold.json`, `results/atlas_results.csv`).

Note: some script docstrings reference the authors' research notes
(`p-prime/...`); those notes are not part of this public artifact set and are
available from the authors on request. In this public copy, `high_moments.py`
writes its markdown probe report under `results/` instead.
