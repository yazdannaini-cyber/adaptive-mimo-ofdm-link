# V15.2.1 — Soft-Decision FEC

V15.2.1 preserves V15.2's statistical and packet-level experiment but upgrades the coded receiver from hard-decision to soft-decision Viterbi.

- Max-log soft QAM LLRs for Gray 4/16/64-QAM.
- Per-symbol LLR variance is derived from MMSE or Alamouti post-equalization SINR.
- Terminated `(7,5)` rate-1/2 Viterbi uses soft branch metrics.
- CRC-16 is checked after decoding.
- Adds `FECDecoder` metadata and a separate result output folder.

## Run tests
```matlab
clear functions
addpath(genpath(pwd));
which run_tests_v15_2_1 -all
run_tests_v15_2_1
```
Expected: 141 tests.

## Run full experiment
```matlab
addpath('examples');
results = pilot_mapped_adaptive_validation_run_v15_2_1;
summary = pilot_mapped_adaptive_results_summary_v15_2_1(results);
```
Expected: 216 summary rows, 21,600 frame-log rows, and 10,800 burst-decision rows at the default 100 frames/condition. Output is written to `pilot_mapped_results_v15_2_1`.

The max-log LLR receiver assumes a scalar Gaussian residual-disturbance model per equalized stream. This is more informative than hard decisions but remains approximate under strongly correlated interference/ICI.
