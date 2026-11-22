# V8 Validation Status

- V7.1 user-reported tests: 35/35 passed.
- V7.1 spatial-correlation study: 30 rows completed (5 profiles x 2 channels x 3 SNRs), 500 frames/run and 1,000 bootstrap resamples.
- V8 source adds joint spatial-frequency LMMSE, posterior uncertainty, paired estimator comparison and tests 36-42.
- V8 has not been run in MATLAB yet. Do not report V8 test or performance success until `run_tests` and `spatial_aware_lmmse_comparison_run` have been executed in MATLAB.
