# V15.1.1 Validation Status

- Static inspection confirms the ICI innovation estimator now supplies an
  `[N x Nr x 1 x Ns]` prior slice matching each isolated Tx observation.
- Test 123 covers the exact two-Tx full-prior condition reported by the user.
- Expected suite total: 123 tests.
- MATLAB/Octave execution is not available in the build environment; user-run
  MATLAB validation is required before the suite can be marked passed.

Run from the extracted project root:

```matlab
clear functions
addpath(genpath(pwd));
which run_tests_v15_1_1 -all
run_tests_v15_1_1
```

After all tests pass, run the existing V15.1 adaptive validation:

```matlab
addpath('examples');
results = pilot_mapped_adaptive_validation_run_v15_1;
summary = pilot_mapped_adaptive_results_summary_v15_1(results);
```
