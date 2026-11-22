# V15.3 Validation Status - Final Portfolio Release

## User-run MATLAB evidence

### Regression suite

- The project owner reported that the V15.3 test run completed successfully.
- The bundled test entry point is `run_tests_v15_3` and the test source contains 160 test assertions, ending with `All 160 tests passed successfully.`
- The full console transcript of that test run was not included in the final evidence supplied here; the successful status is recorded as user-reported.

### 30-frame-per-condition confirmation run

The project owner supplied the completed MATLAB console summary with these dimensions:

- Frames per condition: 30.
- Conditions: 36 (2 channels x 2 spatial profiles x 3 Doppler values x 3 pilot/data SNR values).
- Summary rows: 252.
- Candidate frame-log rows: 7,560.
- Joint burst-decision rows: 3,240.
- Reference burst-prediction rows: 3,240.
- Bootstrap resamples: 200 (from the supplied short-run invocation).

The aggregate results are included in `results/v15_3_final_summary/`. They were transcribed from the user-provided console output; the raw per-frame MAT/CSV files were not uploaded.

## Headline measured results

| Channel | Joint CRC-goodput | Best retained reference | Relative increase | Joint PER | PER-target-feasible decisions |
|---|---:|---:|---:|---:|---:|
| Rayleigh | 0.66 bit/s/Hz | 0.58 bit/s/Hz | +13.8% | 34.01% | 7.3% |
| Rician | 1.13 bit/s/Hz | 0.95 bit/s/Hz | +18.9% | 12.35% | 30.5% |

The best retained reference was the adaptive convolutional-FEC candidate. Joint CRC-goodput improved, but the aggregate PER was higher than the reference in the supplied run. No universal PER-target guarantee is claimed.

## Statistical limitations

- The 30-frame experiment is a confirmation run, not publication-grade validation.
- Worst-condition upper bounds are broad: Rayleigh payload-BER U95 = 0.2136 and PER U95 = 1.0; Rician payload-BER U95 = 0.00620 and PER U95 = 0.5667.
- The controller's one-sided Wilson-style bound is based on a fractional prior-plus-count representation and is a decision statistic, not a pure frequentist interval from raw Bernoulli observations.
- The link-model PER prior assumes independent bit errors within a packet; burst fading and residual ICI can violate this approximation.
- The FEC coding-gain curve is an initialization heuristic. Measured CRC outcomes are the feedback evidence.

## Not run / not supplied

- The default 300-frame-per-condition V15.3 experiment was not supplied; only the 30-frame confirmation run is documented here.
- MATLAB release/version and runtime were not recorded in the supplied console output.
- Raw per-frame result files were not uploaded, so the included CSVs are curated aggregate summaries only.

## Reproduction commands

Regression suite:

```matlab
clear functions
addpath(genpath(pwd));
run_tests_v15_3
```

30-frame confirmation run:

```matlab
addpath('examples');
opts = struct('numFrames',30,'numBootstrap',200);
results = pilot_mapped_adaptive_validation_run_v15_3(opts);
summary = pilot_mapped_adaptive_results_summary_v15_3(results);
```

Default 300-frame run:

```matlab
results = pilot_mapped_adaptive_validation_run_v15_3;
summary = pilot_mapped_adaptive_results_summary_v15_3(results);
```
