# Covariance-Aware Detection — Run Instructions

From the project root in MATLAB:

```matlab
run_tests
```

V10 expects 54 tests. Then run the paired end-to-end comparison:

```matlab
addpath('examples');
results = spatial_aware_lmmse_comparison_run;
```

The joint estimator arm now passes its full per-tone posterior receive covariance to both the linear-MMSE spatial detector and the Alamouti combiner. The per-link estimator arm retains the legacy scalar uncertainty path. Both arms retain matched deterministic channel/pilot frame seeds.

Output is written to `spatial_aware_lmmse_results/` as CSV, MAT, and plots. The results should be interpreted as a comparison of complete receiver configurations, not estimator-only MSE. For isolated estimator quality, run:

```matlab
results = channel_estimator_mse_validation_run;
```
