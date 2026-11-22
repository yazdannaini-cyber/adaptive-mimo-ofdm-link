# Version 2 — Changes after baseline CSV review

## What changed

1. **Channel-model-matched LMMSE priors**
   - AWGN: deterministic identity channel, zero diffuse covariance.
   - Rayleigh: zero-mean diffuse channel with the configured PDP.
   - Rician: deterministic LOS mean and diffuse covariance scaled by `1/(K+1)`.

2. **Posterior channel-estimation uncertainty**
   - Added `channel_estimation_error_variance.m` based on the LMMSE posterior covariance.
   - Adaptive detection uses per-subcarrier effective noise `noiseVariance + nTx * posteriorChannelVariance`.
   - `mimo_equalize` accepts either scalar or per-subcarrier noise variance.

3. **BER-prediction-based link adaptation**
   - Added `qam_ber_approx.m` for Gray-coded square QAM.
   - Default policy evaluates predicted BER over the post-detection SINR values and selects the highest modulation meeting `cfg.adaptive.targetBER`.
   - When no modulation meets the target, QPSK remains the minimum-service fallback.
   - The previous percentile/threshold policy remains available with `cfg.adaptive.selectionMethod = 'sinr-threshold'`.

4. **Reporting and validation**
   - Adaptive CSV now includes target BER, selected-mode predicted BER, target-met frame ratio, and mean posterior channel-estimation error variance.
   - Adaptive plots show measured BER, predicted BER, target BER, and spectral efficiency.
   - Test suite expanded to 13 deterministic checks.
   - Original user-supplied CSVs are retained as `results/baseline_*.csv` and summarized in `BASELINE_DATA_ANALYSIS.md`.

## Revalidation procedure

Run from the extracted project root in MATLAB:

```matlab
run_tests
```

All 13 tests should pass. Then run the reduced smoke test:

```matlab
addpath('examples');
quick_run
```

Finally run the full experiment:

```matlab
[fixedResults, adaptiveResults] = run_project;
```

Review the new adaptive columns and confirm that the measured BER is consistent with (or at least reasonably close to) the predicted BER. A predicted BER below target is not a mathematical guarantee in fading channels; measured BER, BLER, confidence intervals, and target-met frame ratio remain the acceptance criteria.

## Validation status

The source files, function references, ZIP archive, and key estimation/BER equations were statically reviewed. Independent numerical checks confirmed the noiseless Rician LMMSE projection, posterior-variance trends, and QAM BER approximation. MATLAB itself was not available in the build environment, so the 13 MATLAB tests and new simulations remain pending execution on a MATLAB installation.
