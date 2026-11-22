# V9: Direct paired channel-estimator validation

V9 adds an estimator-only experiment so that channel-estimation quality is not
inferred indirectly from final adaptive-link BER. Every frame generates one
MIMO channel and one noisy orthogonal-pilot observation; per-link LMMSE and joint
spatial-frequency LMMSE process exactly that same observation.

Run the regression suite (48 tests):

```matlab
run_tests
```

Run the direct estimator study:

```matlab
addpath('examples');
results = channel_estimator_mse_validation_run;
```

The study uses 5 spatial-correlation profiles, Rayleigh/Rician channels, pilot
SNR = [0 8 16 24 32] dB, 500 frames per condition, and 1,000 bootstrap
resamples. It produces 50 rows with direct per-link/joint NMSE, empirical MSE,
paired MSE-difference confidence intervals, and posterior-MSE calibration
ratios. CSV/MAT and figures are saved in `estimator_mse_results/`.

This source archive has not been executed in MATLAB. Treat the expected test
count and experiment design as planned until the user runs them.
