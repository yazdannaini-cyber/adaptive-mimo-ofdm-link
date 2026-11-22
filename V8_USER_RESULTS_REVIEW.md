# V8 User-Run Results Review

## What was validated

The supplied MATLAB console output reports that tests 1-42 all passed and that the
paired spatial-aware estimator study completed 60 rows (5 spatial profiles x 2
channel models x 3 SNR points x 2 estimators), with 500 frames per run.

## Interpretation of the adaptive-link results

The per-link and joint-estimator runs have identical or very similar end-to-end
BER, spectral efficiency, and mode-selection results in most printed conditions.
This does **not** prove that the two channel estimators have the same estimation
accuracy. The adaptive simulation's final BER is also affected by the selected
MIMO architecture and QAM order. The console output does not print the realized
channel-estimation NMSE for each estimator.

The V8 comparison table contains a posterior channel-error-variance field, which
is a model-based uncertainty prediction, not a direct measurement of
`H_est - H_true`. A direct estimator-only validation is therefore needed before
claiming that joint spatial-frequency LMMSE improves empirical channel estimation
in the tested settings.

## V9 response

V9 adds a paired estimator-only experiment. For every frame, both estimators use
the same generated channel and the same noisy orthogonal-pilot observation. The
experiment directly measures:

- Per-frame and mean channel-estimation NMSE.
- Empirical MSE per channel coefficient.
- Percent NMSE improvement of joint LMMSE over per-link LMMSE.
- Frame-bootstrap confidence intervals for each estimator's mean NMSE.
- A paired bootstrap interval for the per-link-minus-joint empirical MSE.
- Model-predicted posterior MSE versus observed empirical MSE (calibration ratio).

The study sweeps Rayleigh/Rician channels, five TX/RX spatial-correlation
profiles, and pilot SNR values from 0 to 32 dB. It uses 500 frames and 1,000
bootstrap resamples per condition, for 50 result rows.

No V9 MATLAB execution result is available yet. The 48-test suite and the new
experiment must be run in MATLAB before V9 performance claims are made.
