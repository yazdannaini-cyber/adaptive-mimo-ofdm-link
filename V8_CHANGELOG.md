# V8 Changelog

- Added exact joint spatial-frequency LMMSE for the separable Kronecker channel prior.
- Added a precompute/apply split so covariance eigensystems are prepared once per SNR, not once per frame.
- Added posterior per-tone channel-estimation uncertainty using a conservative worst-link variance bound.
- Added opt-in `cfg.estimation.spatialAwareLMMSE`; default is false for V6/V7 backward compatibility.
- Integrated the new estimator and its posterior uncertainty into adaptive MIMO-mode/QAM selection.
- Added a paired comparison run with frame-level reseeding so channel and pilot-noise realizations match between per-link and joint estimators even if prior mode choices differ.
- Added explicit correlation-profile labels to the original V7 correlation-sweep console output.
- Added seven tests (36-42): zero-correlation equivalence, correlated-channel MSE gain, noiseless Rician recovery, posterior uncertainty behavior, covariance validation, end-to-end integration, and paired-run reproducibility.
- Total intended tests: 42.

MATLAB/Octave is not available in the build environment. New tests and V8 Monte Carlo comparisons require execution in MATLAB before any performance claim is treated as validated.
