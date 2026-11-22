# V9 Changelog

- Added `posteriorMeanErrorVarianceByTone` to the joint LMMSE model while
  preserving the existing conservative worst-link variance output.
- Added frame-cluster percentile bootstrap confidence intervals for means.
- Added a paired estimator-only validation function that reuses the exact same
  channel and noisy orthogonal-pilot observation for both estimators.
- Added direct frame-level NMSE, empirical coefficient MSE, paired MSE-difference
  confidence intervals, and posterior-MSE calibration ratios.
- Added a 50-row study over five TX/RX correlation profiles, Rayleigh/Rician
  channels, and five pilot-SNR points.
- Added plots for per-link versus joint LMMSE NMSE and relative improvement.
- Added six tests (43-48): deterministic bootstrap, posterior mean/bound
  outputs, direct-estimator smoke test, paired reproducibility, input validation,
  and IID estimator equivalence.
- The V8 user run passed tests 1-42 and completed the 60-row adaptive paired
  experiment. The printed end-to-end BER results alone do not establish an
  empirical NMSE advantage; V9 measures that directly.
- MATLAB/Octave is unavailable in the build environment. V9's 48 tests and
  estimator-MSE study still require execution in MATLAB.
