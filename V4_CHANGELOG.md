# V4 changelog — frame-cluster statistical validation

## Why this version exists

The V3 quick run showed that bit-level Wilson confidence intervals are useful, but bits within the same fading frame share one channel realization and are not fully independent. V4 therefore adds a frame-cluster bootstrap interval in addition to (not instead of) the Wilson interval.

## Changes

- Added `cluster_bootstrap_ber_ci`: resamples complete frames and estimates a percentile interval for aggregate BER without Statistics Toolbox dependencies.
- Adaptive simulation now stores per-frame bit errors and payload bits and reports frame-cluster bootstrap 95% intervals, separate Wilson/cluster compliance flags, and a combined flag requiring both upper bounds below target.
- Fixed-modulation simulation reports frame-cluster intervals for each receiver variant.
- Added bootstrap settings to `default_config.m` and validation.
- Expanded deterministic test suite from 15 to 19 tests.
- Added `examples/adaptive_validation_run.m`: a 1,000-frame-per-point Monte-Carlo study for Rayleigh and Rician channels over 12–32 dB, with 2,000 bootstrap resamples.
- Kept the existing BER point estimate and Wilson interval unchanged; cluster-bootstrap intervals are additional evidence, not a replacement for the existing metric.

## Interpretation cautions

- A frame-cluster bootstrap interval can be degenerate at zero when every simulated frame has zero errors. The bit-IID Wilson interval is retained because it gives a nonzero upper bound in this case.
- The bootstrap measures empirical between-frame variability for the configured channel model and number of frames; it does not establish a universal guarantee.
- Run `run_tests` and then `adaptive_validation_run` in MATLAB before using V4 results in a CV, report, or interview.
