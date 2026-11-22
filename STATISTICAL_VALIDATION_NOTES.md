# Statistical validation notes

## Why two confidence intervals are reported

The standard Wilson interval in `bit_error_metrics` treats each transmitted bit as a Bernoulli observation. In a block-fading simulation, bits transmitted within a frame share the same channel realization, so treating all bits as independent can understate uncertainty.

V4 additionally resamples whole frames with replacement. For each bootstrap replicate, it calculates total errors divided by total payload bits across the sampled frames. Percentile quantiles provide a frame-cluster interval. This retains within-frame dependence while reflecting frame-to-frame channel variability.

The point estimate remains the pooled BER, not the average of per-frame BERs. The Wilson interval remains in the output for continuity and to avoid a misleading zero-width interval when no errors occur. Inspect both intervals.

## Recommended validation procedure

1. Run `run_tests` and require all 19 tests to pass.
2. Run `addpath('examples'); quick_run` as a short smoke test.
3. Run `addpath('examples'); adaptive_validation_run` for 1,000 independent frames per channel/SNR point.
4. For a target BER of 1e-3, only claim statistically supported compliance when the measured BER is below target and the upper confidence bound is also below target. Inspect the Wilson and frame-cluster intervals; if they disagree materially, increase the number of independent frames.
5. Treat predicted-feasible frame ratio, QPSK fallback ratio, measured BER, and confidence-bound compliance as separate metrics.

## Current model limitations

The channel is block-static during a frame and independently regenerated between frames. Results are conditional on the configured 2x2 MIMO topology, PDP, pilot structure, Rician K-factor, synchronization assumptions, and linear receiver. The bootstrap does not correct model mismatch, imperfect CSI assumptions outside the implemented estimator, or unmodelled RF impairments.
