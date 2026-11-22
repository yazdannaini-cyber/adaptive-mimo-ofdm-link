# V10.1 — User Validation Review

## Scope of the supplied run

The supplied `spatial_aware_lmmse_comparison_run` log contains 60 rows:

- 5 spatial-correlation profiles: IID, TX-correlated, RX-correlated, joint-moderate, joint-high.
- 2 channel models: Rayleigh and Rician.
- 2 estimator/receiver configurations: per-link LMMSE with scalar effective noise, and joint spatial-frequency LMMSE with covariance-aware detection.
- 3 SNR points: 16, 24 and 32 dB.
- 500 frames per run.

The user-provided log confirms completion of all 60 rows. It does not include the `run_tests` summary, so this report does not claim a new test-suite pass based solely on this log.

## Main observations

1. **IID profile:** per-link and joint configurations produce matching BER and spectral-efficiency results. This is expected because the joint spatial prior reduces to independent spatial links in the IID case.
2. **Moderate spatial correlation:** direct estimator-only validation from V9 showed measurable NMSE gains for joint LMMSE, particularly at low pilot SNR. In the end-to-end adaptive receiver study, however, BER and spectral-efficiency changes are small.
3. **High spatial correlation:** the largest end-to-end differences in the supplied run occur at low SNR. For example, high-correlation Rayleigh at 16 dB changes BER from 3.12e-4 to 2.95e-4, while spectral efficiency remains 1.72 bit/s/Hz. High-correlation Rician at 16 dB changes BER from 2.34e-4 to 2.29e-4 and spectral efficiency from 3.45 to 3.50 bit/s/Hz.
4. Some points show a very small regression or no measurable change. Therefore the current evidence supports an estimator-MSE improvement under spatial correlation, but does not support a general claim that covariance-aware detection always improves end-to-end BER or throughput.

## Interpretation

The estimator-only experiment and the link-level experiment answer different questions:

- Estimator-only paired Monte Carlo measures whether the channel estimate is closer to the true channel (NMSE/MSE).
- Adaptive link Monte Carlo measures the combined effect of channel estimation, detector, modulation choice, MIMO mode selection and finite-frame errors.

A better channel estimate need not change the selected mode or modulation. It may therefore improve NMSE while leaving BER and spectral efficiency almost unchanged. The receiver's model assumptions and the current adaptation policy can also dominate the link-level result.

## Recommended next step

Validate temporal channel variation before adding more receiver complexity. The current channel generator is block-static within a frame. The next release adds a sample-varying Doppler channel and measures temporal autocorrelation and OFDM ICI. A later release can add pilot tracking and covariance-aware prediction across OFDM symbols.
