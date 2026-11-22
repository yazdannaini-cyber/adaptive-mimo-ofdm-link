# Baseline Simulation Data Review (Original Implementation)

## Dataset scope

The supplied baseline CSV files contain 87 rows in total:

- Fixed modulation: 63 rows = 3 channel models × 3 QAM orders × 7 Eb/N0 points.
- Adaptive modulation: 24 rows = 3 channel models × 8 Es/N0 points.
- Each operating point used 80 Monte-Carlo frames.

These files were generated with the original implementation, before the model-aware estimator and BER-prediction update. They are retained as baseline evidence and must not be mixed with new-version output files.

## Fixed-modulation findings

Average channel-estimation NMSE across the 21 operating points for each channel:

| Channel | Mean LS NMSE | Mean LMMSE NMSE | LS/LMMSE ratio |
|---|---:|---:|---:|
| AWGN | 0.145514 | 0.010535 | 13.81× |
| Rayleigh | 0.071758 | 0.005401 | 13.29× |
| Rician | 0.073225 | 0.005413 | 13.53× |

Example at Rayleigh, 4-QAM, Eb/N0 = 24 dB:

| Receiver | BER |
|---|---:|
| LS + MMSE | 1.904e-3 |
| LMMSE + MMSE | 5.980e-4 |
| Perfect CSI + MMSE | 5.550e-4 |
| LMMSE + ZF | 9.340e-4 |

This example shows the LMMSE estimate approaching the perfect-CSI reference and the ZF detector exhibiting a higher BER than MMSE in this operating point.

## Adaptive-modulation findings

At the highest simulated Es/N0 point for each channel:

| Channel | Es/N0 | BER | BLER | Spectral efficiency (bit/s/Hz) | Goodput (bit/s/Hz) | QPSK / 16-QAM / 64-QAM frame mix |
|---|---:|---:|---:|---:|---:|---:|
| AWGN | 28 dB | 0 | 0 | 7.680 | 7.680 | 0% / 0% / 100% |
| Rayleigh | 28 dB | 2.000e-3 | 0.8375 | 4.352 | 4.343 | 37.5% / 55% / 7.5% |
| Rician | 28 dB | 1.670e-4 | 0.1625 | 7.488 | 7.487 | 0% / 7.5% / 92.5% |

The original fixed SINR-threshold policy did not meet the configured target BER of 1e-3 in Rayleigh at 20, 24, or 28 dB. The original adaptation used a 10th-percentile post-detection SINR and thresholds based on approximate uncoded-QAM values. Those values were optimistic in the frequency-selective Rayleigh case.

## Interpretation and limitations

1. The initial LMMSE implementation assumed zero-mean fading with the configured PDP for every channel. That prior is mismatched to the deterministic identity AWGN reference and to the nonzero LOS mean in the Rician channel.
2. The initial adaptive selector used estimated post-detection SINR but did not include posterior channel-estimation uncertainty in its effective noise model.
3. A finite-simulation BER of zero is not a true zero-error guarantee. The 95% Wilson confidence bounds in the CSV must be reviewed; more frames are needed for stronger high-SNR claims.
4. A BLER near one can coexist with a modest BER when each frame contains many bits. For example, even a small bit-error probability can make it likely that a long frame contains at least one bit error.
5. The adaptive target is an empirical design objective, not an unconditional guarantee. The updated predicted-BER selector must be validated against new measured BER and BLER results.

## Changes in this revision

- Added channel-model-matched LMMSE priors: deterministic identity for AWGN, zero-mean diffuse prior for Rayleigh, and nonzero LOS mean plus scaled diffuse covariance for Rician.
- Added posterior channel-estimation error variance.
- Extended the MIMO equalizer to accept per-subcarrier effective noise variance.
- Added an analytical Gray-coded square-QAM BER approximation and made predicted-BER link adaptation the default policy.
- Retained the previous SINR-threshold policy as a configurable baseline option.
- Expanded the deterministic test suite from 8 to 12 tests.

The updated MATLAB code must be run on the user's MATLAB installation before claiming that the new 12-test suite or new simulations passed.
