# MIMO-OFDM Simulation Report

> Fill this report only after running the MATLAB project. Do not estimate or invent measured values.

## 1. Objective

Describe the communication-system design problem and the algorithms under comparison.

## 2. Configuration

| Parameter | Value |
|---|---:|
| FFT size | 64 |
| Cyclic prefix | 16 samples |
| Transmit / receive antennas | 2 / 2 |
| Channel taps | 5 |
| Data OFDM symbols per frame | 8 |
| Monte Carlo frames per point | 80 (default) |
| Modulation orders | QPSK, 16-QAM, 64-QAM |
| Channel models | AWGN, Rayleigh, Rician |

Record the MATLAB release, operating system, actual configuration, and run date.

## 3. Channel estimation results

Insert `NMSE_*.png`. Discuss LS versus model-aware LMMSE. Explain the prior mean and diffuse covariance used for AWGN, Rayleigh, and Rician, and the effect of pilot noise.

## 4. Detection results

Insert `BER_Fixed_*.png`. Compare LS+MMSE, LMMSE+MMSE, perfect CSI+MMSE, and LMMSE+ZF. Discuss noise enhancement and conditioning.

## 5. Adaptive modulation results

Insert `Adaptive_Performance_*.png`. Report measured BER/BLER and its 95% interval, selected-mode ratios, payload-bit-weighted predicted BER, predicted-target-feasible frame ratio, QPSK fallback frame ratio, measured-BER target status, posterior channel-estimation error variance, spectral efficiency, and goodput. Do not describe predicted-feasible frame ratio as measured BER compliance. State clearly that the adaptive x-axis is Es/N0 and compare payload-weighted predicted BER against measured Monte-Carlo BER.

## 6. Statistical reliability

Report total transmitted bits and 95% Wilson confidence intervals. Distinguish zero observed errors from zero true BER. Increase Monte-Carlo frames where the confidence interval is too wide.

## 7. Interpretation

Explain which receiver performed best under each channel model and why. Distinguish measured observations from theoretical expectations.

## 8. Limitations and future work

Discuss synchronization assumptions, idealized Rician geometry, lack of coding, Doppler, antenna correlation, and possible extensions.


## 9. MIMO transmission-mode comparison (V5)

Insert the BER and spectral-efficiency figures from `comparison_results/`. Compare 2x2 spatial multiplexing with linear MMSE detection against rate-1 Alamouti STBC with coherent combining. State the QAM order, channel model, SNR, frame count, confidence intervals, FER, and spectral efficiency for each claim. The two modes have different raw spectral efficiencies; interpret reliability together with rate rather than comparing BER alone.

## 10. Frame reliability versus bit reliability

The current uncoded simulator defines frame error rate as the fraction of data frames with at least one payload bit error. The legacy adaptive result table calls this `BLER`, but it is not post-FEC block error rate. Current goodput is bit-level goodput and counts correctly received bits in partially erroneous frames. Do not call it CRC-accepted packet throughput.

## Joint Spatial-Frequency LMMSE (V8)

Document the configured `R_TX`, `R_RX`, PDP covariance, Rician LOS prior,
posterior uncertainty model, and paired comparison against the original per-link
LMMSE estimator. Include the CSV-generated BER/CI and spectral-efficiency plots.
Do not claim an improvement until the paired MATLAB study has been run and
validated. Report whether the BER confidence-interval upper bounds remain below
the configured target for each profile.

## V10 covariance-aware receiver experiment

Report the exact MATLAB release, frame count, random-seed schedule, channel model, TX/RX correlation coefficients, SNR, and BER target. Compare the per-link LMMSE + scalar-uncertainty receiver against the joint spatial-frequency LMMSE + covariance-aware receiver.

Include measured BER and FER with Wilson and frame-cluster bootstrap intervals, spectral efficiency, bit goodput, and selected spatial-multiplexing/Alamouti ratios. Keep the estimator-only NMSE results from `channel_estimator_mse_validation_run` in a separate subsection. Explicitly state that the V10 end-to-end comparison changes both the channel estimator and how posterior uncertainty is propagated into detection.

## V14 — Full ICI Equalization

- State clearly that `perfect-csi` is a diagonal-channel oracle.
- State clearly that `perfect-csi-ici` uses the exact time-varying tap sequence
  to construct the full OFDM frequency-coupling matrix and is an oracle
  diagnostic, not a pilot-estimated implementation.
- Compare BER for the diagonal perfect-CSI and full-ICI oracle by channel,
  spatial profile, Doppler and pilot spacing.
- Report the computational cost of the full block solve relative to one-tap
  per-subcarrier MMSE.
- Verify that the static-channel full matrix collapses to independent
  per-subcarrier MIMO channels.
- Do not claim a deployable ICI-aware receiver until the full coupling matrix
  can be estimated from a practical pilot design.
