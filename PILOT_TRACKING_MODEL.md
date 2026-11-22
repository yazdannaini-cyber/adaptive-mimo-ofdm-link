# V12 — Pilot-Based Temporal Channel Tracking under Doppler

## Objective

V11 established that Doppler makes the MIMO channel evolve over time and creates within-symbol ICI. V12 adds a separate temporal-estimation experiment: given noisy pilot-based channel observations at selected OFDM-symbol indices, estimate the channel on intervening data symbols.

## Symbol-rate channel snapshots

The time-varying tap generator is sampled once per OFDM-symbol period

\[
T_{sym}=\frac{N_{FFT}+N_{CP}}{N_{FFT}\Delta f}.
\]

The tap snapshot at each symbol index is transformed to a frequency response by an FFT over the tap dimension. The temporal process therefore follows the configured Doppler model between OFDM-symbol snapshots. The channel is treated as block-static within each snapshot for this tracking experiment; the within-symbol ICI effect remains separately measured by V11.

## Pilot observation abstraction

For pilot symbol indices \(\mathcal P\), the experiment forms

\[
\widehat{H}_{p}[k,r,t]=H[k,r,t,p]+v_p,
\]

where \(v_p\) is circular complex Gaussian observation noise. The observation is interpreted as a full-grid, orthogonal per-link LS pilot estimate. This deliberately isolates temporal tracking from pilot waveform/resource mapping. It does not claim to model the exact resource-element overhead of a practical pilot pattern.

The Rician deterministic LOS trajectory is supplied as a known prior mean to all three trackers. This is a model-aided benchmark and is optimistic relative to a receiver that must also learn the LOS steering vector and LOS Doppler.

## Compared estimators

1. **Hold-last (causal):** use the most recent pilot's diffuse residual for subsequent data symbols.
2. **Linear interpolation/extrapolation (offline reference):** linearly interpolate between bracketing pilots and extrapolate from the last two pilots at the end of a frame. This method uses future pilots between update times and is therefore a noncausal reference, not a real-time receiver.
3. **Jakes temporal LMMSE (causal):** use only pilot observations at or before the target symbol. The diffuse temporal covariance is

\[
R_h(i,j)=\sigma_h^2J_0(2\pi f_D|i-j|T_{sym}).
\]

For a target time \(t\) and its available past pilot times \(\mathcal P_t\), the estimate of the diffuse component is

\[
\widehat{h}_d(t)=R_{t,p}(R_{p,p}+\sigma_v^2I)^{-1}
                 (y_p-\mu_p),
\]

and the known prior mean \(\mu(t)\) is restored. The posterior variance is

\[
\sigma^2_{post}(t)=\sigma_h^2-R_{t,p}(R_{p,p}+\sigma_v^2I)^{-1}R_{p,t}.
\]

The implementation precomputes prediction weights for each temporal configuration and reuses them across Monte-Carlo frames.

## Metrics

NMSE and MSE are evaluated on data-symbol indices only; pilot symbols are excluded from the score. Frame-level values are resampled as clusters to obtain percentile-bootstrap confidence intervals. For the Jakes-LMMSE predictor, empirical data-symbol MSE is also compared with the mean model-predicted posterior variance.

## Default experiment

- Channels: Rayleigh and Rician.
- Doppler: 0, 30, 120 and 300 Hz.
- Spatial profiles: IID, joint-moderate (TX/RX rho 0.7), and joint-high (TX/RX rho 0.9).
- Pilot update period: 2, 4, 8 and 16 OFDM symbols.
- Pilot observation SNR: 0, 10 and 20 dB.
- 64 OFDM-symbol snapshots per frame; 40 independent frames per condition.
- 1,000 frame-bootstrap resamples at 95% confidence.

The expected channel power used to set pilot observation-noise variance is one for Rayleigh. For Rician it follows the configured diffuse variance plus the configured first-tap LOS power.

## Limitations and interpretation

- The pilot measurements are an abstract full-grid per-link LS observation, not a waveform-level pilot resource mapper.
- The temporal tracker works on per-OFDM-symbol channel snapshots and does not cancel within-symbol ICI. V11 measures ICI separately.
- The Jakes-LMMSE predictor assumes the maximum Doppler frequency and diffuse covariance are known.
- The Rician prior assumes the LOS trajectory is known exactly.
- Linear interpolation is noncausal and should be treated as an offline benchmark.
- No BER or packet-goodput claim is made by this estimator-only experiment. Its principal endpoints are data-symbol NMSE/MSE and posterior-variance calibration.
