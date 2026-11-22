# V11 — Sample-Varying MIMO Doppler Channel

## Purpose

The existing channel generator creates a block-static frequency-selective MIMO channel. That model is useful for frequency-domain equalization, but it cannot represent channel evolution within an OFDM symbol or the resulting inter-carrier interference (ICI).

V11 adds a sample-varying FIR channel:

\[
y_r[n] = \sum_{t=1}^{N_t}\sum_{\ell=0}^{L-1}
 h_{r,t,\ell}[n]x_t[n-\ell] + w_r[n].
\]

The channel taps may change for every time sample. The time-varying channel therefore produces ICI when the channel changes within an OFDM useful symbol.

## Diffuse fading model

Each diffuse link and tap is generated using a finite sum-of-sinusoids approximation:

\[
h(t)=\frac{1}{\sqrt{N_{osc}}}\sum_{m=1}^{N_{osc}}
 e^{j(2\pi f_D\cos\theta_m\,t+\phi_m)},
\]

where the arrival angles \(\theta_m\) and phases \(\phi_m\) are independently drawn. Its target isotropic-scattering autocorrelation is

\[
R_h(\tau)=J_0(2\pi f_D\tau).
\]

The finite oscillator count approximates this reference; it is not an exact continuous Jakes spectrum.

## Spatial correlation

The diffuse MIMO tap matrix is shaped sample-by-sample using the existing separable Kronecker factors:

\[
H_\ell[n]=L_{RX}W_\ell[n]L_{TX}^{H}.
\]

The unit-diagonal correlation matrices preserve each individual link's average tap power while introducing cross-link correlation.

## Rician component

Rician channels combine the diffuse process with a deterministic first-tap LOS component. The diffuse and LOS components are scaled according to the configured K-factor. The LOS Doppler shift is configured as a fraction of the maximum Doppler frequency.

## Sampling and OFDM numerology

The validation example uses 15 kHz subcarrier spacing and

\[
f_s=N_{FFT}\Delta f.
\]

The approximate coherence-time metadata uses \(T_c\approx0.423/f_D\). This is a conventional rule-of-thumb, not a universal exact definition of coherence time.

## New functions

- `src/channel/generate_time_varying_mimo_channel.m`
- `src/channel/apply_time_varying_mimo_ofdm_channel.m`
- `src/metrics/doppler_temporal_metrics.m`
- `examples/doppler_channel_validation_run.m`

## Validation outputs

The Doppler study exports:

- `doppler_channel_validation.csv`
- `doppler_temporal_acf.csv`
- `doppler_channel_validation.mat`
- per-channel plots comparing ICI leakage and estimated/theoretical temporal ACF

## Current limitation

V11 validates sample-varying channel generation and its impact on OFDM ICI. It does not yet implement pilot-based temporal channel tracking or a Doppler-aware receiver. Those are planned follow-up work and must be validated separately.


## V12 relationship

V12 uses the same Doppler channel generator but samples one channel snapshot per OFDM-symbol period to isolate temporal pilot tracking between symbols. It does not replace V11's sample-by-sample FIR/ICI validation; within-symbol ICI remains a separate V11 metric. See `PILOT_TRACKING_MODEL.md`.
