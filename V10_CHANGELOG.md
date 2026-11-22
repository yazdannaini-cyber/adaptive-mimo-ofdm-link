# V10 Changelog — Covariance-Aware MIMO Detection

## Motivation

V9 directly measured the benefit of joint spatial-frequency LMMSE. It showed up to 34.3% Rayleigh and 42.1% Rician empirical NMSE reduction for rhoTx=rhoRx=0.9 at 0 dB pilot SNR. However, end-to-end BER and spectral efficiency changed very little. The receiver previously compressed posterior uncertainty into a scalar per-tone variance, discarding off-diagonal receive-space covariance.

## Changes

1. `channel_estimation_joint_lmmse_setup` now also returns the full posterior covariance among spatial links for every subcarrier.
2. Added `channel_error_effective_noise_covariance`, which forms

   `R_eff[k] = sigma_n^2 I + sum_t Cov(e_t[k])`

   for independent, unit-power transmit streams.
3. `mimo_equalize` now accepts scalar, per-tone scalar, or per-tone receive covariance noise. The covariance form implements generalized linear MMSE weighting and computes output noise power using the full covariance.
4. `alamouti_decode` now accepts the same covariance form and whitens the receive vector and channel before orthogonal combining.
5. `run_adaptive_mimo_mode` uses the full covariance-aware effective-noise model when `cfg.estimation.spatialAwareLMMSE=true`. The legacy per-link path and scalar/vector noise interface remain supported.
6. Added six tests, bringing the suite from 48 to 54 tests.

## New test coverage

- Posterior link covariance Hermitian/PSD and diagonal consistency.
- IID reduction to the legacy isotropic effective-noise model.
- Generalized MMSE equalization versus explicit noise whitening.
- Covariance-aware Alamouti combining versus explicit noise whitening.
- Rejection of invalid noise covariance matrices.
- Scalar/vector detector-noise backward compatibility.

## Important scope note

The V10 per-link-versus-joint end-to-end experiment compares two complete receiver configurations: per-link LMMSE with the legacy scalar uncertainty approximation, against joint LMMSE with covariance-aware detection. The V9 paired estimator-only NMSE experiment remains the clean estimator-only comparison. Do not describe the V10 end-to-end comparison as isolating estimator quality alone.
