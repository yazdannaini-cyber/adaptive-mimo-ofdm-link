# Spatially Correlated MIMO Channel Model (V7)

## Model

V7 optionally correlates the diffuse multipath component using the separable Kronecker model:

\[
H_\ell^{(c)} = R_{RX}^{1/2} H_\ell^{(i.i.d.)} (R_{TX}^{1/2})^H,
\]

where each entry of the independent matrix is circular complex Gaussian and

\[
[R(\rho)]_{i,j}=\rho^{|i-j|}, \qquad 0\le \rho\le 0.99.
\]

The correlation matrices have unit diagonal, so spatial correlation changes
cross-antenna dependence while preserving each link's average tap power and the
configured power-delay profile.

## Rician interpretation

For Rician channels, the deterministic LOS component is preserved and only the
diffuse/scattered component is spatially correlated. The configured K-factor
therefore retains its meaning as the LOS-to-diffuse power ratio under this model.

## Important estimator limitation

The existing V6 channel estimator performs frequency-domain LMMSE smoothing
independently for each TX/RX link. It does **not** use the cross-antenna covariance
matrices. Consequently, V7 is a robustness/sensitivity study of the existing
adaptive controller under correlated channels, not a correlation-aware optimal
LMMSE receiver. A joint spatial-frequency LMMSE estimator is a later enhancement.

## Validation sweep

`examples/spatial_correlation_validation_run.m` evaluates:

| Profile | TX rho | RX rho |
|---|---:|---:|
| iid | 0.0 | 0.0 |
| TX-only | 0.7 | 0.0 |
| RX-only | 0.0 | 0.7 |
| Joint moderate | 0.7 | 0.7 |
| Joint high | 0.9 | 0.9 |

It tests Rayleigh and Rician channels at 16, 24 and 32 dB, with 500 frames per
channel/SNR/profile and 1,000 frame-cluster bootstrap resamples. These are
simulation settings, not results; run the MATLAB script to generate empirical
CSV, MAT and figures.
