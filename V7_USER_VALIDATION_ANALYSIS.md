# V7.1 User Validation Analysis

Source: MATLAB console output supplied by the user on 2026-11-21, mapped using
the exact loop order in `examples/spatial_correlation_validation_run.m`.

The script iterates profiles in this order, then Rayleigh followed by Rician,
with SNR values 16, 24 and 32 dB:

1. `iid`: rhoTx=0.0, rhoRx=0.0
2. `tx-correlated`: rhoTx=0.7, rhoRx=0.0
3. `rx-correlated`: rhoTx=0.0, rhoRx=0.7
4. `joint-moderate`: rhoTx=0.7, rhoRx=0.7
5. `joint-high`: rhoTx=0.9, rhoRx=0.9

The supplied console log does not print the profile label before each six-row
block. Profile attribution below is derived from the source-controlled loop
order. The CSV written by the MATLAB script includes `CorrelationProfile`,
`TxCorrelationCoefficient`, and `RxCorrelationCoefficient` columns and is the
preferred machine-readable source for final reporting.

## Measured spectral efficiency (bit/s/Hz)

| Profile | Rayleigh 16 dB | Rayleigh 24 dB | Rayleigh 32 dB | Rician 16 dB | Rician 24 dB | Rician 32 dB |
|---|---:|---:|---:|---:|---:|---:|
| iid | 2.48 | 4.01 | 5.33 | 3.25 | 6.67 | 7.64 |
| TX-only 0.7 | 2.32 | 3.87 | 4.75 | 3.29 | 6.68 | 7.58 |
| RX-only 0.7 | 2.37 | 3.84 | 4.77 | 3.26 | 6.77 | 7.59 |
| Joint 0.7/0.7 | 2.20 | 3.75 | 4.28 | 3.37 | 6.78 | 7.54 |
| Joint 0.9/0.9 | 1.71 | 3.37 | 3.85 | 3.59 | 6.65 | 7.52 |

## BER and observed confidence bounds

All 30 reported rows had measured BER below the target `1e-3`; all Wilson and
frame-bootstrap upper bounds were also below `1e-3`. The tightest margins were
in Rayleigh at 16 dB under high correlation: BER=0.000308, Wilson upper=0.00035
and frame-bootstrap upper=0.00037.

## Interpretation

- In Rayleigh, stronger joint correlation is associated with lower spectral
efficiency, particularly at 32 dB: 5.33 bit/s/Hz for iid versus 3.85 bit/s/Hz
for rhoTx=rhoRx=0.9. The mode mix shifts toward Alamouti to retain reliability.
- In Rician, spectral efficiency is relatively stable at 32 dB (7.52-7.64
bit/s/Hz), because the configured deterministic LOS component is preserved and
only the diffuse component is spatially correlated.
- Lower BER at high correlation must not be interpreted as correlation improving
the physical channel. The adaptive controller changes its architecture and QAM
order; the associated spectral-efficiency reduction is the reliability cost.
- The V7 estimator remains per-link and ignores cross-antenna covariance. V8
introduces a joint spatial-frequency LMMSE estimator to address this model
mismatch. No V8 performance improvement is claimed until its paired MATLAB
experiment has completed.

The V7.1 sweep used 500 frames per channel/SNR/profile and 1,000 frame-cluster
bootstrap resamples. Keep the generated CSV as the reproducible numeric record;
summary values above are rounded console values.
