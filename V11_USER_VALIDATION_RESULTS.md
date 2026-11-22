# V11.1 User-Reported Validation Results

The user executed the V11.1 project in MATLAB and reported:

- All 69 tests passed.
- Doppler validation completed for Rayleigh and Rician at 0, 30, 120 and 300 Hz.
- The Rician first-tap ensemble K-factor estimate was 3.948 versus configured 3.981.

## Doppler measurements

| Channel | Doppler (Hz) | ACF MAE | ICI leakage ratio | Approx. coherence time (s) |
|---|---:|---:|---:|---:|
| Rayleigh | 0 | 0 | 0 | Inf |
| Rayleigh | 30 | 0.03126 | 6.181e-06 | 0.0141 |
| Rayleigh | 120 | 0.1239 | 1.412e-04 | 0.003525 |
| Rayleigh | 300 | 0.06436 | 5.514e-04 | 0.00141 |
| Rician | 0 | 0 | 0 | Inf |
| Rician | 30 | 0.06552 | 9.541e-06 | 0.0141 |
| Rician | 120 | 0.03870 | 1.742e-04 | 0.003525 |
| Rician | 300 | 0.08862 | 1.166e-03 | 0.00141 |

ICI leakage rises with Doppler in both channel models. ACF MAE is not monotonic because it is a finite-sample/finite-oscillator estimate of the Jakes reference; it should not be interpreted as a monotonic measure of Doppler severity. The coherence-time values use the stated rule-of-thumb \(T_c≈0.423/f_D\), not an exact universal boundary.

At 15 kHz subcarrier spacing and 64-point FFT with 16-sample CP, the OFDM-symbol interval is approximately 83.33 microseconds. The approximate coherence time therefore corresponds to about 169, 42.3 and 16.9 OFDM-symbol intervals at 30, 120 and 300 Hz respectively.
