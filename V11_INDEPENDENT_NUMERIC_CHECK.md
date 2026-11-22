# V11 Independent Numerical Sanity Check

This check was performed independently in Python using the same finite sum-of-sinusoids model and the sample-varying FIR equation. It is a mathematical/algorithmic sanity check, not a MATLAB execution.

## Temporal autocorrelation

For 128 oscillators, `fD = 150 Hz`, `fs = 2 kHz`, 30,000 samples and four independent link processes, the estimated real ACF at lags `[0, 1, 2, 4, 6]` samples was approximately:

`[1.0000, 0.9470, 0.7965, 0.3103, -0.1704]`

The corresponding `J0(2*pi*fD*tau)` values were approximately:

`[1.0000, 0.9452, 0.7900, 0.2906, -0.1962]`

Maximum absolute difference: approximately `0.0258`.

## OFDM ICI trend

A single-subcarrier, 2x2 OFDM probe with 64-point FFT, 16-sample CP and 15 kHz subcarrier spacing was applied to a sample-varying channel. A separate Python check observed ICI leakage ratios close to:

| Maximum Doppler | ICI leakage ratio |
|---:|---:|
| 0 Hz | 0 |
| 30 Hz | 4.5e-6 |
| 120 Hz | 5.0e-5 |
| 300 Hz | 6.98e-4 |

These values depend on the realization and are not target performance results. They show the expected monotonic trend for the chosen probe: higher intra-symbol Doppler variation produces more off-subcarrier energy.

## Limitations

The Python check does not replace the MATLAB test suite. Exact stochastic values will differ because the random-number generators and implementation details differ.
