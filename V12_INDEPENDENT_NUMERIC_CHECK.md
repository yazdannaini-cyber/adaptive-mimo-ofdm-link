# V12 Independent Numeric Check

A separate Python/SciPy calculation was used to check the temporal-LMMSE covariance equations independently of MATLAB. It used a single complex Jakes fading process with 64 symbols, maximum Doppler 120 Hz, symbol period 1/12000 s, pilot period 8 symbols, and complex pilot-observation noise variance 0.1 (10 dB relative to unit diffuse variance). Across 1,000 independent 64-oscillator sum-of-sinusoids realizations, the data-symbol MSE values were:

| Tracker | Mean data-symbol MSE |
|---|---:|
| Hold-last | 0.13840 |
| Linear interpolation (noncausal reference) | 0.08894 |
| Causal Jakes-LMMSE | 0.09868 |

The Jakes-LMMSE predictor outperformed causal hold-last in this synthetic check. The linear method performed better in this particular case because it uses future pilot observations and is an offline smoother/reference, not a causal receiver. This Python check verifies the mathematical direction and implementation structure only; it does not validate the MATLAB functions or the full MIMO simulation.

The analytical posterior variance from the Jakes covariance model was positive and decreased when pilot-noise variance was reduced from 0.1 to 0.01:

- Mean posterior variance at 10 dB pilot SNR: 0.09432
- Mean posterior variance at 20 dB pilot SNR: 0.01956

MATLAB execution remains required for V12 tests 70–78 and the Monte-Carlo runner.
