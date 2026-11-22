# V12 Validation Status

## Inherited, user-reported validation

- V11.1: all 69 tests passed in the user's MATLAB environment.
- V11.1: Doppler study completed for Rayleigh and Rician at 0, 30, 120 and 300 Hz.
- V10.1: the user previously reported all 69 tests passing; V11.1 extends and revalidates that inherited suite.

## New V12 work

- Added tests 70–78 for temporal tracking primitives and validation behavior.
- Expected total test count: 78.
- Expected output rows: 864 (2 channels × 3 spatial profiles × 4 Doppler values × 4 pilot periods × 3 pilot SNR values × 3 estimators).
- A separate Python/SciPy numerical check supports the Jakes-LMMSE covariance equations; it is not a MATLAB test.

MATLAB execution of V12 has not yet been performed. Do not describe V12's new tests or pilot-tracking simulation results as passed until run in MATLAB.
