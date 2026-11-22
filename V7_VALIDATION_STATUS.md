# V7 Validation Status

## Inherited and user-confirmed V6 status

The user supplied MATLAB console output confirming that all 29 V6 tests passed and that the 1,000-frame adaptive MIMO simulation completed for 12 Rayleigh/Rician operating points. The complete transcribed measurements are in `V6_USER_VALIDATION_RESULTS.md`.

## V7 additions

- Optional separable Kronecker spatial correlation for Rayleigh/Rician diffuse fading.
- Independent TX/RX correlation coefficients using exponential Toeplitz matrices.
- Default remains disabled to preserve V6 independent-fading behavior.
- New sweep: five correlation profiles x two channel models x three SNR values.
- Six new tests, intended total 35.

## Not yet executed

The V7 test suite and spatial-correlation Monte Carlo sweep have not been run in this environment because MATLAB and GNU Octave are unavailable. Do not describe the new V7 tests as passed until the user runs `run_tests` in MATLAB. The same applies to the new V7 numerical sweep.
