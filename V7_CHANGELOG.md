# V7 Changelog

- Added optional Kronecker spatial correlation for diffuse Rayleigh/Rician fading.
- Supports separate TX and RX exponential Toeplitz correlation coefficients.
- Keeps the feature disabled by default, preserving prior V6 independent-fading runs.
- Preserves per-link average PDP because the correlation matrices have unit diagonal.
- Preserves the deterministic Rician LOS component and correlates only diffuse scattering.
- Added a 5-profile, 2-channel, 3-SNR adaptive robustness validation script.
- Added BER and spectral-efficiency plots versus spatial-correlation profile.
- Added six tests (30-35) for backward compatibility, TX/RX correlation calibration,
  per-link power preservation, Rician LOS preservation and invalid-input rejection.
- Total intended test count: 35.

## Validation status

The package was created from the MATLAB V6 source and inspected statically. MATLAB
and GNU Octave are not installed in this build environment; tests 30-35 and the V7
Monte Carlo sweep have **not** been executed here. Run `run_tests` in MATLAB first.
