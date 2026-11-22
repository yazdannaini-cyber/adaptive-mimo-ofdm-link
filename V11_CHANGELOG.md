# V11 Changelog

## Added

- Sample-varying MIMO FIR channel with per-sample tap coefficients.
- Finite sum-of-sinusoids diffuse fading approximation targeting the Jakes `J0` temporal autocorrelation.
- Rician first-tap LOS component with configurable LOS Doppler fraction.
- Compatible TX/RX Kronecker spatial correlation for time-varying diffuse taps.
- Time-varying channel application that preserves OFDM-symbol chronology and allows intra-symbol channel variation/ICI.
- Temporal ACF estimator and Doppler validation runner.
- ICI leakage study for 0, 30, 120 and 300 Hz maximum Doppler.
- Doppler configuration fields and validation.
- 14 new deterministic/unit/integration tests (test IDs 56–69).

## Preserved

- All V10/V10.1 channel-estimation, covariance-aware detection, adaptive MIMO, metrics and existing tests.
- Existing block-static channel functions and simulation runners remain available and unchanged for backward compatibility.

## Validation status

The package has been statically inspected and an independent Python numerical sanity check has been run for the sum-of-sinusoids ACF and Doppler-induced ICI trend. MATLAB/Octave is not installed in the build environment, so the 69 MATLAB tests and the Doppler validation runner still require execution in the user's MATLAB environment.


## V11.1 hotfix

- Hardened Test 61 (Rician first-tap K-factor) to estimate the ensemble power ratio across eight independent seeded realizations.
- Reduced this test to a one-tap PDP so the extra independent realizations do not increase the test's computational load materially.
- Added an exact deterministic LOS-power check and retained a 20% tolerance for the finite-realization diffuse-power estimate.
- No channel-generation equations, Doppler process, or receiver algorithms were changed.
