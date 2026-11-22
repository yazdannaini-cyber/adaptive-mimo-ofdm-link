# V12 Changelog

- Added causal hold-last temporal channel tracking.
- Added linear interpolation/extrapolation as a noncausal offline reference.
- Added a causal Jakes-covariance temporal LMMSE predictor.
- Added precomputed temporal prediction weights and posterior error-variance output.
- Added a pilot-tracking Monte-Carlo experiment across Rayleigh/Rician channels, four Doppler values, three spatial profiles, four pilot periods and three pilot SNR values.
- Added frame-cluster bootstrap confidence intervals for data-symbol NMSE.
- Added empirical-versus-posterior MSE calibration reporting for Jakes-LMMSE.
- Added nine deterministic tests (Tests 70–78).
- Preserved V11 sample-varying FIR/ICI experiment without changing its channel-generation equations.

## Validation status

The V12 files are authored and statically reviewed; an independent Python/SciPy check of the Jakes-LMMSE covariance equations was also completed. MATLAB/Octave is not available in the build environment, so the 78-test suite and pilot-tracking Monte-Carlo runner must be executed in MATLAB before their results are considered validated.
