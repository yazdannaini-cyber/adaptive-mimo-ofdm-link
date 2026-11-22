# V10 Validation Status

- V9 user-reported baseline: 48/48 MATLAB tests passed; all 50 direct estimator-validation rows completed.
- V10 implementation: covariance propagation added to the joint estimator, MIMO MMSE detector, Alamouti combiner, and adaptive link runner.
- Added tests: 49–54; expected total: 54 tests.
- Independent numerical checks performed in Python: full posterior per-tone link covariance agrees with direct full-matrix Gaussian posterior covariance to approximately 1e-15 in a representative complex Hermitian case; generalized colored-noise MMSE agrees with explicit whitening to approximately 5e-15 for tested random matrices.
- MATLAB execution status: **not yet executed by the user**. Do not report the 54-test suite or V10 simulations as passed until the user runs them in MATLAB.
