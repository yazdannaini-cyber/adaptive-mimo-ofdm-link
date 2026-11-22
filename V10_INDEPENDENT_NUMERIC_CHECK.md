# V10 Independent Numerical Checks

These checks are implementation cross-checks, not MATLAB execution or Monte Carlo validation.

## 1. Posterior covariance slice

For a representative complex Hermitian positive-semidefinite frequency covariance and link covariance, a Python/NumPy calculation compared:

- the full Gaussian posterior covariance `C - C (C + sigma^2 I)^(-1) C`, with frequency as the fastest vectorized dimension; and
- the V10 separable-eigen expression for each per-subcarrier link-covariance slice.

The maximum absolute entrywise difference across tested subcarriers was approximately `1.01e-15`.

## 2. Generalized MMSE weighting

For 100 representative random 2x2 MIMO channel/noise-covariance cases, the covariance-weighted filter was compared with the equivalent explicitly whitened receiver. Maximum differences were approximately:

- Equalized symbol output: `4.56e-15`.
- Post-detection SINR: `3.11e-15`.

## 3. Alamouti covariance path

The covariance-aware implementation explicitly applies the Hermitian covariance whitening transform before the usual Alamouti equations. The MATLAB regression test (Test 52) compares this path against a separately whitened scalar-noise reference. This test still requires MATLAB execution.

## Scope

These checks validate the key matrix transformations in isolation. They do not validate MATLAB syntax/runtime compatibility, the full link simulation, or whether V10 improves measured BER/FER. The user must run the V10 MATLAB regression and end-to-end comparison.
