# V8 Independent Numerical Check (Not a MATLAB Test Run)

The joint estimator's eigen-domain transform was checked independently in a
NumPy reference calculation against a dense LMMSE covariance solution for a
complex Hermitian separable covariance model.

- Maximum absolute difference between the efficient Kronecker/eigen solution
  and the dense `C*(C + sigma2*I)^(-1)` reference: approximately `7.6e-14`.
- Maximum absolute difference between the factored posterior diagonal and the
  diagonal of the dense posterior covariance: approximately `1.1e-14`.
- A separate 250-realization Rayleigh reference simulation with rhoTx=rhoRx=0.8
  and LS-noise variance 0.4 yielded normalized MSE about `0.0539` for per-link
  LMMSE and `0.0434` for joint spatial-frequency LMMSE (ratio about `0.805`).

These are independent numerical checks of the mathematical implementation, not
MATLAB execution results and not the V8 paired Monte-Carlo study. The actual
MATLAB tests and experiment must still be run before making project performance
claims.
