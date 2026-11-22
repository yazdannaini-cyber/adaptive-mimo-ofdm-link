# V9 User Validation — Direct Paired Channel-Estimation Results

## Execution status

The user executed the V9 MATLAB project and reported:

- All 48 tests passed.
- The direct estimator-only validation completed all 50 rows.
- Each row used 500 paired channel/pilot-noise realizations and 1,000 frame-bootstrap resamples.
- The two estimators saw the same channel and the same pilot-noise realization in each paired frame.

## Main conclusion

The joint spatial-frequency LMMSE estimator is equivalent to per-link LMMSE in the IID spatial case, as expected. Its benefit grows with spatial correlation and is largest at low pilot SNR. At high pilot SNR the gains become small because both estimators already have low estimation error.

The experiment measures channel-estimator MSE directly. It does not establish that BER or spectral efficiency improves by the same percentage; end-to-end link metrics also depend on detection, mode selection, modulation, and the receiver's use of posterior uncertainty.

## Representative measured improvements

The reported improvement is the percentage reduction in empirical NMSE of joint LMMSE relative to per-link LMMSE.

| Spatial profile | Channel | Pilot SNR | Per-link NMSE | Joint NMSE | Improvement |
|---|---|---:|---:|---:|---:|
| IID (0, 0) | Rayleigh | 0 dB | 0.07226 | 0.07226 | 0.0% |
| TX 0.7 | Rayleigh | 0 dB | 0.07465 | 0.06913 | 7.4% |
| RX 0.7 | Rayleigh | 0 dB | 0.07486 | 0.06934 | 7.4% |
| Joint 0.7/0.7 | Rayleigh | 0 dB | 0.07878 | 0.06768 | 14.1% |
| Joint 0.9/0.9 | Rayleigh | 0 dB | 0.09028 | 0.05930 | 34.3% |
| Joint 0.7/0.7 | Rician | 0 dB | 0.04941 | 0.03860 | 21.9% |
| Joint 0.9/0.9 | Rician | 0 dB | 0.04899 | 0.02835 | 42.1% |
| Joint 0.9/0.9 | Rayleigh | 8 dB | 0.01550 | 0.01257 | 18.9% |
| Joint 0.9/0.9 | Rician | 8 dB | 0.01135 | 0.007708 | 32.1% |
| Joint 0.9/0.9 | Rayleigh | 16 dB | 0.002515 | 0.002336 | 7.1% |
| Joint 0.9/0.9 | Rician | 16 dB | 0.002073 | 0.001702 | 17.9% |
| Joint 0.9/0.9 | Rayleigh | 32 dB | 6.513e-5 | 6.489e-5 | 0.4% |
| Joint 0.9/0.9 | Rician | 32 dB | 5.197e-5 | 5.118e-5 | 1.5% |

## Paired uncertainty interpretation

The paired difference is defined in the implementation as

`frameMsePerLink - frameMseJoint`.

Therefore a strictly positive 95% paired-bootstrap interval supports lower MSE for joint LMMSE. The strongest evidence occurs in the highly correlated profile. For example:

- Rayleigh, rhoTx=rhoRx=0.9, pilot SNR 0 dB: paired MSE-difference CI [0.0230, 0.0249].
- Rician, rhoTx=rhoRx=0.9, pilot SNR 0 dB: paired MSE-difference CI [0.0192, 0.0209].

At several high-pilot-SNR points the paired interval includes zero, consistent with a negligible estimator advantage there. In the IID profile, paired differences are numerical zero, which is the intended equivalence baseline.

## Engineering interpretation

The V9 end-to-end comparison did not show material BER/spectral-efficiency changes despite these direct NMSE gains. Source inspection showed why the detector could underuse the joint estimator: the receiver reduced posterior estimation uncertainty to a scalar, per-tone conservative variance. This discarded off-diagonal receive-space error covariance.

V10 therefore propagates the full posterior spatial covariance into a colored effective-noise covariance for both linear-MMSE spatial detection and Alamouti combining. This is a receiver-model improvement, not a claim that the V9 estimator-only results themselves changed.

## Limitations

- Results are specific to the project's 2x2 Kronecker-correlated channel model and configured PDP.
- Correlation is assumed known to the estimator; covariance mismatch is not covered by this experiment.
- The study uses block-static channel frames and does not yet include Doppler/time variation.
- V9 results validate estimation MSE, not coded packet reliability or hardware performance.
