# V8: Joint Spatial-Frequency LMMSE Channel Estimation

## Motivation

V7 added a Kronecker-correlated diffuse MIMO channel but retained a per-link LMMSE
estimator. That estimator uses the frequency-domain PDP prior independently for
each transmit/receive link and therefore ignores cross-antenna covariance.

V8 adds a joint LMMSE posterior-mean estimator for the separable covariance
model used by the channel generator.

## Covariance model

With `Hls` stored as `[nFFT x nRx x nTx]`, and MATLAB column-major vectorization
of `reshape(Hls,nFFT,nRx*nTx)`, the diffuse-channel covariance is

\[
C_h = (R_{TX}\otimes R_{RX})\otimes R_f,
\]

where `R_f` is induced by the configured PDP and `R_TX`, `R_RX` are the spatial
correlation matrices. The pilot LS observation is modeled as

\[
\hat h_{LS}=h+n_{LS},\qquad
\mathrm{Cov}(n_{LS})=\frac{\sigma_n^2}{P_p}I.
\]

For prior mean `m`, the LMMSE estimate is

\[
\hat h=m+C_h(C_h+\sigma_{LS}^2I)^{-1}(\hat h_{LS}-m).
\]

The implementation does not form or invert the full covariance matrix. It uses
eigendecompositions of the frequency and link covariance factors, exploiting
Kronecker separability. The resulting estimator is equivalent to the full joint
LMMSE solution under the stated model and white orthogonal-pilot LS noise.

## Estimation uncertainty

The posterior diagonal variance is computed in the same eigenbasis. For use by
the existing detector, V8 reports the maximum posterior channel-coefficient
variance over spatial links at each subcarrier. This is a conservative common
per-tone bound when the detector accepts one effective-noise value per tone.

## Rician model

The configured deterministic LOS response remains the prior mean. The joint
spatial-frequency covariance describes only the diffuse component, scaled by
`1/(K+1)`. The diffuse component receives the configured TX/RX correlation;
the LOS mean is not altered.

## Backward compatibility

`cfg.estimation.spatialAwareLMMSE` defaults to `false`; existing V6/V7 behavior is
preserved. The new estimator is opt-in. Set it to `true` when running
`run_adaptive_mimo_mode` to use the joint estimator and its posterior uncertainty.

## Paired validation

Run:

```matlab
addpath('examples');
results = spatial_aware_lmmse_comparison_run;
```

The validation compares the original per-link estimator with joint spatial-
frequency LMMSE using frame-level matched seeds for each channel/profile pair. Channel and pilot-noise realizations are exactly matched; payload/data-noise streams may diverge if the estimators select different modes or QAM orders. The
result CSV explicitly identifies channel, profile, TX/RX rho, estimator, SNR,
BER, confidence intervals, spectral efficiency, mode mix, and feasibility.

The paired study contains 5 spatial profiles, 2 channel models, 3 SNR points,
2 estimators and 500 frames per run. It produces 60 rows. The frame-seed schedule
restarts each frame from a deterministic seed so an earlier different QAM/mode
choice cannot desynchronize later channel and pilot realizations. Results must be run in
MATLAB; this source package does not claim pre-executed MATLAB outcomes.

V9 additionally exposes `posteriorMeanErrorVarianceByTone`, the mean posterior
variance over spatial links at each subcarrier. This complements the existing
`posteriorErrorVarianceByTone`, which remains the worst-link value used as a
conservative common bound by the detector. V9 uses the mean field only for
posterior-MSE calibration; the conservative detector bound remains unchanged.

## V10: full posterior spatial covariance at the detector

V10 extends the joint model with `posteriorLinkCovarianceByTone`, a full
`[nRx*nTx x nRx*nTx x nFFT]` posterior covariance array. This retains the
receive-space covariance structure that is not visible in a scalar worst-link
variance.

For independent, zero-mean, unit-power streams, channel-estimation error
contributes the following effective receive-noise covariance on each tone:

\[
R_{\mathrm{eff}}[k] = \sigma_n^2 I_{N_r} +
\sum_{t=1}^{N_t} C_{e_t}[k].
\]

Here `C_e_t[k]` is the `nRx` by `nRx` diagonal transmit-link block of the
posterior covariance associated with transmit antenna `t`. Cross-TX error
covariance does not enter this second moment when independent data streams have
zero cross-correlation. The helper
`channel_error_effective_noise_covariance` implements this construction.

The spatial-MMSE equalizer accepts the resulting per-tone covariance and uses

\[
W[k] = \left(H^H[k]R_{\mathrm{eff}}^{-1}[k]H[k]+I\right)^{-1}
H^H[k]R_{\mathrm{eff}}^{-1}[k].
\]

The Alamouti combiner whitens `H[k]` and each received slot using the same
receive covariance before applying the orthogonal combining equations. The
existing scalar and per-tone scalar noise interfaces remain supported for
legacy callers.

The new receiver path is only enabled when
`cfg.estimation.spatialAwareLMMSE = true`. It should be evaluated against the
per-link/scalar receiver using matched frame seeds and measured BER/FER,
not merely inferred from estimator NMSE gains. The effective covariance relies
on the assumed channel prior, known spatial correlations, independent
unit-power streams, and block-static channel during the STBC pair.
