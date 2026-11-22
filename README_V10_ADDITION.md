# V10 — Covariance-Aware Detection

V10 closes the gap between spatially-aware channel estimation and the detector. The joint LMMSE setup now retains the full posterior covariance of channel-estimation error across spatial links at each subcarrier. For independent, unit-power transmit streams, this induces the receive-side effective-noise covariance

\[
R_{\mathrm{eff}}[k] = \sigma_n^2 I_{N_r} + \sum_{t=1}^{N_t} C_{e_t}[k].
\]

The linear MMSE detector uses the generalized form

\[
W[k] = \left(H^H[k]R_{\mathrm{eff}}^{-1}[k]H[k]+I\right)^{-1}
H^H[k]R_{\mathrm{eff}}^{-1}[k].
\]

The Alamouti receiver whitens the channel and received vectors using the same per-tone covariance before orthogonal combining. Existing scalar and per-tone scalar noise inputs remain supported.

## Run

From the project root:

```matlab
run_tests
addpath('examples');
results = spatial_aware_lmmse_comparison_run;
```

The comparison script runs paired per-link versus joint-estimator conditions over IID, TX-correlated, RX-correlated, moderately correlated, and highly correlated spatial profiles in Rayleigh and Rician channels. The joint-estimator arm now uses covariance-aware detection. The experiment is therefore an end-to-end receiver comparison; use `channel_estimator_mse_validation_run` for the isolated estimator-only NMSE study.
