# Adaptive MIMO-OFDM Communication System

**A reproducible MATLAB project for wireless-system design and digital signal-processing algorithm evaluation.**

## 1. Project objective

This project implements and evaluates a frequency-selective MIMO-OFDM link. It is designed as a technical portfolio project for roles involving communication-system design, digital signal processing, channel estimation, and algorithm development.

The implementation includes:

- A configurable OFDM transmitter and receiver with IFFT/FFT and cyclic-prefix processing.
- 2x2 MIMO spatial multiplexing with orthogonal time-multiplexed training pilots.
- Rate-1 Alamouti STBC encoder/combiner and a matched spatial-multiplexing versus transmit-diversity comparison.
- AWGN-reference, Rayleigh multipath, and Rician multipath channel models.
- Gray-coded 4-QAM/QPSK, 16-QAM, and 64-QAM modulation/demodulation.
- Least-squares (LS) and frequency-domain linear-MMSE (LMMSE) channel estimation.
- Zero-forcing (ZF) and linear-MMSE MIMO detection.
- Fixed-modulation BER, frame-error-rate (FER; legacy output labels it BLER), and channel-estimation NMSE analysis.
- Channel-aware adaptive modulation using post-detection SINR and an analytical QAM-BER predictor, with explicit fallback and target-compliance reporting.
- CSV/MAT outputs, plots, deterministic tests, and both bit-level Wilson and frame-cluster bootstrap confidence intervals for BER.

The code uses base MATLAB language functions and does not require Communications Toolbox.

## 2. System model

For each OFDM subcarrier, the received MIMO signal is modelled as

\[
\mathbf y[k] = \mathbf H[k]\mathbf x[k] + \mathbf n[k],
\]

where \(\mathbf H[k]\) is the frequency response of the multipath MIMO channel. The default configuration is 2 transmit antennas, 2 receive antennas, 64 subcarriers, and a 16-sample cyclic prefix. The channel is assumed constant during one simulated frame.

### Channel estimation

Orthogonal pilot OFDM symbols are transmitted in time: pilot symbol \(t\) activates only transmit antenna \(t\). The LS estimate is

\[
\widehat H_{LS}[k] = Y_p[k]/\sqrt{P_p}.
\]

The LMMSE estimator uses a model-aware prior mean and frequency-domain diffuse covariance:

\[
\widehat{\mathbf h}_{LMMSE} = \boldsymbol\mu_h + \mathbf R_h
  (\mathbf R_h + \sigma^2/P_p\,\mathbf I)^{-1}
  (\widehat{\mathbf h}_{LS}-\boldsymbol\mu_h).
\]

The prior matches each configured channel: the AWGN reference has a deterministic identity response, Rayleigh has zero mean and the configured PDP, and Rician includes its deterministic LOS response plus diffuse covariance scaled by \(1/(K+1)\). The posterior channel-estimation error variance is also calculated for link adaptation.

### MIMO detection

The ZF detector uses the Moore-Penrose pseudoinverse. The linear MMSE detector uses

\[
\mathbf W_{MMSE} = (\mathbf H^H\mathbf H + \sigma_n^2\mathbf I)^{-1}\mathbf H^H.
\]

A diagonal gain correction is applied before hard QAM decisions to compensate for the shrinkage introduced by the MMSE filter.

## 3. Simulation experiments

### Experiment A — fixed modulation

For each channel model, QAM order, and Eb/N0 point, the same generated payload and received data are evaluated with four receiver configurations:

1. LS channel estimate + MMSE detector.
2. LMMSE channel estimate + MMSE detector.
3. Perfect channel knowledge + MMSE detector (ideal-CSI reference).
4. LMMSE channel estimate + ZF detector.

The results include BER, frame-error ratio (FER; exposed under legacy BLER column names), bit-level Wilson 95% intervals, frame-cluster bootstrap 95% intervals, and channel-estimation NMSE. The cluster bootstrap resamples complete frames to account for bits sharing a channel realization within each fading frame.

For fixed modulation, the complex frequency-domain noise variance is set to

\[
\sigma_n^2 = 1/(\log_2(M)\,10^{E_b/N_0/10}),
\]

assuming unit average QAM-symbol energy and an uncoded link. OFDM time-domain noise is scaled by \(1/N_{FFT}\) so that the FFT-domain noise variance matches the configured value.

### Experiment B — adaptive modulation

The receiver estimates post-detection SINR using the model-aware LMMSE channel estimate. Posterior channel-estimation uncertainty is added to the effective noise variance used by the MMSE detector. The default link-adaptation policy estimates the expected Gray-coded square-QAM BER over all subcarrier/stream SINRs using an analytical AWGN approximation, then selects the highest modulation whose predicted BER is no greater than the configured target. If no mode meets the target, QPSK is retained as a minimum-service fallback.

A legacy policy based on a lower post-SINR percentile and configurable thresholds remains available by setting `cfg.adaptive.selectionMethod = 'sinr-threshold'`. The adaptive experiment is indexed by Es/N0, not Eb/N0, because the selected number of bits per QAM symbol changes from frame to frame. The predicted BER is a link-selection estimate, not a guarantee. The output distinguishes (i) the fraction of frames for which a target-compliant mode was predicted to exist, (ii) the QPSK fallback fraction when no mode met the predicted target, and (iii) whether the measured aggregate BER and the Wilson and frame-cluster confidence upper bounds meet the target. The predicted aggregate BER is also weighted by each frame's actual payload-bit count before comparison with measured BER.

Net spectral efficiency accounts for pilot-symbol and cyclic-prefix overhead:

\[
\eta = N_t\,\bar k\,\frac{N_s}{N_s+N_t}\,\frac{N_{FFT}}{N_{FFT}+N_{CP}},
\]

where \(\bar k\) is the average number of bits per selected QAM symbol and \(N_s\) is the number of data OFDM symbols.

## 4. Repository structure

```text
MIMO_OFDM_Project_v12/
├── config/
│   └── default_config.m
├── src/
│   ├── channel/       # AWGN, Rayleigh, Rician channel generation and propagation
│   ├── detection/     # ZF/MMSE MIMO detectors and Alamouti combining
│   ├── estimation/   # LS/LMMSE estimators, channel priors, posterior error variance
│   ├── metrics/       # BER, BLER, confidence intervals and adaptive-link summaries
│   ├── phy/           # OFDM modulation/demodulation and Alamouti encoding
│   ├── plots/         # Results figures
│   ├── simulation/    # Fixed/adaptive experiments and MIMO transmission-mode comparison
│   └── utils/         # QAM modulator/demodulator and QAM BER approximation
├── tests/
│   └── run_all_tests.m
├── results/           # New run outputs plus preserved baseline CSV data
├── run_project.m
├── run_tests.m
└── README.md
```

## 5. How to run

1. Open MATLAB and set the current directory to this project folder.
2. Run the deterministic test suite:

   ```matlab
   run_tests
   ```

3. For a quick smoke test with only five frames per point, run these commands from the project root. The output is written to `quick_results/`:

   ```matlab
   addpath('examples');
   quick_run
   ```

4. Run the longer adaptive statistical validation (1,000 frames per point for Rayleigh and Rician):

   ```matlab
   addpath('examples');
   results = adaptive_validation_run;
   ```

   This focused run writes `validation_results/adaptive_validation_results.csv` and `.mat`.

5. Run the complete simulation:

   ```matlab
   [fixedResults, adaptiveResults] = run_project;
   ```

6. Inspect the generated files in `results/`:
   - `fixed_modulation_results.csv`
   - `adaptive_modulation_results.csv`
   - `complete_results.mat`
   - `BER_Fixed_*.png`
   - `NMSE_*.png`
   - `Adaptive_Performance_*.png`
   - `Adaptive_Feasibility_*.png`

All simulation settings are centralized in `config/default_config.m`. Reduce `cfg.simulation.numFrames` for a quick smoke test or increase it for more reliable Monte Carlo estimates.

## 6. Validation and interpretation

The 42-test suite checks QAM mapping, OFDM IFFT/FFT and cyclic-prefix processing, channel propagation, LS/LMMSE estimation, matched AWGN and Rician priors, posterior uncertainty, QAM BER prediction, ZF/MMSE detection, BER/frame-error calculations, frame-cluster bootstrap intervals, Alamouti STBC, spatial correlation, and the joint spatial-frequency LMMSE estimator. The suite also checks the joint spatial-frequency LMMSE estimator against the per-link estimator at zero correlation, validates its correlated-channel MSE behavior and Rician prior, and exercises its end-to-end adaptive integration. The simulation provides a perfect-CSI receiver reference to separate channel-estimation loss from detection loss.

A BER value of zero in a finite simulation does not mean the true BER is exactly zero. Inspect both Wilson and frame-cluster intervals and increase the number of independent frames before making high-SNR performance claims. If every frame has zero errors, the nonparametric cluster bootstrap is degenerate at zero; use the nonzero Wilson upper bound and collect more samples.

## 7. Assumptions and current limitations

- Perfect timing and carrier-frequency synchronization; no CFO, phase noise, or sampling-frequency offset.
- A block-static channel over each OFDM frame; independent channel realizations between frames.
- Cyclic prefix is long enough to cover the configured channel memory.
- No channel coding or retransmission protocol. Optional Kronecker spatial correlation is supported for the diffuse Rayleigh/Rician channel component; when configured, the V8 joint spatial-frequency LMMSE estimator can exploit that prior.
- The cluster bootstrap assumes independently generated frames; its validity is conditional on the configured frame-generation process and is not a guarantee under an unmodelled channel.
- Pilot symbols are orthogonal in time and are assumed known at the receiver.
- The Rician LOS component is a deterministic phase matrix on the first tap; this is a controllable abstract model, not a geometry-calibrated propagation model. Its known mean is included in the LMMSE prior.
- The AWGN option is a deterministic identity MIMO channel plus noise, included as a sanity-check baseline; its matched prior is therefore exact.

These assumptions are explicit so the simulation results are interpretable. Suitable future extensions include correlated MIMO channels, pilot contamination/overhead studies, Doppler time variation, CFO compensation, soft detection, coded OFDM, and ML-assisted channel estimation.

## 8. Suggested interview explanation

“I implemented a modular MIMO-OFDM baseband link and evaluated how channel estimation and detection algorithms affect reliability and spectral efficiency. Orthogonal pilots provide LS channel estimates; the model-aware LMMSE estimator exploits the channel PDP and includes the Rician LOS mean. At the receiver, I compared ZF and MMSE detection, including a perfect-CSI reference. I evaluated BER, BLER, NMSE, and confidence intervals under AWGN, Rayleigh, and Rician channels, then added an estimation-uncertainty-aware link-adaptation layer that predicts uncoded QAM BER from post-detection SINR and selects the highest-rate mode meeting a target. The simulation is reproducible and includes unit tests for the main DSP blocks.”

Only describe experiments as completed after running the MATLAB project and reviewing the generated output on your machine.


## 9. V5 experiment — spatial multiplexing versus Alamouti STBC

V5 adds a second MIMO transmission strategy to expose the reliability/rate trade-off:

- **Spatial multiplexing:** two independent QAM streams per channel use, recovered with linear MMSE detection.
- **Alamouti STBC:** two transmit antennas send a rate-1 orthogonal space-time code over adjacent OFDM symbols, recovered with a coherent Alamouti combiner. It trades the two-stream multiplexing gain for transmit diversity.

Both modes use the same per-antenna symbol energy, the same channel realization for each paired Monte-Carlo frame, and the same orthogonal channel-estimation pilots. The channel is block-static over each OFDM frame and must remain sufficiently constant over each Alamouti symbol pair.

Run the regression tests (expected: 29 tests after V6):

```matlab
run_tests
```

Run the transmission-mode comparison:

```matlab
addpath('examples');
comparisonResults = mimo_mode_comparison_run;
```

The default sweep covers AWGN identity, Rayleigh, and Rician channels; 4/16/64-QAM; Es/N0 = 0:8:32 dB; and 200 frames per configuration. CSV/MAT outputs and BER/FER/spectral-efficiency figures are saved in `comparison_results/`.

See `MIMO_MODE_COMPARISON_NOTES.md` for equations and interpretation. See `V4_VALIDATION_RESULTS.md` and `results/adaptive_validation_v4_user_results.csv` for the complete 1,000-frame-per-point adaptive results provided from the previous MATLAB run.

### FER and goodput terminology

The legacy adaptive result column called `BLER` is calculated as the fraction of simulated data frames with at least one bit error. Because the current link has no FEC decoder, CRC, or ARQ, this is more precisely a frame error rate (FER), not post-FEC block error rate. The current goodput is bit-level goodput: it counts correctly received bits even when a frame also contains errors. It is not CRC-accepted packet goodput.

The V5 comparison reports this metric as `FER` to make the definition explicit. Neither spatial multiplexing nor Alamouti comparison includes channel coding, CRC, or retransmissions.


## 10. V6 experiment — joint MIMO-mode and modulation adaptation

V6 adds a frame-by-frame link adaptation policy that selects both the transmission architecture and QAM order from pilot-based channel estimates:

- **2x2 spatial multiplexing + linear MMSE detection** for higher multiplexing gain.
- **2x2 Alamouti STBC + coherent combining** for transmit diversity.

For every frame, the receiver obtains one channel estimate and evaluates the predicted uncoded Gray-QAM BER for all six mode/modulation candidates. The policy first keeps candidates whose predicted BER is at or below the configured target. It then chooses the candidate with the highest net spectral efficiency after pilot and cyclic-prefix overhead. If no candidate meets the target, it chooses the minimum-predicted-BER candidate and explicitly marks that frame as target-infeasible.

Run the 29-test regression suite for the historical V6 snapshot:

```matlab
run_tests
```

Run the joint adaptation validation (1,000 frames per channel/SNR point, 2,000 frame-bootstrap resamples):

```matlab
addpath('examples');
results = adaptive_mimo_mode_run;
```

The default study covers Rayleigh and Rician channels at Es/N0 = 12, 16, 20, 24, 28, and 32 dB, with QPSK, 16-QAM, and 64-QAM candidates. CSV/MAT results and BER, spectral-efficiency, and mode-selection figures are written to `adaptive_mimo_results/`.

See `examples/README_ADAPTIVE_MIMO_MODE.md` for the algorithm and interpretation. The new V6 code and tests are provided for MATLAB validation; do not describe the 29-test suite or adaptive study as executed until they have run successfully in MATLAB.


## V7 spatial-correlation extension

See `README_V7_ADDITION.md` and `SPATIAL_CORRELATION_MODEL.md`. V7 adds an optional
Kronecker spatial-correlation stress test for Rayleigh/Rician diffuse fading.
The new feature is disabled by default for backward compatibility.

## V8: Joint spatial-frequency LMMSE

V8 adds an optional joint LMMSE estimator using the Kronecker TX/RX spatial
covariance together with the frequency covariance induced by the PDP. It is
backward-compatible and disabled by default:

```matlab
cfg = default_config();
cfg.channel.spatialCorrelation.enabled = true;
cfg.channel.spatialCorrelation.txRho = 0.7;
cfg.channel.spatialCorrelation.rxRho = 0.7;
cfg.estimation.spatialAwareLMMSE = true;
```

For a paired robustness comparison of the per-link estimator and the new joint
estimator, run:

```matlab
addpath('examples');
results = spatial_aware_lmmse_comparison_run;
```

See `SPATIAL_AWARE_LMMSE.md`, `V8_CHANGELOG.md`, and
`V8_INDEPENDENT_NUMERIC_CHECK.md` for the model, limitations, and independent
numerical cross-check. The latter is not a substitute for MATLAB execution. Run `run_tests` before using the new estimator.


### V8 validation and paired estimator comparison

Run the full regression suite (42 tests):

```matlab
run_tests
```

Compare the original per-link LMMSE estimator with the joint spatial-frequency
LMMSE estimator using matched per-frame channel and pilot-noise seeds:

```matlab
addpath('examples');
results = spatial_aware_lmmse_comparison_run;
```

The comparison evaluates five TX/RX correlation profiles, Rayleigh and Rician
channels, 16/24/32 dB, and 500 frames per estimator/profile/channel/SNR run.
The CSV has explicit estimator and correlation-profile columns. V8 MATLAB
results are not included in this source archive; run the experiment before
claiming an empirical advantage for the joint estimator.

## V9: direct paired channel-estimator validation

The V8 user-run console output showed that per-link and joint LMMSE often lead
to very similar end-to-end adaptive BER and mode choices. Since final BER can
hide estimator differences behind adaptive mode selection, V9 adds a direct
estimator-only comparison.

Run the current regression suite (48 tests):

```matlab
run_tests
```

Run the direct paired NMSE experiment:

```matlab
addpath('examples');
results = channel_estimator_mse_validation_run;
```

For every frame, both estimators see the same channel and same noisy orthogonal
pilots. The output includes frame-bootstrap confidence intervals for NMSE, the
paired empirical MSE difference, and a comparison between empirical MSE and the
model-predicted posterior MSE. The 50-row sweep covers Rayleigh/Rician channels,
five spatial-correlation profiles, and pilot SNR = [0 8 16 24 32] dB.
See `V8_USER_RESULTS_REVIEW.md` and `README_V9_ADDITION.md` for interpretation.
The V9 tests and study have not yet been run in MATLAB.

## V10: covariance-aware detection

V9's estimator-only paired experiment established that joint spatial-frequency
LMMSE reduces empirical channel-estimation MSE under correlated fading, while
the earlier adaptive link comparison showed little end-to-end change. The V10
receiver update addresses an information-loss point: the previous link adapter
reduced posterior uncertainty to a scalar upper bound and discarded the
receive-space cross-covariance.

The joint estimator now provides the full posterior covariance of spatial
channel links on each subcarrier. For independent unit-power transmit streams,
the induced receive-side effective-noise covariance is

\[
R_{\mathrm{eff}}[k] = \sigma_n^2 I_{N_r} +
\sum_{t=1}^{N_t} C_{e_t}[k],
\]

where `C_e_t[k]` is the posterior receive-vector covariance for channel errors
in transmit column `t`. The spatial MMSE detector uses the generalized
covariance-weighted filter, and the Alamouti receiver first whitens the received
vectors and channel by this covariance before combining.

The original scalar/vector noise inputs are retained. V10's covariance path is
opt-in through `cfg.estimation.spatialAwareLMMSE = true`; the default per-link
configuration remains backward-compatible.

### V10 regression and end-to-end comparison

Run the expanded regression suite (55 tests):

```matlab
run_tests
```

Then rerun the per-link versus joint estimator comparison:

```matlab
addpath('examples');
results = spatial_aware_lmmse_comparison_run;
```

Important distinction: in V10, the per-link arm uses the existing scalar
uncertainty approximation, while the joint arm uses the full posterior spatial
covariance in both spatial-MMSE detection and Alamouti combining. This is an
end-to-end receiver comparison. Use `channel_estimator_mse_validation_run` for
the separate estimator-only NMSE comparison. V10 MATLAB tests and simulations
must be run before claiming a measured receiver-level gain.

See `V9_USER_RESULTS_ANALYSIS.md`, `V10_CHANGELOG.md`,
`V10_VALIDATION_STATUS.md`, and `README_V10_ADDITION.md`.

## V10.1 Hotfix

The V10.1 hotfix corrects MATLAB's trailing-singleton `ndims` behavior for one-symbol detector inputs. See `V10_1_HOTFIX.md`. The test suite now includes Test 55 for this regression.


## V11 — Sample-Varying Doppler Channels

V11 extends the channel layer beyond block-static fading. A sample-varying MIMO FIR model now supports Rayleigh and Rician diffuse fading with a finite sum-of-sinusoids approximation to the Jakes temporal autocorrelation. The same TX/RX Kronecker spatial-correlation model is retained. Since taps may vary inside an OFDM symbol, the model produces inter-carrier interference (ICI).

Run the expanded test suite (69 tests):

```matlab
run_tests
```

Run the Doppler/ICI validation study:

```matlab
addpath('examples');
[results, acfResults] = doppler_channel_validation_run;
```

The validation sweeps maximum Doppler frequencies of 0, 30, 120 and 300 Hz, compares the diffuse temporal ACF with the Jakes `J0` reference, and measures off-subcarrier leakage from a single-tone OFDM probe. Results are saved in `doppler_validation_results/`.

V11 validates the time-varying channel model and ICI mechanism; it does not yet implement pilot-based temporal tracking. See `DOPPLER_CHANNEL_MODEL.md`, `V10_1_USER_RESULTS_ANALYSIS.md`, `V11_CHANGELOG.md` and `V11_VALIDATION_STATUS.md`.


## V11.1 test hotfix

V11.1 makes the Rician K-factor unit test statistically robust by averaging independent seeded realizations. See `V11_1_HOTFIX.md`. The channel model itself is unchanged.


## V12 — Pilot tracking under Doppler

V12 adds three temporal tracking methods for channel estimates sampled at OFDM-symbol rate: causal hold-last, noncausal linear interpolation/extrapolation (offline reference), and a causal Jakes-covariance temporal LMMSE predictor. The experiment sweeps Rayleigh/Rician channels, IID/moderate/high spatial correlation, four Doppler values, four pilot update periods, and three pilot-observation SNR values.

Run the new experiment with:

```matlab
addpath('examples');
results = pilot_tracking_validation_run;
```

See `PILOT_TRACKING_MODEL.md` and `README_V12_ADDITION.md` for equations, assumptions and limitations. V12 adds tests 70–78; the expected inherited suite total is 78 tests. Run `run_tests` in MATLAB before treating the new V12 code as validated.


## V12.1 Hotfix

V12.1 fixes MATLAB name shadowing in `temporal_channel_tracking_apply` and adds uppercase-method regression tests. See `V12_1_HOTFIX.md`. Run `run_tests`; expected test count: 81.

## V13 — Waveform-Level Pilot Mapping and Time-Varying OFDM

V12 used abstract full-grid orthogonal LS observations to isolate temporal
tracking. V13 adds an explicit 2x2 time-division pilot resource mapper: each
pilot occasion reserves two consecutive OFDM symbols, with one transmit
antenna active per pilot symbol and known full-band QPSK pilots. The remaining
OFDM symbols carry 16-QAM payload on both transmit antennas. The staggered
per-transmitter pilot timestamps are passed to separate causal temporal
trackers.

The waveform is modulated with IFFT/CP, propagated through the existing
sample-varying MIMO FIR channel, and demodulated with CP removal/FFT. V13
computes a useful-symbol-averaged diagonal channel reference, extracts LS
pilot observations from the actual received grid, and measures the residual
between the waveform output and its diagonal-channel prediction as ICI.

The compared receivers are hold-last, linear interpolation (offline
noncausal reference), Jakes-LMMSE and perfect diagonal-CSI. The current
per-subcarrier MMSE detector does **not** cancel the full ICI coupling matrix;
ICI is measured and remains part of the BER impairment. The Jakes tracker is
provided thermal pilot-noise variance, while pilot contamination including
ICI is reported separately. This makes the model mismatch visible rather
than silently treating ICI as thermal noise.

Run V13 with:

```matlab
run_tests
addpath('examples');
results = pilot_mapped_ofdm_validation_run;
summary = pilot_mapped_results_summary(results);
```

The default screening experiment sweeps Rayleigh/Rician channels, IID and
joint-high spatial correlation, Doppler 0/120/300 Hz, pilot spacing 4/8/16,
and 15 dB pilot SNR. It uses 8 frames per condition and 500 frame-bootstrap
resamples and writes 144 estimator/condition rows to `pilot_mapped_results/`.
This is an initial screening study; increase `numFrames` in the example
before making narrow statistical claims. V13 adds Tests 82–91; the expected
suite total is 91. The user subsequently confirmed that the V13 tests passed
and supplied the completed 144-row summary. See `V13_VALIDATION_STATUS.md`
and `V13_USER_RESULTS_ANALYSIS.md`.

See `README_V13_ADDITION.md`, `PILOT_TRACKING_MODEL.md`,
`DOPPLER_CHANNEL_MODEL.md`, `V12_USER_RESULTS_ANALYSIS.md`,
`V13_CHANGELOG.md` and `V13_VALIDATION_STATUS.md` for assumptions and details.

## V14 — Full ICI-Matrix MMSE Reference

V13 established an actual waveform-level OFDM pilot path and measured the
ICI residual, but its `perfect-csi` receiver still used a diagonal,
per-subcarrier MMSE detector. V14 adds the exact full frequency-coupling
matrix computed from the sample-varying channel taps and a full-band block
MMSE detector called `perfect-csi-ici`.

The new method is deliberately an oracle reference: it uses the true channel
tap sequence to construct every diagonal and off-diagonal subcarrier
coefficient. The existing TDM pilot pattern does not yet estimate the complete
ICI matrix. Comparing `perfect-csi` with `perfect-csi-ici` quantifies how much
of the V13 BER is caused by unmodelled ICI versus thermal noise and MIMO
spatial conditioning.

Run V14 with:

```matlab
run_tests
addpath('examples');
results = pilot_mapped_ofdm_validation_run;
summary = pilot_mapped_results_summary(results);
```

The V14.1 runner evaluates five references, including `perfect-csi-ici`, and
writes 180 rows (36 operating conditions x 5 receiver references). Section E
prints 12 diagonal-perfect-CSI versus full-ICI BER comparisons. V14 adds Tests
92–98; the expected suite total is 98 tests. The initial V14 runner omitted
the fifth method; V14.1 corrects this packaging defect (see
`V14_1_VALIDATION_RUNNER_FIX.md`).
The full matrix adds a larger block solve per data OFDM symbol, so the
simulation is more computationally demanding than the V13 diagonal receiver.

See `README_V14_ADDITION.md`, `ICI_EQUALIZATION_MODEL.md`,
`V13_USER_RESULTS_ANALYSIS.md`, `V14_CHANGELOG.md` and
`V14_VALIDATION_STATUS.md` and `V14_1_VALIDATION_RUNNER_FIX.md` for derivation,
assumptions and validation status.

## V15: pilot-driven adaptive MIMO-OFDM

The latest V15 extension integrates frame-level selection between 2-stream Spatial Multiplexing and rate-one Alamouti STBC, jointly with QPSK/16-QAM selection, into the actual mapped-pilot time-varying MIMO-OFDM waveform. The transmitter makes its decision from the first orthogonal pilot occasion; later pilots feed a causal Jakes-LMMSE tracker. A two-slot Alamouti MMSE detector uses separate channel estimates at each payload time.

Run `run_tests` (107 expected tests), then:

```matlab
addpath('examples');
results = pilot_mapped_adaptive_validation_run;
summary = pilot_mapped_adaptive_results_summary(results);
```

The screening run is expected to produce 108 strategy/condition rows and 288 per-frame log rows. It compares the adaptive link with fixed SM-16QAM and fixed Alamouti-QPSK using paired channel, pilot and noise realizations. Read `README_V15_ADDITION.md` and `ADAPTIVE_PILOT_MAPPED_MIMO_MODEL.md` for the full model and limitations. V15 MATLAB execution is not yet confirmed; see `V15_VALIDATION_STATUS.md`.


## V15.1: burst-adaptive pilot-mapped MIMO-OFDM

V15.1 supersedes the V15 frame-level adaptation experiment. It makes a new
Spatial-Multiplexing/Alamouti and QPSK/16-QAM decision after each orthogonal
TDM pilot occasion. Each decision uses only pilot observations available at
that time, so later pilots are not used to select an earlier data burst.

The predictor now includes temporal-LMMSE posterior channel-error variance
and a causal pilot-innovation estimate of excess pilot disturbance. The
innovation is measured against the Jakes prediction error plus thermal LS
noise; unexplained excess is treated as pilot contamination/ICI/model-mismatch
variance. At the first occasion, where no innovation history exists, the code
uses a small-normalized-Doppler analytical ICI floor. Per-pilot measurement
noise variances are supported by the temporal tracker, and the same estimated
impairment levels are passed to the actual MMSE detectors. The receiver
also uses each burst’s own causal channel estimate, so future pilots do not
leak into detection of previously transmitted data.

Run the uniquely named V15.1 suite to avoid stale V15 functions remaining on
the MATLAB path:

```matlab
clear functions
addpath(genpath(pwd));
run_tests_v15_1
```

Then run the burst-adaptive study:

```matlab
addpath('examples');
results = pilot_mapped_adaptive_validation_run_v15_1;
summary = pilot_mapped_adaptive_results_summary_v15_1(results);
```

The screening experiment is expected to produce 108 strategy/condition rows,
288 frame rows and 864 per-burst decision rows. Detailed CSV/MAT/PNG outputs
are written to `pilot_mapped_results_v15_1/`.

The uncertainty model is an equivalent scalar receive-disturbance
approximation: temporal posterior variance is summed across transmit links,
and the pilot-innovation excess is treated as additional receiver noise. It
does not perform full ICI-matrix cancellation or model all cross-link error
covariances. The predicted BER remains an uncoded Gray-QAM AWGN approximation
applied to post-detection SINR. The measured BER and its confidence intervals
must be used to judge target compliance. The default 8 frames per condition
is a screening run, not publication-grade statistical evidence.

See `README_V15_1_ADDITION.md`, `BURST_ADAPTIVE_MCS_MODEL.md`,
`V15_1_CHANGELOG.md` and `V15_1_VALIDATION_STATUS.md` for details.


## V15.1.1 Hotfix

V15.1.1 fixes the per-transmitter LOS-prior dimension mismatch in causal ICI
innovation estimation. Run `run_tests_v15_1_1` from the extracted project root;
123 tests are expected. See `V15_1_1_HOTFIX.md` and `V15_1_1_VALIDATION_STATUS.md`.

---

## V15.2 — Statistical validation and packet-level reliability

V15.2 adds a 100-frame default validation campaign, CRC-16/CCITT-FALSE packet
checking per pilot/data burst, a terminated rate-1/2 (7,5) convolutional FEC
profile with hard-decision Viterbi decoding, and separate payload BER, coded BER,
packet error rate, frame error rate, net payload spectral efficiency, bit-goodput
and CRC-accepted packet-goodput metrics.

Run the version-specific test suite from the project root:

```matlab
clear functions
addpath(genpath(pwd));
run_tests_v15_2
```

Run the default full validation (100 frames per condition):

```matlab
addpath('examples');
results = pilot_mapped_adaptive_validation_run_v15_2;
summary = pilot_mapped_adaptive_results_summary_v15_2(results);
```

See `V15_2_MODEL.md` for coding and metric definitions, and
`V15_2_VALIDATION_STATUS.md` for the expected row counts and validation caveats.

---

## V15.2.1 — Soft-decision FEC receiver update

V15.2.1 preserves the 36-condition, 100-frame/condition statistical experiment and changes only the coded receiver's demapping/decoding path. The V15.2 hard QAM decisions are replaced by max-log soft LLRs derived from each MMSE/Alamouti stream's post-equalization SINR. The rate-1/2 terminated (7,5) convolutional decoder now uses soft-decision Viterbi branch metrics, followed by the same CRC-16 check. Uncoded packets use the LLR signs for hard CRC checking.

New entry points:

```matlab
clear functions
addpath(genpath(pwd));
run_tests_v15_2_1

addpath('examples');
results = pilot_mapped_adaptive_validation_run_v15_2_1;
summary = pilot_mapped_adaptive_results_summary_v15_2_1(results);
```

The new suite expects 141 tests. Full validation expects 216 summary rows, 21,600 frame-log rows, and 10,800 burst-decision rows at the default 100 frames/condition. Results are written to `pilot_mapped_results_v15_2_1` so V15.2 files are preserved.

The soft demapper uses max-log LLRs with a scalar Gaussian approximation derived from post-SINR. It improves on discarding reliability through hard decisions, but it is not a full joint MAP receiver for correlated residual interference or ICI. MATLAB tests have not been executed in the build environment; confirm the 141-test suite locally before relying on the results.

---

## V15.3 — Joint adaptive MIMO, modulation and FEC

V15.3 adds a single online controller that selects MIMO mode (Spatial
Multiplexing or Alamouti), modulation (QPSK or 16-QAM), and coding (uncoded
CRC-16 or rate-1/2 convolutional FEC + CRC-16) for every pilot-to-pilot data
burst. The controller maximizes posterior expected CRC-accepted goodput under
a configurable PER target using causal CRC ACK/NACK feedback. If no candidate
is feasible, it uses an explicitly marked maximum-goodput fallback.

The six V15.2.1 adaptive/fixed references are retained, and a seventh
`joint-adaptive` candidate is appended. The default validation now uses 300
frames per condition over 36 channel/profile/Doppler/SNR conditions.

Run the version-specific suite:

```matlab
clear functions
addpath(genpath(pwd));
run_tests_v15_3
```

Expected suite size: 160 tests.

Run the default 300-frame/condition experiment:

```matlab
addpath('examples');
results = pilot_mapped_adaptive_validation_run_v15_3;
summary = pilot_mapped_adaptive_results_summary_v15_3(results);
```

Expected: 252 summary rows, 75,600 frame-log rows, 32,400 joint burst
decisions and 32,400 reference burst predictions. Results are written to
`pilot_mapped_results_v15_3/`. A short 30-frame/condition functional run is
possible with `struct('numFrames',30,'numBootstrap',200)`, but it is not
sufficient for stringent low-PER claims.

The PER prior and FEC coding-gain curve are link-adaptation heuristics, not
guaranteed bounds. The controller updates its per-action estimates with
previously completed CRC ACK/NACK outcomes. Review measured BER/PER and
frame-cluster confidence intervals before making compliance claims. The
build environment did not include MATLAB or GNU Octave, so the 160-test suite
must be run locally before calling this release MATLAB-validated.

See `README_V15_3_ADDITION.md`, `V15_3_MODEL.md`, `V15_3_CHANGELOG.md` and
`V15_3_VALIDATION_STATUS.md`.
