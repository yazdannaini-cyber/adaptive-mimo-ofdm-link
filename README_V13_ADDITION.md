# V13 — Waveform-Level Pilot Mapping and Time-Varying OFDM Validation

## Why V13 was added

The V12 pilot-tracking benchmark deliberately supplied abstract full-grid per-link LS observations at selected OFDM-symbol indices. It isolated temporal prediction, but it did not map pilots into a real OFDM resource grid or expose the pilot overhead and intra-symbol ICI seen by a waveform-level receiver.

V13 adds a practical, explicit pilot waveform to close that modelling gap.

## Pilot resource mapping

For a 2x2 MIMO link, each pilot occasion reserves two consecutive OFDM symbols:

- First pilot symbol: TX1 sends a known full-band QPSK pilot; TX2 is silent.
- Second pilot symbol: TX2 sends a known full-band QPSK pilot; TX1 is silent.
- Every remaining OFDM symbol carries independent QAM payload on both transmit antennas.

The two per-transmitter pilot timelines are staggered by one OFDM symbol. The temporal tracker is invoked independently for each transmit antenna using its own actual pilot timestamps. `pilotSpacing` is the period between occasions for the same transmitter and must be at least `nTx` to prevent pilot-slot overlap.

This is a time-division orthogonal training design. It is simple and identifiable, but costs two pilot OFDM symbols per occasion. Pilot overhead is included in the reported net spectral efficiency.

## Waveform-level channel and observation model

The grid is passed through OFDM IFFT and cyclic-prefix insertion. The complete serialized waveform is propagated through the V11 sample-varying MIMO FIR channel. The receiver removes the CP and applies the FFT before extracting LS pilot observations by dividing by the known full-band pilot sequence.

For each OFDM symbol, the reference diagonal frequency response is defined by averaging the time-varying taps over the useful FFT interval and taking an FFT over the tap dimension. The received FFT output is also compared with the diagonal-channel prediction. Their difference measures off-diagonal inter-carrier leakage due to within-symbol channel variation:

\[
Y_{\mathrm{ICI}}[k,r,s]=Y_{\mathrm{noiseless}}[k,r,s]-
\sum_t H_{\mathrm{diag}}[k,r,t,s]X[k,t,s].
\]

`ICI_Ratio` is the total power of this residual divided by the total diagonal-reference receive power over data symbols. The pilot-observation contamination variance similarly includes both thermal noise and ICI at the mapped pilot resources.

## Compared receivers

1. `hold-last`: causal last-pilot tracking per transmit antenna.
2. `linear`: noncausal interpolation/extrapolation reference.
3. `jakes-lmmse`: causal Jakes temporal LMMSE predictor using the configured Doppler and thermal pilot-noise variance.
4. `perfect-csi`: oracle diagonal-channel reference, used to isolate channel-estimation loss from unmodelled ICI and thermal noise.

The Rician deterministic LOS trajectory is supplied as a known model prior, consistent with the previous V12 benchmark. The tracker uses the thermal LS-observation noise variance; the additional ICI contamination is measured and reported but is not yet included in its observation-noise model. This deliberately exposes model mismatch under mobility.

## Receiver and metrics

The receiver currently performs per-subcarrier linear MMSE detection using the estimated diagonal channel and thermal-noise variance. It does **not** jointly equalize the full inter-carrier coupling matrix. As a result, measured BER includes the effects of residual ICI and channel-estimation error. V13 measures ICI; it does not claim ICI cancellation.

Per method and condition, the experiment reports:

- data-symbol channel NMSE and frame-bootstrap interval;
- uncoded BER, bit-level Wilson interval and frame-bootstrap interval;
- frame error ratio (FER; any payload bit error in a frame);
- pilot overhead and net spectral efficiency including CP;
- uncoded goodput approximation;
- measured ICI power ratio;
- pilot LS observation contamination variance;
- Jakes-LMMSE posterior-MSE calibration ratio.

## Default screening experiment

- Channels: Rayleigh and Rician.
- Spatial profiles: IID and joint-high (TX/RX correlation 0.9).
- Maximum Doppler: 0, 120 and 300 Hz.
- Same-transmitter pilot spacing: 4, 8 and 16 OFDM symbols.
- Pilot SNR: 15 dB.
- Payload: 16-QAM.
- Frame length: 24 OFDM symbols.
- Frames per condition: 8; 500 frame-bootstrap resamples.
- Four receiver references; 144 result rows.

This is a screening validation, not a publication-grade Monte Carlo campaign. Increase `numFrames` in `examples/pilot_mapped_ofdm_validation_run.m` (e.g. 50–100 or more) and rerun before making narrow confidence-interval or deployment claims.

## Run

```matlab
run_tests
addpath('examples');
results = pilot_mapped_ofdm_validation_run;
summary = pilot_mapped_results_summary(results);
```

Outputs are written to `pilot_mapped_results/`:

- `pilot_mapped_ofdm_validation.csv`
- `pilot_mapped_ofdm_validation.mat`
- Rayleigh/Rician NMSE, BER and ICI plots
- `pilot_mapped_results_summary` prints compact console summaries so the CSV does not need to be uploaded

## Current limitations

- TDM pilots use dedicated full-band pilot OFDM symbols; no comb-pilot or 2D time-frequency pilot optimisation is attempted.
- The Rician LOS trajectory and Doppler statistics are known to the tracking prior.
- Jakes-LMMSE currently models thermal LS noise but not the measured ICI contamination variance.
- Per-subcarrier MMSE detection treats ICI as residual impairment; there is no full ICI-matrix equalizer in V13.
- Perfect CSI means perfect knowledge of the diagonal channel response, not of the complete time-varying ICI coupling operator.
- No FEC, CRC, ARQ, CFO, phase noise or synchronization error is included.
