# V15 addition: pilot-driven adaptive MIMO-OFDM

V15 integrates adaptive spatial-mode and QAM selection into the waveform-level pilot-mapped 2x2 MIMO-OFDM simulation.

## What's new

- Frame-level mode selection between 2-stream Spatial Multiplexing and rate-one Alamouti STBC.
- Joint QPSK / 16-QAM selection from pilot-derived post-equalization SINR.
- First-pilot-only causal link decision; subsequent pilots support Jakes-LMMSE channel tracking.
- A generalized Alamouti MMSE receiver that uses separate channel estimates for the two data slots.
- Paired fixed-link references: SM-16QAM and Alamouti-QPSK, sharing channel taps, pilots and noise with the adaptive link.
- Measured BER, FER, Wilson and frame-cluster bootstrap intervals, mode mix, modulation mix, Net Spectral Efficiency and uncoded bit-goodput.
- Nine new regression/integration tests (V14.1's 98 tests become 107).

## Run in MATLAB

From the project root:

```matlab
run_tests
```

Expected test count: **107**.

Then run the pilot-mapped adaptive experiment:

```matlab
addpath('examples');
results = pilot_mapped_adaptive_validation_run;
summary = pilot_mapped_adaptive_results_summary(results);
```

Expected result table: **108 rows** (36 channel/profile/Doppler/SNR conditions × 3 strategies), plus **288 frame-log rows** (36 conditions × 8 frames).

Output files are written under `pilot_mapped_results_v15/`:

- `adaptive_link_results.csv`
- `adaptive_link_frame_log.csv`
- `adaptive_link_results.mat`
- per-channel BER and mode-selection PNG plots

## Important interpretation

The MCS selector uses a post-SINR-based analytical BER approximation. Its target-feasible flag records the prediction only. Target-BER compliance must be judged from measured BER plus Wilson and frame-cluster upper confidence limits.

The Alamouti detector models channel variation between the two codeword OFDM slots, but it is still a per-subcarrier detector. It does not cancel the full within-symbol ICI matrix introduced by the V11/V13 time-varying channel.

The default eight frames per condition are for initial screening. Increase `numFrames` in `examples/pilot_mapped_adaptive_validation_run.m` (e.g., 100 or more) before reporting precise low-BER performance.

The summary helper reports predicted BER/target-feasibility alongside measured BER and maximum condition-level Wilson/frame-bootstrap upper bounds. These maxima are not pooled confidence intervals.
