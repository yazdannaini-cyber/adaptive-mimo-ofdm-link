# README V15.1 Addition

V15.1 upgrades V15 frame-level link adaptation to **pilot-driven burst-level
link adaptation**. The MIMO mode (Spatial Multiplexing or Alamouti STBC) and
QAM order (QPSK or 16-QAM) can be re-selected after every orthogonal pilot
occasion.

New improvements:
- Strictly causal per-burst decisions using only pilots observed so far.
- Per-pilot measurement-noise variance support in the Jakes-LMMSE tracker.
- Posterior channel-error variance included in predicted and actual receiver
  effective noise.
- Causal pilot-innovation excess estimate for ICI/pilot contamination, with
  a small-Doppler floor for the first occasion. Thermal noise is converted
  to LS variance by dividing by pilot power; the innovation-derived ICI
  variance is already in LS channel units and is added without re-scaling.
- Burst-specific prediction uses H estimates at the actual OFDM data slots.
- Two-slot Alamouti MMSE supports different noise levels in its two slots.
- Frame, burst-decision and summary CSV logs are exported independently.
- The actual receiver uses the CSI estimate available at each burst boundary;
  later pilot observations do not leak into detection of earlier symbols.

Run:

```matlab
clear functions
addpath(genpath(pwd));
run_tests_v15_1  % 122 deterministic tests
addpath('examples');
results = pilot_mapped_adaptive_validation_run_v15_1;
summary = pilot_mapped_adaptive_results_summary_v15_1(results);
```

Expected screening output: 108 summary rows, 288 frame records and 864 burst
decisions. See `BURST_ADAPTIVE_MCS_MODEL.md` for formulas, assumptions and
limitations. The package also includes a unit-conversion regression test.
V15.1 MATLAB execution remains unverified until the user runs the
new test suite and validation script in MATLAB.
