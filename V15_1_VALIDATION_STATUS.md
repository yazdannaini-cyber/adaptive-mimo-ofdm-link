# V15.1 Validation Status

## Implemented

- Burst-adaptive link mode and QAM selection at every pilot occasion.
- Causal Jakes-LMMSE prediction with per-pilot observation-noise variance.
- Posterior channel-error variance and pilot-innovation ICI excess included in
  predicted and detected effective noise.
- Two-slot Alamouti prediction and detection with slot-dependent disturbance.
- Burst-causal receiver CSI assembly; no future pilot estimates leak into
  detection of earlier data bursts.
- Expected summary/frame/burst logs and independent versioned entry points.

## Static review completed

- Reviewed function signatures, result-field consistency and output dimensions.
- Preserved scalar-noise compatibility for existing temporal tracking and
  Alamouti detector calls.
- Added 15 new tests, bringing the expected test suite total to 122.
- Added a unit-consistent pilot LS noise profile: thermal variance is scaled
  by pilot power, whereas the already-normalized pilot-innovation ICI estimate
  is added directly in LS channel units.
- Archive must be tested by the user in MATLAB; this environment has neither
  MATLAB nor GNU Octave, so no MATLAB execution is claimed.

## User validation required

1. Run `clear functions; addpath(genpath(pwd)); run_tests_v15_1` from the
   extracted `MIMO_OFDM_Project_v15_1` root. Expected: 122 tests pass.
2. Run `addpath('examples'); results = pilot_mapped_adaptive_validation_run_v15_1;`
   followed by `summary = pilot_mapped_adaptive_results_summary_v15_1(results);`.
3. Expected rows for the default screening configuration: 108 summary rows,
   288 frame-log rows and 864 burst-decision rows.

The 8-frame default is a screening run. BER-target compliance must be
established from measured BER confidence intervals, not predicted feasibility.
