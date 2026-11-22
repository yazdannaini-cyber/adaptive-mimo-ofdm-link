# V15.1 Changelog

- Replaced one-shot frame-level adaptation with causal adaptation at every
  pilot/data burst boundary.
- Added per-burst predicted channel arrays and temporal-LMMSE posterior
  variance propagation into candidate post-SINR and receiver MMSE detection.
- Added causal pilot-innovation excess estimation for residual pilot ICI/model
  mismatch, including a small-Doppler analytical floor at the first occasion.
- Extended `temporal_channel_tracker_setup` to accept scalar or per-pilot
  observation-noise variances while preserving scalar backward compatibility.
- Extended the generalized Alamouti detector to support different noise
  variances for its two OFDM slots.
- Replaced retrospective all-frame receiver CSI with per-burst CSI assembly,
  preventing future-pilot leakage in detection of earlier data symbols.
- Added burst decision and frame-level diagnostic logs and summaries.
- Added Tests 108–122 (15 additional tests; 122 expected total).
- Corrected the LS-noise unit conversion: thermal noise is divided by pilot
  power, while pilot-innovation ICI variance is already in LS channel units
  and is not divided by pilot power again.
- Added unique V15.1 test and validation entry-point names to avoid MATLAB
  path collisions with older V15 files.
