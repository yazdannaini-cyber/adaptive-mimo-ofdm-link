# V14 changelog

## Full ICI matrix and block MMSE oracle

- Added `build_time_varying_ofdm_ici_matrix` to construct the exact
  frequency-coupling matrix from the sample-varying MIMO FIR tap sequence.
- Added `mimo_ici_equalize` for full-band linear MMSE/ZF detection across all
  subcarriers and transmit streams in one OFDM symbol.
- Added the `perfect-csi-ici` diagnostic receiver to the pilot-mapped OFDM
  experiment. It uses the true full ICI matrix and is explicitly an oracle
  reference; practical pilot estimation of this matrix is not implemented.
- Kept `perfect-csi` as the diagonal-channel oracle using the existing one-tap
  MMSE detector, enabling a direct ICI-cancellation comparison.
- Expanded the pilot-mapped experiment from 144 to 180 rows.
- Expanded `pilot_mapped_results_summary` with a direct diagonal-perfect-CSI
  versus full-ICI-oracle BER comparison.
- Added Tests 92–98 for static-channel reduction, exact waveform reconstruction,
  full/per-tone MMSE equivalence, noiseless Doppler recovery and validation.
- Updated the pilot-mapped frame smoke test to expect five receiver methods.

## Scope

No FEC, CRC, ARQ, CFO, phase-noise, timing-recovery or pilot-estimated full ICI
matrix is added in V14. The full-ICI receiver is a simulation oracle used to
measure the upper-bound value of ICI-aware detection.


## V14.1 runner fix

- Fixed the pilot-mapped experiment driver to include the already implemented
  `perfect-csi-ici` receiver.
- The run now evaluates all five methods and should produce 180 rows; the
  section-E BER comparison should contain 12 operating-point comparisons.
- Added `V14_1_VALIDATION_RUNNER_FIX.md` documenting the correction.
