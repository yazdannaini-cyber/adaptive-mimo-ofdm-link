# V14.1 — Pilot-mapped validation runner correction

## Why this patch exists

The initial V14 package added the `perfect-csi-ici` receiver to
`simulate_pilot_mapped_ofdm_frame` and added the section-E comparison to
`pilot_mapped_results_summary`, but the experiment driver still listed only
the four V13 methods. Consequently, that driver generated 144 rows and left
section E empty, exactly as in the V13 baseline. This was a packaging defect
in V14, not a problem with the user's MATLAB run.

## Correction

- Added `perfect-csi-ici` to the method list in
  `examples/pilot_mapped_ofdm_validation_run.m`.
- Updated run metadata and comments to describe the oracle reference.
- Updated summary comments to identify V14 result support.
- Retained all V14 receiver code and tests.

## Expected output

There are 36 operating conditions (2 channels x 2 spatial profiles x 3 Doppler
values x 3 pilot spacings) and five receiver methods, so the run should produce
180 rows. Section E should print 12 comparisons (2 channels x 2 profiles x 3
Doppler values) for pilot spacing 8.

## MATLAB verification

MATLAB execution is pending user-side. Do not call the V14.1 experiment
MATLAB-validated until the 98-test suite completes and the pilot-mapped run
reports 180 rows with non-empty section E.
