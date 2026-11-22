# V14 validation status

- Inherited V13 baseline: user reported that the V13 test suite was successful
  and supplied the 144-row `pilot_mapped_results_summary` output.
- New V14 tests added: 92–98 (7 tests).
- Expected deterministic MATLAB test count: 98.
- New receiver reference: `perfect-csi-ici` (full ICI oracle).
- New experiment output: 180 rows when the V14 pilot-mapped run completes.
- Static source review: completed.
- Independent NumPy numeric checks:
  - Full frequency-coupling matrix reproduced the time-domain waveform with
    maximum complex errors between 5.0e-16 and 9.1e-16 for tested examples.
  - Static-channel full matrix reduced to diagonal per-tone channel blocks
    with maximum mismatch approximately 3.6e-15.
  - Full-block MMSE output matched independent per-tone MMSE to approximately
    3.1e-15 for a static block-diagonal OFDM channel.
- MATLAB execution of the new V14 tests and simulation: pending user run.

Do not describe V14 as MATLAB-validated until `run_tests` and the V14 pilot-
mapped experiment have completed successfully in MATLAB.


## V14.1 runner correction

- The initial V14 archive accidentally retained a four-method list in the
  pilot-mapped validation driver, although the frame simulator already
  supported `perfect-csi-ici`. This caused a 144-row run and an empty section E.
- V14.1 corrects the driver method list to include all five methods.
- Expected output: 180 rows and 12 diagonal-versus-full-ICI comparisons.
- MATLAB run remains pending user-side confirmation.
