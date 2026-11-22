# V15.1.1 Changelog

- Fixed dimension mismatch in pilot-innovation ICI estimation with multiple
  transmit antennas and a full-MIMO Rician LOS prior.
- Slice the corresponding Tx prior before calling the single-link temporal
  Jakes-LMMSE tracker, while retaining support for a broadcast single-link prior.
- Validate prior tensor dimensions and finite values up front.
- Added regression Test 123 for two-transmitter prior slicing.
- Added version-specific entry point `run_tests_v15_1_1`.
