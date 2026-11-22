# V15 validation status

- [x] V14.1 project source recovered and used as the parent version.
- [x] Frame-level adaptive MIMO/QAM selector implemented for pilot-mapped OFDM.
- [x] First-pilot-only decision path implemented; later pilots are used by the causal Jakes-LMMSE tracker.
- [x] Time-varying two-slot Alamouti MMSE decoder implemented.
- [x] Fixed SM-16QAM and fixed Alamouti-QPSK comparison links added with shared channel/pilot/noise realizations.
- [x] Net spectral-efficiency formula includes TDM pilot and CP overhead.
- [x] Added tests 99–107; expected total is 107.
- [x] Independent NumPy algebra check: noiseless time-varying Alamouti recovery maximum error approximately 4.3e-13.
- [x] Independent NumPy static-channel comparison: generalized versus conventional Alamouti maximum difference approximately 8.9e-16.
- [ ] Run `run_tests` in MATLAB and confirm all 107 tests pass.
- [ ] Run `pilot_mapped_adaptive_validation_run` and confirm 108 summary rows / 288 frame-log rows.
- [ ] Review measured BER versus predicted BER and both confidence intervals.
- [ ] Increase frames per condition before any low-BER or publication-grade claim.

**Important:** MATLAB and GNU Octave were unavailable in the build environment. The new MATLAB files have not been executed by a MATLAB runtime. Do not treat the expected test count or the expected output row count as confirmed until the user runs them.
