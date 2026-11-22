# V15 changelog

## Adaptive transmission on mapped-pilot OFDM

- Added `simulate_pilot_mapped_adaptive_frame.m` to perform a frame-level decision using the first orthogonal pilot occasion, then transmit the selected payload on the sample-varying time-domain MIMO-OFDM channel.
- Added `predict_pilot_mapped_mcs.m` to score Spatial Multiplexing / Alamouti and QPSK / 16-QAM candidates from estimated post-detection SINR, then maximize net spectral efficiency subject to a predicted BER target.
- Added `alamouti_decode_timevarying.m`, a two-slot linear MMSE detector that uses separate channel estimates at the two payload times and gain-corrects the detected symbols.
- Added `pilot_mapped_adaptive_validation_run.m`, with paired comparisons against fixed SM-16QAM and fixed Alamouti-QPSK baselines.
- Added `pilot_mapped_adaptive_results_summary.m` and `plot_pilot_mapped_adaptive_results.m`.
- Added V15 mathematical model, scope and limitations documentation.
- Added tests 99–107 for time-varying Alamouti decoding, static-channel compatibility, MCS prediction and selection, end-to-end reproducibility, spectral-efficiency accounting and input validation.

## Validation status

The source and test additions are authored. The current build environment does not contain MATLAB or GNU Octave, so the new MATLAB test suite and V15 Monte-Carlo runner have **not** been executed here. Independent NumPy checks confirmed the two-slot time-varying Alamouti model recovers noiseless payload symbols to numerical precision and reduces to the standard static-channel Alamouti combiner. MATLAB execution is still required to confirm the 107-test suite and 108-row validation run.
