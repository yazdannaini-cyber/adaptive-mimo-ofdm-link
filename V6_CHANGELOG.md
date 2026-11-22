# V6 changelog — joint MIMO-mode and QAM adaptation

## Added

- `select_mimo_link_mode.m`: joint candidate selector for two-stream spatial multiplexing and rate-1 Alamouti STBC across configured QAM orders.
- `mimo_mode_spectral_efficiency.m`: centralized accounting for stream count, pilot overhead, and cyclic-prefix overhead.
- `run_adaptive_mimo_mode.m`: frame-by-frame adaptive link simulation using one pilot-derived channel estimate to score all candidate mode/modulation pairs.
- `adaptive_mimo_mode_run.m`: reproducible 1,000-frame/channel/SNR study for Rayleigh and Rician fading with 2,000 frame-cluster bootstrap resamples.
- `plot_adaptive_mimo_mode_results.m`: measured BER, net spectral efficiency, and selected mode-ratio figures.
- `examples/README_ADAPTIVE_MIMO_MODE.md`: algorithm description, run instructions, and limitations.
- Five additional tests (Tests 25–29), extending the intended deterministic suite from 24 to 29 tests.

## Selection policy

1. Predict uncoded QAM BER for each combination of spatial multiplexing/Alamouti and QPSK/16-QAM/64-QAM.
2. Retain candidates with predicted BER no greater than the configured target.
3. Select the highest net spectral-efficiency candidate among those feasible options.
4. If none is feasible, select the candidate with minimum predicted BER and record a target-infeasible fallback.

## Validation status

The V5 MATLAB run supplied by the user confirms that the 24-test baseline and the fixed V5 comparison completed. The V6 additions have been created and statically reviewed in this build environment. MATLAB/Octave is not installed here, so Tests 25–29 and `adaptive_mimo_mode_run` remain to be executed on the user's MATLAB installation. No V6 performance result is claimed yet.
