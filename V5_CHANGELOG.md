# V5 changelog

## Added

- `alamouti_encode.m`: two-transmit-antenna, rate-1 Alamouti OFDM encoder.
- `alamouti_decode.m`: coherent orthogonal combiner supporting multiple receive antennas and per-subcarrier effective-noise variance.
- `run_mimo_mode_comparison.m`: matched-channel Monte-Carlo comparison of 2x2 spatial multiplexing + MMSE detection versus Alamouti STBC.
- `mimo_mode_comparison_run.m`: executable comparison sweep for AWGN, Rayleigh, and Rician channels.
- `plot_mimo_mode_comparison.m`: BER, FER, and spectral-efficiency plot generation.
- Five new tests for Alamouti encoding, noiseless multipath recovery, noisy combining, and input validation. Total test count: 24.
- Frame-level Wilson 95% confidence intervals for FER in the comparison table.
- `MIMO_MODE_COMPARISON_NOTES.md`: mathematical model, metric definitions, run instructions, and limitations.

## Methodology notes

- Both modes use identical channel realizations within each Monte-Carlo frame and identical channel-estimation pilots.
- Both modes use the same per-antenna QAM-symbol energy and Es/N0 convention.
- Spatial multiplexing carries two QAM streams per channel use; Alamouti is rate one and carries one QAM symbol per channel use for two transmit antennas.
- `FER` explicitly means the fraction of data frames containing at least one bit error. No FEC, CRC, or ARQ is included.

## Validation status

The existing V4 baseline was run by the user and passed 19/19 tests. V5 adds five MATLAB tests (24 total) and new simulation code. MATLAB/Octave is not available in the build environment, so V5's new tests and Monte-Carlo study must be executed in the user's MATLAB installation before the new results are treated as validated.
