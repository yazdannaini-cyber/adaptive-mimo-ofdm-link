# V13 changelog

## Added

- Explicit 2x2 time-division orthogonal pilot scheduling with non-overlap validation.
- Full-band known QPSK pilot and QAM payload mapping into the OFDM frequency-domain grid.
- Per-transmitter pilot extraction from the received FFT grid using known pilot sequences.
- Per-transmitter temporal tracking at staggered real pilot timestamps.
- Sample-varying waveform propagation through the existing Doppler MIMO FIR channel.
- Per-symbol diagonal frequency-response reference computed from the time-averaged useful-interval taps.
- Measured within-symbol ICI residual and pilot-observation contamination variance.
- End-to-end data-symbol NMSE, uncoded BER, FER, net spectral efficiency and goodput metrics for hold-last, linear, Jakes-LMMSE and perfect-diagonal-CSI references.
- Frame-bootstrap intervals for NMSE and BER.
- CSV/MAT exports and Rayleigh/Rician plots.
- Tests 82–91 for pilot scheduling, grid mapping, LS extraction, static-channel response, dimensions, end-to-end waveform flow, zero-Doppler ICI, Doppler-induced ICI and pilot-overhead accounting.

## Methodological notes

- The V12 abstract pilot-observation benchmark is preserved. V13 adds a separate waveform-level experiment rather than silently changing the meaning of the existing V12 results.
- TDM pilot occasions reserve `nTx` consecutive OFDM symbols; the per-TX temporal trackers use their own staggered pilot indices.
- ICI is measured, not cancelled. The V13 per-subcarrier MMSE detector uses the diagonal channel and thermal-noise variance, so the reported BER includes residual ICI.
- The Jakes tracker is supplied thermal pilot-noise variance; pilot ICI contamination is reported separately. This makes the effect of observation-model mismatch visible.
- `perfect-csi` is a diagonal-channel oracle benchmark, not a full ICI-aware oracle detector.

## Validation status

V13 code and test cases were added and statically reviewed in the build environment. MATLAB/Octave is not available in the build environment, so the new tests and Monte Carlo experiment have **not** yet been executed here. Run `run_tests` in MATLAB; expected suite size is 91 tests.
