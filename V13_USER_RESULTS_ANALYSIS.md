# V13 User-Run Results — Analysis

The user confirmed the V13 tests were successful and supplied the 144-row
pilot-mapped OFDM summary. The default study used 8 frames per condition,
16-QAM payload, 15 dB pilot SNR, two channel models, two spatial profiles,
three Doppler values and three pilot spacings.

## Overall channel/estimator summary

| Channel | Hold-last NMSE | Linear NMSE | Jakes-LMMSE NMSE | Perfect diagonal-CSI BER |
|---|---:|---:|---:|---:|
| Rayleigh | 0.1615 | 0.06777 | 0.1122 | 0.1304 |
| Rician | 0.07396 | 0.04240 | 0.05062 | 0.01355 |

The displayed values are arithmetic means across the tested operating
conditions, not bit-count-pooled statistics.

## Findings

- Jakes-LMMSE reduces aggregate NMSE relative to the causal hold-last method
  by about 30.5% in Rayleigh and 31.6% in Rician, using the displayed mean
  NMSE values.
- Linear interpolation has lower average NMSE than Jakes-LMMSE in this
  waveform-level test. Linear uses future pilot observations and is therefore
  an offline/noncausal reference, not a like-for-like online receiver.
- The Jakes tracker is given thermal pilot-noise variance only. Pilot
  observations are contaminated by both thermal noise and within-symbol ICI;
  the latter is measured in V13 but is not included in the tracker noise
  covariance. This model mismatch is a plausible contributor to the gap
  between Jakes-LMMSE and the noncausal linear reference.
- Mean ICI ratios were 0.000258 for Rayleigh and 0.000426 for Rician. In the
  high-correlation profile at 300 Hz, the reported values were 0.000635 and
  0.00107, respectively.
- All methods reported FER=1. With long uncoded frames, even modest BER makes
  it overwhelmingly likely that every frame contains at least one bit error.
  FER=1 is therefore not evidence that the estimators are identical. FEC/CRC
  are not present yet.
- The `perfect-csi` method knows the exact diagonal channel response but still
  uses one-tap-per-subcarrier MMSE detection. It does not know/cancel the full
  off-diagonal ICI coupling matrix, so it is not a full ICI oracle.
- The reported net spectral efficiency is determined by pilot overhead and CP
  and is the same for all estimator methods at a fixed pilot spacing. The
  spacing results show the expected trade-off: spacing 4 yields 3.20 bit/s/Hz
  at 50% pilot overhead; spacing 8 yields 4.80 bit/s/Hz at 25% overhead; and
  spacing 16 yields 5.33 bit/s/Hz at 16.7% overhead.

## Why V14 adds a full ICI oracle

The diagonal perfect-CSI BER (0.1304 Rayleigh and 0.01355 Rician when averaged
across the V13 conditions) shows that the remaining receiver error is not
only a channel-estimation issue. V14 adds an oracle full frequency-coupling
matrix and block MMSE equalizer. The comparison against diagonal perfect CSI
will isolate whether unmodelled ICI is a material contributor or whether
spatial multiplexing conditioning and thermal noise dominate.

The full-ICI method uses true channel taps and is a diagnostic upper
reference, not a deployable pilot-estimated receiver.
