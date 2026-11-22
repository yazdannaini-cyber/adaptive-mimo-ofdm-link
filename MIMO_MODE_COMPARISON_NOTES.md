# MIMO transmission-mode comparison (V5)

## Purpose

V5 adds a controlled comparison between:

1. **2x2 spatial multiplexing** with linear MMSE detection.
2. **2x2 Alamouti space-time block coding (STBC)** with coherent orthogonal combining.

Both modes use the same channel realization for each Monte-Carlo frame, the same two-symbol orthogonal pilot structure, the same channel estimate, the same per-antenna average QAM-symbol energy, and the same per-subcarrier noise variance. Noise samples and payload bits are independently generated for the two transmission modes.

## Signal models

For spatial multiplexing, two independent symbols are transmitted during each OFDM data symbol:

\[
\mathbf y[k,m]=\mathbf H[k]\mathbf x[k,m]+\mathbf n[k,m].
\]

For Alamouti STBC, two information symbols are sent across two adjacent OFDM symbols:

\[
\mathbf X_{\rm A}[k]=
\begin{bmatrix}
 s_1[k] & -s_2^*[k]\\
 s_2[k] & s_1^*[k]
\end{bmatrix}.
\]

The code has rate one: it sends two complex information symbols over two channel uses. With two transmit antennas, spatial multiplexing sends two independent complex symbols per channel use, so its raw spectral efficiency is twice that of Alamouti at the same QAM order. Alamouti trades that multiplexing gain for transmit diversity and simpler orthogonal detection.

For each receive antenna, the decoder combines the two slots as

\[
\hat s_1=\frac{\sum_r h_{1,r}^*y_{1,r}+h_{2,r}y_{2,r}^*}{\sum_r(|h_{1,r}|^2+|h_{2,r}|^2)},
\]

\[
\hat s_2=\frac{\sum_r h_{2,r}^*y_{1,r}-h_{1,r}y_{2,r}^*}{\sum_r(|h_{1,r}|^2+|h_{2,r}|^2)}.
\]

The derivation assumes the channel remains approximately constant over each adjacent OFDM-symbol pair. This is satisfied by the current block-static channel simulator; Doppler/time variation is not yet modeled.

## Metrics and interpretation

- **BER:** aggregate payload bit errors divided by transmitted payload bits.
- **FER:** fraction of simulated data frames with at least one payload bit error; the comparison also reports a Wilson interval based on the number of frames.
- **Wilson CI:** bit-level interval under an approximate independent-bit model.
- **Frame-bootstrap CI:** nonparametric interval obtained by resampling complete frames.
- **Spectral efficiency:** theoretical payload-bit rate normalized by occupied time-frequency resources, including the current pilot and cyclic-prefix overhead.
- **Goodput:** correctly received payload bits normalized by the same simulated resource count. It is a bit-level goodput, not CRC-accepted packet goodput.

The legacy adaptive output column named `BLER` is calculated as `number of frames with >=1 bit error / number of frames`. In the current uncoded simulator this is best interpreted as FER (frame error rate). No FEC decoder, CRC, ARQ, or packet retransmission is included; do not claim post-FEC BLER or packet-delivery reliability from these results.

## Run

From the project root in MATLAB:

```matlab
addpath('examples');
comparisonResults = mimo_mode_comparison_run;
```

Default sweep:

- Channels: AWGN identity reference, Rayleigh, Rician
- Modulations: 4-QAM, 16-QAM, 64-QAM
- Es/N0: 0, 8, 16, 24, 32 dB
- Frames per configuration: 200
- Bootstrap resamples: 1000

This produces a CSV table, MAT file, BER, FER, and spectral-efficiency plots in `comparison_results/`.

## Scope and limitations

- Quasi-static frequency-selective fading per frame.
- Perfect synchronization and no carrier-frequency/timing offset.
- No Doppler within an Alamouti pair.
- No channel coding, CRC, or retransmission.
- Linear MMSE detection for spatial multiplexing; orthogonal Alamouti combining for STBC.
- Same per-antenna symbol energy rather than forcing equal spectral efficiency. Rate and reliability must therefore be interpreted together.
- AWGN is a deterministic identity-MIMO reference, not an additive-noise-only fading channel with a random matrix.

The comparison is a simulation study, not an over-the-air measurement.
