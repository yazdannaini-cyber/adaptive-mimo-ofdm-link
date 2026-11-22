# V15: Pilot-Driven Adaptive MIMO-OFDM

## Scope

V15 carries the earlier link-adaptation idea into the actual pilot-mapped, sample-varying 2x2 MIMO-OFDM waveform. It jointly selects:

- Spatial Multiplexing (two independent spatial streams), or rate-one Alamouti STBC;
- QPSK or 16-QAM (the selector can also score further square-QAM orders if included).

The first orthogonal TDM pilot occasion is reserved for a frame-level transmission decision. Each transmit antenna has one dedicated full-band pilot symbol. The decision is therefore made from pilot observations only; it does not inspect the payload or the true channel.

## Causal decision and channel tracking

For each transmit channel column, the first pilot observation is processed with the configured Rician LOS prior (zero mean for Rayleigh) and a one-observation Jakes-LMMSE predictor to estimate the channel at the first data OFDM symbol. Candidate post-detection SINR is then calculated for:

1. two-stream MMSE Spatial Multiplexing;
2. generalized two-slot Alamouti detection.

The uncoded square-QAM BER approximation is evaluated on the candidate post-SINR values. Among candidates whose predicted BER does not exceed `targetBER`, the selector chooses the largest net spectral efficiency. If none is predicted feasible, it chooses the candidate with the smallest predicted BER and marks the frame as a fallback.

The selected mode and modulation remain fixed for the remainder of the frame. Later orthogonal pilots feed the causal Jakes-LMMSE tracker for data detection. The estimator uses past pilot observations only; it does not use future pilots to decode an earlier data symbol.

## Time-varying Alamouti receiver

The Alamouti encoder maps two QAM information symbols per subcarrier over two data OFDM slots:

\[
\mathbf{x}_1=[s_1,s_2]^T,\qquad
\mathbf{x}_2=[-s_2^*,s_1^*]^T.
\]

The two data slots need not be adjacent in the OFDM frame; pilot symbols may lie between payload slots. For each pair, the receiver uses its channel estimates at both actual data times. The augmented observation is

\[
\begin{bmatrix}\mathbf{y}_1\\\mathbf{y}_2^*\end{bmatrix}
=
\underbrace{\begin{bmatrix}
\mathbf{h}_1(t_1)&\mathbf{h}_2(t_1)\\
\mathbf{h}_2^*(t_2)&-\mathbf{h}_1^*(t_2)
\end{bmatrix}}_{\mathbf{A}}
\begin{bmatrix}s_1\\s_2\end{bmatrix}+\mathbf{n}.
\]

It applies a two-symbol linear MMSE estimate

\[
\widehat{\mathbf{s}}=(\mathbf{A}^H\mathbf{A}+\sigma^2\mathbf{I})^{-1}
\mathbf{A}^H\mathbf{y},
\]

followed by diagonal gain correction before QAM slicing. Unlike a textbook orthogonal combiner that assumes a constant channel over both slots, this receiver accommodates different channel estimates at `t1` and `t2`. Intra-symbol ICI is not jointly cancelled in this V15 detector; it remains part of the effective disturbance seen by the per-subcarrier receiver.

## Net spectral efficiency

For `D` payload OFDM symbols in a frame of `S` total OFDM symbols, FFT length `N`, CP length `Ncp`, and modulation order `M`:

\[
\eta_{SM}=2\log_2(M)\frac{D}{S}\frac{N}{N+N_{cp}},
\qquad
\eta_{Alamouti}=\log_2(M)\frac{D}{S}\frac{N}{N+N_{cp}}.
\]

Thus both pilot-symbol and cyclic-prefix overhead are included. Alamouti conveys one information stream per payload channel use while Spatial Multiplexing conveys two.

`UncodedGoodput_bpsHz = NetSpectralEfficiency_bpsHz * (1 - measured BER)` is a bit-level approximation. It is not CRC-accepted packet goodput or post-FEC throughput.

## Validation comparisons

The V15 screening runner compares three strategies on the same channel taps, pilot symbols and FFT-domain thermal-noise realization:

- `adaptive`: the pilot-driven selected mode and modulation;
- `fixed-sm-16qam`: a conventional fixed two-stream 16-QAM reference;
- `fixed-alamouti-qpsk`: a robust fixed diversity-mode reference.

Each candidate gets an independently seeded payload, but all candidates share the same channel and noise realization. This is a paired physical-channel comparison, not a bit-for-bit common-payload comparison.

## Limitations

- The selection is frame-level, based on the first pilot occasion. It does not adapt per subcarrier or within a frame.
- The Jakes temporal statistics, maximum Doppler, Rician K-factor and configured LOS mean are assumed known to the receiver model.
- The BER predictor is an AWGN-based uncoded approximation of post-equalization SINR using point channel estimates and the thermal-noise variance. It does not explicitly fold the tracker posterior variance or pilot ICI contamination into the selection metric. It must be compared with measured BER and its confidence intervals; predicted feasibility is not a guarantee.
- The current per-subcarrier detector does not perform full ICI matrix cancellation. ICI from within-symbol channel variation remains.
- The current comparison uses Rayleigh and Rician fading, two spatial profiles, Doppler 0/120/300 Hz, SNR 8/14/20 dB, pilot spacing 8, and eight frames per condition. This is a screening run; increase the frame count substantially before making low-BER or publication-grade claims.
- Alamouti requires an even number of payload OFDM symbols after pilot reservation. The simulator rejects configurations that leave an odd count.
