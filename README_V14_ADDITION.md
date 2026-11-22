# V14 — Full ICI-Matrix MMSE Equalization (Oracle Reference)

## Motivation from the V13 pilot-mapped results

The V13 waveform-level study reported substantial uncoded BER even for the
`perfect-csi` method. That method knows the exact per-subcarrier diagonal
channel response, but it still uses a one-tap-per-subcarrier detector. It does
not know or cancel the off-diagonal coupling between subcarriers created by a
channel that changes within an OFDM symbol.

V14 adds a full frequency-coupling matrix and an ICI-aware block MMSE detector.
The new method is called `perfect-csi-ici`. It is an oracle diagnostic reference:
it uses the exact sample-varying channel taps to construct the full matrix. The
current mapped pilot design does not yet estimate this full matrix, so the new
method must not be described as a deployable pilot-estimated receiver.

## Exact channel matrix

For one OFDM symbol, let `m` index useful time samples, `l` the channel tap,
`q` the input subcarrier and `k` the output subcarrier. With a CP covering the
channel memory, the frequency-domain coupling coefficient is

\[
H_{kq}^{(r,t)} = \frac{1}{N}\sum_{m=0}^{N-1}\sum_{l=0}^{L-1}
 h_{r,t,l}[n_m] e^{-j2\pi ql/N}e^{j2\pi(q-k)m/N}.
\]

For a block-static channel, this matrix is diagonal in subcarrier index and
reduces to the usual FFT of the FIR taps. For a time-varying channel,
off-diagonal terms represent ICI.

`build_time_varying_ofdm_ici_matrix` returns a matrix with dimensions
`[N*nRx x N*nTx]`, using MATLAB column-major stream ordering:

```matlab
xVec = reshape(Xfd, [], 1);  % all subcarriers of TX1, then TX2, ...
yVec = reshape(Yfd, [], 1);  % all subcarriers of RX1, then RX2, ...
yVec = Hblock*xVec;
```

The matrix is calculated from the true tap sequence at the useful FFT
sampling instants. It includes all same-subcarrier and cross-subcarrier terms.

## Full-band MMSE detector

`mimo_ici_equalize` computes the block linear MMSE filter

\[
W=(H^H H+\sigma_n^2 I)^{-1}H^H,
\]

then applies a per-stream diagonal-gain correction before QAM hard decisions.
The receiver operates jointly across all subcarriers and both transmit streams
for one OFDM symbol. It is more computationally expensive than the existing
per-tone equalizer.

## V14 receiver references

The pilot-mapped experiment now compares five methods:

1. `hold-last` — causal per-link pilot hold.
2. `linear` — noncausal interpolation/extrapolation reference.
3. `jakes-lmmse` — causal temporal predictor based on the Jakes covariance.
4. `perfect-csi` — exact diagonal channel response with one-tap MMSE detection.
5. `perfect-csi-ici` — exact full frequency-coupling matrix with block MMSE detection.

The key diagnostic is the BER difference between `perfect-csi` and
`perfect-csi-ici`. A large reduction indicates that the omitted off-diagonal
ICI coupling is a material part of the V13 receiver loss. A small reduction
indicates that spatial multiplexing conditioning, thermal noise, or other
receiver assumptions dominate instead.

## Run

```matlab
run_tests
addpath('examples');
results = pilot_mapped_ofdm_validation_run;
summary = pilot_mapped_results_summary(results);
```

The default experiment retains the V13 conditions and adds one receiver
reference, producing 180 result rows (36 channel/profile/Doppler/pilot-spacing
conditions x 5 methods). The full-matrix oracle is evaluated only for data
symbols and does not change pilot mapping or the three practical tracking
methods.

## Important limitations

- `perfect-csi-ici` uses the true channel tap sequence. It is an oracle upper
  reference, not an estimated receiver.
- The actual pilot pattern supplies a full-band LS diagonal observation for
  each transmit antenna at its own pilot time. It does not provide enough
  independent pilot equations to identify the complete `N x N` ICI coupling
  matrix for every MIMO link. A new orthogonal time/frequency training design
  would be needed for practical full-matrix estimation.
- The Jakes tracker still uses the thermal LS pilot-noise variance; the
  measured pilot contamination including ICI is reported but not yet fed into
  its covariance model.
- No FEC, CRC, ARQ, CFO, phase noise, timing error or synchronization loop is
  included.
- The screening experiment retains 8 frames per condition. Increase frame
  count before making narrow statistical or deployment claims.
