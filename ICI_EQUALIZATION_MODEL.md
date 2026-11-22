# Full ICI Matrix and MMSE Equalization — V14

## Frequency-coupling derivation

For an OFDM symbol with useful time index `m`, the cyclically extended input
sample on transmit antenna `t` is

\[
x_t[(m-l)\bmod N] = \frac{1}{N}\sum_{q=0}^{N-1}
X_t[q]e^{j2\pi q(m-l)/N}.
\]

For sample-varying FIR taps, the useful received sample is

\[
y_r[m]=\sum_t\sum_l h_{r,t,l}[n_m]x_t[(m-l)\bmod N].
\]

Taking the N-point FFT gives

\[
Y_r[k]=\sum_t\sum_q H_{kq}^{(r,t)}X_t[q],
\]

where

\[
H_{kq}^{(r,t)}=\frac{1}{N}\sum_l e^{-j2\pi ql/N}
\sum_m h_{r,t,l}[n_m]e^{-j2\pi(k-q)m/N}.
\]

The inner summation is an N-point FFT of the tap trajectory during the
useful interval, evaluated at index `(k-q) mod N`. This allows the matrix to
be assembled efficiently with vectorized FFT indexing.

## Static-channel check

If each tap is constant during the OFDM symbol, the FFT of its time sequence
is zero at all non-DC Doppler bins. Therefore `H(k,q)=0` for `k != q`, and

\[
H_{kk}^{(r,t)}=\sum_l h_{r,t,l}e^{-j2\pi kl/N},
\]

which is the conventional per-subcarrier MIMO channel matrix.

## MMSE receiver

With unit-energy independent QAM symbols and white FFT-domain noise variance
`σ²`, the full linear MMSE filter is

\[
W=(H^H H+\sigma^2 I)^{-1}H^H.
\]

A diagonal gain correction is applied to the filtered estimates to remove the
linear-MMSE shrinkage bias before nearest-constellation hard decisions.

## Interpretation

The full matrix assumes perfect knowledge of the sample-varying channel taps.
It is used to quantify the performance gap between a diagonal perfect-CSI
receiver and an ICI-aware receiver with oracle knowledge. It is not an estimate
obtained from the V13 pilot observations.
