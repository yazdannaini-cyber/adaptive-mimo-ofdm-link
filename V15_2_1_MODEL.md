# V15.2.1 Soft-Decision Receiver Model

## QAM bit LLRs
For an equalized symbol `z` and effective complex-noise variance `sigma2`, the max-log LLR for bit position `b` is

`LLR_b = ( min_{s in S_b=1} |z-s|^2 - min_{s in S_b=0} |z-s|^2 ) / sigma2`.

Positive LLR favors bit zero. The QAM constellation is generated through the project's own Gray mapper, so the emitted LLR order is identical to the bit order expected by `qam_demod` and packet serialization.

For each stream and subcarrier, the detector's post-equalization SINR is converted to an approximate variance using `sigma2_eff = 1 / postSINR`. In the linear MMSE and time-varying Alamouti detectors this post-SINR includes modeled residual inter-stream interference plus thermal, channel-estimation, and estimated ICI disturbance terms.

## Soft Viterbi
The terminated rate-1/2 K=3 convolutional code with generator polynomials `(7,5)_8` is unchanged. Viterbi branch costs use the two coded-bit LLRs rather than counting hard bit mismatches. The trellis starts and ends in state zero, matching the encoder's two zero tail bits.

## Packet processing
- `uncoded-crc`: max-log LLR signs are hard-decided and the CRC-16/CCITT-FALSE is verified.
- `conv-crc-1/2`: soft LLRs enter the Viterbi decoder; two termination bits are removed; then CRC is verified.
- CodedBER remains the hard-sign BER before FEC, preserving comparability with V15.2.
- PayloadBER is measured after Viterbi for the coded profile.

## Approximation boundary
The residual disturbance after MMSE/Alamouti equalization is summarized by a per-stream scalar Gaussian variance. Max-log LLRs are therefore approximate when interference is strongly non-Gaussian, correlated across streams/subcarriers, or dominated by unmodeled ICI. This revision adds reliability information but does not claim full joint MAP detection.
