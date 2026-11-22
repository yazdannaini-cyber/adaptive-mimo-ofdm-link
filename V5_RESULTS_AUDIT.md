# V5 comparison result audit

The user-provided MATLAB output reports 24/24 tests passed and a 90-row mode-comparison sweep: three channel references, three QAM orders, five Es/N0 points, two transmission modes, and 200 frames per configuration.

## Main observations

- Rayleigh QPSK at 16 dB: spatial multiplexing BER 0.00958 and FER 0.94; Alamouti had zero observed bit/frame errors in 200 frames. The observed zero is not proof of zero true BER; the FER Wilson upper bound was about 0.019.
- Rayleigh 64-QAM at 24 dB: spatial multiplexing BER 0.0264 and FER 1.0; Alamouti BER 4.07e-5 and FER 0.075. Net spectral efficiency was 7.68 bit/s/Hz for spatial multiplexing and 3.84 bit/s/Hz for Alamouti.
- Rician 16-QAM at 16 dB: spatial multiplexing BER 0.00443 and FER 0.925; Alamouti BER 1.22e-5 and FER 0.02.

## Fairness and scope audit from source

- Both candidate modes use the same generated channel taps and the same pilot-based channel estimate within a Monte-Carlo frame.
- The per-antenna QAM symbol energy is equal. Both modes use two active transmit antennas in each data OFDM symbol; spatial multiplexing carries two independent QAM symbols per channel use, while Alamouti rate-1 STBC carries one. This is a rate-versus-reliability comparison, not an equal-spectral-efficiency comparison.
- Payload bits and receiver noise are independently generated for the two modes. This does not bias the expected comparison, but paired common-random-number noise could reduce Monte-Carlo variance in a future refinement.
- The Alamouti orthogonal combiner assumes the channel remains constant over each adjacent OFDM-symbol pair. The simulator is block-static over the frame; Doppler is not included.
- `Goodput_bpsHz` is bit-level correctly received throughput. It is not CRC-accepted packet goodput. FER counts frames containing at least one bit error.
- The AWGN channel is a deterministic identity-MIMO reference plus noise, not a random fading channel.

## Interpretation

The results demonstrate the expected trade-off: Alamouti substantially improves uncoded error performance under frequency-selective Rayleigh/Rician fading, while two-stream spatial multiplexing offers twice the raw payload spectral efficiency at the same QAM order. The finite 200-frame sweep is useful for a mode comparison but should not be used to claim zero BER; inspect confidence bounds and run longer trials for final low-error-rate claims.
