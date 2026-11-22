# Joint MIMO-mode and modulation adaptation (V6)

V6 evaluates two transmission choices on each quasi-static MIMO-OFDM frame:

- **2x2 spatial multiplexing + linear MMSE detection**: two independent QAM streams per OFDM data symbol.
- **2x2 Alamouti STBC + coherent combining**: one QAM information symbol per channel use, with transmit diversity.

The receiver uses orthogonal pilots to estimate the channel, computes posterior channel-estimation uncertainty, and predicts uncoded Gray-coded square-QAM BER from the post-detection SINR for every candidate mode/modulation pair.

## Decision rule

1. Enumerate all combinations of transmission mode and QAM order.
2. Keep candidates whose predicted uncoded BER is no greater than the configured target.
3. Among feasible candidates, choose the maximum net spectral efficiency after pilot and cyclic-prefix overhead.
4. If no candidate is feasible, choose the minimum predicted-BER candidate and flag the frame as target-infeasible.

This is a prediction-based link adaptation policy, not a guarantee. The measured BER and both Wilson and frame-cluster confidence intervals are reported separately.

## Run

From the project root in MATLAB:

```matlab
run_tests
addpath('examples');
results = adaptive_mimo_mode_run;
```

Default settings: Rayleigh and Rician channels, Es/N0 = 12:4:32 dB, QPSK/16-QAM/64-QAM, 1,000 frames per channel/SNR point, 2,000 frame-cluster bootstrap resamples, and target BER = 1e-3.

Results are saved in `adaptive_mimo_results/` as CSV/MAT and plots.

## Important interpretation

- The Alamouti rate is one QAM symbol per channel use; spatial multiplexing transmits two. With equal QAM order, Alamouti net spectral efficiency is half that of 2x2 spatial multiplexing under the same pilot and CP overhead.
- Both candidate modes are scored using the same channel realization and the same pilot-based channel estimate within a frame. Only the selected mode transmits payload in the adaptive experiment.
- Payload bits and receiver noise are generated for the selected mode; there is no simultaneous transmission of both candidate modes.
- The channel is assumed constant over the complete frame and over each adjacent Alamouti symbol pair. Doppler, CFO, synchronization error, FEC, CRC, and ARQ are not included.
- `BitGoodput_bpsHz` counts correctly received bits, even when other bits in the same frame are wrong. It is not CRC-accepted packet goodput.
- This is a baseband simulation, not an over-the-air measurement.
