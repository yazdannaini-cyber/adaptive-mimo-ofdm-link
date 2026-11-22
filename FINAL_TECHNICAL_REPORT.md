# Final Engineering Report - Adaptive 2x2 MIMO-OFDM (V15.3)

**Release:** V15.3 joint per-burst adaptive MIMO / modulation / FEC  
**Validation run:** 30 frames per condition, 36 channel conditions  
**Implementation:** MATLAB  
**Result status:** User-executed MATLAB run; aggregate console results archived in this release

## Executive summary

This project implements a pilot-mapped 2x2 MIMO-OFDM link simulator and an online per-burst controller that jointly selects spatial transmission mode, modulation order and forward-error-correction profile. The final V15.3 controller evaluates eight action combinations: Spatial Multiplexing or Alamouti STBC, QPSK or 16-QAM, and uncoded-plus-CRC or rate-1/2 convolutional FEC-plus-CRC. It uses causal channel predictions and CRC ACK/NACK feedback from previously completed bursts. When an action cannot be certified against the configured packet-error target, the controller uses an explicitly marked maximum-posterior-goodput fallback.

In the 30-frame confirmation run, joint adaptation achieved aggregate CRC-accepted goodput of **0.66 bit/s/Hz in Rayleigh** and **1.13 bit/s/Hz in Rician**. The strongest retained reference candidate by CRC-goodput was adaptive convolutional FEC, at 0.58 and 0.95 bit/s/Hz respectively. The joint controller therefore improved the aggregate metric by approximately **13.8% in Rayleigh** and **18.9% in Rician**.

These are promising prototype results, not a claim that the 10% PER target is met in all conditions. The measured joint-controller PER was 34.01% in Rayleigh and 12.35% in Rician, while only 7.3% and 30.5% of data decisions, respectively, were classified as PER-target-feasible by the controller's one-sided posterior bound. Worst-condition uncertainty remained broad.

## 1. Project objectives and scope

The objective was to build an end-to-end communications engineering project that connects waveform-level MIMO-OFDM simulation to channel estimation, receiver detection, adaptive link selection, packet coding, CRC feedback and statistical performance reporting.

The final validation grid contained 36 conditions:

- Channels: Rayleigh and Rician.
- Spatial profiles: IID and joint-high spatial correlation (rhoTx = rhoRx = 0.9 for the joint-high profile).
- Maximum Doppler: 0, 120 and 300 Hz.
- Pilot/data SNR: 8, 14 and 20 dB.
- OFDM frame: 24 OFDM symbols, with pilot spacing of 8 symbols.
- Frames per condition: 30; frame-cluster bootstrap setting: 200 resamples.

## 2. System architecture

The simulation chain includes:

1. Payload generation and packet formatting.
2. CRC-16 attachment and optional terminated rate-1/2 convolutional encoding.
3. Gray-coded QPSK or 16-QAM mapping.
4. MIMO pilot/data resource mapping and OFDM modulation with cyclic prefix.
5. Time-varying multipath 2x2 MIMO channel with Doppler and spatial correlation.
6. CP removal, FFT, pilot extraction and causal Jakes-LMMSE channel tracking.
7. Per-subcarrier MIMO equalization or time-varying Alamouti detection.
8. Max-log soft-bit reliability generation, soft-decision Viterbi decoding and CRC verification.
9. Per-burst controller feedback, logging and performance aggregation.

The joint controller has eight possible actions. It estimates a packet-error prior from the current link prediction, then updates per-action success/failure counts only after a burst has completed and its CRC result is available. The current burst's CRC result cannot influence its own action selection.

## 3. Controller decision rule

For each candidate action, an initial packet-error proxy is formed from the predicted payload BER under an independent-bit approximation. The proxy is combined with a prior weight and that action's earlier CRC ACK/NACK history. A one-sided 95% Wilson-style upper bound is used to classify an action as target-feasible. Among feasible actions, the controller selects the maximum posterior expected CRC-goodput. If no action passes the bound, it selects the maximum expected-goodput fallback and retains an explicit infeasibility flag.

The initial FEC coding-gain curve is a link-adaptation heuristic. It is not a physical guarantee or a substitute for measured coded-packet outcomes. CRC feedback is what updates the controller's per-action history.

## 4. Main results

### 4.1 Aggregate joint-controller performance

| Channel | Payload BER | Coded BER | PER | FER | Net SE (bit/s/Hz) | Bit goodput | CRC goodput |
|---|---:|---:|---:|---:|---:|---:|---:|
| Rayleigh | 0.8992% | 2.025% | 34.01% | 57.22% | 0.98 | 0.97 | **0.66** |
| Rician | 0.05787% | 0.1294% | 12.35% | 25.93% | 1.33 | 1.33 | **1.13** |

### 4.2 CRC-goodput against retained reference candidates

| Channel | Best reference | Best-reference CRC goodput | Joint adaptive | Relative change |
|---|---|---:|---:|---:|
| Rayleigh | Adaptive, conv-crc-1/2 | 0.58 | **0.66** | **+13.8%** |
| Rician | Adaptive, conv-crc-1/2 | 0.95 | **1.13** | **+18.9%** |

The reference set also includes fixed Spatial Multiplexing with 16-QAM and fixed Alamouti with QPSK, each evaluated with and without convolutional FEC. The joint controller is not uniformly superior on every metric: its CRC-goodput is higher, but its aggregate PER is higher than the strongest reliability-oriented reference, especially in Rayleigh.

### 4.3 Per-SNR action behavior

| Channel | SNR | PER | CRC goodput | SM selected | 16-QAM selected | FEC selected | PER feasible |
|---|---:|---:|---:|---:|---:|---:|---:|
| Rayleigh | 8 dB | 41.67% | 0.34 | 1% | 0% | 96% | 7% |
| Rayleigh | 14 dB | 27.22% | 0.67 | 0% | 1% | 46% | 11% |
| Rayleigh | 20 dB | 33.15% | 0.98 | 1% | 34% | 25% | 4% |
| Rician | 8 dB | 13.70% | 0.62 | 0% | 0% | 66% | 32% |
| Rician | 14 dB | 9.63% | 1.11 | 5% | 1% | 0% | 37% |
| Rician | 20 dB | 13.70% | 1.67 | 19% | 46% | 1% | 22% |

The action mix is qualitatively sensible: Rayleigh predominantly selects Alamouti and frequently invokes FEC, while the stronger Rician channel allows more spatial multiplexing and 16-QAM as SNR rises. The controller does not use FEC uniformly; it trades coding overhead against observed packet reliability.

### 4.4 Ten-frame and thirty-frame run consistency

| Channel | 10 frames/condition | 30 frames/condition |
|---|---:|---:|
| Rayleigh joint CRC goodput | 0.683 | 0.66 |
| Rician joint CRC goodput | 1.154 | 1.13 |

The aggregate values are close across the two runs. This is a useful stability check, but it is not a formal confidence interval for the difference and does not replace larger Monte Carlo studies.

## 5. Statistical interpretation and limitations

- The final confirmation run used 30 frames per condition. This is stronger than a 10-frame screening run but remains modest for low-PER claims.
- The worst-condition one-sided upper bounds remained high: payload-BER upper bound 0.2136 and PER upper bound 1.0 in Rayleigh; payload-BER upper bound 0.00620 and PER upper bound 0.5667 in Rician.
- Only 7.3% of Rayleigh and 30.5% of Rician joint-controller data decisions were certified as PER-target-feasible by the controller's posterior decision bound. Therefore, do not claim that the configured 10% PER target is guaranteed or globally achieved.
- The packet-error prior assumes independent bit errors within a packet. Burst fading and residual ICI can violate this approximation.
- The controller's initial convolutional-FEC coding-gain curve is heuristic. Measured CRC and packet outcomes, rather than the heuristic prior, are the acceptance evidence.
- The one-sided Wilson-style bound is formed from a fractional prior-plus-count representation. It is a controller decision statistic, not a purely frequentist confidence interval from raw Bernoulli counts.
- The 30-frame summary CSVs in this release were transcribed from the user's MATLAB console output. Raw per-frame CSV/MAT files were not uploaded and are not included.
- MATLAB release/version and exact elapsed runtime were not captured in the supplied output.

## 6. Reproducibility

From the extracted `MIMO_OFDM_Project_v15_3` folder in MATLAB:

```matlab
clear functions
addpath(genpath(pwd));
run_tests_v15_3
```

The project test source contains 160 tests. The project owner reported that the V15.3 test run completed successfully.

A shorter validation run matching the supplied 30-frame confirmation setting is:

```matlab
addpath('examples');
opts = struct('numFrames',30,'numBootstrap',200);
results = pilot_mapped_adaptive_validation_run_v15_3(opts);
summary = pilot_mapped_adaptive_results_summary_v15_3(results);
```

The default experiment uses 300 frames per condition and 2000 bootstrap resamples:

```matlab
results = pilot_mapped_adaptive_validation_run_v15_3;
summary = pilot_mapped_adaptive_results_summary_v15_3(results);
```

## 7. Conclusion

V15.3 demonstrates a complete pilot-mapped MIMO-OFDM simulation with causal, per-burst joint adaptation across MIMO mode, modulation and FEC, supported by soft-decision Viterbi decoding and CRC feedback. In the supplied 30-frame validation, it improved aggregate CRC-accepted goodput over the strongest retained references by 13.8% in Rayleigh and 18.9% in Rician. The project is suitable as an engineering portfolio prototype and interview discussion case. Claims about stringent PER guarantees or publication-grade statistical confidence require larger experiments and tighter worst-condition uncertainty bounds.
