# Adaptive 2x2 MIMO-OFDM Link Simulation (MATLAB)

**Final implemented release: V15.3 - Joint per-burst MIMO, modulation and FEC adaptation**

This repository contains an end-to-end pilot-mapped 2x2 MIMO-OFDM simulation. Its joint link controller chooses one transmission action per data burst by combining causal channel-quality predictions with CRC ACK/NACK feedback from completed bursts.

## What is implemented

- 2x2 MIMO-OFDM waveform simulation with cyclic prefix and pilot/data mapping.
- Rayleigh and Rician fading, multipath PDP, Doppler/time variation and spatial correlation.
- LS/LMMSE channel estimation, joint spatial-frequency estimation and causal Jakes-LMMSE tracking.
- Spatial Multiplexing and Alamouti STBC receivers, including time-varying Alamouti detection.
- Gray-coded QPSK and 16-QAM, max-log soft LLRs, rate-1/2 convolutional coding, soft-decision Viterbi and CRC-16.
- Full-band ICI matrix construction and an oracle full-ICI detector benchmark.
- Per-burst joint selection among 8 actions: SM/Alamouti x QPSK/16-QAM x uncoded/FEC.
- Causal per-action ACK/NACK statistics, PER-target feasibility checks, explicit fallback behavior and packet-level metrics.
- BER, coded BER, PER, FER, net spectral efficiency, bit goodput, CRC-accepted goodput, Wilson-style bounds and frame-cluster bootstrap summaries.

## Final confirmation result (30 frames per condition)

The user-executed V15.3 confirmation run covered 36 conditions: 2 channels x 2 spatial profiles x 3 Doppler points x 3 pilot/data SNR points. It reported 252 summary rows, 7,560 frame-log rows, 3,240 joint burst decisions and 3,240 reference burst records.

| Channel | Joint payload BER | Joint PER | Joint CRC goodput | Best reference CRC goodput | Gain |
|---|---:|---:|---:|---:|---:|
| Rayleigh | 0.8992% | 34.01% | **0.66 bit/s/Hz** | 0.58 bit/s/Hz | **+13.8%** |
| Rician | 0.05787% | 12.35% | **1.13 bit/s/Hz** | 0.95 bit/s/Hz | **+18.9%** |

The best reference in both channels was the retained adaptive convolutional-FEC candidate. The joint controller improves CRC-accepted goodput but does not provide a universal PER guarantee: only 7.3% of Rayleigh and 30.5% of Rician data decisions were classified as PER-target-feasible by its one-sided posterior decision bound in this run. Worst-condition bounds remain broad.

Curated aggregate CSV files are in `results/v15_3_final_summary/`. They were transcribed from the supplied MATLAB console output; raw per-frame CSV/MAT files were not provided.

## Quick start

Open this folder as the MATLAB working directory:

```matlab
clear functions
addpath(genpath(pwd));
run_tests_v15_3
```

The test source contains 160 tests. The project owner reported a successful V15.3 test run.

Shorter confirmation experiment (30 frames/condition, 200 bootstrap resamples):

```matlab
addpath('examples');
opts = struct('numFrames',30,'numBootstrap',200);
results = pilot_mapped_adaptive_validation_run_v15_3(opts);
summary = pilot_mapped_adaptive_results_summary_v15_3(results);
```

Default experiment (300 frames/condition, 2000 bootstrap resamples):

```matlab
results = pilot_mapped_adaptive_validation_run_v15_3;
summary = pilot_mapped_adaptive_results_summary_v15_3(results);
```

Generated files are written under `pilot_mapped_results_v15_3/` and are excluded from Git by `.gitignore`.

## Important model notes

- The controller's initial packet-error model assumes independent bit errors within a packet; burst fading and residual ICI can violate that approximation.
- Its FEC coding-gain curve is an initialization heuristic, not a guaranteed coding gain. Measured CRC feedback updates the per-action history.
- Its one-sided Wilson-style bound uses a fractional prior-plus-count representation. It is a decision statistic, not a purely frequentist interval from raw Bernoulli samples.
- When no action meets the PER bound, the controller uses a maximum posterior expected CRC-goodput fallback and marks the selected action infeasible.
- The 30-frame result is a confirmation run, not publication-grade validation. Increase the frame count and inspect worst-condition frame-cluster intervals before making stringent reliability claims.
- MATLAB release/version and runtime were not captured in the supplied run log.

## Documentation

- `FINAL_TECHNICAL_REPORT.md` - complete final engineering report and measured results.
- `docs/V15_3_FINAL_REPORT.docx` / `.pdf` - formatted report.
- `docs/FINAL_INTERVIEW_BRIEF_FA.md` - Persian technical interview preparation notes.
- `docs/CV_PROJECT_DESCRIPTION_EN_FA.md` - concise CV/project descriptions in English and Persian.
- `V15_3_MODEL.md` - control model and mathematical details.
- `V15_3_VALIDATION_STATUS.md` - release-specific validation notes.
- `GITHUB_RELEASE_CHECKLIST.md` - publishing checklist and evidence status.

## Repository structure

- `config/` - default configuration and validation.
- `src/` - channel, PHY, estimation, detection, coding, simulation and metrics.
- `examples/` - executable study runners and summary functions.
- `tests/` - regression and feature tests.
- `results/v15_3_final_summary/` - curated aggregate results from the supplied 30-frame run.
- `docs/` - final report, interview brief, CV text and figures.

## Scope

This is a communications engineering portfolio/research prototype for simulation and comparison. It is not a certified modem implementation or a guarantee of field performance.
