# V15.3 — Joint adaptive MIMO, modulation and FEC

V15.3 adds an online per-burst link controller that jointly chooses:

- Spatial Multiplexing or Alamouti STBC;
- QPSK or 16-QAM;
- uncoded + CRC-16 or rate-1/2 convolutional FEC + CRC-16.

The selector maximizes posterior expected CRC-goodput among actions whose one-sided 95% PER decision bound satisfies the configured target (default 10%). If no action is feasible, it uses a maximum-goodput fallback and records the infeasibility. CRC ACK/NACK feedback is strictly causal and is retained across frames within each channel condition.

The six V15.2.1 references remain in the output; the seventh candidate is `joint-adaptive`.

## Run tests

```matlab
clear functions
addpath(genpath(pwd));
run_tests_v15_3
```

Expected: 160 tests.

## Run the default experiment

```matlab
addpath('examples');
results = pilot_mapped_adaptive_validation_run_v15_3;
summary = pilot_mapped_adaptive_results_summary_v15_3(results);
```

The default run uses 300 frames per condition, 36 conditions, and 2000 frame-cluster bootstrap resamples. Expected: 252 summary rows, 75,600 candidate frame-log rows, 32,400 joint-controller burst-decision rows and 32,400 reference burst-prediction rows. Files are written to `pilot_mapped_results_v15_3/`.

A short functional check is available with `struct('numFrames',30,'numBootstrap',200)`, but it is not intended for low-PER statistical claims.

## Caveats

The link-model PER and FEC coding gain are initialization heuristics; measured CRC results update each action's history. The Wilson-style upper bound uses a fractional prior-plus-feedback count. It is a decision statistic, not a hard guarantee. Inspect measured BER/PER and frame-cluster bootstrap intervals before making performance claims.

See `V15_3_MODEL.md`, `V15_3_CHANGELOG.md`, `V15_3_VALIDATION_STATUS.md` and `V15_3_INDEPENDENT_NUMERIC_CHECK.md`.
