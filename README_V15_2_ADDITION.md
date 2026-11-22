# V15.2 — Statistical validation + packet-level CRC/FEC

V15.2 builds on V15.1.1 and adds 100-frame default validation, burst packet
CRC-16, a rate-1/2 convolutional FEC option with hard-decision Viterbi decoding,
separate payload/coded BER, CRC PER, FER, net payload SE, bit-goodput and
CRC-accepted packet-goodput metrics.

## Run tests

```matlab
clear functions
addpath(genpath(pwd));
which run_tests_v15_2 -all
run_tests_v15_2
```

Expected: 134 tests pass. MATLAB/Octave execution is required to verify this.

## Run full validation

```matlab
addpath('examples');
results = pilot_mapped_adaptive_validation_run_v15_2;
summary = pilot_mapped_adaptive_results_summary_v15_2(results);
```

Default is 100 frames/condition. For a short diagnostic run:

```matlab
results = pilot_mapped_adaptive_validation_run_v15_2(struct('numFrames',10));
summary = pilot_mapped_adaptive_results_summary_v15_2(results);
```

A run below 30 frames is not intended for strong statistical conclusions.

See `V15_2_MODEL.md`, `V15_2_CHANGELOG.md`, and `V15_2_VALIDATION_STATUS.md`
for definitions, model assumptions and validation status.
