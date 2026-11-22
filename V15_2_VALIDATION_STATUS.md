# V15.2 Validation Status

## Static and independent checks completed

- CRC algorithm checked against the CRC-16/CCITT-FALSE standard ASCII vector
  `123456789` -> `0x29B1` using an independent implementation.
- Convolutional encoder/Viterbi trellis checked independently for noiseless
  terminated recovery and deterministic single-bit error correction.
- The vectorized time-varying FIR expression was numerically compared against
  the original scalar sample/tap/antenna loop; maximum observed difference in
  the independent numeric check was below `2e-15`.
- Source-level checks confirm the runner produces six strategy/coding profiles:
  three transmission strategies x two packet-coding modes.

## MATLAB validation still required

This build environment does not have MATLAB or GNU Octave, so the 134-test
suite has not been executed here. Do not treat the expected test count as a
passed result until it is confirmed in MATLAB.

From the extracted V15.2 project root:

```matlab
clear functions
addpath(genpath(pwd));
which run_tests_v15_2 -all
run_tests_v15_2
```

Expected terminal footer: `All 134 tests passed successfully.`

## Full statistical run

Default: 100 frames per condition. There are 36 conditions, 3 transmission
strategies, and 2 coding profiles. Expected outputs:

- 216 summary rows
- 21,600 frame-log rows
- 10,800 burst-decision rows

Run:

```matlab
addpath('examples');
results = pilot_mapped_adaptive_validation_run_v15_2;
summary = pilot_mapped_adaptive_results_summary_v15_2(results);
```

For a diagnostic smoke run only, use e.g. `struct('numFrames',5)`. Runs below
30 frames are warned and should not be used as confidence-grade evidence.
