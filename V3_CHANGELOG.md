# Version 3 — Clear adaptive-link acceptance metrics

## Why this version exists
The V2 console output printed `target met` for a frame-level prediction condition. That label could be mistaken for the measured aggregate BER meeting the target. V3 separates these concepts and corrects the averaging used to compare predicted BER with measured BER.

## Changes

1. **Explicit adaptive feasibility and fallback reporting**
   - `PredictedTargetMetFrameRatio`: fraction of frames for which the selected mode's predicted BER is at or below the target.
   - `FallbackQPSKFrameRatio`: fraction of frames where no modulation mode met the prediction target and QPSK was used as minimum-service fallback.
   - `ObservedBER_TargetSatisfied`: whether the measured aggregate Monte-Carlo BER is at or below the target.
   - `ObservedBER_CI95UpperBelowTarget`: stronger statistical check; the upper Wilson 95% confidence bound is also below the target.

2. **Payload-bit-weighted predicted BER**
   - Added `PredictedBER_PayloadWeighted`, weighting each frame's predicted BER by that frame's actual number of payload bits. This is the appropriate aggregate prediction to compare against aggregate measured BER when QAM order changes across frames.
   - Retained `MeanPredictedBER_FrameAverage` as a diagnostic only.

3. **Better plot interpretation**
   - Adaptive performance plots now use payload-bit-weighted predicted BER.
   - Added `Adaptive_Feasibility_*.png` to distinguish target-feasible frames from QPSK fallback frames.

4. **Validation**
   - Expanded deterministic test suite from 13 to 15 tests.
   - Added tests for highest-rate target-compliant mode selection, explicit fallback status, payload-bit-weighted BER aggregation, and separation of prediction feasibility from measured BER target compliance.

## Reading the adaptive output correctly

A QPSK fallback still transmits data, but it does **not** mean the target BER is met. For example, if `PredictedTargetMetFrameRatio = 0.60` and `FallbackQPSKFrameRatio = 0.40`, then 60% of simulated frames had at least one predicted target-compliant QAM mode, while 40% required fallback. Whether the overall link actually meets the target must be read from `BER`, `BER_CI95_Upper`, `ObservedBER_TargetSatisfied`, and `ObservedBER_CI95UpperBelowTarget`.

## Revalidation

```matlab
run_tests
addpath('examples');
quick_run
[fixedResults, adaptiveResults] = run_project;
```

MATLAB execution is required to validate the new test cases and regenerate results. Do not claim the new MATLAB tests or V3 simulations have passed until run on MATLAB.
