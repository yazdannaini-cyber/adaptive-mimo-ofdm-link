# Adaptive MIMO-OFDM metric definitions (V3)

| Column | Meaning | What it does **not** mean |
|---|---|---|
| `BER` | Monte-Carlo measured aggregate bit error rate over all payload bits at this SNR | Not a theoretical BER prediction |
| `BER_CI95_Lower`, `BER_CI95_Upper` | 95% Wilson interval for measured aggregate BER | Not a confidence interval for model mismatch |
| `MeanPredictedBER_FrameAverage` | Arithmetic mean of the per-frame predicted BER values | Not weighted for varying frame payload sizes |
| `PredictedBER_PayloadWeighted` | Sum of predicted BER times payload bits divided by total payload bits | Not a guarantee; it remains an analytical predictor |
| `PredictedTargetMetFrameRatio` | Fraction of frames where the selected mode's predicted BER is at/below target | Not the fraction of measured bits meeting target |
| `FallbackQPSKFrameRatio` | Fraction of frames where no candidate mode met predicted target and QPSK fallback was used | Not a dropped-frame/outage ratio; the simulation still sends QPSK payload |
| `ObservedBER_TargetSatisfied` | Measured aggregate BER is at/below configured target | Not sufficient statistical evidence by itself |
| `ObservedBER_CI95UpperBelowTarget` | The upper 95% Wilson bound for measured BER is at/below target | Stronger evidence, still conditional on the simulation model |

## Recommended acceptance language

- **Prediction feasibility:** report `PredictedTargetMetFrameRatio` and `FallbackQPSKFrameRatio`.
- **Measured target check:** report `BER`, `BER_CI95_Upper`, and both measured-target flags.
- **Model calibration:** compare `PredictedBER_PayloadWeighted` with `BER`. A persistent gap suggests predictor mismatch or simulation/model assumptions requiring review.
- **Statistical confidence:** increase the number of Monte-Carlo frames when the confidence interval is too wide to support the target claim.
