# V15.3 — Joint per-burst MIMO, modulation and FEC adaptation

## Objective

V15.3 adds one joint link controller to the V15.2.1 pilot-mapped 2x2 MIMO-OFDM simulator. At every pilot-to-pilot data burst it chooses one of eight actions:

- Spatial Multiplexing or Alamouti STBC;
- QPSK or 16-QAM;
- uncoded + CRC-16 or terminated rate-1/2 convolutional FEC + CRC-16.

The six V15.2.1 adaptive/fixed reference candidates remain in the experiment. A seventh candidate, `joint-adaptive`, is evaluated alongside them.

## Causal control loop

For each burst, the transmitter-side controller receives only the current burst's causal link prediction and CRC outcomes from previously completed bursts. It evaluates every supported action, chooses an action, the simulator sends and decodes that burst, and only then is the CRC ACK/NACK entered into the controller state. The current burst's CRC result cannot affect its own selection.

The controller state is retained across frames within the same channel/profile/Doppler/SNR condition and reset when the condition changes. Each action has its own success/failure counters indexed by MIMO mode, modulation and coding mode.

## Candidate utility and PER constraint

For an action with payload length `L` and predicted payload BER `p_b`, the initial independent-bit PER proxy is

`p_model = 1 - (1 - p_b)^L`.

The proxy is blended with previous CRC feedback using a prior weight `w` (default 4):

`p_post = (w*p_model + failures)/(w + successes + failures)`.

A one-sided 95% Wilson-style upper bound is evaluated from the same fractional pseudo-count representation. An action is considered PER-feasible only when this upper bound is at or below the configured target PER (default 0.10). Among feasible actions, the controller maximizes

`expected CRC-goodput = payload spectral efficiency * (1 - p_post)`.

If no action satisfies the target, it selects the action with the largest posterior expected CRC-goodput and explicitly marks the selection as target-infeasible. No target-compliance guarantee is implied by the fallback.

## FEC link abstraction

For the controller's initial model only, the rate-1/2 convolutional profile applies a bounded effective-SNR gain to the post-equalization SINR before the QAM BER approximation is calculated. Defaults are 2 dB below raw BER 0.05, 1.5 dB from 0.05 to 0.12, and zero gain at or above 0.12. These values are a link-adaptation prior, not a physical analytical bound or a guaranteed coding gain. Actual packet outcomes from soft-decision Viterbi and CRC update the per-action feedback statistics. Measured BER/PER and their confidence intervals remain the acceptance criteria.

## Rate accounting

The controller accounts for payload bits after CRC and FEC overhead, the chosen spatial mode, OFDM cyclic prefix and the already allocated pilot/data schedule. The objective is CRC-accepted payload goodput, not uncoded bit-goodput. The simulation separately reports payload BER, coded BER, packet error rate, frame error rate, gross and net spectral efficiency, bit-goodput and CRC-accepted packet-goodput.

## Default experiment

- Channels: Rayleigh and Rician.
- Spatial profiles: IID and joint-high correlation.
- Doppler: 0, 120 and 300 Hz.
- Pilot/data SNR: 8, 14 and 20 dB.
- Conditions: 36.
- Frames per condition: 300 by default.
- Pilot spacing: 8 OFDM symbols; 24 OFDM symbols per frame.
- Target PER: 0.10.
- Frame-cluster bootstrap resamples: 2000.

Expected default output: 252 summary rows (36 conditions x 7 candidates), 75,600 candidate frame-log rows, 32,400 joint burst-decision rows (36 x 300 x 3 bursts), and 32,400 reference burst-prediction rows.

## Interpretation and limitations

1. The initial PER model assumes independent bit errors within a packet. Burst fading and residual ICI violate that assumption; causal CRC feedback is intended to correct the prior as observations accumulate.
2. The controller's FEC coding-gain curve is a configurable heuristic used only for initial ranking. It must not be presented as a proven convolutional-code performance curve.
3. The one-sided Wilson-style bound is computed on a fractional prior-plus-count representation. It is a conservative decision statistic, not a frequentist confidence interval based on raw Bernoulli observations alone.
4. The joint controller is online and causal, but it is not yet a full reinforcement-learning policy. It uses an interpretable per-action ACK/NACK estimator.
5. Results with fewer than 30 frames per condition are explicitly screening-only. Even 300 frames can be insufficient to tightly bound very low packet-error rates; inspect the frame-cluster intervals and increase the frame count for claims near a stringent target.
6. MATLAB execution is required to validate the MATLAB implementation. The build environment used to package this release did not contain MATLAB or GNU Octave.
