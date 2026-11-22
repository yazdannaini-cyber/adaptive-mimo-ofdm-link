# V15.1 Burst-Adaptive Link Adaptation Model

## 1. Pilot-to-data burst structure

The orthogonal TDM pilot schedule reserves one OFDM symbol per transmit antenna
at each pilot occasion. Immediately after the second antenna's pilot, the
following data symbols form one adaptation burst. A new decision is made at
each later pilot occasion. In the current 24-symbol, spacing-8 configuration,
the pilot pairs are at symbols (1,2), (9,10), and (17,18), followed by data
bursts (3:8), (11:16), and (19:24).

For burst `b`, transmitter `t`'s temporal tracker receives only pilot samples
1 through `b`. The decision process therefore excludes future pilots. The
current burst's predicted channel and posterior variance are evaluated at
each data-symbol time in the burst.

## 2. Posterior-error-aware disturbance

For each predicted data symbol `s`, the tracker supplies a posterior diffuse
channel-error variance `v_t(s)` for transmit-link column `t`. The scalar
equivalent receiver-noise variance used by the link predictor is

` sigma_eff^2(s) = sigma_thermal^2 + sigma_ICI^2(s) + sum_t v_t(s) `.

This treats channel-estimation error as independent additive disturbance
across the receive antennas and independent transmit streams. It is a
first-order approximation; it does not preserve the full spatial covariance
of channel-estimation error.

Spatial-Multiplexing post-SINR is evaluated with a per-subcarrier MMSE
receiver using this effective disturbance. For Alamouti, the two-slot
augmented observation uses the predicted channel matrices at the actual first
and second data slots. Its augmented noise covariance is block diagonal, with
a possibly different effective noise variance in each slot.

## 3. Causal pilot-innovation ICI estimate

For pilot occasion `b > 1`, a one-step Jakes-LMMSE prediction is formed at the
current pilot time from earlier pilot observations only. For each transmit
link, the innovation power is compared with the expected thermal LS noise
plus the previous-pilot-conditioned posterior prediction variance. Positive
excess is treated as unmodelled pilot impairment (ICI, residual model error or
other pilot contamination). The per-link excess terms are summed to form an
equivalent scalar receive disturbance for the data burst. This is an
approximation over receive branches and subcarriers, not a full covariance
reconstruction.

At the first pilot occasion no temporal innovation exists. The model uses the
small-normalized-Doppler approximation

` eta_ICI ~= min(0.1, (pi^2/3) * (f_D / Delta_f)^2 ) `

multiplied by the estimated receive signal power as a conservative floor.
This approximation is only a floor; later pilot innovations update the
estimate from measured observations.

The per-pilot temporal tracker accepts an individual observation-noise
variance for each pilot occasion. Its LS-domain noise profile is

` sigma_LS^2(i) = sigma_thermal^2 / P_pilot + sigma_ICI,LS^2(i) `.

The thermal term is divided by pilot power because it is present before LS
pilot division. The innovation-derived ICI term is already normalized by the
known pilot and is therefore added directly, without a second pilot-power
scaling. Earlier ICI estimates inflate later channel-estimation covariance
calculations.

## 4. Candidate prediction and selection

For every data burst, the predictor gathers post-detection SINRs across all
subcarriers and symbols for both 2-stream Spatial Multiplexing and generalized
two-slot Alamouti. It maps SINR to uncoded Gray-coded square-QAM BER using the
project's AWGN approximation for QPSK and 16-QAM. The selector chooses the
highest net spectral-efficiency candidate predicted to meet the configured
BER target. If no candidate is predicted feasible, it chooses the candidate
with the lowest predicted BER and marks the burst infeasible.

The same effective-noise model is passed to the actual MMSE detectors. The
transmitter mode and modulation are fixed within a burst and may change at the
next pilot occasion.

## 5. Metrics and limitations

The frame's net spectral efficiency is the sum of each burst's payload bits
per OFDM symbol, weighted by the number of data symbols and including pilot
and cyclic-prefix overhead. Uncoded goodput is reported as `NetSE*(1-BER)`.

The pilot-innovation excess is not a direct estimate of the complete
frequency-coupling matrix. It is a scalar disturbance approximation. The
receiver does not cancel the full ICI matrix. The predictor does not model all
cross-link channel-estimation error correlations, decision errors or the full
non-Gaussian distribution of equalizer output. Therefore predicted BER and
predicted feasibility remain advisory link-adaptation metrics. Target
compliance must be judged from measured BER and confidence intervals.

The default 8 frames per condition is a screening run. Increase the number of
frames and bootstrap resamples before making narrow statistical claims.
