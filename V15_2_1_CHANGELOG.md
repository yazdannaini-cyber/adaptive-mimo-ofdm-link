# V15.2.1 Changelog — Soft-Decision FEC

## Why this revision
V15.2 used hard QAM decisions before the convolutional Viterbi decoder. This discarded the reliability information available from the MMSE and time-varying Alamouti detectors. The V15.2 results showed that rate-1/2 coding could fail to improve payload BER in difficult 16-QAM spatial-multiplexing conditions.

## Changes
- Added `qam_soft_demod_llr.m`: Gray-QAM max-log bit LLR generation, positive LLR favors bit 0; output order matches `qam_demod`.
- Added `viterbi_decode_rate12_soft.m`: terminated K=3 (7,5) soft-decision Viterbi with max-log branch metrics.
- Added `packet_soft_decode_v15_2.m`: soft QAM LLRs, optional soft Viterbi decoding, and CRC verification.
- Updated the pilot-mapped V15.2 frame simulator to pass per-symbol `1/postSINR` equivalent noise variances from both MMSE MIMO and time-varying Alamouti detectors into the soft QAM demapper.
- Uncoded packets use the hard sign of the LLRs for CRC checking; coded packets use soft-decision Viterbi.
- Added `FECDecoder` metadata to the summary result table.
- Added a separate validation runner and results directory to preserve V15.2 outputs.
- Added seven tests (135–141) covering QAM LLR mapping, soft Viterbi, weak wrong evidence, CRC packet decoding, validation, LLR scaling, and end-to-end integration.

## Preserved experiment design
The six strategy/coding configurations, 36 channel/spatial/Doppler/SNR conditions, default 100 frames per condition, 2,000 bootstrap replicates, pilot layout, adaptation policy, and gross/net payload-rate definitions are unchanged. This allows a controlled comparison against V15.2.

## Important limitation
The soft demapper uses a max-log approximation and a scalar Gaussian residual model per equalized stream. This is an approximation to the post-equalizer disturbance distribution; it is not an exact joint detector LLR for correlated residual interference or ICI.
