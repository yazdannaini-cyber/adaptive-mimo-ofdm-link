# V15.2 Changelog

## Statistical validation and packet-level reliability

- Increased the default validation count from 8 to 100 independent frame/channel
  realizations per condition (36 conditions, 3,600 frame realizations).
- Added a configurable runner option `struct('numFrames',N)` for smaller
  diagnostic runs; runs below 30 frames are explicitly warned as screening runs.
- Each pilot-to-pilot data burst is now a transport packet protected by
  CRC-16/CCITT-FALSE (poly `0x1021`, init `0xFFFF`).
- Added an uncoded + CRC profile and a terminated rate-1/2 convolutional FEC
  profile using constraint length 3 and generator polynomials (7,5) octal.
- Added a hard-decision Viterbi decoder. The two zero tail bits are stripped
  after decoding; CRC is checked on the reconstructed packet payload.
- Compared three transmission strategies (adaptive, fixed SM-16QAM, fixed
  Alamouti-QPSK) under both packet coding profiles.
- Added distinct measures for post-decoder payload BER, pre-decoder coded-bit
  BER, packet error rate (CRC failures), frame error rate, gross SE, net payload
  SE, bit-level goodput, and CRC-accepted packet goodput.
- Added Wilson intervals and frame-cluster bootstrap intervals for payload
  BER, coded BER, packet error rate, and CRC goodput.
- Added full frame-level and burst-decision logs. Expected outputs for the
  default configuration are 216 summary rows, 21,600 frame-log rows, and
  10,800 burst-decision rows.
- Vectorized the sample dimension in the time-varying FIR channel application.
  A deterministic scalar-reference test checks numerical equivalence.

## Interpretation notes

- The adaptive MCS selector still uses the existing pilot-driven pre-FEC QAM
  BER estimate. The same selected MCS is evaluated under both coding profiles;
  V15.2 does not yet optimize a separate FEC-aware MCS policy.
- `conv-crc-1/2` uses hard-decision Viterbi, not soft-decision decoding.
- Each data burst is one packet spanning all active spatial layers in that burst.
- CRC-goodput counts only payload bits from CRC-accepted packets. This is the
  packet-level service metric; bit-goodput is retained as a separate estimate.
- FEC is a rate-1/2 convolutional reference profile, not LDPC/Polar coding.
