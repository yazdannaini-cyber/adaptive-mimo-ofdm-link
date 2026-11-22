# V15.2 Packet and Statistical Model

## Packet formation per data burst

The data symbols between two orthogonal pilot occasions form one transport
packet. A single CRC covers the concatenated bit stream over all active spatial
layers for that burst. The CRC is CRC-16/CCITT-FALSE:

- Polynomial: `0x1021`
- Initial register: `0xFFFF`
- Input/CRC bit order: most-significant bit first
- Final XOR: zero

Two profiles are compared using the same channel, pilots, temporal CSI estimates,
receiver thermal-noise realization, MCS selection, and burst boundaries.

### `uncoded-crc`

The complete QAM bit capacity of a burst is filled with `payload || CRC16`.
The payload length is `capacity - 16` bits.

### `conv-crc-1/2`

The payload is followed by CRC-16 and encoded by a terminated rate-1/2
convolutional code with constraint length 3 and generators `(7,5)` octal.
Two zero tail bits terminate the trellis. The payload length is
`capacity/2 - 16 - 2`, so the complete terminated codeword fills the QAM bit
capacity exactly. The receiver uses hard-decision Viterbi decoding, strips the
tail, and checks the CRC.

## Metric definitions

- **PayloadBER:** BER on the original payload bits after decoding (before FEC
  for `uncoded-crc`, after Viterbi for `conv-crc-1/2`).
- **CodedBER:** hard-demapper BER on the transmitted bit positions before FEC
  decoding; this includes CRC and trellis bits where applicable.
- **PER:** fraction of burst packets whose CRC fails.
- **FER:** fraction of frames containing at least one CRC-failed burst packet.
- **GrossSE:** physical MCS spectral efficiency after pilot and CP overhead,
  before CRC/FEC payload overhead.
- **NetPayloadSE:** original payload bits per normalized time-frequency resource,
  after CRC and FEC overhead but before packet rejection.
- **BitGoodput:** `NetPayloadSE * (1 - PayloadBER)`, an expected correct-bit rate.
- **CRCGoodput:** payload bits from CRC-accepted packets divided by the frame's
  normalized time-frequency resources. This is the measured packet-delivery rate.

## Confidence intervals

Wilson intervals are supplied as a bit-IID or packet-IID reference. Frame-cluster
bootstrap intervals resample complete frames to preserve within-frame fading and
packet dependence. CRC-goodput confidence intervals use a frame bootstrap over
per-frame accepted goodput. For low error counts, inspect both Wilson and cluster
intervals; the nonparametric cluster interval can be degenerate if no frame
contains an error.

## Limits

The MCS predictor estimates pre-FEC QAM BER and the same MCS policy is used for
both coding profiles. V15.2 therefore compares coding gain under a fixed
pilot-driven policy; it does not claim that the MCS policy is optimal after FEC.
The Viterbi decoder is hard-decision. The convolutional code is a reference FEC
profile and does not replace an LDPC/Polar implementation for a modern coded
system.
