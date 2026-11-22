# Long-run adaptive validation — V4 user run

## Run configuration

- Channel models: Rayleigh and Rician
- SNR axis: Es/N0 per subcarrier, 12–32 dB in 4 dB steps
- Frames per channel/SNR point: 1,000
- Frame-cluster bootstrap resamples: 2,000
- Target BER: 1e-3
- Total result points: 12

The full user-reported result table is preserved in `results/adaptive_validation_v4_user_results.csv`.

## Statistical target status

A point is counted as statistically confirmed only when **both** the upper 95% Wilson bound and the upper 95% frame-cluster bootstrap bound are below 1e-3.

| Channel | Es/N0 | Measured BER | Wilson 95% upper | Frame-bootstrap 95% upper | Confirmed |
|---|---:|---:|---:|---:|---|
| Rayleigh | 12 dB | 2.1601e-2 | 2.1801e-2 | 2.2450e-2 | No |
| Rayleigh | 16 dB | 9.2534e-3 | 9.3855e-3 | 9.7678e-3 | No |
| Rayleigh | 20 dB | 3.9314e-3 | 4.0169e-3 | 4.2334e-3 | No |
| Rayleigh | 24 dB | 1.4075e-3 | 1.4564e-3 | 1.5329e-3 | No |
| Rayleigh | 28 dB | 5.7070e-4 | 5.9816e-4 | 6.3672e-4 | Yes |
| Rayleigh | 32 dB | 3.9013e-4 | 4.0991e-4 | 4.2258e-4 | Yes |
| Rician | 12 dB | 1.4106e-3 | 1.4630e-3 | 1.6104e-3 | No |
| Rician | 16 dB | 3.5108e-4 | 3.7486e-4 | 4.1065e-4 | Yes |
| Rician | 20 dB | 1.7401e-4 | 1.8775e-4 | 1.9790e-4 | Yes |
| Rician | 24 dB | 2.4061e-4 | 2.5408e-4 | 2.6294e-4 | Yes |
| Rician | 28 dB | 6.4750e-5 | 7.1503e-5 | 7.5900e-5 | Yes |
| Rician | 32 dB | 1.7872e-5 | 2.1557e-5 | 2.4570e-5 | Yes |

**Result:** the configured uncoded adaptive link meets the BER target with both statistical upper bounds at 2/6 Rayleigh points (28 and 32 dB) and 5/6 Rician points (16–32 dB).

## Spectral efficiency

- Rayleigh: 2.56 bit/s/Hz at 12–16 dB, rising to 5.03 bit/s/Hz at 32 dB.
- Rician: 2.56 bit/s/Hz at 12 dB, rising to 7.62 bit/s/Hz at 32 dB.
- At 24 dB, Rician reaches 6.72 bit/s/Hz while its measured BER is 2.4061e-4 and both confidence upper bounds remain below target.

## BER prediction calibration

Payload-weighted predicted BER tracks the measured BER reasonably closely across the 12 reported points. Examples:

- Rayleigh 24 dB: predicted 1.4119e-3; measured 1.4075e-3.
- Rayleigh 28 dB: predicted 5.5988e-4; measured 5.7070e-4.
- Rayleigh 32 dB: predicted 3.9059e-4; measured 3.9013e-4.
- Rician 24 dB: predicted 2.5471e-4; measured 2.4061e-4.
- Rician 32 dB: predicted 1.6363e-5; measured 1.7872e-5.

This is evidence of useful calibration under the tested simulation assumptions, not a guarantee for an unmodelled physical channel.

## Important interpretation: BLER/FER

In the V4 source, the reported `BLER` counter is incremented whenever a data frame has **one or more payload bit errors** (`errors > 0`), then divided by the number of simulated frames. There is no FEC decoder in this link. For clarity, this is best called **frame error rate (FER)** for the current uncoded simulation; it must not be represented as post-FEC block error rate.

The result illustrates why BER and frame reliability are distinct:

- Rayleigh 32 dB: BER = 3.9013e-4, but FER = 0.554.
- Rayleigh 28 dB: BER = 5.7070e-4, but FER = 0.524.
- Rician 24 dB: BER = 2.4061e-4, but FER = 0.471.
- Rician 32 dB: BER = 1.7872e-5, FER = 0.045.

A longer adaptive payload and a larger QAM order can change the chance that a frame contains at least one error. A BER target alone does not ensure error-free frame delivery. The current `Goodput_bpsHz` is bit-level goodput (it counts correct bits even in a frame containing errors); it is **not** CRC-accepted packet goodput.

## Engineering conclusion

V4's adaptive modulation and BER predictor show statistically supported BER-target regions for both fading models. However, frame-level delivery reliability remains weaker than bit-level BER suggests, particularly in Rayleigh fading. The next meaningful system-level extension is to compare spatial multiplexing against Alamouti transmit diversity, then consider channel coding/CRC and retransmission as separate later layers.

## Reproducibility and scope

This file records the numbers supplied from the user's MATLAB run. It is not a fresh MATLAB execution performed in the build environment. The new V5 code is a separate extension and must pass its own MATLAB tests before its results are cited.
