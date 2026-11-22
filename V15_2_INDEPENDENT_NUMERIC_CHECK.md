# V15.2 Independent Numeric Checks

These checks were run with an independent Python implementation of the core
bit-level algorithms and FIR equation. They are not a substitute for the
MATLAB test suite.

- CRC-16/CCITT-FALSE standard check string `123456789` produced `0x29B1`.
- The terminated K=3 rate-1/2 convolutional code with generators (7,5) octal
  was independently checked for noiseless Viterbi recovery, single coded-bit
  correction, and the deterministic multi-bit test pattern used by the MATLAB
  regression test.
- The vectorized sample-time FIR expression was compared with the direct
  scalar `n / tap / rx / tx` equation. Maximum absolute complex difference was
  `1.776e-15` for the deterministic 7-sample-symbol, 3-symbol, 2x2 MIMO test.
- Rate-1/2 packet sizing was checked for encoded capacities of 256, 768, 1536,
  and 3072 bits; all codewords fill their assigned QAM bit capacity exactly.

## MATLAB status

The 134-test MATLAB suite and the 100-frame full validation have not been run
inside this build environment because MATLAB and GNU Octave are unavailable.
The user must run the version-specific MATLAB test command before treating the
suite or the 3,600-frame statistical study as passed.
