# V15.2.1 Independent Numeric Checks

Static checks performed outside MATLAB:

1. A separate reference implementation of the terminated `(7,5)` convolutional encoder and max-log soft Viterbi confirmed exact recovery in a seeded weak-error / low-confidence test. The hard-decision decoder failed the same constructed observation while soft Viterbi recovered all payload and tail bits.
2. Independent Gray-QAM constellation checks for 4-, 16-, and 64-QAM confirmed that the max-log LLR sign selects the transmitted bit for every exact constellation point.
3. Source inspection confirms the MIMO and time-varying Alamouti detectors provide post-SINR values, which are converted to per-symbol equivalent variance before max-log demapping.

MATLAB/Octave is not installed in the build environment. The 141-test suite and the full 3,600-frame MATLAB validation have not been executed here; they must be run by the user in MATLAB.
