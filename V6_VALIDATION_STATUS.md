# V6 validation status

## Confirmed from the supplied V5 run

- The user ran the V5 regression suite: all 24 tests passed.
- The user ran `mimo_mode_comparison_run`: 90 result rows were produced from the configured 3 channels x 3 QAM orders x 5 SNR points x 2 transmission modes, with 200 frames per configuration.

## V6 additions in this package

- Joint candidate selection across spatial multiplexing/Alamouti and QPSK/16-QAM/64-QAM.
- Explicit predicted-target-feasible and fallback status per frame.
- Spectral-efficiency-aware choice among target-feasible candidates.
- Minimum-predicted-BER fallback when no candidate is predicted to meet the target.
- Wilson and frame-cluster confidence intervals for aggregate BER, plus FER Wilson intervals.
- Five added tests, for an intended total of 29 tests.

## Not yet validated on MATLAB

MATLAB and GNU Octave are not installed in the package-build environment. Therefore the new V6 test cases and the 1,000-frame adaptive validation are not represented as executed. The V6 test suite and simulation must be run on a MATLAB installation before the new feature is claimed as validated.
