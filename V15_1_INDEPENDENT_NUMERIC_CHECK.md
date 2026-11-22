# V15.1 Independent Numeric Checks

These checks were executed in Python/NumPy/SciPy as independent algebraic
checks. They are not MATLAB execution and do not replace the 120-test MATLAB
suite.

## Pilot LS-noise units

For thermal receive-noise variance 0.08, pilot power 4, and already
pilot-normalized ICI variances [0.01, 0.03], the correct LS noise profile is

`0.08/4 + [0.01, 0.03] = [0.03, 0.05]`.

This verifies that ICI variance is not divided by pilot power twice.

## Causal pilot-innovation test case

For a zero-Doppler synthetic channel whose first pilot LS observation is 1.0
and whose next observation is 1.3, with LS thermal variance 0.001, the
per-transmitter excess innovation variance is approximately 0.08860. Summing
the two link excess terms yields approximately 0.17720, satisfying the
purpose-built regression test threshold (>0.1).

## Time-varying Alamouti augmented model

A separate numerical experiment generated random complex 2x2 channels at two
different time slots, encoded the two-symbol Alamouti codeword, built the
augmented observation `[y1; conj(y2)]`, and applied the generalized two-slot
MMSE matrix followed by diagonal-gain correction. With near-zero additive
noise, the maximum recovered-symbol error was approximately `1.74e-10`.

## Execution boundary

MATLAB and GNU Octave are not installed in the build environment. Therefore
MATLAB syntax/runtime behavior, toolbox compatibility and the full V15.1
smoke tests must still be confirmed by running `run_tests_v15_1` in MATLAB.
