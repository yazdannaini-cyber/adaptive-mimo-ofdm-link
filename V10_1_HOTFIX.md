# V10.1 Hotfix — MATLAB Singleton-Dimension Validation

## Reported issue

MATLAB reported `Yfd must be a finite numeric 3-D array` in Test 7 when `mimo_equalize` received `Y0` created as `complex(zeros(N,Nt,1))`. MATLAB suppresses trailing singleton dimensions in `ndims`, so an array logically shaped `[N x Nt x 1]` reports `ndims(Y0) == 2`.

## Fix

- Changed detector input validation to reject arrays with more than three dimensions (`ndims(A) > 3`) rather than requiring `ndims(A) == 3`.
- Applied the same convention to `alamouti_decode` for consistency.
- Updated joint-LMMSE application validation to check dimensions with `size(A,dim)`, which also handles a singleton transmit-link dimension correctly.
- Added Test 55 to compare a single-symbol matrix input against the equivalent slice of a multi-symbol call.

## Scope

No signal-processing equations, estimator parameters, channel models, or mode-selection logic were changed. This is an input-shape validation fix and a regression test.

## Validation status

The source and archive were updated and inspected in this environment. MATLAB/Octave are not installed here, so the 55-test suite still needs to be executed in the user's MATLAB installation.
