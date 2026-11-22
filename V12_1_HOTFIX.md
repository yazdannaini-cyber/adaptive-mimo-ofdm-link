# V12.1 Hotfix — Temporal tracking method-name parsing

## Root cause

`temporal_channel_tracking_apply.m` assigned a local variable named `lower` in
its linear-interpolation branch. MATLAB resolves function and variable names
within a function scope, so that local variable shadowed the built-in `lower`
function used earlier to normalize the method string. The first call therefore
failed with `Unrecognized function or variable 'lower'`.

## Fix

- Renamed the interpolation index to `lowerPilotPos`.
- Replaced lower-case conversion with explicit case-insensitive canonicalization
  using `strcmpi` for `hold-last`, `linear`, and `jakes-lmmse`.
- Added three regression tests for uppercase method names.
- No channel, Doppler, pilot-tracking, or estimator equations were changed.

## Validation status

The package and source changes have been inspected. The MATLAB test suite has
not been executed in this build environment. Run `run_tests` in MATLAB; the
expected suite count is 81 tests.
