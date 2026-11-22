# V15.1.1 Hotfix — Full-MIMO LOS Prior Slicing

## Reported failure

The V15.1 adaptive frame passed the full MIMO LOS prior tensor to a temporal
tracker that processed one isolated Tx pilot stream at a time. Legacy callers
may also provide a single-link prior that is shared across isolated links. The tracker
therefore received observations shaped `[N x Nr x 1 x nPilots]` but a prior
shaped `[N x Nr x Nt x Ns]`, causing `priorMean` dimension validation to fail
when `Nt > 1`.

## Fix

`estimate_pilot_ici_variance` now accepts either a single-link prior (broadcast
across isolated links for legacy/test use) or a full MIMO prior. For the full
MIMO case it slices the matching Tx before temporal prediction.
Each Tx innovation is therefore evaluated against its own LOS prior, with the
original symbol timeline and causal pilot history preserved.

## Regression coverage

Test 123 exercises the ICI innovation estimator with a two-Tx full-MIMO prior
and verifies finite nonnegative outputs for both transmit links. The versioned
runner is `run_tests_v15_1_1` and the suite expects 123 tests.

## Validation status

Static code and archive checks are performed in this environment. MATLAB/Octave
are not installed here; the 123 MATLAB tests must be executed by the user.
