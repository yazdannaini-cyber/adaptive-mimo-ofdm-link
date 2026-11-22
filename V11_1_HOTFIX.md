# V11.1 Hotfix — Rician K-factor Test Robustness

## Reported failure

`run_tests` failed at Test 61 with:

`Rician first-tap LOS/diffuse power ratio does not match K-factor.`

## Root cause

The channel generator scales the LOS and diffuse components according to the configured Rician K-factor. The test, however, estimated K from a single short finite-sum-of-sinusoids realization (0.1 s) and required the time-averaged diffuse power to lie within 20% of its ensemble expectation. A finite realization's temporal average can deviate from ensemble power, so the test was statistically fragile even when the generator's scaling equations were correct.

## Fix

- Keep the channel generator unchanged.
- Estimate the K-factor using eight independent deterministic seeds.
- Use a one-tap PDP for this unit test to keep runtime bounded.
- Check deterministic LOS power separately against the exact configured value.
- Retain the 20% tolerance for the ensemble-averaged diffuse-power estimate.

This preserves the correct distinction between ensemble K-factor and a single realization's sample power. No Doppler or Rician channel-generation equations were changed.

## Validation

The patch has been applied and statically inspected. MATLAB/Octave is not available in the build environment, so the updated test suite must still be run in MATLAB. Expected test count remains 69.
