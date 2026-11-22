# V11 Validation Status

- Inherited MATLAB tests: 55 (from the V10.1 baseline).
- Added MATLAB tests: 14.
- Expected total: 69.
- MATLAB test execution: pending user-side MATLAB run.
- Doppler study: pending user-side MATLAB run.
- Independent numerical sanity check: completed with Python for the Jakes ACF approximation and ICI trend; see `V11_INDEPENDENT_NUMERIC_CHECK.md`.

Do not describe V11 as MATLAB-validated until `run_tests` and `doppler_channel_validation_run` have completed successfully in MATLAB.


## V11.1 hotfix note

The original Test 61 used one short finite-sum-of-sinusoids realization to estimate the Rician K-factor. Because K is an ensemble-power ratio, the time-averaged diffuse power of one finite realization can fluctuate enough to fail a strict 20% test threshold. V11.1 averages eight independent seeded one-tap realizations and separately checks the deterministic LOS power. MATLAB execution of the hotfix remains pending until the updated `run_tests` and Doppler validation runner are run by the user.
