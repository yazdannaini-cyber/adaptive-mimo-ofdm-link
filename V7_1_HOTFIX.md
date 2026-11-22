# V7.1 Hotfix

This maintenance update removes an accidental standalone backslash (`\\`) between Test 29 and Test 30 in `tests/run_all_tests.m`.

The character caused MATLAB's parser error:

`Invalid use of operator.`

No simulation algorithms, channel models, metrics, or test assertions were changed. The intended test suite remains 35 tests.

After extracting this package, run from the project root:

```matlab
run_tests
```

The source-level syntax defect has been removed. MATLAB execution still needs to be confirmed in the user's MATLAB installation.
