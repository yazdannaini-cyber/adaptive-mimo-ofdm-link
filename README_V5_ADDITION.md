# V5: Spatial multiplexing vs Alamouti STBC

V5 adds a second transmission architecture to the OFDM system so reliability/throughput trade-offs can be studied instead of relying only on QAM-order adaptation.

Run the existing regression suite:

```matlab
run_tests
```

Expected result: 24 tests pass.

Run the new mode comparison:

```matlab
addpath('examples');
comparisonResults = mimo_mode_comparison_run;
```

Results are saved in `comparison_results/` as CSV, MAT, BER figures, and spectral-efficiency figures. See `MIMO_MODE_COMPARISON_NOTES.md` for equations, assumptions, metric definitions, and limitations.
