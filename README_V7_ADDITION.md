# V7: Spatial Correlation Stress Test

V7 adds optional TX/RX spatial correlation to the diffuse MIMO channel component.
The new feature is disabled by default, so existing simulations remain comparable
with V6 until explicitly enabled.

## Run regression tests

```matlab
run_tests
```

Expected suite size: 35 tests.

## Run the spatial-correlation experiment

```matlab
addpath('examples');
results = spatial_correlation_validation_run;
```

The script sweeps five TX/RX correlation profiles, Rayleigh and Rician fading,
three SNR points, 500 frames per point and 1,000 bootstrap resamples. Results are
written to `spatial_correlation_results/`.

See `SPATIAL_CORRELATION_MODEL.md` for the equations, assumptions and estimator
limitation.
