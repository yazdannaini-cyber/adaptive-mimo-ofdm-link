# V11 — Doppler and Time-Varying MIMO Channel

V11 adds a sample-varying multipath MIMO channel model using a finite sum-of-sinusoids approximation to isotropic-scattering Doppler fading. It supports Rayleigh and Rician channels, the existing transmit/receive spatial-correlation model, reproducible seeds, and sample-by-sample time-varying FIR filtering.

The validation example compares estimated diffuse temporal autocorrelation with the Jakes `J0` reference and measures OFDM off-subcarrier leakage (ICI) for maximum Doppler values of 0, 30, 120 and 300 Hz.

## Run

```matlab
run_tests
addpath('examples');
[results, acfResults] = doppler_channel_validation_run;
```

See `DOPPLER_CHANNEL_MODEL.md` for equations, assumptions and limitations.
