# V12.1 user-result analysis

The supplied 864-row V12.1 output shows:

- Overall mean NMSE: Rayleigh 0.50358 for hold-last, 0.16231 for Jakes-LMMSE, and 0.38242 for linear.
- Overall mean NMSE: Rician 0.41001 for hold-last, 0.083617 for Jakes-LMMSE, and 0.32695 for linear.
- The table's mean condition-level reported improvement of Jakes-LMMSE over hold-last is 60.892% for Rayleigh and 65.482% for Rician.
- By Doppler, Jakes-LMMSE remains lower-NMSE than hold-last and linear for all reported Rayleigh and Rician rows.
- The Jakes posterior-MSE calibration ratio is close to one: overall 1.001 for Rayleigh and 1.0006 for Rician; the Doppler-specific values are also close to one.
- As pilot spacing increases from 2 to 16 symbols, Jakes-LMMSE NMSE increases: 0.075737 to 0.28754 for Rayleigh and 0.04730 to 0.12867 for Rician.
- In the supplied paired table, positive Jakes-LMMSE improvement is reported in 100% of conditions for each channel/Doppler group and each of the three displayed spatial profiles.

Interpretation caveats:

- The 60.892% and 65.482% values are averages of reported condition-level improvements; they are not the same statistic as computing one ratio from the overall mean NMSE values.
- Linear interpolation is a noncausal offline reference because it uses future pilot observations between pilot times.
- The provided spatial-profile summary displays IID, joint-moderate and joint-high only; it does not provide separate TX-only and RX-only summaries in that table.
- V12 uses abstract full-grid orthogonal per-link LS pilot observations and excludes waveform-level pilot resource mapping and intra-symbol ICI from that tracking experiment.

These source-derived results motivate V13's waveform-level pilot mapping and explicit ICI measurement.
