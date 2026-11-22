# V15.3 Independent Numeric / Static Review

This note records checks made during package preparation. It is not a substitute for the MATLAB regression suite.

- Reviewed the eight-action selector dimensions: 2 MIMO modes x 2 modulation orders x 2 coding modes.
- Checked the packet payload capacities against CRC-16 and the terminated rate-1/2 convolutional packet wrapper formulas.
- Checked the posterior PER update and one-sided Wilson-style upper-bound implementation for valid fractional error/observation counts.
- Added regression tests for all-feasible maximum CRC-goodput selection, the no-feasible-action fallback, CRC ACK/NACK updates, state-counter consistency, persistent state across frames, fixed-seed repeatability and reference-candidate slot weighting.
- The six V15.2.1 candidates are retained; the joint controller is appended as candidate seven.
- The runner dimensions are 36 conditions, 300 frames/condition, seven candidates and three pilot-to-pilot data bursts per frame. Hence expected table dimensions are 252, 75,600, 32,400 and 32,400 rows, respectively.
- Reviewed result semantics: joint PER-target feasibility is stored separately from the baseline BER-target feasibility field. The selector fallback remains explicitly marked infeasible.

MATLAB and GNU Octave are unavailable in the build environment. No MATLAB execution result is claimed here.
