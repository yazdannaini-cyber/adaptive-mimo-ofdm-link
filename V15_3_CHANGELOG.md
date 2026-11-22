# V15.3 Changelog

## Joint adaptive link controller

- Added an eight-action controller that jointly selects SM/Alamouti, QPSK/16-QAM and uncoded/rate-1/2 convolutional FEC per pilot-to-pilot burst.
- Added per-action CRC ACK/NACK state that persists across frames within one simulated channel condition.
- Added independent success/failure counters for each mode/modulation/coding combination.
- Added link-model PER prior, feedback-updated posterior PER and one-sided 95% Wilson-style upper decision bound.
- Added maximum expected CRC-goodput selection subject to a PER target, with an explicit infeasible fallback if no action meets the constraint.
- Preserved the six V15.2.1 baseline candidates as comparators and appended a seventh `joint-adaptive` result.
- Added burst-level controller decision logs and reference burst prediction logs.
- Increased the default experiment to 300 frames per condition, with configurable frame/bootstrap counts.
- Corrected reference-candidate slot-weight orientation and separated BER-target feasibility from joint-controller PER-target feasibility in result logs.
- Added selector fallback, feedback consistency and reference weighting regression tests. Expected suite size: 160 tests.

## Important model note

The FEC coding gain in the selector is an initial link-adaptation prior, not a guaranteed analytical convolutional-code gain. The controller updates its estimate using observed CRC outcomes. See `V15_3_MODEL.md` for assumptions and limitations.

## Validation status

MATLAB and GNU Octave were not available in the packaging environment. The 160 MATLAB tests and the 300-frame/condition study must be run locally before describing V15.3 as MATLAB-validated.

## Final portfolio release update

- Added the user's 30-frame-per-condition MATLAB summary as curated aggregate CSVs with provenance notes.
- Added final technical report, figures, Persian interview brief and bilingual CV project description.
- Recorded user-reported successful regression run; the bundled suite contains 160 tests.
- Recorded the 30-frame evidence and explicitly marked the default 300-frame run as not supplied.
- Added the measured CRC-goodput improvements (+13.8% Rayleigh, +18.9% Rician) while retaining PER and uncertainty caveats.
