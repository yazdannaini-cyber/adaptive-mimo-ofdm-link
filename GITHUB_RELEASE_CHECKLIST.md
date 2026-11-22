# Final GitHub Release Checklist - V15.3

## Source and package integrity
- [x] MATLAB source, config, examples and tests included in the release tree.
- [x] README rewritten as the single GitHub landing page.
- [x] Full historical README preserved as `README_LEGACY_HISTORY.md`.
- [x] Final engineering report, interview brief and CV description added.
- [x] Aggregate results CSVs and figure assets included with provenance notes.
- [x] Generated result folders and personal machine paths are excluded from the source release.

## User-run validation evidence
- [x] Project owner reported the V15.3 test run completed successfully; the bundled test source contains 160 tests. The detailed MATLAB pass log was not included in the supplied output.
- [x] V15.3 30-frame-per-condition validation completed on the user's Windows MATLAB setup.
- [x] Console reported 36 conditions, 252 summary rows, 7,560 frame-log rows, 3,240 joint burst decisions and 3,240 reference burst records.
- [x] Ten-frame and thirty-frame joint CRC-goodput summaries were compared as a limited stability check.
- [ ] Default 300-frame-per-condition experiment has not been supplied / validated in this release record.
- [ ] Exact MATLAB release, runtime and random seed details should be added by the project owner before a formal paper submission.
- [ ] Raw per-frame CSV/MAT output was not provided; the included aggregate CSVs were transcribed from the supplied console summaries.

## Claims and caveats
- [x] CRC-accepted goodput gain is reported relative to the best retained reference candidate in the same 30-frame run.
- [x] PER is reported alongside CRC-goodput; no universal PER-target guarantee is claimed.
- [x] Only 7.3% of Rayleigh and 30.5% of Rician joint-controller data decisions were certified PER-feasible in the final run.
- [x] The controller's prior model, FEC gain curve and fractional-count Wilson-style bound are documented as heuristics/decision statistics.
- [x] The 30-frame run is identified as a confirmation run, not publication-grade validation.

## Optional next release steps
- [ ] Run 300 or more frames per condition and retain raw MAT/CSV output.
- [ ] Record MATLAB release and host details.
- [ ] Review all generated figures from the full run.
- [ ] If preparing a paper, run sensitivity studies for PER target, prior weight, FEC gain heuristic and Doppler/spatial-correlation grid.
