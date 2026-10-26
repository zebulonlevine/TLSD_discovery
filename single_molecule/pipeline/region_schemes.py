"""Subunit-specific diffusion-region schemes for the paper SPT figures.

All boundaries are set EMPIRICALLY BY EYE from the posterior-density overlays
(`PaperFigures_updated/overlay_60S_datasets.py`,
`PaperFigures_updated/build_bin_selection_overlay.py`) and are reported as
"manually selected". An automated MLE-trough + posterior-GMM procedure
(`agreed_cutoffs`) exists in the codebase and was explored, but it is NOT used
to set any boundary in the paper -- do not describe these cutoffs as
data-derived.

The selection principle is consistency across experiments: a boundary is
preferred when it sits in a density minimum in more than one dataset, so the
same scheme can be carried across the 60S experiments rather than tuned per
figure.

60S -- [0.004, 0.025, 0.35, 3.0], five bins (frozen 2026-09-16).
  Claims 1, 2 and 3 all image a 60S tag at 15.091 ms / 0.11 um, so their
  diffusion axes are directly comparable and share this scheme.
    0.004  minimum in all three datasets (observed 0.0029-0.0037).
    0.025  NOT a density minimum anywhere -- no 60S dataset has a minimum
           between 0.004 and 0.19. Retained from the previous scheme, and
           kept for resolution: homoharringtonine depletes occupancy below
           ~0.05 and refills above it, so without a boundary in that window
           the run-off redistribution cancels within a single bin and ~37% of
           the claim-2 effect goes unmeasured. Any boundary in 0.025-0.1
           recovers 98-102% of it, so the result does not depend on the exact
           position.
    0.35   minimum in claims 1 and 3 (0.3551 in both); claim 2's nearest is
           0.2816.
    3.0    minimum in claims 1 and 3 (3.09, 2.86); separates the ~0.8 and
           ~4-5 peaks, which are resolved in all three datasets.
  Dropped from the previous six-bin scheme: 0.1, which sat on the CREST of
  the ~0.08 shoulder in claims 1 and 2 rather than in a minimum.

40S -- [0.004, 0.04, 0.3, 2.0], five bins. Chosen for 40S's own broader
  structure; deliberately NOT the 60S boundaries.

60S is used for the 60S contrasts (claim 1 60S, claim 2 S3, claim 3 ER=L29).
40S is used for the 40S contrasts (claim 1 40S, claim 2 L29).
"""

SCHEME_60S = {
    'immobile':     (0.0,   0.004),
    'very_slow':    (0.004, 0.025),
    'slow_medium':  (0.025, 0.35),
    'free_60S':     (0.35,  3.0),
    'free_protein': (3.0,   30.0),
}

SCHEME_40S = {
    'immobile':    (0.0,   0.004),
    'very_slow':   (0.004, 0.04),
    'slow_medium': (0.04,  0.3),
    'free':        (0.3,   2.0),
    'fast':        (2.0,   30.0),
}
