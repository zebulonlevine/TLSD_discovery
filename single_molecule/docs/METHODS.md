# Methods — single-molecule sections

Draft text for the three currently-empty Methods headings in
`ZebRibosomeTurnover_draft2 260912.pdf`: **Single-molecule microscopy**,
**Single-molecule analysis**, and **ER co-labeling and analysis**. Every
numeric value below is read from the deployed `settings.pkl` of the
corresponding dataset or from the ND2 metadata, not asserted from memory; the
machine-readable dump is `verification/saspt_equivalence_report.json` (test T4).

---

## Single-molecule microscopy

> Cells were imaged on a Nikon Ti2 inverted microscope with an Ilas2
> illuminator in HiLo (highly-inclined thin illumination) mode and a dual
> Photometrics Prime 95B sCMOS camera, through a 100× oil objective (effective
> pixel size 0.11 µm). Ribosomes were labelled sub-stoichiometrically so that
> individual particles could be resolved, and the Perfect Focus System was
> engaged throughout. HaloTag-labelled ribosomal proteins were imaged in the
> far-red channel (`em640`; `JFX650` in the single-colour experiments) and
> SNAP-Sec61β in the orange channel (`em549`). Movie length differed between
> experiments (table below): 1000 frames at a 15.091 ms frame interval (15.1 s)
> for the 60S/40S UBA5 comparison and for the ribosome phase of the ER
> experiment, 500 frames at ~15 ms (7.5 s) for the drug experiment, and 2000
> frames.

Channel assignment is read off the acquisition itself, not assumed: in the ER
ND2s the `em549` frames are bright and photobleaching over frames 0–99 and then
absent, while `em640` is absent over 0–99 and present from 100 to 1099 — so
`em549` is the briefly-imaged ER marker and `em640` is the long single-molecule
movie. The excitation lines recorded as active are 561 and 640 nm; the metadata
does not bind them to channels in a machine-readable way, so **confirm the
ex/em pairing and the dye names before submission**.

Confirmed per-dataset acquisition (report as a supplementary table):

| Figure panels | Dataset | SPT frames | Frame interval | Movie length | Pixel size | Frame (px) | Field of view | FOV / condition |
|---|---|---|---|---|---|---|---|---|
| 5C, S7A–C | `2025_06_06_zeb` | 1000 | 15.091 ms | 15.1 s | 0.11 µm | 1024 × 704 | 112.6 × 77.4 µm (8720 µm²) | 10 |
| S7D, S7E | `2026_07_08` | **500** | 14.75–15.09 ms | 7.4–7.5 s | 0.11 µm | 532 × 348 | 58.5 × 38.3 µm (2240 µm²) | 20 (19 for sgNT/emetine) |
| 5E, 5F, S7H–J | `2025_09_15_zeb_ER` | 1000 (after 100 ER frames) | 15.091 ms | 16.6 s total | 0.11 µm | 800 × 700 | 88.0 × 77.0 µm (6780 µm²) | 10 DMSO, 14 PT-179 |

Frame dimensions are the acquired sensor region; no further cropping is applied
in analysis (no `crop` is set in any `settings.pkl`), so the analysed area
equals the acquired area. The imaged area differs about six-fold between
experiments, largely because a shorter sensor region buys a faster frame rate on
a rolling-shutter sCMOS. One consequence worth stating: a "field of view" is the
resampling unit for every confidence interval and permutation null here, and it
is a constant amount of material *within* an experiment but not *between* them —
10 fields of view in the 60S/40S experiment cover ~4× the area of 10 in the
drug experiment. All resampling is within-experiment, so this does not affect
any interval; it does mean field-of-view counts should not be read as comparable
sample sizes across figures.

The ER mask area fractions in Fig S7I are fractions of the full 800 × 700 field
(the generator asserts that the on-structure, near-ER and off-ER fractions sum
to 1), so the ~29.7% broad-ER coverage corresponds to ~2000 µm² per field of
view.

Not used in any figure: the same 60S/40S experiment also contains a 5 ms arm
acquired over a shorter 1024 × 240 sensor region (2970 µm²); only the 15 ms
files feed Figures 5C and S7A–C.

Frame counts are from the ND2 acquisition metadata and confirmed against the
frame range actually present in the tracked trajectories (0–999, 0–499, 0–999
and 0–1999 respectively).

**The ER experiment is a two-phase acquisition, not simultaneous two-colour.**
Each 1100-frame ND2 holds the ER channel in frames 0–99 and the ribosome
channel in frames 100–1099, the other channel being written as gap frames in
each phase. The ER is therefore illuminated for **100 frames (1.51 s)** and the
ribosome SPT for **1000 frames (15.1 s)**; the ER mask is a median projection
over those 100 ER frames, and tracking sees only the 1000 SPT frames,
re-indexed to 0–999. Note the draft legend for Fig 5E says the ER channel was
imaged "for 150 ms prior" to the ribosome channel; the metadata gives 100 frames
at 15.091 ms ≈ 1.5 s, so that number needs correcting (it cannot be a per-frame
exposure either, since the sampling interval is 15.091 ms).

Frame intervals in the drug experiment are the per-movie values measured by the
acquisition software and span 14.75–15.09 ms, with untreated movies at the long
end and emetine at the short end; a single nominal 15.091 ms was used for all of
them. Since D scales as 1/Δt this biases the emetine conditions' diffusion
coefficients low by at most 2.3%, i.e. 0.010 decades — under a third of one grid
state, and ~1.7% of the narrowest bin's width. Recorded for completeness rather
than because it affects a conclusion.

> Because the recovered diffusion-coefficient distributions depend on the
> acquisition regime, all comparisons are made strictly within a dataset. We do
> not compare absolute diffusion coefficients or state occupancies between
> figures, and diffusion-state boundaries are defined per dataset.

---

## Single-molecule analysis

### Localization and tracking

> Spots were detected and localized with `quot`
> (github.com/alecheckert/quot): log-likelihood-ratio detection (window 15 px,
> threshold 15.0, k = 2.0) followed by least-squares fitting of an
> integrated 2-D Gaussian (window 15 px, PSF σ = 2.5 px, camera offset 100 ADU).
> Localizations were linked by Euclidean nearest-neighbour assignment with a
> 1.2 µm search radius and no gap-closing (`max_blinks = 0`). Trajectories whose
> mean spot intensity fell below 100 ADU above background were discarded.
> Single-frame images of sparse labels are close to background (real spots sit
> only ~80–130 ADU above a ~148 ADU background), so displayed stills are
> maximum- or median-intensity projections; projection is never applied to the
> data used for tracking.

### Diffusion-state inference

> Population diffusion spectra were estimated with state-array SPT (saSPT
> v0.4.0; Heckert et al.), which places a Dirichlet prior over occupancies on a
> fixed grid of diffusion states and infers the posterior by variational Bayes
> expectation–maximization. We used the regular-Brownian-motion-with-error
> (RBME) likelihood on a grid of 150 diffusion coefficients log-spaced from
> 10⁻³ to 10² µm² s⁻¹ crossed with 36 localization-error states from 0 to
> 0.070 µm in 0.002 µm steps — saSPT's own default error grid — with
> trajectories split into sub-trajectories of at most 10
> displacements (`splitsize = 10`), a prior concentration of 1.0 pseudocount per
> state, and 200 VBEM iterations — saSPT's own defaults for the last two. No
> sub-sampling was applied: every sub-trajectory of a condition entered its fit.
> Each sub-trajectory is weighted by its number of displacements, which is
> saSPT's definition of occupancy.

> Occupancies were corrected for defocalization exactly as implemented in
> saSPT: the occupancy of each diffusion state is divided by the fraction of
> particles of that diffusion coefficient still within the focal slab one frame
> later, computed by numerical propagation of a uniformly-filled slab through
> the corresponding Brownian kernel with absorbing boundaries. A focal depth of
> 0.7 µm was used throughout (see *Focal depth* below).

### Localization error

`quot`'s detection kernel (`k = 2.0`) and fitted PSF width (`sigma = 2.5`) are
specified in **pixels**. The per-localization fit uncertainty is roughly
constant in pixels across the three experiments (0.27–0.36 px median), and at a
common 0.110 µm pixel size therefore also in microns:

| dataset | pixel size | median fit uncertainty (px) | (µm) |
|---|---|---|---|
| `2025_06_06_zeb` | 0.110 | 0.289 | 0.032 |
| `2026_07_08` | 0.110 | 0.362 | 0.040 |
| `2025_09_15_zeb_ER` | 0.110 | 0.273 | 0.030 |

**What `x_err`/`y_err` are, and what they are not.** These are the per-spot
standard errors returned by the Gaussian fit and written to the trajectory
tables. Neither stage of the analysis consumes them: `quot`'s linker never
references them (no occurrence in `quot/track.py`), and saSPT reads only the
`trajectory`, `frame`, `y` and `x` columns (`saspt/constants.py:39,42,48`;
`saspt/utils.py:125,142`). Localization error enters saSPT solely as the
`loc_errors` grid — a model parameter marginalized over, not per-localization
data. Pixel size likewise enters at exactly one point, the conversion of
positions to microns (`saspt/utils.py:328`).

So the fit uncertainties below are an *independent estimate of what the grid
ought to span*, not a quantity the pipeline was wrongly ignoring. They are a
reasonable but imperfect proxy: the effective error the grid represents also
absorbs stage drift, pixelation and registration, so agreement is not expected
to be exact. That is why the grid question was settled by refitting rather than
by comparing numbers.

**All four datasets are fit on saSPT's default 36-state error grid, 0 to
0.070 µm in 0.002 µm steps.** That grid spans every measured value in the table
above with room at both ends, so no dataset sits outside it.

This replaced a 5-state grid spanning only 0.025–0.032 µm, and the reason for
the change is worth recording, because it is what the robustness work below was
originally written about. Reading the localization-error posterior off the
cached likelihoods on that grid
(`verification/loc_error_grid_check.csv`) showed it was essentially **flat** in
every dataset — ~0.20 of the mass on each of the five states, posterior mean
0.0285 µm throughout, where a uniform prior gives exactly 0.20. The grid was not
piling up against an edge; it was too narrow (a 28% range) for the data to
discriminate within, and in practice acted as a **fixed localization error of
≈0.0285 µm** for every experiment rather than as something inferred. On the
intensity-filtered set saSPT actually sees, the measured fit uncertainty is
0.029 / 0.035 / 0.028 µm for claims 1–3 — above that grid for the drug
experiment.

The consequence was measured before the change, by refitting one headline
contrast per dataset on a wide grid (0.010–0.055 µm, 10 states) spanning every
measured value, against the narrow grid, on identical subsampled data
(`verification/loc_error_refit_comparison.csv`):

The cutoff-free statistic was stable; the per-bin calls were not uniformly so:

| contrast | W₁ narrow grid | W₁ wide grid | bin calls flipping |
|---|---|---|---|
| claim 1, sgNT vs sgUBA5 | 0.271 | 0.240 | 2 of 6 (immobile, medium) |
| claim 2, none vs harringtonine | 0.465 | 0.522 | 0 of 6 |
| claim 3, noDegron vs Degron | 0.133 | 0.148 | 0 of 6 |

Five of 24 bin calls moved. Two points on how to read that. First, W₁ did not
weaken on the wide grid in the two datasets whose error sat outside the narrow
one — it *increased* (claim 2, 0.465→0.522) — so the distributional
difference between conditions was never an artefact of the grid.
Second, the flips are boundary mis-registration, and the mechanism is measured
(`verification/grid_shift_mechanism.csv`). saSPT's jump variance along one axis
is 2(D·Δt + σ²), so D and σ trade off directly at fixed data:
D·Δt = var/2 − σ². On the narrow grid the localization-error posterior is flat
and the fit effectively uses the midpoint, 0.0285 µm. On the wide grid it is no
longer flat and saSPT infers **0.0320 µm** — *higher*, not lower, and close to
the narrow grid's upper edge. (Note this also means the per-spot fit
uncertainty, 0.0223 µm for this dataset, *underestimates* the effective error,
as anticipated: the effective σ absorbs drift, pixelation and registration.)

Raising σ from 0.0285 to 0.0320 µm forces D down by
(0.0285² − 0.0320²)/0.01 = **−0.021 µm² s⁻¹**, and the observed median log₁₀D
moves −1.104 → −1.174, i.e. −0.070 decades. The shift is *absolute* in D, so on
a log axis it is wildly non-uniform:

| D (µm² s⁻¹) | 0.03 | 0.2 | 0.8 | 2 | 10 |
|---|---|---|---|---|---|
| shift (decades) | **−0.53** | −0.048 | −0.012 | −0.005 | −0.001 |

With boundaries held fixed, only the slow end can move — and the measured
occupancies do exactly that: mass leaves `slow` (0.320→0.284) for `very_slow`
(0.203→0.237) and `immobile` (0.133→0.144), while `medium`, `free_60S` and
`free_protein` change by ≤0.005. Those are precisely the three bins that flip.
So the effect is not lost; it is redistributed across boundaries that were
registered to a mis-specified σ.

Both halves of that problem have since been addressed, and the panels in this
work carry neither. The error grid is now the 36-state default, which spans
every measured value, so no dataset is fit at a mis-specified σ; and the bin
boundaries were re-placed on each dataset's own density minima and frozen
(2026-09-16).

### Localization error is set by spot brightness

Across the four experiments the fit uncertainty tracks spot brightness exactly
(Spearman ρ = −1.00 against median I0), and a 1/√N scaling predicts its
magnitude to within 10% in every dataset. The drug experiment has the largest
uncertainty because its spots are ~1.8× dimmer than the others'
(`verification/loc_error_origin.csv`) — not because of lower illumination: it
delivered the most 640 nm light of the four. The ND2 records the illumination
path in two separate blocks — `NIDAQ, MultiLaser(AOTF)` (the modulator) and
`NIDAQ, MultiLaser(Genesis Heads)` (the sources), plus the `Shutter(AOTF)`,
`Shutter(AOTF Blank)` and `Shutter(Genesis Heads)` states — so both have to be
read together:

| dataset | AOTF 640 | Genesis 640 | Ilas2 preset | median I0 |
|---|---|---|---|---|
| claim 1 | 100% | 35% | TIRF 561 | 598 |
| claim 2 | 100% | **50%** | **TIRF 640** | **339** |
| claim 3 | 100% | 30% | TIRF 561 | 615 |

Brightness does not track delivered power: claim 2 receives the most 640 nm
light of the three and is the dimmest. (With three datasets a rank correlation
is not worth quoting.) The leading
remaining candidate is the illumination angle — claim 2 is the only experiment
whose Ilas2 preset is `TIRF 640` rather than `TIRF 561`, and a steeper,
more TIRF-like sheet brightly illuminates only molecules close to the coverslip,
so molecules sampled deeper in the cell would appear dimmer despite high power.
That would also imply a shallower effective focal slab for this dataset than for
the others. It is a hypothesis the metadata cannot settle; dye identity and
labelling chemistry are the other candidates.

### Focal depth

> The defocalization correction used a focal depth of 0.7 µm, saSPT's default,
> for every dataset (`focal_depth` in each `settings_36state.pkl`). The axial
> detection range was also estimated from the data, as a check on that value.
> Because fast particles leave the focal slab, the
> observed distribution of trajectory length falls with diffusion coefficient,
> while photobleaching and blinking do not depend on it; taking the ratio of
> the length-survival curve in each diffusion bin to that of the slowest bin
> therefore isolates the defocalization term, which was fit for focal depth by
> least squares. This gave 1.30 µm for the 60S/UBA5 dataset and 1.37 µm for the ER dataset, against a diffraction-limited depth of field of 0.51 µm for this objective (Berek's formula at NA 1.49). The bootstrap interval on the fit is narrow (1.22–1.31 µm over 192,701 sub-trajectories) but reflects sampling error only; the systematic uncertainty is much larger and is what the sweep below addresses. Because linking failure and
> motion blur also shorten fast trajectories and are absorbed into the same
> decay, this estimate is a lower bound on the true axial detection range. All
> reported contrasts were additionally recomputed across focal depths of 0.1,
> 0.3, 0.5, 0.7, 1.0 and 1.5 µm and with the correction disabled
> (`verification/focal_depth_calibration.json`); the Wasserstein statistic for the sgNT-vs-sgUBA5 60S contrast varied only between 0.227 and 0.249 across that range, so no reported conclusion depends on the value chosen.

> Note that the per-field-of-view occupancies on which all statistics below are
> computed are built from the posterior responsibilities of individual
> sub-trajectories. saSPT applies the defocalization correction after the
> variational loop, to the population-level sufficient statistic only, so those
> responsibilities — and hence every per-bin panel in this work — are
> independent of the focal-depth choice by construction. Only the population
> posterior curves and the Wasserstein statistic depend on it.

### What the y-axes show

saSPT's own vocabulary for the quantity is **occupation**: its plotting module
labels these axes "Posterior occupation", "Naive occupation" and "Marginal
posterior occupation" (`saspt/plot.py:231, 294–295`). The word "density" never
appears as an axis label in the package.

The posterior is a discrete Dirichlet mean over a fixed grid of diffusion
states — here 150 states log-spaced from 10⁻³ to 10² µm² s⁻¹ — so the native
quantity is a probability **mass** pₖ with Σpₖ = 1. Two representations of it
appear in these figures, and they carry different units:

- **Overlay panels** (5C, 5F left, S7B, S7J) plot pₖ divided by the grid
  spacing Δlog₁₀D = 5/149 ≈ 0.0336 decade, i.e. a probability **density with
  respect to log₁₀D**, which integrates to 1 over the plotted axis. Axis label:
  *Probability density*, with the measure given in the legend (below).
- **Bin panels** (5F right, S7C, S7D) plot pₖ summed over a bin: a
  dimensionless **fraction** of posterior occupancy — not a density, so no
  qualifier would make "density" correct here. Axis label: *ΔProbability
  fraction*.

**Legend sentence** (paste into any panel using the density axis):

> "Probability density is per decade of diffusion coefficient (i.e. per unit
> log₁₀ D), so area under the curve between two diffusion coefficients is that
> interval's occupancy; read by area, not by peak height."

Note the axis vocabulary is a style choice, not a correctness one. saSPT's own
plotting module calls this quantity "occupation" (`saspt/plot.py:231, 294–295`)
and never labels an axis "density"; "probability density" is equally accurate
provided the measure is stated, which is what the legend sentence does. What is
**not** interchangeable is density vs. fraction between the two panel types.

Why the density rescaling for the overlays, when saSPT itself plots the raw
mass: the height of a raw-mass curve depends on how many grid states were
chosen. The same simulated peak reaches 0.057 on saSPT's default 100-state
grid, 0.038 on the 150-state grid used here, and 0.019 on a 300-state grid;
after dividing by the spacing it is 1.137, 1.139 and 1.140. The density is the
grid-invariant representation, and it is the one under which **area** under an
interval equals that interval's occupancy — which is exactly what the bin
panels report, so the two panel types become readable as the same number.

**Y-axis range.** The density overlays are drawn on a fixed 0–1.0 axis rather
than auto-scaled. The floor-state pile-up is several-fold taller than any
resolved peak, so auto-scaling flattens the part of the curve the figure is
about; saSPT's own `plot_posterior` carries a `truncate_y_axis` flag for the
same reason. At 1.0 the only features clipped in the reported panels are
floor-state spikes at D = 10⁻³ (Fig 5F left, S7J), while
the tallest resolved feature — the free-protein peak of Fig 5C at
D = 7.2 µm² s⁻¹, height 0.96 — is left intact. Each build logs any curve it
clips and by how much, so a truncated peak is never silent.

**Decade.** A decade is a factor of 10 in D — one unit on the log₁₀D axis.
The state grid spans 10⁻³ to 10² µm² s⁻¹, i.e. 5 decades, over 150 states, so
adjacent states differ by 10^0.0336 = 1.08 (8% in D) and there are 29.8 states
per decade. "Probability per decade" therefore means probability mass per
factor-of-10 interval of D, and the five diffusion bins are **not** equal widths
on that scale:

| bin | D range (µm² s⁻¹) | width (decades) |
|---|---|---|
| I immobile | ≤ 0.004 (grid floor 10⁻³) | 0.60 |
| II very_slow | 0.004 – 0.025 | 0.80 |
| III slow_medium | 0.025 – 0.35 | 1.15 |
| IV free_60S | 0.35 – 3 | 0.93 |
| V free_protein | 3 – 30 | 1.00 |

Because the widths differ by up to 1.9×, **peak height is not proportional to
bin occupancy**: in the sgNT 60S curve bins III and IV hold almost the same
mass (0.290 and 0.299), yet bin III's mean density is 21% lower (0.253 versus
0.321 decade⁻¹) purely because it is 1.23× wider. Read the overlay by area, not
by height (`verification/fig5C_decade_units_check.csv`).

Two further points for reading these axes. First, the measure: this is a density
in log₁₀D, not in D; a density in D differs by the Jacobian D·ln10. Second, the
spike at the left edge of every overlay is the grid floor, not a resolved state
— sub-trajectories slower than the localization-error limit pile onto the lowest
grid point, so bin I's apparent peak height (0.66 decade⁻¹ for sgNT) is an edge
artefact and only its integrated mass is meaningful. Mass above the top bin
(> 30 µm² s⁻¹) is 0.0006, so the five bins effectively partition the posterior.
See `PaperFigures/METHODS_posterior_density_note.md` for the derivation.

### Diffusion-state boundaries

> Boundaries between diffusion states were chosen manually by inspection of the
> pooled posterior density overlays, as round values stated in advance of the
> contrasts rather than fit to them, and applied identically to every condition
> within a dataset. They were not derived by any automated procedure. Placement
> was guided by two criteria: a boundary was preferred where the density showed
> a minimum in more than one independent experiment, so that a single scheme
> could be carried across datasets rather than tuned per figure; and a boundary
> was retained where it was required to resolve a redistribution of occupancy
> that would otherwise cancel within a single bin. The three 60S experiments
> (Figs 5C, 5F, S7C, S7D, S7J) were acquired under identical conditions
> (15.091 ms, 0.11 µm pixels) and therefore share one scheme.

| Scheme | Boundaries (µm² s⁻¹) | Bins | Used for |
|---|---|---|---|
| 60S | 0.004, 0.025, 0.35, 3 | 5 | Fig 5C, 5F, S7C, S7D, S7J |
| 40S | 0.004, 0.04, 0.3, 2 | 5 | Fig S7B |

Of the four 60S boundaries, 0.004, 0.35 and 3 each sit in a density minimum in
at least two of the three experiments (observed at 0.0029–0.0037, 0.355 in two
datasets, and 2.86/3.09). The 0.025 boundary does not: no 60S dataset shows a
minimum between 0.004 and 0.19. It is retained on the second criterion —
homoharringtonine depletes occupancy below ≈0.05 µm² s⁻¹ and refills it above,
so with no boundary in that window the redistribution cancels and 37% of the
measured drug effect is lost. Any boundary between 0.025 and 0.1 recovers
98–102% of it, so the reported effect does not depend on its exact position. A
boundary at 0.1, used in an earlier version of this scheme, was dropped because
it falls on a local density *maximum* in two of the three datasets.
Selection is reproduced by `overlay_60S_datasets.py` and
`build_bin_selection_overlay.py` (`verification/bins_*`).

### Translation-state assignment (Fig S7E)

> Diffusion bins were assigned to translation states a priori, from the
> diffusion regime alone, and identically for both genotypes: bins I and II
> (D ≤ 0.025 µm² s⁻¹) as polysome-associated, bin IV (0.35–3 µm² s⁻¹) as free
> 60S subunit, and bin III (0.025–0.35 µm² s⁻¹) as ambiguous. Bin V
> (> 3 µm² s⁻¹) is faster than a free 60S subunit can diffuse and is not
> assigned a translation state. The two drug treatments provide an independent
> test of these assignments rather than their basis. Harringtonine blocks
> initiation and permits run-off, so a polysome-associated species is expected
> to decrease and free subunits to accumulate; emetine freezes elongating
> ribosomes, so a polysome-associated species is expected to be stabilised and
> free subunits depleted. Both predictions were borne out in both genotypes
> for every bin on which they make a directional claim.

The assignment is stated rather than inferred deliberately. An earlier version
called each bin from the sign of its own two drug responses, which made the
genotypes disagree — sgNT identified bin II as the polysome and sgUBA5 bin III
— because the label was being read off the same data it was then used to
interpret. Fixing the assignment removes that circularity and turns the drug
responses into a falsifiable prediction. The sgUBA5 bin III responses do trend
the way a polysome would (harringtonine −0.030, emetine +0.034); bin III is
nonetheless reported as ambiguous, which is the conservative call, and that
trend is consistent with — not evidence against — an intermediate-mobility
translating species accumulating under UBA5 knockdown.

### Statistics

> For every contrast, each group's posterior was fit independently over that
> group's sub-trajectories, so the compared posteriors carry no cross-group
> coupling. The unit of analysis is the field of view: each sub-trajectory
> contributes its full posterior over D, these are averaged within a field of
> view, and fields of view are averaged within a condition. Sub-trajectories
> cluster within a cell and a field of view, so resampling them individually
> would be anticonservative.

> All uncertainty was estimated by resampling fields of view and **re-running
> the full variational inference on every resample**, rather than resampling a
> single fixed posterior. Confidence intervals are percentile intervals from
> 500 bootstrap resamples drawn with replacement within each condition;
> p-values are from 500 permutations of the field-of-view labels between
> conditions, with group sizes preserved. Each draw therefore requires two
> independent variational fits, and the reported intervals propagate
> uncertainty in the fitted posterior as well as sampling variability.
> Intervals obtained by resampling a single fixed posterior are systematically
> narrower — by 1.2× to 2.7× across the bins of Fig. S7C, median 2.3× — and
> are not reported.

> Whether two conditions' diffusion distributions differed overall was tested
> with the cutoff-free Wasserstein-1 distance between their posteriors on the
> log₁₀D axis, compared against the same field-of-view permutation null.
> Per-state effects are the between-condition difference in field-of-view
> occupancy of each diffusion bin. Per-bin quantities are integrals of the
> stored per-draw posterior differences, so the state boundaries can be revised
> without repeating the inference.

> Permutation p-values use the add-one estimator (1 + k)/(1 + B). **With
> B = 500 the smallest attainable value is 0.002**, so a p-value reported at
> 0.002 means no permutation reached the observed statistic and is not
> comparable to a smaller value obtained at higher B; such values are flagged
> in the per-panel statistics tables. Per-bin p-values are corrected across
> the bins of each panel by the Benjamini–Hochberg procedure; the panel is the
> correction family because its bins partition the diffusion axis and "which
> state changed" is a selection over that closed set.

> The observed statistic and its resampled counterparts are computed by the
> same code path, so the two cannot diverge. Resampling for each contrast is
> seeded independently, so statistics formed across contrasts are not
> correlated through the random number stream.

> **Zone contrasts (Fig. 5F).** Every field of view contains both zones, so
> relabelling whole fields of view cannot generate a zone contrast; the null
> instead permutes zone labels among trajectories within each field of view,
> with all sub-trajectories of a trajectory moving together and each field of
> view's zone composition preserved. Because the two zones are fit jointly over
> their pooled sub-trajectories and then partitioned, this relabelling leaves
> the fit unchanged, so the permutation is exact with respect to the inference
> and was run at B = 2000; its p-values are not subject to the floor above. The
> bootstrap does re-fit, since resampling fields of view changes the pooled set.

> **Zone × treatment interaction (Fig. 5F).** The zone contrast is taken within
> each field of view, which removes the zone main effect before any resampling
> and reduces the response to one value per field of view; treatment labels
> were then permuted across fields of view (B = 10,000), which is exact under
> the hypothesis that treatment does not change the zone difference. The
> interval was formed from the two treatments' independent bootstraps: because
> the two are independent samples the pairing between their draws is arbitrary,
> so all B × B differences were used, giving the empirical convolution rather
> than a B-sample approximation of it. The draws are not mutually independent —
> each appears B times — so this sharpens the quantiles without extending tail
> resolution, and the p-value remains floored at the per-arm resolution.

### A note on 2 × 2 interaction nulls

The two-factor interaction analysis was developed for a TLSD1 dose series that
is not part of this figure set. It is retained here because it constrains how
Fig. 5F's interaction may be tested.

Relabelling either factor of a 2 × 2 design pools units that differ by that
factor's main effect, so the resulting null is inflated by it, and the two
single-factor nulls disagree — measured at 54% in null width for the bin
carrying the largest main effect, moving p from 0.0005 to 0.013. The standard
remedy, Freedman–Lane residual permutation, cannot be combined with refitting:
its null draws are synthetic per-field-of-view occupancies (fitted values plus
permuted residuals) with no sub-trajectories beneath them, so the variational
inference cannot be re-run on them, and applied to already-fitted occupancies
it misses the uncertainty in the shared localization-error posterior — the term
responsible for the interval widening quoted above.

Resampling within each cell avoids the problem entirely: nothing is relabelled,
so both main effects are preserved by construction, and every draw is a real
fit of real data. Fig. 5F's interaction is additionally paired within fields of
view, which removes the zone main effect before resampling and admits an exact
treatment permutation; that is the route used there.

### Software and reproducibility

> Analysis used saSPT v0.4.0 and `quot` (commit 1b9051e) under Python 3.8.
> Because the per-field-of-view resampling requires slicing a fixed
> sub-trajectory likelihood many times, the variational loop was re-expressed
> in a form that avoids materializing the full responsibility array each
> iteration. This re-expression was verified to reproduce
> `saspt.StateArray.infer_posterior` to machine precision on the same
> likelihoods (maximum absolute occupancy difference 2.0 × 10⁻¹⁵), with the
> likelihood, state grid and defocalization correction taken unmodified from
> the saSPT package; the full verification suite is provided
> (`verification/saspt_equivalence_report.json`).

---

## ER co-labeling and analysis

> The ER was imaged in a separate channel from the ribosome: SNAP-Sec61β was
> labelled to completion and acquired immediately before the
> sub-stoichiometrically labelled 60S acquisition, as the first 100 frames of a
> dual-channel acquisition whose remaining 1000 frames are the SPT movie.

> **ER reference image.** The two channels were separated and the ER channel
> reduced to a single reference image by a **median projection over its 100
> frames**; the median rather than a mean so that any residual single-molecule
> signal bleeding into the ER channel does not survive into the mask. The SPT
> channel was written out from frame 100 onward and passed to `quot` unmodified —
> no projection, filtering or denoising is ever applied to the data used for
> tracking.

> **Illumination correction.** The reference image was corrected for uneven
> field illumination by dividing by a heavily blurred copy of itself (Gaussian,
> σ = 100 px), with the blurred estimate floored at 1 to avoid division by
> near-zero, and the quotient multiplied by the input mean to restore
> approximately the original intensity scale. Correction
> output was kept in the image's natural intensity units; no per-image min–max
> rescaling was applied, since that would make intensities incomparable between
> fields of view and the mask threshold below is partly a global one.

> **Segmentation pre-processing.** Before thresholding, the corrected image was
> passed through a second, identical σ = 100 px flat-field division to remove
> residual large-scale structure, then deconvolved by Richardson–Lucy against a
> Gaussian PSF (σ = 0.74 px, 2 iterations), again rescaled to the input mean. The deconvolution is deliberately shallow: it is there to
> sharpen tubule edges enough for a local threshold to follow them, not to
> reconstruct structure, and at 2 iterations it does not amplify noise into
> spurious objects.

> **Mask construction.** The ER was segmented into three classes by adaptive
> local thresholding at two spatial scales. A pixel passes a given scale if it
> exceeds the local Gaussian-weighted mean of a surrounding block by more than
> 1.5%, computed at a **fine** block of 25 px (individual tubules) and a **broad**
> block of 101 px (the region the ER network occupies). Each thresholded mask was
> morphologically closed with a 2 px disk and objects smaller than 500 px
> removed. Because a purely local criterion will find "structure" in an empty,
> dark part of the field, a global floor was applied on top: a pixel must also
> exceed the **60th percentile of the whole processed image** to be eligible for
> either mask. `on_structure` is then the fine-scale mask, `near_ER` the broad
> mask excluding `on_structure`, and `off_ER` everything in neither. The broad
> mask (`on_structure` ∪ `near_ER`) is the "on-ER" zone used for the contrasts.

> Mask quality was checked per field of view by the
> intensity profile across the mask boundary (Fig S7H) and by the total masked
> area (Fig S7I), which did not differ between DMSO and 24 h PT-179.

> Each trajectory was tested individually for enrichment on or off the mask
> rather than assigned by a single-point overlap. `on_ER` denotes the broad
> on-ER zone (mask-proximal trajectories); `off_ER` is defined positively as
> ≥95% spatial purity outside the broad mask, not as a residual category — a
> catch-all residual would introduce a trajectory-length floor that confounds
> the diffusion comparison. Trajectories crossing the mask boundary are
> classified as intermediate and enter neither zone. Because the binomial test
> has limited power in fields of view where the mask covers a large fraction of
> the cell, the smallest testable trajectory length was computed per field of
> view.

> For zone contrasts the two zones were fit jointly, so that on-ER and off-ER
> sub-trajectories share one posterior and are directly comparable. Because
> every field of view contains both zones, the field-of-view relabelling null
> used elsewhere cannot generate a zone contrast; the null here is instead a
> within-field-of-view, trajectory-level permutation of zone labels that
> preserves each field of view's zone counts and moves all of a trajectory's
> sub-trajectories together. Permuting at the sub-trajectory level is
> anticonservative and was shown by a same-vs-same negative control to produce
> false positives.

**Deployed parameters and provenance.** The numbers above are the deployed
invocation, not the function defaults: `rezone_array.sh` runs
`rerun_er_zoning.py` against `all_data/2025_09_15_zeb_ER/preprocessed_corrected`
with `--pct-floor 60 --require-broad false --min-frames 2`, leaving
`pct_above = 1.5`, `block_fine = 25`, `block_broad = 101`,
`min_object_size = 500`, `pct_above_strict = None` and
`enrichment_alpha = 0.05` at their defaults
(`ER_SPT/2026_05_26/er_analysis_lib.py`, `three_class_mask` and
`process_er_for_segmentation`).

Two of those differ from the module defaults and are worth stating, because
both change what `on_structure` means. `--pct-floor 60` raises the global floor
from 50; `--require-broad false` removes the requirement that an `on_structure`
pixel also pass the broad threshold, so the global floor rather than a broad
spatial gate is what keeps isolated bright pixels in dark regions out of the
mask. With `require_broad = True` the mask reads "on the ER has to also be near
the ER", which is the more conservative reading but excludes genuine peripheral
tubules.

A Noise2Void denoising model was trained on these ER images
(`all_data/2025_09_15_zeb_ER/n2v_model_er_hilo/`, outputs in `denoised_er/`)
but is **not** part of this pipeline: no deployed script references either path,
and the masks are built from `preprocessed_corrected/`. It should not be
described in the methods.

---

## Verification summary

`verification/saspt_equivalence_report.json` — 12 tests, all passing, on two
datasets against the installed saspt 0.4.0 (verified pristine against its
dist-info RECORD hashes; a normal pip wheel, not an editable or forked install).

| Check | Result |
|---|---|
| Variational loop vs `StateArray.infer_posterior`, same likelihood | max abs occupancy difference **2.0e-15** |
| saSPT's own distance from its fixed point at its default `max_iter=200` | 5.5e-02 — i.e. ~838× larger than any difference between the two implementations |
| Defocalization correction vs an independent `saspt.defoc.defoc_corr` call | **bitwise identical** |
| Focal-slab survival curve, likelihood rebuilt from cache vs inside a `StateArray` | **bitwise identical** |
| Correction leaves the converged `exp_log_tau` untouched | **bitwise identical** — this is what makes every per-bin panel focal-depth-independent |
| Localization-error marginalization vs saSPT's `marginalize_on_diff_coef` | 1.7e-14 |
| End to end: stock saSPT `posterior_occs` vs the route used here | W₁ = 4.1e-05 |

For scale on the last row: the biological effect in Fig 5C is W₁ ≈ 0.22, so the
implementation difference is four orders of magnitude below the signal.

`verification/weighting_verdict_table.csv` — 96 bin-level comparisons between
saSPT's displacement weighting and the equal weighting used previously, across
7 panels. Four calls change, three of them in the dose-series interaction.

One number worth keeping in view: at the inherited `focal_depth = 0.1 µm` the
defocalization correction moves the population posterior by W₁ = 0.65 —
larger than any biological effect in these figures. That is the reason the
focal depth had to be pinned down rather than left at a placeholder.

---

## Points the draft text needs reconciled

Ordered by how much they change the text. Items 1–4 arise from re-running every
statistic with a variational refit per resample; the point estimates did not
move (they reproduce the previous values to <1e-5) but the intervals did.

1. **"a significant shift in the diffusion rates of the 60S ribosome, but not
   the 40S" (Fig 5C / S7A–B) is not supported.** The 40S shift is significant
   against its own field-of-view permutation null (W₁ = 0.181, p at the 0.002
   floor) and has two significant bins of its own (`immobile` −0.065,
   `free` +0.046). Against W₁(60S) = 0.219 the ratio is **1.2-fold**, not the
   ~3-fold the specificity-control framing implies. What the data support is
   that both subunits shift and that the *shapes* differ — the 60S moves mass
   into the intermediate band, the 40S out of the immobile floor. Recommend
   rewording to that contrast and dropping the null claim for the 40S.

2. **Fig S7C: the immobile bin no longer reaches significance.** −0.028
   (BH p = 0.11, interval spanning zero), previously ** on a fixed-posterior
   bootstrap. The `slow_medium` (+0.130) and `free_60S` (−0.121) transfer that
   carries the claim is unaffected, and `very_slow` (+0.036) remains
   significant. Only the immobile bin needs removing from the text.

3. **Fig S7J: the total distribution shift is no longer significant by W₁.**
   (S7J itself shows no statistic, and the refit behind these numbers is not
   part of this repository.)
   0.121 against a null 95th percentile of 0.130, p = 0.074 (previously 0.032).
   The localized transfer is still significant per bin (`slow_medium` +0.094,
   `free_60S` −0.058), so the defensible statement is that the change is
   confined to specific states rather than a global shift detectable at
   n = 10/14 fields of view. This is consistent with the draft's existing
   hedged wording on this panel.

4. **Fig S7E: emetine corroborates fewer bins than stated.** All six testable
   bin–genotype pairs agree in direction, but only four have *both* drug
   responses individually significant: every harringtonine prediction holds,
   while emetine's does not reach significance for the immobile bin in either
   genotype or for `very_slow` in sgUBA5. "Both drugs confirm each assignment"
   is too strong for bins I and II; the assignment remains corroborated.

5. **Fig S7D supports a stronger claim than currently made.** The UBA5
   genotype effect is *abolished by harringtonine*: the mass transfer out of
   `immobile` into `slow_medium` is present untreated (+0.133) and under
   emetine (+0.115) but absent under harringtonine (−0.038, n.s.). Emetine
   traps elongating ribosomes and preserves the effect; harringtonine runs them
   off and erases it. This is direct evidence the UBA5 phenotype requires
   active translation rather than being a property of free subunits.

6. **Resolved (2026-10-25): Fig S7A and S7C use independent fits.** The placed
   panels come from the refit analysis, which fits each condition independently
   (W₁ = 0.219 for the 60S, 0.181 for the 40S), as every other panel does. An
   earlier no-refit S7A was built from a *joint* fit of the two conditions,
   which gives 0.161 and 0.115 (a 60S/40S ratio of 1.40× rather than 1.21×):
   a joint fit shares one Dirichlet prior across conditions and so shrinks the
   two posteriors toward each other. That joint-fit code is not part of this
   repository.

7. **Fig S7D legend says "n = 20 FOV per condition."** sgNT/emetine has **19**.

8. **The TLSD1 dose-series panels (old 7E, S8E, S8F) are not part of this
   figure set**, and their code is not included in this repository.
