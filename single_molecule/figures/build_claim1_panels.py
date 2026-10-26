"""
build_claim1_panels.py -- Figure 5C and Supplementary Figure S7B.

Dataset `2025_06_06_zeb`, tier `q=15p0_splitsize10`, 15.091 ms / 0.11 um.
RPL29-Halo (60S) and Halo-RPS24 (40S), sgNT vs sgUBA5, 10 FOV per condition.
The 40S arm is a pre-specified specificity control, not a co-equal test.

Panels
------
  5C   60S posterior density overlay, sgNT vs sgUBA5, bin boundaries drawn
  S7B  40S posterior density overlay

Each condition is fitted independently, at saSPT's own settings (displacement
weighting, defocalization on) -- see panel_core. The statistics panels of the
same contrasts, S7A (W1 against the FOV-relabelling null) and S7C (per-bin
60S difference), refit saSPT on every resampling draw and are drawn by the
refit chain: run_panel_refit_batch.py -> aggregate_panel_refit.py ->
replot_refit_panels.py.

Run:  sbatch call_build_panels.slurm build_claim1_panels.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_core as pc
import paperfig_style as ps
from region_schemes import SCHEME_60S, SCHEME_40S

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from common import paths as tlsd_paths  # noqa: E402  repo root; data/output roots
TIER = tlsd_paths.spt('all_data', '2025_06_06_zeb', 'tracking_output_q=15p0_splitsize10')
OUT = tlsd_paths.outp('PaperFigures_updated')

SUBUNITS = {
    '60S': dict(cache='occupancy_compare_cache/60S_bothgeno_likelihood_cache.npz',
                a='sgNT_60S', b='sgUBA5_60S', scheme=SCHEME_60S,
                label='60S', construct='RPL29-HaloTag'),
    '40S': dict(cache='occupancy_compare_cache/40S_bothgeno_likelihood_cache.npz',
                a='sgNT_40S', b='sgUBA5_40S', scheme=SCHEME_40S,
                label='40S', construct='HaloTag-RPS24'),
}


def fit_subunit(key):
    cfg = SUBUNITS[key]
    # condition_ctx hides the layout difference: a merged cache with a
    # source_member filter by default, one cache per condition under a variant.
    ctx, ga = pc.condition_ctx(TIER, cfg['a'], cfg['cache'])
    ctx_b, gb = pc.condition_ctx(TIER, cfg['b'], cfg['cache'])
    pmf_a, ns_a = pc.population_pmf(ctx, ga)
    pmf_b, ns_b = pc.population_pmf(ctx_b, gb)
    return dict(ctx=ctx, cfg=cfg, pmf_a=pmf_a, pmf_b=pmf_b,
                n_sub_a=ns_a, n_sub_b=ns_b)


def overlay(fit, out_pdf):
    """Density overlay, both conditions filled in their own colour.

    Previously only the perturbed curve was filled, on the argument that it made
    the two separable in greyscale. It also made the panel read as unfinished,
    and it undercut the y-axis: this is a density to be read by area, and a
    single fill draws the eye to peak heights instead. Greyscale separation is
    carried by the strokes, which differ in value (CONTROL #6D6E71 vs PERTURBED
    #D96227), not only in hue."""
    c = fit['cfg']
    pc.plot_overlay(
        [dict(pmf=fit['pmf_a'], label=f"sgNT (n={fit['n_sub_a']:,})",
              color=ps.CONTROL, fill=True),
         dict(pmf=fit['pmf_b'], label=f"sgUBA5 (n={fit['n_sub_b']:,})",
              color=ps.PERTURBED, fill=True)],
        fit['ctx']['diff_coefs'], c['scheme'], out_pdf,
        xlabel=f"Estimated Diffusion Coefficient {c['label']} "
               r"(µm$^2$/sec)")
    return pc.tag(out_pdf)


def main():
    ps.apply()
    manifest = {}
    for key in ('60S', '40S'):
        print(f'fitting {key} ...', flush=True)
        fit = fit_subunit(key)
        if key == '60S':
            manifest['fig5C'] = overlay(fit, f'{OUT}/figure5/C/fig5C_60S_posterior_overlay.pdf')
        else:
            manifest['figS7B'] = overlay(fit, f'{OUT}/figureS7/B/figS7B_40S_posterior_overlay.pdf')
    for k, v in manifest.items():
        print(f'  {k}: {v}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
