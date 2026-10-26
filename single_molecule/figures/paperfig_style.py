"""
paperfig_style.py -- shared look-and-feel for the SPT panels of Figures 5 and S7.

Everything here is measured from `ZebRibosomeTurnover_draft2 260912.pdf`, not
invented: the palette is sampled from the placed figure images, and the panel
sizes are fractions of the 468 pt (6.5 in) text width every figure is placed at.

Why this module exists
----------------------
The pre-existing panels were emitted at three different figure sizes
(8x4.5, 7x4, 8.6x4.5 in) with per-script colors, so they could not be dropped
into the layout without rescaling -- which changes effective font size and
line weight panel-by-panel. Building every panel through `figure()` here makes
them mutually consistent and drop-in at 100%.

Illustrator compatibility
-------------------------
`pdf.fonttype = 42` embeds TrueType outlines with the text still live, so type
can be edited and restyled after placement. (The matplotlib default, Type 3,
rasterizes glyph programs and Illustrator cannot edit it.) Nothing is
rasterized and `bbox_inches='tight'` is deliberately NOT used, so the saved
media box equals the requested figsize exactly and panels line up on a grid.
"""

import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Palette -- sampled at 300 dpi from the placed figures in draft2 260912
# ---------------------------------------------------------------------------
CONTROL = '#6D6E71'       # sgNT / DMSO / untreated (Fig 5C, 5F)
CONTROL_LIGHT = '#A5A7AA'  # second control series (e.g. DMSO off-ER)
PERTURBED = '#D96227'     # sgUBA5 / PT-179 (Fig 5C, 5F-right)
PERTURBED_LIGHT = '#F7941D'  # second perturbed series (e.g. PT-179 off-ER)
PERTURBED_DARK = '#8C3A14'   # third perturbed shade, so the sgUBA5 family has
                             # three steps to match the grey family's
                             # BLACK / CONTROL / CONTROL_LIGHT (Fig S7E)
BLUE = '#23AAE1'          # emetine (Fig S7D drug panels)
BLUE_DARK = '#0F75BC'
BLACK = '#231F20'

# Semantic aliases, so panel scripts never hard-code a hex value.
COLORS = {
    'sgNT': CONTROL, 'control': CONTROL, 'DMSO': CONTROL, 'untreated': CONTROL,
    'sgUBA5': PERTURBED, 'PT-179': PERTURBED, 'perturbed': PERTURBED,
    'sgTLSD1': CONTROL_LIGHT, 'oeTLSD1': PERTURBED,
    'sgTLSD1_sgUBA5': BLUE,
}

# Zone styling for Fig 5F-left: zone is encoded by LINE STYLE, condition by
# color, so the two factors stay visually separable.
ZONE_LINESTYLE = {'on_ER_broad': '-', 'on_structure': '-', 'off_ER': '--'}

ROMAN = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII']

# ---------------------------------------------------------------------------
# Geometry -- fractions of the 6.5 in (468 pt) figure width used in the draft
# ---------------------------------------------------------------------------
FIG_WIDTH_IN = 6.5

SIZES = {
    # density-vs-D overlay, occupies half the figure width (Fig 5C, 5F-left,
    # S7B, S7E, S7J). Aspect ~2.3, close to the ~2.6 the draft places 5C at.
    'overlay': (3.25, 1.42),
    # per-bin effect panel with error bars (Fig 5F-right, S7C, S7D)
    'effect': (2.95, 1.95),
    # small single-statistic panel (Fig S7A)
    'small': (2.10, 1.60),
    # full-width variants when a panel spans the whole text block
    'overlay_wide': (6.50, 1.90),
    'effect_wide': (6.50, 2.10),
    # as effect_wide but taller: a 4-series panel needs a title plus a legend
    # row, which effect_wide's 0.07 bottom margin clips
    'effect_wide_tall': (6.50, 2.60),
}

# Fixed axes rectangles, as (left, right, bottom, top) figure fractions.
#
# Deliberately NOT tight_layout/constrained_layout: those size the axes to fit
# whatever decorations a particular panel happens to carry, so two panels built
# at the same nominal size end up with different plot boxes and will not line
# up when placed side by side. Fixing the margins per size class means every
# panel of a class has a byte-identical axes rectangle, and the bin-label row
# above the axes always has the same reserved headroom.
MARGINS = {
    'overlay': (0.175, 0.985, 0.255, 0.865),
    'effect': (0.230, 0.985, 0.150, 0.845),
    'small': (0.235, 0.985, 0.110, 0.845),
    'overlay_wide': (0.090, 0.992, 0.215, 0.885),
    'effect_wide': (0.110, 0.992, 0.070, 0.885),
    'effect_wide_tall': (0.110, 0.992, 0.150, 0.855),
}

# Type sizes are set for placement at 100% inside a 6.5 in figure.
RC = {
    'pdf.fonttype': 42,          # editable text in Illustrator (NOT Type 3)
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'font.size': 7,
    'axes.titlesize': 7,
    'axes.labelsize': 7,
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5,
    'legend.fontsize': 6,
    'legend.frameon': False,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'lines.linewidth': 1.0,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.transparent': False,
    'image.composite_image': False,   # keep vector art vector
}


def apply():
    """Install the paper rcParams globally. Idempotent."""
    plt.rcParams.update(RC)


def figure(size='overlay', **kw):
    """A figure at one of the standard panel sizes, with a fixed axes rectangle.

    `size` may be a SIZES key or an explicit (w, h) tuple in inches; an explicit
    tuple falls back to the 'overlay' margins.
    """
    apply()
    key = size if isinstance(size, str) else 'overlay'
    w, h = SIZES[size] if isinstance(size, str) else size
    fig, ax = plt.subplots(figsize=(w, h), **kw)
    left, right, bottom, top = MARGINS[key]
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
    return fig, ax


def save(fig, path, **kw):
    """Save as vector PDF with live text and an exact media box.

    `bbox_inches='tight'` is intentionally not passed: it would shrink the
    media box to the drawn content, so two panels built at the same nominal
    size would land at different sizes in the layout.
    """
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    fig.savefig(path, format='pdf', **kw)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Diffusion-bin decoration (the "Diffusion Bin: I II III ..." header used on
# every density and effect panel in the draft)
# ---------------------------------------------------------------------------

def _bin_row_label(ax, label, y):
    """Left-hand caption for the bin-numeral row.

    Drawn in FIGURE coordinates, flush to the left edge, rather than in axes
    coordinates: the axes box is narrow enough that an axes-relative offset
    large enough to clear the first numeral pushes the caption off the canvas.
    """
    if not label:
        return
    fig = ax.figure
    pos = ax.get_position()
    yb = pos.y1 + (y - 1.0) * pos.height
    # Right-aligned into the y-label gutter, at the height of the numeral row.
    # The y-label is vertically centred, so the top of the gutter is free and
    # the two never collide however long the y-label is.
    fig.text(pos.x0 - 0.012, yb, label, ha='right', va='bottom', fontsize=6.5)


def bin_edges(scheme):
    """Interior boundaries of a region scheme, in ascending D."""
    lo = sorted(v[0] for v in scheme.values())
    return [e for e in lo if e > 0]


def annotate_bins_log(ax, scheme, label='Bin:', y=1.04,
                      line_kw=None, text_kw=None):
    """Dotted vertical boundaries plus roman-numeral bin labels above a
    log-D axis. Matches Fig 5C / 5F-left.

    The caption is "Bin:" rather than the draft's "Diffusion Bin:" purely for
    fit -- the longer string overruns the left edge of a 3.25 in panel -- and it
    matches the bin-effect panels. Pass label='Diffusion Bin:' with a wider
    left margin if the full wording is wanted.

    Call AFTER the x-limits are final -- label positions are the geometric
    midpoints of each bin as clipped to the current view.
    """
    lk = dict(color='0.35', lw=0.5, ls=':', zorder=0)
    lk.update(line_kw or {})
    tk = dict(ha='center', va='bottom', fontsize=6.5,
              transform=ax.get_xaxis_transform())
    tk.update(text_kw or {})

    edges = bin_edges(scheme)
    for e in edges:
        ax.axvline(e, **lk)

    x0, x1 = ax.get_xlim()
    bounds = [x0] + [e for e in edges if x0 < e < x1] + [x1]
    for i in range(len(bounds) - 1):
        mid = np.sqrt(bounds[i] * bounds[i + 1])
        ax.text(mid, y, ROMAN[i], **tk)
    _bin_row_label(ax, label, y)


def annotate_bins_categorical(ax, n_bins, label='Bin:', y=1.04, sep=True):
    """Roman-numeral headers over a categorical per-bin axis (Fig 5F-right,
    S7C, S7D), with light separators between bins."""
    for i in range(n_bins):
        ax.text(i, y, ROMAN[i], ha='center', va='bottom', fontsize=6.5,
                transform=ax.get_xaxis_transform())
    if sep:
        for i in range(n_bins - 1):
            ax.axvline(i + 0.5, color='0.80', lw=0.5, zorder=0)
    _bin_row_label(ax, label, y)


def stars(p):
    """Draft convention: *** p<0.001, ** p<0.01, * p<0.05, else blank."""
    if p is None or not np.isfinite(p):
        return ''
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''


def sig_bracket(ax, x0, x1, y, text, dy=0.012, lw=0.6, color='0.15'):
    """A significance bracket spanning two x positions in data coordinates,
    as used between the two series of a bin in Fig 5F-right."""
    ax.plot([x0, x0, x1, x1], [y, y + dy, y + dy, y], lw=lw, color=color,
            clip_on=False, solid_joinstyle='miter')
    ax.text((x0 + x1) / 2, y + dy, text, ha='center', va='bottom',
            fontsize=7, color=color, clip_on=False)


XLABEL_D = r'Estimated Diffusion Coefficient (µm$^2$/sec)'

# Labels follow the draft's own vocabulary ("Probability Density" on the
# overlays, "Probability Fraction" on the bin panels), with the measure stated
# in the legend and Methods rather than on the axis. saSPT itself calls the
# quantity "occupation" (saspt/plot.py:231, 294-295) and never labels an axis
# "density"; either vocabulary is defensible, but a bare "probability density"
# does need the measure spelled out somewhere, because the x-axis is labelled in
# D while the density is with respect to log10 D -- the two differ by D*ln(10).
# See LEGEND_DENSITY_NOTE below for the sentence to paste.
#
# Two representations of the SAME quantity, and they need different labels:
#
#   overlay panels  the occupation mass p_k divided by the grid spacing
#                   Delta log10 D, i.e. a probability DENSITY with respect to
#                   log10 D. Integrates to 1 over the log axis, so area under an
#                   interval is that interval's occupancy. Unlike the raw mass,
#                   its height does not depend on how many grid states were
#                   chosen -- on this 150-state grid a peak of height 1.14
#                   decade^-1 would read as 0.038 raw, and as 0.057 on saSPT's
#                   default 100-state grid, for identical biology.
#
#   bin panels      the same mass summed over a bin: a dimensionless FRACTION of
#                   posterior occupancy, NOT a density. Labelling these
#                   "probability density" is simply wrong.
#
# The measure matters: this is a density in log10 D, not in D. A density in D
# would differ by the Jacobian D*ln(10). See
# PaperFigures/METHODS_posterior_density_note.md.
YLABEL_DENSITY = 'Probability density\n(per unit log$_{10}$D)'
# Two lines, not one: the y-label is rotated, so its length is bounded by the
# 1.42 in panel HEIGHT, not the width. Every single-line form carrying the
# measure -- 'per unit log10D', ', (unit log10D)^-1' -- renders clipped at the
# canvas edge; the two-line form fits with room to spare. Verified by
# rasterizing each candidate and testing for ink on the border.
# The qualifier is on the axis rather than left to the legend because a bare
# 'probability density' invites the reasonable objection that it should sum to
# 1 -- it integrates to 1 over log10 D, and sums to 1/dlog10D = 29.8.
# The bin panels are a FRACTION, not a density -- no qualifier rescues "density"
# for an integrated mass. This is the draft's own Fig 5F-right wording.
YLABEL_DELTA = r'$\Delta$Probability fraction'

# Paste into any legend whose panel uses YLABEL_DENSITY.
LEGEND_DENSITY_NOTE = (
    'Density is per unit log10 D, so the area under the curve between two '
    "diffusion coefficients is that interval's occupancy; read by area, not by "
    'peak height.')
