"""
Shared, data-source-agnostic plotting utilities for the ROC plotting scripts.

Everything here is used identically (or near-identically) by both the real-data
pipeline (roc_real_data.py) and the simulated-data pipeline (roc_sim_data.py),
and by any script that draws from either — currently plot_main_figures.py, and
in the future plot_supplemental_figures.py. Nothing in this module knows about
call-file formats, sample-name parsing, or metadata CSVs — that lives in
roc_real_data.py / roc_sim_data.py instead.

Extracted from plot_roc.py / plot_roc_simulated.py, where this code was
duplicated near-verbatim across both files (see those scripts' git history /
CLAUDE.md for the original derivation notes this carries forward).
"""

import math
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

METHODS = ["SA", "PR", "LL", "PP"]

GENOME_LENGTH = 29903  # SARS-CoV-2 reference length; denominator for call_rate
                        # over all positions, real and simulated data alike.

# Expressed as linear probabilities (used for display/labeling via _fmt_threshold).
# PP_val in the call files is log10(probability) as of the "_precise" file format
# (see THRESHOLDS_LOG10 below) — comparisons against PP_val must use that list, not
# this one, since deltas below ~1e-15 cannot be represented as distinct linear floats.
# Geometric grid in delta = 1-t: 1x/5x steps per decade from delta=1e-15 (the
# float64/log10 precision floor) down through delta=1e-1, plus t=0.0 kept as an
# explicit sentinel ("always passes" — see find_pp_sa_crossover / _log10_threshold
# in roc_real_data.py / roc_sim_data.py). Density falls off naturally toward t=0
# since only one decade of delta separates t=0.9 from t=0.0, vs. ten decades
# between t=0.999999999999999 and t=0.99999.
THRESHOLDS = []
for _n in range(15, 0, -1):
    THRESHOLDS.append(1 - 1 * 10 ** (-_n))
    THRESHOLDS.append(1 - 5 * 10 ** (-_n))
THRESHOLDS.append(0.0)
N_THRESHOLDS = len(THRESHOLDS)


def _log10_threshold(t):
    """Precise log10(t) for a threshold expressed as t = 1 - delta, where delta may
    be far smaller than float64 can resolve near 1.0 (log1p avoids the cancellation
    that a plain math.log10(t) would suffer for t this close to 1). t == 0.0 (the
    "always passes" sentinel) maps to -inf."""
    if t <= 0.0:
        return float("-inf")
    return math.log1p(t - 1.0) / math.log(10.0)


# Log10-space equivalents of THRESHOLDS, for comparison against PP_val (which is
# already log10(probability) in the "_precise" call files, real or simulated).
# Computed via log1p so the extra precision the log10 storage format provides
# isn't thrown away by round-tripping through a linear float first.
THRESHOLDS_LOG10 = [_log10_threshold(t) for t in THRESHOLDS]

_PP_MARKER = "D"  # PP threshold points are always diamonds, everywhere.

# Fixed real-world (inch) size for gradient-bar legends (used by the sa_pp
# family in both real and simulated pipelines), so they render identically
# small regardless of the enclosing figure's own size — this is what lets
# e.g. sa_pp_1v5.pdf (4in tall) and sim_sa_pp_1v5.pdf (~9.5in tall, two
# stacked error_rate blocks) end up with matching gradient bars despite
# their very different overall figure heights.
_GRADIENT_BAR_WIDTH_IN  = 0.22
_GRADIENT_BAR_HEIGHT_IN = 0.55
_GRADIENT_BAR_GAP_IN    = 0.25


# ---------------------------------------------------------------------------
# Accumulator helpers
# ---------------------------------------------------------------------------

def _new_method_acc():
    return [0.0, 0.0, 0]


def _new_threshold_acc():
    return [[0.0, 0.0, 0] for _ in range(N_THRESHOLDS)]


def _add_score(bucket, acc, cr):
    bucket[0] += acc
    bucket[1] += cr
    bucket[2] += 1


def _safe_mean(values):
    return sum(values) / len(values) if values else float("nan")


# ---------------------------------------------------------------------------
# Formatting / color helpers
# ---------------------------------------------------------------------------

def _fmt_threshold(t):
    """Format a threshold value for diamond labels.
    delta = 1-t < 0.001 → '1-Xe-Y' notation; otherwise plain '{t:g}'."""
    delta = 1.0 - t
    if delta < 0.001:
        s = f"{delta:.2e}"
        coeff, exp = s.split("e")
        coeff = coeff.rstrip("0").rstrip(".")
        return f"1-{coeff}e{int(exp)}"
    return f"{t:g}"


def _lerp_to_black(color, fraction):
    """Interpolate an RGB color toward black. fraction=0 → color, fraction=1 → black."""
    r, g, b = color[:3]
    t = 1.0 - fraction
    return (r * t, g * t, b * t)


# ---------------------------------------------------------------------------
# Gradient-bar legend (sa_pp family)
# ---------------------------------------------------------------------------

def _draw_pp_gradient_bars(make_axes, condition_colors, x0, y_top, bar_w, bar_h, gap, lfs):
    """Draw one vertical gradient colorbar-style strip per condition, stacked
    downward from y_top: top = condition color (least stringent PP
    threshold, i.e. threshold=0), bottom = black (most stringent, i.e.
    threshold=THRESHOLDS[0]) — the same color-to-black gradient already used
    for the PP threshold points themselves, so each condition's bar
    auto-matches its curve/point color rather than needing a hand-labeled
    "low"/"high" proxy legend entry. Tick labels give the actual threshold
    values via _fmt_threshold, not a "least/most stringent" description.
    `make_axes(rect)` creates each bar's Axes — pass fig.add_axes for a
    figure-level placement (rect in figure fractions) or ax.inset_axes for a
    placement anchored to one panel (rect in axes fractions), so this one
    drawing routine serves both the shared-legend and per-panel-legend cases.
    bar_w/bar_h/gap are already-converted fractions (see
    _place_pp_gradient_bars, which computes them from the fixed inch sizes
    above)."""
    top_label = f"threshold = {_fmt_threshold(THRESHOLDS[-1])}"
    bottom_label = f"threshold = {_fmt_threshold(THRESHOLDS[0])}"
    y = y_top
    for condition, color in condition_colors.items():
        cax = make_axes([x0, y - bar_h, bar_w, bar_h])
        cmap = LinearSegmentedColormap.from_list(f"{condition}_pp_grad", [color, "black"])
        gradient = np.linspace(0, 1, 256).reshape(-1, 1)
        cax.imshow(gradient, aspect="auto", cmap=cmap, origin="upper")
        cax.set_xticks([])
        cax.yaxis.tick_right()
        cax.set_yticks([0, 255])
        cax.set_yticklabels([top_label, bottom_label], fontsize=max(lfs - 3, 5))
        cax.set_title(condition, fontsize=max(lfs - 2, 6), pad=3)
        y -= bar_h + gap
    return y


def _place_pp_gradient_bars(fig, legend, condition_colors, lfs, panel_axes, ref_ax=None):
    """Draw the sa_pp-family gradient bars (see _draw_pp_gradient_bars),
    horizontally centered under `legend`'s color swatches (the SA/PP marker
    and line handles, not the legend box's own left edge — the box also
    spans the text labels to the right of the swatches, so aligning to its
    raw left edge left the bars sitting a bit left of where the swatches
    actually are), vertically centered on the midpoint of the topmost panel
    as a whole. `panel_axes` is the
    sequence of Axes making up that one panel — for a broken-y-axis panel
    (ax_top, ax_bottom) that's both halves, whose *union* bounding box is
    what "the panel" means (centering on ax_top alone would land the
    gradient a full half-panel too high, at the midpoint of just the top
    half); for an unbroken single-axis panel it's a single-element sequence
    containing that one axis. Rather than positioned relative to the
    legend's own height, which made the offset hard to predict since it
    moved whenever the legend's entry count changed. Bar size is fixed in
    real inches (_GRADIENT_BAR_*_IN) and converted to the target's own
    fraction units using its actual rendered size, so the bars come out the
    same absolute size whether anchored to the whole figure (ref_ax=None) or
    to one panel Axes (ref_ax=that panel) — and the same absolute size
    across figures despite very different overall figure heights."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    panel_bboxes_px = [ax.get_window_extent(renderer) for ax in panel_axes]
    panel_y0_px = min(b.y0 for b in panel_bboxes_px)
    panel_y1_px = max(b.y1 for b in panel_bboxes_px)
    panel_y_center_px = (panel_y0_px + panel_y1_px) / 2

    # Center under the legend's color swatches specifically (the first
    # handle with a valid rendered box) -- matplotlib aligns every entry's
    # handle to the same left column, so any one handle's x-extent gives
    # that column's true position regardless of which entry (SA or a PP
    # condition) it belongs to. A scatter-based handle's legend proxy
    # (e.g. SA's) reports an inf/-inf bbox rather than a real one, so it's
    # skipped in favor of the first handle (typically a PP entry's Line2D)
    # that reports a finite box.
    handles = getattr(legend, "legend_handles", None) or legend.legendHandles
    legend_bbox_px = legend.get_window_extent(renderer)
    swatch_bbox_px = next(
        (b for h in handles for b in [h.get_window_extent(renderer)] if np.isfinite(b.x0)),
        None,
    )
    if swatch_bbox_px is not None:
        x0_px = (swatch_bbox_px.x0 + swatch_bbox_px.x1) / 2 - (_GRADIENT_BAR_WIDTH_IN * fig.dpi) / 2
    else:
        x0_px = legend_bbox_px.x0

    if ref_ax is None:
        fig_bbox_px = fig.get_window_extent(renderer)
        ref_w_in = fig_bbox_px.width / fig.dpi
        ref_h_in = fig_bbox_px.height / fig.dpi
        x0 = fig.transFigure.inverted().transform((x0_px, 0))[0]
        y_center = fig.transFigure.inverted().transform((0, panel_y_center_px))[1]
        make_axes = fig.add_axes
    else:
        ax_bbox_px = ref_ax.get_window_extent(renderer)
        ref_w_in = ax_bbox_px.width / fig.dpi
        ref_h_in = ax_bbox_px.height / fig.dpi
        x0 = ref_ax.transAxes.inverted().transform((x0_px, 0))[0]
        y_center = ref_ax.transAxes.inverted().transform((0, panel_y_center_px))[1]
        make_axes = ref_ax.inset_axes

    bar_w = _GRADIENT_BAR_WIDTH_IN / ref_w_in
    bar_h = _GRADIENT_BAR_HEIGHT_IN / ref_h_in
    gap = _GRADIENT_BAR_GAP_IN / ref_h_in

    n_bars = len(condition_colors)
    block_height = n_bars * bar_h + (n_bars - 1) * gap
    y_top = y_center + block_height / 2

    _draw_pp_gradient_bars(make_axes, condition_colors,
                            x0=x0, y_top=y_top,
                            bar_w=bar_w, bar_h=bar_h, gap=gap, lfs=lfs)


# ---------------------------------------------------------------------------
# Broken-y-axis limit helpers
# ---------------------------------------------------------------------------

def _find_ylim_break(values, pad_frac=0.08, min_gap_frac=0.15, break_buffer_frac=0.2):
    """Find the largest gap in a list of y-values and return
    ((bottom_lo, bottom_hi), (top_lo, top_hi)) for a broken y-axis, or None if
    no gap is large enough (>= min_gap_frac of the full value range) to be worth
    breaking on."""
    vals = sorted(v for v in values if v == v)  # drop NaN
    if len(vals) < 2:
        return None
    total_range = vals[-1] - vals[0]
    if total_range <= 0:
        return None
    gap, idx = max((vals[i + 1] - vals[i], i) for i in range(len(vals) - 1))
    if gap / total_range < min_gap_frac:
        return None
    low_max, high_min = vals[idx], vals[idx + 1]
    bottom_pad = max((low_max - vals[0]) * pad_frac, gap * 0.05)
    top_pad = max((vals[-1] - high_min) * pad_frac, gap * 0.05)
    break_buffer = gap * break_buffer_frac
    return ((vals[0] - bottom_pad, low_max + break_buffer),
            (high_min - break_buffer, vals[-1] + top_pad))


def _find_ylim_by_group(bottom_values, top_values, pad_frac=0.03, break_buffer_frac=0.08):
    """Like _find_ylim_break, but the bottom/top split is given directly (e.g. by
    condition/method identity) rather than auto-detected from the largest gap —
    guarantees every bottom_values point lands in the bottom half and every
    top_values point lands in the top half, regardless of how close the two
    groups' ranges get. Returns None if the groups overlap (bottom's max >=
    top's min), since no clean break exists between them."""
    bottom_vals = sorted(v for v in bottom_values if v == v)
    top_vals = sorted(v for v in top_values if v == v)
    if not bottom_vals or not top_vals:
        return None
    bottom_min, bottom_max = bottom_vals[0], bottom_vals[-1]
    top_min, top_max = top_vals[0], top_vals[-1]
    if bottom_max >= top_min:
        return None
    gap = top_min - bottom_max
    bottom_pad = max((bottom_max - bottom_min) * pad_frac, gap * 0.05)
    top_pad = max((top_max - top_min) * pad_frac, gap * 0.05)
    break_buffer = gap * break_buffer_frac
    return ((bottom_min - bottom_pad, bottom_max + break_buffer),
            (top_min - break_buffer, top_max + top_pad))
