"""
Main-figure ROC plots: the subset of plot_roc.py / plot_roc_simulated.py's
figures currently in active use, combined into one script. Draws from the
shared roc_real_data.py / roc_sim_data.py (data pipelines) and plot_common.py
(generic plotting utilities) modules rather than duplicating them — those
modules also back plot_supp_figures.py, which additionally imports this
module directly to reuse its _make_row_grid/_add_row_labels layout helpers
and its two pickle caches (REAL_CACHE_PATH, SIMS_CACHE_PATH).

Produces 6 plot files plus one CSV, all written to OUTPUT_DIR
("main_figures/"): sa_pp_1v5.pdf, sim_sa_pp_1v5.pdf,
sample_metadata_1v5.pdf, quality_flag_1v5.pdf, gene_region_1v5.pdf,
replicate_acc_diff.pdf, pp_sa_crossover_callrate.csv (the call-rate-based
crossover -- see roc_real_data.find_pp_sa_crossover; its accuracy-based
counterpart, pp_sa_crossover_accuracy.csv, is produced by
plot_supp_figures.py instead).

See plot_roc.py / plot_roc_simulated.py (left untouched) for the full history
and every other plot function these were extracted from — including
gene_acc_diff_kde, dropped from this script since replicate_acc_diff now
covers essentially the same structural/nonstructural-by-center comparison;
gene_acc_diff_kde is still fully defined and documented in plot_roc.py if
needed again.

Caching: three independent pickle caches, one per data source, each with its
own rebuild flag — since each data source feeds a disjoint subset of the
figures above and reading the underlying call files is slow:
  --rebuild-real : re-reads the three main subsampled call files
                   (CONDITIONS in roc_real_data.py) -> sa_pp_1v5,
                   sample_metadata_1v5, quality_flag_1v5, gene_region_1v5,
                   and (via mle_mode_t_idx/all_genes) replicate_acc_diff.
  --rebuild-reps : re-reads the sub20_only replicate call files
                   -> replicate_acc_diff.
  --rebuild-sims : re-reads the three simulated call files
                   (CONDITIONS in roc_sim_data.py) -> sim_sa_pp_1v5.
  --rebuild      : shorthand for all three of the above at once.
Pass any combination of the three individual flags to rebuild just those.
Metadata CSVs (sample/filter/gene) are always reloaded fresh — they're small.
"""

import os
import pickle
import random
import re
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.lines import Line2D

import roc_real_data as real
import roc_sim_data as sim
from plot_common import (
    N_THRESHOLDS,
    _PP_MARKER,
    _lerp_to_black,
    _place_pp_gradient_bars,
    _find_ylim_break,
    _find_ylim_by_group,
)

# ---------------------------------------------------------------------------
# Visual-encoding constants specific to this figure set
# ---------------------------------------------------------------------------

# Shared sa_pp-family panel sizing, used by both sa_pp_1v5 (real) and
# sim_sa_pp_1v5, so their panels stay the same absolute size without either
# script guessing at the other's literal values. Square (width == height) —
# for sa_pp_1v5 "the panel" is the union of its broken-axis top+bottom
# halves, matching sim_sa_pp_1v5's single unbroken panel.
_SA_PP_PANEL_WIDTH  = 3.0
_SA_PP_PANEL_HEIGHT = 3.0

# Explicit inch-based margins for the sa_pp family (sa_pp_1v5, sim_sa_pp_1v5),
# shared so both figures reserve the exact same fraction of their canvas for
# the legend/title/xlabel — the only way to guarantee their rendered panels
# come out the same absolute size. Plain tight_layout() cannot guarantee
# this: it recomputes margins per-figure from that figure's own content, so
# two figures with nominally-matching canvas dimensions can still render
# their actual plot areas at different sizes (as sa_pp_1v5 used to, before
# adopting this same explicit-margin technique from sim_sa_pp_1v5).
_SA_PP_GRIDSPEC_LEFT   = 0.08
_SA_PP_GRIDSPEC_RIGHT  = 0.80
_SA_PP_TOP_MARGIN_IN   = 0.4   # above the (top) row's own column titles
_SA_PP_BOTTOM_MARGIN_IN = 0.45  # below the (bottom) row's xlabel
_SA_PP_BLOCK_GAP_IN    = 0.7   # between stacked blocks (sim_sa_pp_1v5 only)
_SA_PP_LEGEND_ANCHOR_X = 0.83  # just past _SA_PP_GRIDSPEC_RIGHT
_SA_PP_WSPACE          = 0.24  # fraction of panel width, gap between columns

# replicate_acc_diff: marker shapes cycled across replicate IDs (rep1, rep2,
# ...), independent of CenterName's color encoding. 10 distinct shapes
# covers the sub20_only replicate count seen in practice; if more replicates
# ever appear, shapes repeat (color still disambiguates by CenterName).
_REPLICATE_MARKERS = ["o", "s", "^", "v", "D", "P", "X", "*", "h", "<", ">", "p"]

REAL_CACHE_PATH = "plot_main_real_cache.pkl"  # --rebuild-real to force re-read
REPS_CACHE_PATH = "plot_main_reps_cache.pkl"  # --rebuild-reps to force re-read
SIMS_CACHE_PATH = "plot_main_sims_cache.pkl"  # --rebuild-sims to force re-read

OUTPUT_DIR = "main_figures"  # all 7 plot PDFs are written here


def _out(filename):
    """Path for a plot output file, under OUTPUT_DIR (created on first use)."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    return os.path.join(OUTPUT_DIR, filename)


# ---------------------------------------------------------------------------
# sa_pp_1v5 (real data)
# ---------------------------------------------------------------------------

def _draw_sa_pp_figure(mean_stats, mean_thr_stats, subsamples, filename, wspace=_SA_PP_WSPACE):
    """Shared drawing routine for the sa_pp family (sa_pp_1v5 here;
    sim_sa_pp_1v5 mirrors this layout for simulated data). SA + PP curves
    excluding pr_scaling=0.001, no PR/LL. One broken-y-axis panel per
    subsample in `subsamples`. All panels share one identity-based break (SA
    + PP(pr_scaling=MLE) always top, PP(pr_scaling=0.0) always bottom — see
    _find_ylim_by_group) rather than an auto-detected gap, so the grouping is
    guaranteed even where one panel's clusters sit closer together than
    another's. Equal top/bottom height.
    Every threshold plotted as a point, gradient condition color (least
    stringent) → black (most stringent), threaded together by a connecting
    line in the condition's base color (plain line, not gradient — the
    gradient lives on the points) — the line makes the curve shape legible at
    a glance, while the diamonds still mark each individual threshold's exact
    position. SA gets its own color (tab10 index 1, otherwise unused here
    since pr_scaling=0.001 is excluded). The legend's PP entries are a
    zero-length line+diamond proxy in the condition's color (so the swatch
    shows both the line and the marker at once), and SA's legend entry is
    just its circle — every panel has exactly one subsample, so no
    subsample-distinguishing marker/style is needed and every panel produces
    identical legend labels, letting one legend (and one set of gradient
    bars below it) cover the whole figure. Rather than a hand-labeled "PP
    threshold (low/high)" proxy legend entry, each condition gets its own
    vertical gradient colorbar-style strip (see plot_common), left-aligned
    under the SA/PP legend but vertically centered on the topmost panel's
    midpoint rather than tied to the legend's own height, labeled with the
    actual threshold values (via _fmt_threshold) rather than "least/most
    stringent" text, and colored to match that condition's points (top =
    condition color/threshold=0, bottom = black/most stringent threshold).

    Sizing uses the same explicit inch-based margins as sim_sa_pp_1v5
    (_SA_PP_GRIDSPEC_LEFT/RIGHT, _SA_PP_TOP/BOTTOM_MARGIN_IN — see those
    constants), rather than tight_layout()'s automatic per-figure spacing,
    specifically so the two figures' rendered panels come out the same
    absolute size — tight_layout cannot guarantee that on its own, since it
    recomputes margins independently per figure from that figure's own
    content."""
    conditions = {k: v for k, v in real.CONDITIONS.items() if k != "pr_scaling=0.001"}
    condition_colors = {c: real.ALL_CONDITION_COLORS[c] for c in conditions}
    first_cond = next(iter(conditions))
    sa_color = cm.tab10.colors[1]  # unused by any plotted condition here

    sc, dia, lfs = 60, 40, 8

    # Precompute per-subsample data once, and a single break shared across every
    # panel (rather than one break per panel), so all panels share the same y-axis.
    # The break is identity-based, not gap-based: pr_scaling=0.0 always goes in
    # the bottom half, SA + pr_scaling=MLE always go in the top half.
    sa_by_subsample = {}
    curves_by_subsample = {}
    bottom_group_values = []  # pr_scaling=0.0
    top_group_values = []     # SA + pr_scaling=MLE
    for subsample in subsamples:
        sa_acc, sa_cr = mean_stats[first_cond][subsample].get("SA", (float("nan"), float("nan")))
        sa_by_subsample[subsample] = (sa_cr, sa_acc)
        if sa_acc == sa_acc:
            top_group_values.append(sa_acc)
        subsample_curves = []
        for condition, color in condition_colors.items():
            curve_cr, curve_acc = [], []
            group_values = bottom_group_values if condition == "pr_scaling=0.0" else top_group_values
            for t_idx in range(N_THRESHOLDS):
                acc, cr = mean_thr_stats[condition][subsample].get(
                    t_idx, (float("nan"), float("nan"))
                )
                curve_cr.append(cr)
                curve_acc.append(acc)
                if acc == acc:
                    group_values.append(acc)
            subsample_curves.append((condition, color, curve_cr, curve_acc))
        curves_by_subsample[subsample] = subsample_curves

    ylims = _find_ylim_by_group(bottom_group_values, top_group_values,
                                 pad_frac=0.03, break_buffer_frac=0.08)
    if ylims is None:
        # pr_scaling=0.0 and SA/MLE overlap somewhere — fall back to a plain
        # largest-gap split, then to one shared range if even that finds nothing.
        all_values = bottom_group_values + top_group_values
        ylims = _find_ylim_break(all_values, pad_frac=0.03, break_buffer_frac=0.08)
    if ylims is None:
        all_values = bottom_group_values + top_group_values
        lo, hi = min(all_values), max(all_values)
        pad = (hi - lo) * 0.03 or 0.01
        bottom_ylim = top_ylim = (lo - pad, hi + pad)
    else:
        bottom_ylim, top_ylim = ylims

    # fig_width/fig_height scale from the fixed per-panel dimensions and the
    # explicit margins (see docstring), so every sa_pp-family figure keeps
    # the same rendered panel size regardless of panel count.
    n_panels = len(subsamples)
    usable_w_frac = _SA_PP_GRIDSPEC_RIGHT - _SA_PP_GRIDSPEC_LEFT
    fig_width = _SA_PP_PANEL_WIDTH * (n_panels + (n_panels - 1) * wspace) / usable_w_frac
    fig_height = _SA_PP_TOP_MARGIN_IN + _SA_PP_PANEL_HEIGHT + _SA_PP_BOTTOM_MARGIN_IN
    top_frac = (fig_height - _SA_PP_TOP_MARGIN_IN) / fig_height
    bottom_frac = _SA_PP_BOTTOM_MARGIN_IN / fig_height

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(2, n_panels, height_ratios=[1, 1], hspace=0.08, wspace=wspace,
                           left=_SA_PP_GRIDSPEC_LEFT, right=_SA_PP_GRIDSPEC_RIGHT,
                           top=top_frac, bottom=bottom_frac)

    for col, subsample in enumerate(subsamples):
        ax_top = fig.add_subplot(gs[0, col])
        ax_bottom = fig.add_subplot(gs[1, col], sharex=ax_top)

        for ax in (ax_top, ax_bottom):
            sa_cr, sa_acc = sa_by_subsample[subsample]
            ax.scatter(sa_cr, sa_acc, marker="o", color=sa_color, s=sc, zorder=5, label="SA")

            for condition, color, curve_cr, curve_acc in curves_by_subsample[subsample]:
                pp_label = f"PP ({condition})"
                ax.plot(curve_cr, curve_acc, color=color, linewidth=1.5, zorder=4)
                ax.plot([], [], color=color, linewidth=1.5,
                        marker=_PP_MARKER, markersize=6, label=pp_label)
                for t_idx, (cr, acc) in enumerate(zip(curve_cr, curve_acc)):
                    fraction = 1 - t_idx / (N_THRESHOLDS - 1)
                    pt_color = _lerp_to_black(color, fraction)
                    ax.scatter(cr, acc, color=pt_color, marker=_PP_MARKER, s=dia, zorder=6)

        ax_top.set_ylim(*top_ylim)
        ax_bottom.set_ylim(*bottom_ylim)
        ax_top.set_xlim(0, 1.05)
        ax_bottom.set_xlim(0, 1.05)
        ax_top.ticklabel_format(useOffset=False)
        ax_bottom.ticklabel_format(useOffset=False)

        # Keep the boundary spines visible (bottom of top axis, top of bottom
        # axis) so each half is framed right up to the break, making the
        # discontinuity and each cluster's extent more visually explicit.
        ax_top.tick_params(labelbottom=False, bottom=False, labelsize=lfs - 1)
        ax_bottom.tick_params(labelsize=lfs - 1)

        # Diagonal break marks at the shared boundary
        d = 0.5
        break_kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12, linestyle="none",
                             color="k", mec="k", mew=1, clip_on=False)
        ax_top.plot([0, 1], [0, 0], transform=ax_top.transAxes, **break_kwargs)
        ax_bottom.plot([0, 1], [1, 1], transform=ax_bottom.transAxes, **break_kwargs)

        ax_top.set_title(f"mean depth = {subsample[3:]}x", fontsize=lfs + 1)
        ax_bottom.set_xlabel("Mean call rate", fontsize=lfs)
        if col == 0:
            ax_bottom.set_ylabel("Mean accuracy", fontsize=lfs)

    # Every panel produces identical labels here, so pull handles from
    # whichever panel was drawn last and share one legend across all of them.
    handles, labels = ax_bottom.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    legend = fig.legend(by_label.values(), by_label.keys(), fontsize=lfs,
                         loc="upper left", bbox_to_anchor=(_SA_PP_LEGEND_ANCHOR_X, 0.95))
    _place_pp_gradient_bars(fig, legend, condition_colors, lfs,
                             panel_axes=(ax_top, ax_bottom))
    plt.savefig(filename, bbox_inches="tight")
    plt.close(fig)


def sa_pp_1v5(mean_stats, mean_thr_stats):
    """SA + PP curves, sub1 (left panel) vs sub5 (right panel). See
    _draw_sa_pp_figure for the full visual-encoding description.
    Output: sa_pp_1v5.pdf"""
    _draw_sa_pp_figure(mean_stats, mean_thr_stats, ["sub1", "sub5"], _out("sa_pp_1v5.pdf"))


# ---------------------------------------------------------------------------
# sample_metadata_1v5 (real data)
# ---------------------------------------------------------------------------

_SAMPLE_METADATA_ROW_LABELS = ["Center name", "Primer scheme"]

# Reserve less of the canvas at the bottom (no bottom-anchored legend here —
# see _draw_sample_metadata_figure) and more on the right (two per-row
# legends instead), unlike the other three _make_row_grid callers.
_SAMPLE_METADATA_GRIDSPEC_RIGHT = 0.75
_SAMPLE_METADATA_BOTTOM_MARGIN_IN = 0.5
_SAMPLE_METADATA_LEGEND_X = 0.78


def _draw_sample_metadata_figure(stats, threshold_stats, run_metadata,
                                  plot_subsamples, filename):
    """
    Two rows, each stratifying by one sample_metadata field independently
    (collapsing the other field entirely): top row color = CenterName,
    bottom row color = primer_scheme. Split into separate rows — rather than
    a single panel combining both via color=CenterName/linestyle=
    primer_scheme, which read poorly — since each row just aggregates over
    every run sharing that one field's value, ignoring the other field. Only
    pr_scaling=MLE. Columns = plot_subsamples. Sizing comes from
    _make_row_grid (shared with quality_flag_1v5/gene_region_1v5/
    replicate_acc_diff) — every row's panel renders at exactly
    (_ROW_PANEL_WIDTH_IN, _ROW_PANEL_HEIGHT_IN) regardless of row/column
    count — but with its own right/bottom margins (_SAMPLE_METADATA_*):
    since the two rows carry different color encodings, they get two
    independent legends (_place_row_legend, one per row, vertically centered
    on it and anchored just right of the panels) rather than one shared
    bottom-anchored legend, so there's no need to reserve bottom space for
    one and the freed-up canvas goes to the right margin instead. Each row
    is labeled above it ("Center name" / "Primer scheme", via
    _add_row_labels) in the same style as gene_region_1v5's structural/
    nonstructural headers — together with the split legends, the two rows
    read as two concatenated but independent plots rather than one merged
    figure. SA (circle) + PP curve on the same panel; every threshold is
    also plotted as its own diamond point along the curve, at alpha=0.5
    (color already carries the category here, so full-opacity diamonds
    crowded out the SA point — same convention as quality_flag/gene_region).
    Panels within a row share that row's y-axis (accuracy) range (5% pad);
    the two rows keep independent y-ranges since CenterName- and
    primer_scheme-grouped accuracy can sit at different levels. Fixed
    x-axis (call rate) range (0, 1.05).
    """
    lw  = 1.5
    sc  = 60
    dia = 40
    lfs = 8
    lms = 6
    condition = "pr_scaling=MLE"
    n_subs = len(plot_subsamples)

    center_names   = sorted({meta["CenterName"] for meta in run_metadata.values()})
    primer_schemes = sorted({meta["primer_scheme"] for meta in run_metadata.values()})

    center_colors = {
        c: cm.tab10.colors[i % len(cm.tab10.colors)] for i, c in enumerate(center_names)
    }
    scheme_colors = {
        s: cm.tab10.colors[i % len(cm.tab10.colors)] for i, s in enumerate(primer_schemes)
    }

    # Run sets per category value, one field at a time (the other field is
    # collapsed — unlike the combined (center, scheme) grouping this replaced).
    center_runs = {c: set() for c in center_names}
    scheme_runs = {s: set() for s in primer_schemes}
    for run, meta in run_metadata.items():
        center_runs[meta["CenterName"]].add(run)
        scheme_runs[meta["primer_scheme"]].add(run)

    rows = [
        ("CenterName", center_names, center_colors, center_runs),
        ("primer_scheme", primer_schemes, scheme_colors, scheme_runs),
    ]

    fig, axes, _ = _make_row_grid(len(rows), n_subs, right=_SAMPLE_METADATA_GRIDSPEC_RIGHT,
                                   bottom_margin_in=_SAMPLE_METADATA_BOTTOM_MARGIN_IN)

    for row, (field_name, values, colors, run_sets) in enumerate(rows):
        # Aggregate once per category value for this field.
        group_ms  = {}
        group_mts = {}
        for value in values:
            ms, mts = real._aggregate_for_run_subset(
                stats[condition], threshold_stats[condition], run_sets[value]
            )
            group_ms[value]  = ms
            group_mts[value] = mts

        row_acc_values = []
        row_axes = []
        for col, subsample in enumerate(plot_subsamples):
            ax = axes[row, col]
            row_axes.append(ax)

            for value in values:
                color = colors[value]

                sa_acc, sa_cr = group_ms[value][subsample].get(
                    "SA", (float("nan"), float("nan"))
                )
                ax.scatter(sa_cr, sa_acc, color=color, marker="o", s=sc, zorder=5)
                if sa_acc == sa_acc:
                    row_acc_values.append(sa_acc)

                curve_cr, curve_acc = [], []
                for t_idx in range(N_THRESHOLDS):
                    acc, cr = group_mts[value][subsample].get(
                        t_idx, (float("nan"), float("nan"))
                    )
                    curve_cr.append(cr)
                    curve_acc.append(acc)
                    if acc == acc:
                        row_acc_values.append(acc)
                ax.plot(curve_cr, curve_acc, color=color, linewidth=lw, zorder=4)
                for cr, acc in zip(curve_cr, curve_acc):
                    ax.scatter(cr, acc, color=color, marker=_PP_MARKER, s=dia, alpha=0.5, zorder=6)

            if row == 0:
                ax.set_title(real._fmt_subsample(subsample), fontsize=lfs + 1)
            ax.set_xlabel("Mean call rate", fontsize=lfs)
            ax.tick_params(labelsize=lfs - 1)
            ax.ticklabel_format(useOffset=False)

        if row_acc_values:
            lo, hi = min(row_acc_values), max(row_acc_values)
            pad = (hi - lo) * 0.05 or 0.01
            ylim = (lo - pad, hi + pad)
            for ax in row_axes:
                ax.set_ylim(*ylim)
                ax.set_xlim(0, 1.05)

        axes[row, 0].set_ylabel("Mean accuracy", fontsize=lfs)

        row_legend_handles = [
            Line2D([0], [0], color=colors[value], linewidth=lw, label=value)
            for value in values
        ]
        row_legend_handles.append(
            Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
                   markersize=lms, label="SA")
        )
        row_legend_handles.append(
            Line2D([0], [0], color="gray", linestyle="-", marker=_PP_MARKER,
                   markersize=lms, label="PP (all thresholds)")
        )
        _place_row_legend(fig, row_axes, row_legend_handles,
                           [h.get_label() for h in row_legend_handles],
                           lfs, _SAMPLE_METADATA_LEGEND_X)

    _add_row_labels(fig, axes, lfs, _SAMPLE_METADATA_ROW_LABELS)
    plt.savefig(filename, bbox_inches="tight")
    plt.close(fig)


def sample_metadata_1v5(stats, threshold_stats, run_metadata):
    """
    sample_metadata narrowed to sub1 (left column) and sub5 (right column) only.
    Output: sample_metadata_1v5.pdf
    """
    _draw_sample_metadata_figure(stats, threshold_stats, run_metadata,
                                  ["sub1", "sub5"], _out("sample_metadata_1v5.pdf"))


# ---------------------------------------------------------------------------
# quality_flag_1v5 (real data)
# ---------------------------------------------------------------------------

def _draw_quality_flag_figure(mean_filter_stats, mean_filter_thr_stats,
                               plot_subsamples, filename):
    """
    All quality_flag (FILTER) values on one compact figure; color =
    quality_flag value. Only pr_scaling=MLE. Every threshold is plotted as
    its own diamond point in that quality_flag value's color at alpha=0.5
    (unlike the sa_pp family, where the gradient gives each diamond its own
    meaning; here color is already carrying quality_flag, so the diamonds
    are just there to mark the curve's thresholds, and full opacity made the
    SA points hard to pick out among them), threaded together by a
    connecting line in the same color. One panel per subsample in
    plot_subsamples, arranged in a single row. Sizing comes from
    _make_row_grid (shared with sample_metadata_1v5/gene_region_1v5/
    replicate_acc_diff) — every panel renders at exactly
    (_ROW_PANEL_WIDTH_IN, _ROW_PANEL_HEIGHT_IN) regardless of row/column
    count, so this figure's single row matches a row in any of its
    2-row siblings exactly, not just approximately. Every panel shares the
    same y-axis (accuracy) range, computed from every point plotted across
    all panels with a 5% pad, and the same x-axis (call rate) range fixed
    at (0, 1.05).
    """
    lw  = 1.5
    sc  = 60
    dia = 40
    lfs = 8
    lms = 6
    condition   = "pr_scaling=MLE"
    n_subs      = len(plot_subsamples)
    flag_colors = {v: cm.tab10.colors[i] for i, v in enumerate(real.FILTER_VALUES)}
    ms_mle  = mean_filter_stats[condition]
    mts_mle = mean_filter_thr_stats[condition]

    fig, axes, legend_y = _make_row_grid(1, n_subs)
    axes_flat = axes[0]

    all_acc_values = []
    for ax, subsample in zip(axes_flat, plot_subsamples):
        for filt_val in real.FILTER_VALUES:
            color = flag_colors[filt_val]

            sa_acc, sa_cr = ms_mle[subsample][filt_val].get(
                "SA", (float("nan"), float("nan"))
            )
            ax.scatter(sa_cr, sa_acc, color=color, marker="o", s=sc, zorder=5)
            if sa_acc == sa_acc:
                all_acc_values.append(sa_acc)

            curve_cr, curve_acc = [], []
            for t_idx in range(N_THRESHOLDS):
                acc, cr = mts_mle[subsample][filt_val].get(
                    t_idx, (float("nan"), float("nan"))
                )
                curve_cr.append(cr)
                curve_acc.append(acc)
                if acc == acc:
                    all_acc_values.append(acc)
            ax.plot(curve_cr, curve_acc, color=color, linewidth=lw, zorder=4)
            for cr, acc in zip(curve_cr, curve_acc):
                ax.scatter(cr, acc, color=color, marker=_PP_MARKER, s=dia, alpha=0.5, zorder=6)

        ax.set_title(real._fmt_subsample(subsample), fontsize=lfs + 1)
        ax.set_xlabel("Mean call rate", fontsize=lfs)
        ax.tick_params(labelsize=lfs - 1)
        ax.ticklabel_format(useOffset=False)

    if all_acc_values:
        lo, hi = min(all_acc_values), max(all_acc_values)
        pad = (hi - lo) * 0.05 or 0.01
        ylim = (lo - pad, hi + pad)
        for ax in axes_flat:
            ax.set_ylim(*ylim)
            ax.set_xlim(0, 1.05)

    axes_flat[0].set_ylabel("Mean accuracy", fontsize=lfs)

    legend_handles = [
        Line2D([0], [0], color=flag_colors[v], linewidth=lw, label=v)
        for v in real.FILTER_VALUES
    ]
    legend_handles.append(
        Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
               markersize=lms, label="SA")
    )
    legend_handles.append(
        Line2D([0], [0], color="gray", linestyle="-", marker=_PP_MARKER,
               markersize=lms, label="PP (all thresholds)")
    )

    fig.legend(handles=legend_handles, loc="upper center",
               bbox_to_anchor=(0.5, legend_y), fontsize=lfs, ncol=4)
    plt.savefig(filename, bbox_inches="tight")
    plt.close(fig)


def quality_flag_1v5(mean_filter_stats, mean_filter_thr_stats):
    """
    quality_flag narrowed to sub1 (left panel) and sub5 (right panel) only.
    Output: quality_flag_1v5.pdf
    """
    _draw_quality_flag_figure(mean_filter_stats, mean_filter_thr_stats,
                               ["sub1", "sub5"], _out("quality_flag_1v5.pdf"))


# ---------------------------------------------------------------------------
# Shared structural/nonstructural row-grid label helper
# (gene_region_1v5, replicate_acc_diff — both 2-row grids: row 0 = structural
# genes, row 1 = nonstructural, matching real.GENE_REGION_TYPES order)
# ---------------------------------------------------------------------------

_GENE_TYPE_ROW_LABELS = [
    "Structural protein-coding genes",
    "Non-structural protein-coding genes",
]

# Explicit inch-based panel sizing shared by every "categorical-comparison"
# figure (sample_metadata_1v5, quality_flag_1v5, gene_region_1v5,
# replicate_acc_diff), so a row's rendered panel comes out the same absolute
# size in every one of them regardless of how many rows or columns that
# particular figure has. Plain tight_layout() cannot guarantee this: its
# margins are computed per-figure from that figure's own content, and — less
# obviously — a 1-row figure (quality_flag_1v5) ends up with a *taller* row
# than a 2-row figure under the same nominal top/bottom margins, since there
# is no second row competing for the same vertical budget. Explicit margins
# sidestep both problems by fixing each row's height directly, independent
# of row count. Square (width == height).
_ROW_PANEL_WIDTH_IN   = 3.0
_ROW_PANEL_HEIGHT_IN  = 3.0
_ROW_GAP_IN           = 0.75  # between stacked rows (room for a row label, when present)
_ROW_TOP_MARGIN_IN    = 0.35  # above the top row (room for its own column titles)
_ROW_BOTTOM_MARGIN_IN = 1.7   # below the bottom row (room for the bottom-anchored legend)
_ROW_WSPACE           = 0.24  # fraction of panel width, gap between columns
_ROW_GRIDSPEC_LEFT    = 0.06
_ROW_GRIDSPEC_RIGHT   = 0.99


def _make_row_grid(n_rows, n_cols, right=_ROW_GRIDSPEC_RIGHT,
                    bottom_margin_in=_ROW_BOTTOM_MARGIN_IN):
    """Build a figure + (n_rows, n_cols) array of Axes for the categorical-
    comparison family, using the explicit _ROW_* margins above instead of
    tight_layout(), so every row's panel renders at exactly
    (_ROW_PANEL_WIDTH_IN, _ROW_PANEL_HEIGHT_IN) regardless of n_rows/n_cols
    (regardless of `right`/`bottom_margin_in` too — those only change how
    much canvas is reserved outside the panels, e.g. for a right-side legend
    instead of the usual bottom-anchored one, as sample_metadata_1v5 does).
    Returns (fig, axes, legend_y) — legend_y is a figure-fraction y position
    inside the reserved bottom margin, for `fig.legend(..., bbox_to_anchor=
    (0.5, legend_y))`; unused by callers that place their legend(s)
    elsewhere."""
    usable_w_frac = right - _ROW_GRIDSPEC_LEFT
    fig_width = _ROW_PANEL_WIDTH_IN * (n_cols + (n_cols - 1) * _ROW_WSPACE) / usable_w_frac
    fig_height = (_ROW_TOP_MARGIN_IN + n_rows * _ROW_PANEL_HEIGHT_IN
                  + (n_rows - 1) * _ROW_GAP_IN + bottom_margin_in)
    top_frac = (fig_height - _ROW_TOP_MARGIN_IN) / fig_height
    bottom_frac = bottom_margin_in / fig_height
    hspace_frac = _ROW_GAP_IN / _ROW_PANEL_HEIGHT_IN

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(n_rows, n_cols, left=_ROW_GRIDSPEC_LEFT, right=right,
                           top=top_frac, bottom=bottom_frac,
                           wspace=_ROW_WSPACE, hspace=hspace_frac)
    axes = np.empty((n_rows, n_cols), dtype=object)
    for r in range(n_rows):
        for c in range(n_cols):
            axes[r, c] = fig.add_subplot(gs[r, c])

    legend_y = 0.6 * bottom_margin_in / fig_height
    return fig, axes, legend_y


def _place_row_legend(fig, row_axes, handles, labels, lfs, x_frac):
    """Places a legend vertically centered on one row of Axes (the union of
    row_axes' rendered bounding boxes), horizontally anchored at x_frac
    (figure fraction) — same bbox-based vertical-centering technique
    _place_pp_gradient_bars uses for the sa_pp family's gradient bars, here
    for a plain per-row legend instead (sample_metadata_1v5's two rows carry
    different color encodings — CenterName, primer_scheme — so need two
    independent legends rather than one shared across both)."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    y0 = min(ax.get_window_extent(renderer).y0 for ax in row_axes)
    y1 = max(ax.get_window_extent(renderer).y1 for ax in row_axes)
    y_center_frac = fig.transFigure.inverted().transform((0, (y0 + y1) / 2))[1]
    return fig.legend(handles, labels, loc="center left",
                       bbox_to_anchor=(x_frac, y_center_frac), fontsize=lfs)


def _add_row_labels(fig, axes, lfs, labels):
    """Places a section header above each row of a 2-row grid, in the same
    fig.text()-above-the-block style used for sim_sa_pp_1v5's error_rate
    block labels. Positioned from the row's actual rendered top edge — the
    *max* over every column's axes box and (if present) its ax.set_title()
    text, via get_window_extent() after an explicit canvas draw (same
    pixel-based technique _place_pp_gradient_bars uses). Clearing the title
    text specifically matters here: a fixed offset above just the axes box
    (no title check) looks fine with an even column count, but with an odd
    count the row label's centered x position lands exactly on the middle
    column's own title, overlapping it.

    For every row after the first, the label must *also* clear the row
    above's xlabel text (not just its own row's top edge) — the same
    odd-column-count issue, but against the previous row's xlabel instead of
    this row's own title, so a fixed pad above just this row's top edge
    isn't enough on its own: normally the label sits pad_px above this row's
    top (anchored va="bottom" so that's its own bottom edge), but is capped
    from rising any higher than pad_px below the previous row's xlabel
    bottom, whenever that would otherwise be exceeded.

    Horizontally centered on the panels' own span (leftmost to rightmost
    axes box in the top row), not a fixed figure-fraction 0.5 — the two
    coincide only when the figure has no asymmetric off-panel content; a
    figure with e.g. a right-side legend (sample_metadata_1v5) has its
    panels occupying less than the full canvas width, so 0.5 would land the
    label noticeably right of the panels' true center.

    Must be called after the figure's layout is finalized (after
    tight_layout()/subplots_adjust() or an equivalent explicit-margin
    layout, before savefig()) — extra top margin and inter-row spacing must
    already be reserved by the caller so these labels have room."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    pad_px = 6
    text_h_px = (lfs + 2) * fig.dpi / 72.0  # approx rendered text height
    n_rows = axes.shape[0]
    x0_px = min(ax.get_window_extent(renderer).x0 for ax in axes[0])
    x1_px = max(ax.get_window_extent(renderer).x1 for ax in axes[0])
    x_center = fig.transFigure.inverted().transform(((x0_px + x1_px) / 2, 0))[0]
    for row, label in enumerate(labels[:n_rows]):
        row_top_px = max(
            max(ax.get_window_extent(renderer).y1 for ax in axes[row]),
            max((ax.title.get_window_extent(renderer).y1 for ax in axes[row] if ax.get_title()),
                default=0),
        )
        y_bottom_px = row_top_px + pad_px

        if row > 0:
            prev_xlabel_bottoms = [
                ax.xaxis.label.get_window_extent(renderer).y0
                for ax in axes[row - 1] if ax.xaxis.label.get_text()
            ]
            if prev_xlabel_bottoms:
                max_y_bottom = min(prev_xlabel_bottoms) - pad_px - text_h_px
                y_bottom_px = min(y_bottom_px, max_y_bottom)

        y_frac = fig.transFigure.inverted().transform((0, y_bottom_px))[1]
        fig.text(x_center, y_frac, label, ha="center", va="bottom", fontsize=lfs + 2)


# ---------------------------------------------------------------------------
# gene_region_1v5 (real data)
# ---------------------------------------------------------------------------

def _draw_gene_region_figure(mean_gene_stats, mean_gene_thr_stats, all_genes,
                              plot_subsamples, filename):
    """
    Rows are grouped by gene_region_type (structural / nonstructural), each
    labeled with a section header above it via _add_gene_type_row_labels.
    Columns = plot_subsamples, one panel per subsample. Sizing comes from
    _make_row_grid (shared with sample_metadata_1v5/quality_flag_1v5/
    replicate_acc_diff) — every row's panel renders at exactly
    (_ROW_PANEL_WIDTH_IN, _ROW_PANEL_HEIGHT_IN) regardless of row/column
    count. SA + PP
    on same panel. Color = gene region. Only pr_scaling=MLE. "other"
    (positions outside any named gene) is excluded from all_genes here — it's
    uninformative for a gene-by-gene comparison. Every threshold is plotted
    as its own diamond point in that gene's color at alpha=0.5 (unlike the
    sa_pp family, where the gradient gives each diamond its own meaning; here
    color is already carrying gene region, so the diamonds are just there to
    mark the curve's thresholds, and full opacity made the SA points hard to
    pick out among them), threaded together by a connecting line in the same
    color. Panels within a row share that row's y-axis (accuracy) range,
    computed from every point plotted across that row's columns with a 5%
    pad, and a fixed x-axis (call rate) range of (0, 1.05). The two rows
    (structural / nonstructural genes) keep independent y-ranges rather than
    sharing one across the whole figure, since those two gene groups can sit
    at quite different accuracy levels and forcing them onto one scale would
    squash whichever row has the smaller range.
    """
    lw  = 1.5
    sc  = 60
    dia = 40
    lfs = 8
    lms = 6
    condition = "pr_scaling=MLE"
    n_subs    = len(plot_subsamples)

    # cat_colors stays keyed off the full all_genes (below) so a gene's color
    # never shifts regardless of whether "other" is excluded from plotting —
    # only the iteration/legend lists are filtered.
    plot_genes = [g for g in all_genes if g != "other"]

    gene_to_type = {
        gene: ("structural protein genes" if gene in real.STRUCTURAL_GENES else "nonstructural genes")
        for gene in plot_genes
    }
    genes_by_type = {
        t: [g for g in plot_genes if gene_to_type[g] == t] for t in real.GENE_REGION_TYPES
    }

    cmap = cm.tab20 if len(all_genes) > 10 else cm.tab10
    cat_colors = {g: cmap.colors[i % len(cmap.colors)] for i, g in enumerate(all_genes)}

    n_types = len(real.GENE_REGION_TYPES)
    fig, axes, legend_y = _make_row_grid(n_types, n_subs)

    for row, gene_type in enumerate(real.GENE_REGION_TYPES):
        row_acc_values = []
        row_axes = []
        for col, subsample in enumerate(plot_subsamples):
            ax = axes[row, col]
            row_axes.append(ax)
            ms  = mean_gene_stats[condition][subsample]
            mts = mean_gene_thr_stats[condition][subsample]

            for gene in genes_by_type[gene_type]:
                if gene not in ms:
                    continue
                color = cat_colors[gene]

                sa_acc, sa_cr = ms[gene].get("SA", (float("nan"), float("nan")))
                ax.scatter(sa_cr, sa_acc, color=color, marker="o", s=sc, zorder=5)
                if sa_acc == sa_acc:
                    row_acc_values.append(sa_acc)

                curve_cr, curve_acc = [], []
                for t_idx in range(N_THRESHOLDS):
                    acc, cr = mts[gene].get(t_idx, (float("nan"), float("nan")))
                    curve_cr.append(cr)
                    curve_acc.append(acc)
                    if acc == acc:
                        row_acc_values.append(acc)
                ax.plot(curve_cr, curve_acc, color=color, linewidth=lw, zorder=4)
                for cr, acc in zip(curve_cr, curve_acc):
                    ax.scatter(cr, acc, color=color, marker=_PP_MARKER, s=dia, alpha=0.5, zorder=6)

            if row == 0:
                ax.set_title(real._fmt_subsample(subsample), fontsize=lfs + 1)
            ax.set_xlabel("Mean call rate", fontsize=lfs)
            ax.tick_params(labelsize=lfs - 1)
            ax.ticklabel_format(useOffset=False)

        if row_acc_values:
            lo, hi = min(row_acc_values), max(row_acc_values)
            pad = (hi - lo) * 0.05 or 0.01
            ylim = (lo - pad, hi + pad)
            for ax in row_axes:
                ax.set_ylim(*ylim)
                ax.set_xlim(0, 1.05)

        axes[row, 0].set_ylabel("Mean accuracy", fontsize=lfs)

    legend_handles = [
        Line2D([0], [0], color=cat_colors[g], linewidth=lw, marker="o",
               markersize=lms, linestyle="-", label=g)
        for g in plot_genes
    ]
    legend_handles.append(
        Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
               markersize=lms, label="SA")
    )
    legend_handles.append(
        Line2D([0], [0], color="gray", linestyle="-", marker=_PP_MARKER,
               markersize=lms, label="PP (all thresholds)")
    )
    _add_row_labels(fig, axes, lfs, _GENE_TYPE_ROW_LABELS)
    fig.legend(handles=legend_handles, loc="upper center",
               bbox_to_anchor=(0.5, legend_y), fontsize=lfs, ncol=4)
    plt.savefig(filename, bbox_inches="tight")
    plt.close(fig)


def gene_region_1v5(mean_gene_stats, mean_gene_thr_stats, all_genes):
    """
    gene_region narrowed to sub1 (left column) and sub5 (right column) only.
    Output: gene_region_1v5.pdf
    """
    _draw_gene_region_figure(mean_gene_stats, mean_gene_thr_stats, all_genes,
                              ["sub1", "sub5"], _out("gene_region_1v5.pdf"))


# ---------------------------------------------------------------------------
# replicate_acc_diff (real data, sub20_only replicates)
# ---------------------------------------------------------------------------

_GENE_ROW_ORDER = ["Spike", "Nucleocapsid", "Membrane", "Envelope"]


def replicate_acc_diff(gene_rep_stats, gene_rep_thr_stats, all_genes,
                        mle_mode_t_idx, run_metadata):
    """
    4×3 grid: rows = structural gene, in the fixed order _GENE_ROW_ORDER
    (Spike, Nucleocapsid, Membrane, Envelope, top to bottom); columns =
    CenterName (there happen to be 3) — one panel per (gene, CenterName)
    combination. Nonstructural genes are excluded entirely (not just
    "other") — replicate variability is concentrated in the structural
    genes. This used to overlay all 3 centers by color within one panel per
    gene (2×2); splitting centers out into their own columns instead (i)
    cuts per-panel crowding roughly 3x and (ii) lets each (gene, center)
    combination scale independently — see the "independent axes" note below
    — since the point of this figure is to expose exactly how much gene- and
    center-dependent replicate variability there is. pr_scaling=MLE only —
    pr_scaling=0.0 replicate data still exists in gene_rep_stats/
    gene_rep_thr_stats but is not plotted here. Sizing comes from
    _make_row_grid (shared with sample_metadata_1v5/quality_flag_1v5/
    gene_region_1v5) — every panel renders at exactly (_ROW_PANEL_WIDTH_IN,
    _ROW_PANEL_HEIGHT_IN) regardless of row/column count.
    Data comes from the sub20-only replicate call files (REPLICATE_CONDITIONS,
    REPLICATE_DIR in roc_real_data.py), not the main condition files — see
    read_replicate_data().

    Output: replicate_acc_diff.pdf — y-axis: PP_thr minus SA accuracy.
    (A call-rate analogue, replicate_cr_diff.pdf, used to exist here too, but
    was dropped: since find_pp_sa_crossover now picks mle_mode_t_idx as the
    most stringent threshold whose call rate still meets/exceeds SA's, PP's
    call rate at that threshold is by construction always >= SA's, making a
    call-rate diff plot uninformative.)

    Each point is one (Run, replicate, n_reads_cat) observation within that
    (gene, center) panel — replicates are not averaged together, each gets
    its own point, distinguished by marker shape via _REPLICATE_MARKERS,
    since replicate identity turned out to matter: see the sub20_only
    replicate-duplication investigation. Color still encodes CenterName
    (kept for visual continuity/scanning even though CenterName is also now
    a column, at the user's request) via the same tab10 assignment used
    elsewhere; every point within a column is therefore one color, with
    marker shape distinguishing replicates within it. Plotted with
    alpha=0.55 and a small random x jitter (n_reads is bucketed to a handful
    of integer values, so multiple runs/replicates can still share the same
    x — jitter and transparency together make overlapping points visible;
    random.seed(42)).
    x-axis : n_reads value (integer), jittered.
    Color  : CenterName.
    Marker : replicate ID (rep1, rep2, ...).
    Row labels (gene names) come from _add_row_labels (shared with
    gene_region_1v5/sample_metadata_1v5); column titles (CenterName) are set
    on the top row only, matching every other _make_row_grid-based figure's
    convention (e.g. gene_region_1v5's subsample column titles).
    Each panel has its own independent x- and y-axis range — unlike the
    single shared range this replaced, gene- and center-specific variability
    is exactly what this figure is meant to expose, so forcing every panel
    onto a common scale would hide it.
    """
    condition = "pr_scaling=MLE"
    if condition not in gene_rep_stats or condition not in gene_rep_thr_stats:
        print("replicate_acc_diff: no MLE replicate data available, skipping.")
        return

    plot_genes = [
        g for g in _GENE_ROW_ORDER
        if g in real.STRUCTURAL_GENES and g in all_genes and g in gene_rep_stats[condition]
    ]
    center_names = sorted({meta["CenterName"] for meta in run_metadata.values()})
    if not plot_genes or not center_names or mle_mode_t_idx is None:
        print("replicate_acc_diff: no replicate data available, skipping.")
        return

    random.seed(42)
    sc  = 40
    lfs = 8
    x_jitter = 0.15

    center_colors = {
        c: cm.tab10.colors[i % len(cm.tab10.colors)] for i, c in enumerate(center_names)
    }
    color_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=center_colors[c],
               markersize=lfs, label=c)
        for c in center_names
    ]

    def _rep_sort_key(rep):
        m = re.match(r"rep(\d+)$", rep)
        return (0, int(m.group(1))) if m else (1, rep)

    all_replicates = set()
    gs_all = gene_rep_stats[condition]
    for gene in plot_genes:
        for method_dict in gs_all.get(gene, {}).values():
            all_replicates.update(rep for _, rep in method_dict["SA"].keys())
    replicate_labels = sorted(all_replicates, key=_rep_sort_key)
    replicate_markers = {
        rep: _REPLICATE_MARKERS[i % len(_REPLICATE_MARKERS)]
        for i, rep in enumerate(replicate_labels)
    }
    marker_handles = [
        Line2D([0], [0], marker=replicate_markers[rep], color="w",
               markerfacecolor="gray", markeredgecolor="gray",
               markersize=lfs, label=rep)
        for rep in replicate_labels
    ]

    def _collect(gene, center):
        """Return a list of (n_reads_val, diff, replicate) tuples: one point
        per (run, replicate, n_reads_cat) belonging to a run at this
        CenterName, accuracy diff (pr_scaling=MLE). Not averaged across
        replicates -- each replicate is its own point (see marker shape)."""
        points = []
        gs  = gene_rep_stats[condition]
        gts = gene_rep_thr_stats[condition]
        if gene not in gs or gene not in gts:
            return points
        for nr_cat, method_dict in gs[gene].items():
            n_reads_val = int(nr_cat[2:])
            if nr_cat not in gts[gene]:
                continue
            thr_by_rep = gts[gene][nr_cat][mle_mode_t_idx]
            for (run, rep), (acc_sum_sa, cr_sum_sa, n_total) in method_dict["SA"].items():
                if run not in run_metadata or run_metadata[run]["CenterName"] != center:
                    continue
                if (run, rep) not in thr_by_rep:
                    continue
                acc_sum_pp, cr_sum_pp, n_passing = thr_by_rep[(run, rep)]
                sa_val = acc_sum_sa / cr_sum_sa if cr_sum_sa > 0 else float("nan")
                pp_val = acc_sum_pp / cr_sum_pp if cr_sum_pp > 0 else float("nan")
                diff = pp_val - sa_val
                if diff == diff:  # not NaN
                    points.append((n_reads_val, diff, rep))
        return points

    n_rows, n_cols = len(plot_genes), len(center_names)
    fig, axes, legend_y = _make_row_grid(n_rows, n_cols)

    for row, gene in enumerate(plot_genes):
        for col, center in enumerate(center_names):
            ax = axes[row, col]
            color = center_colors[center]
            points = _collect(gene, center)
            by_rep = defaultdict(lambda: ([], []))
            for x, y, rep in points:
                xs_j, ys = by_rep[rep]
                xs_j.append(x + random.uniform(-x_jitter, x_jitter))
                ys.append(y)
            for rep, (xs_j, ys) in by_rep.items():
                ax.scatter(xs_j, ys, color=color, marker=replicate_markers[rep],
                           s=sc, alpha=0.55, linewidths=0)
            ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
            if row == 0:
                ax.set_title(center, fontsize=lfs + 1)
            ax.set_xlabel("n_reads", fontsize=lfs)
            ax.tick_params(labelsize=lfs - 1)
            ax.ticklabel_format(useOffset=False)
            if col == 0:
                ax.set_ylabel("PP_thr − SA accuracy", fontsize=lfs)

    _add_row_labels(fig, axes, lfs, plot_genes)
    # Two independent legends (color=CenterName, marker=replicate) side by
    # side under the panels, rather than one combined legend -- combining
    # both encodings into a single handle set would need one entry per
    # (CenterName, replicate) pair.
    fig.legend(handles=color_handles, loc="upper center",
               bbox_to_anchor=(0.28, legend_y), fontsize=lfs, ncol=1,
               title="CenterName", title_fontsize=lfs)
    fig.legend(handles=marker_handles, loc="upper center",
               bbox_to_anchor=(0.72, legend_y), fontsize=lfs,
               ncol=min(5, len(replicate_labels)),
               title="Replicate", title_fontsize=lfs)
    plt.savefig(_out("replicate_acc_diff.pdf"), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# sim_sa_pp_1v5 (simulated data)
# ---------------------------------------------------------------------------

def sim_sa_pp_1v5(mean_stats, mean_thr_stats):
    """
    SA + PP curves excluding pr_scaling=0.001, no PR/LL — simulated-data
    analogue of sa_pp_1v5.

    2x2 panel grid: rows = error_rate (beta0.001 top block, beta0.007 bottom
    block), columns = mean_depth (sub1 left, sub5 right). Unlike sa_pp_1v5,
    each panel here is a single (unbroken) axis — no y-axis break. Panel
    width and height match sa_pp_1v5's panels exactly (_SA_PP_PANEL_WIDTH /
    _SA_PP_PANEL_HEIGHT, shared with the real-data pipeline — see module
    docstring), and the same wspace (0.3) is used between columns, so the
    overall figure width equals sa_pp_1v5's; the figure is taller only
    because of the second error_rate row (plus the gap between blocks for
    the block label). Each block's y-range is computed independently (not
    shared across both blocks), since the two error rates can sit at very
    different overall accuracy levels — two separate GridSpecs, one per
    block, positioned at different top/bottom figure-fraction ranges
    (converted from inch-based block heights/margins so panel size stays
    fixed in inches regardless of block count).

    Every threshold plotted as a point, gradient condition color (least
    stringent) → black (most stringent), threaded together by a connecting
    line in the condition's base color. SA gets its own color (tab10 index
    1, otherwise unused since pr_scaling=0.001 is excluded). The legend's PP
    entries are a zero-length line+diamond proxy in the condition's color,
    and SA's legend entry is likewise just its circle. Every panel has
    exactly one depth, so no marker/linestyle is needed to disambiguate —
    one shared legend for the whole figure, placed fully outside the plotted
    area. Each condition gets its own vertical gradient colorbar-style strip
    (see plot_common), left-aligned under that SA/PP legend but vertically
    centered on the top block's panel midpoint, labeled with the actual
    threshold values (via _fmt_threshold) — same convention as sa_pp_1v5,
    including matching the bars' absolute (inch) size despite this figure
    being much taller overall.

    Output: sim_sa_pp_1v5.pdf
    """
    error_rate_blocks = ("beta0.001", "beta0.007")
    _ERROR_RATE_LABELS = {"beta0.001": "Low error rate", "beta0.007": "High error rate"}
    depths = ("sub1", "sub5")

    conditions = {k: v for k, v in sim.CONDITIONS.items() if k != "pr_scaling=0.001"}
    all_condition_colors = {
        cond: col for cond, col in zip(sim.CONDITIONS.keys(), cm.tab10.colors[:len(sim.CONDITIONS)])
    }
    condition_colors = {c: all_condition_colors[c] for c in conditions}
    first_cond = next(iter(conditions))
    sa_color = cm.tab10.colors[1]  # unused by any plotted condition here

    sc, dia, lfs = 60, 40, 8
    _SA_MARKER = "o"

    n_cols = len(depths)
    n_blocks = len(error_rate_blocks)
    wspace = _SA_PP_WSPACE

    # Block layout computed in inches (panel width/height/margins shared
    # with sa_pp_1v5 via the _SA_PP_* constants — see that function's
    # docstring), then converted to figure fractions for add_gridspec.
    fig_width = _SA_PP_PANEL_WIDTH * (n_cols + (n_cols - 1) * wspace) / (_SA_PP_GRIDSPEC_RIGHT - _SA_PP_GRIDSPEC_LEFT)
    fig_height = (_SA_PP_TOP_MARGIN_IN + n_blocks * _SA_PP_PANEL_HEIGHT
                  + (n_blocks - 1) * _SA_PP_BLOCK_GAP_IN + _SA_PP_BOTTOM_MARGIN_IN)
    fig = plt.figure(figsize=(fig_width, fig_height))

    last_ax = None
    top_block_ax = None  # any panel from the first (topmost) error_rate block
    for bi, error_rate in enumerate(error_rate_blocks):
        block_top_in = fig_height - _SA_PP_TOP_MARGIN_IN - bi * (_SA_PP_PANEL_HEIGHT + _SA_PP_BLOCK_GAP_IN)
        block_bottom_in = block_top_in - _SA_PP_PANEL_HEIGHT
        block_top = block_top_in / fig_height
        block_bottom = block_bottom_in / fig_height
        gs = fig.add_gridspec(1, n_cols, wspace=wspace,
                               left=_SA_PP_GRIDSPEC_LEFT, right=_SA_PP_GRIDSPEC_RIGHT,
                               top=block_top, bottom=block_bottom)

        # Precompute this block's data and its shared y-range (no break, so
        # SA + all PP conditions share one range per block).
        sa_by_depth = {}
        curves_by_depth = {}
        all_values = []
        for depth in depths:
            sa_acc, sa_cr = mean_stats.get(first_cond, {}).get(error_rate, {}).get(
                depth, {}
            ).get("SA", (float("nan"), float("nan")))
            sa_by_depth[depth] = (sa_cr, sa_acc)
            if sa_acc == sa_acc:
                all_values.append(sa_acc)

            depth_curves = []
            for condition, color in condition_colors.items():
                curve_cr, curve_acc = [], []
                cond_ts = mean_thr_stats.get(condition, {}).get(error_rate, {}).get(depth, {})
                for t_idx in range(N_THRESHOLDS):
                    acc, cr = cond_ts.get(t_idx, (float("nan"), float("nan")))
                    curve_cr.append(cr)
                    curve_acc.append(acc)
                    if acc == acc:
                        all_values.append(acc)
                depth_curves.append((condition, color, curve_cr, curve_acc))
            curves_by_depth[depth] = depth_curves

        if all_values:
            lo, hi = min(all_values), max(all_values)
            pad = (hi - lo) * 0.05 or 0.01
            ylim = (lo - pad, hi + pad)
        else:
            ylim = (0, 1)

        block_axes = []
        for col, depth in enumerate(depths):
            ax = fig.add_subplot(gs[0, col])
            block_axes.append(ax)

            sa_cr, sa_acc = sa_by_depth[depth]
            ax.scatter(sa_cr, sa_acc, marker=_SA_MARKER, color=sa_color, s=sc, zorder=5,
                       label="SA")

            for condition, color, curve_cr, curve_acc in curves_by_depth[depth]:
                pp_label = f"PP ({condition})"
                ax.plot(curve_cr, curve_acc, color=color, linewidth=1.5, zorder=4)
                ax.plot([], [], color=color, linewidth=1.5, marker=_PP_MARKER,
                        markersize=6, label=pp_label)
                for t_idx, (cr, acc) in enumerate(zip(curve_cr, curve_acc)):
                    fraction = 1 - t_idx / (N_THRESHOLDS - 1)
                    pt_color = _lerp_to_black(color, fraction)
                    ax.scatter(cr, acc, color=pt_color, marker=_PP_MARKER, s=dia, zorder=6)

            ax.set_ylim(*ylim)
            ax.set_xlim(0, 1.05)
            ax.ticklabel_format(useOffset=False)
            ax.tick_params(labelsize=lfs - 1)

            ax.set_title(sim._fmt_depth(depth), fontsize=lfs + 1)
            ax.set_xlabel("Mean call rate", fontsize=lfs)
            if col == 0:
                ax.set_ylabel("Mean accuracy", fontsize=lfs)

            last_ax = ax
            if bi == 0:
                top_block_ax = ax

        # Block-level label identifying the error_rate, positioned above this
        # block's own panel titles (not just a fixed offset above the
        # block's nominal top edge, which put it at nearly the same height
        # as "mean depth = Nx" instead of clearly higher) — computed the
        # same pixel-based way _add_row_labels clears a row's own title
        # text, via get_window_extent() after a canvas draw.
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        pad_px = 6
        top_px = max(
            max(ax.get_window_extent(renderer).y1 for ax in block_axes),
            max((ax.title.get_window_extent(renderer).y1 for ax in block_axes if ax.get_title()),
                default=0),
        )
        y_frac = fig.transFigure.inverted().transform((0, top_px + pad_px))[1]
        fig.text(0.44, y_frac, _ERROR_RATE_LABELS[error_rate],
                  ha="center", va="bottom", fontsize=lfs + 2)

    # One shared legend for the whole figure — every panel produced identical
    # labeled handles (SA, PP per condition), so pull from whichever axis was
    # drawn last. Anchored past the panels' right edge (gridspec right=0.80
    # above) so it sits fully outside the plotted area rather than drifting
    # over the rightmost panel. Gradient bars are left-aligned under this
    # legend but vertically centered on the top block's panel midpoint.
    handles, labels = last_ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    legend = fig.legend(by_label.values(), by_label.keys(), fontsize=lfs,
                         loc="upper left", bbox_to_anchor=(_SA_PP_LEGEND_ANCHOR_X, 0.95))
    _place_pp_gradient_bars(fig, legend, condition_colors, lfs, panel_axes=(top_block_ax,))

    filename = _out("sim_sa_pp_1v5.pdf")
    plt.savefig(filename, bbox_inches="tight")
    plt.close(fig)
    print(f"Written {filename}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    rebuild_all = "--rebuild" in sys.argv
    rebuild_real = rebuild_all or "--rebuild-real" in sys.argv
    rebuild_reps = rebuild_all or "--rebuild-reps" in sys.argv
    rebuild_sims = rebuild_all or "--rebuild-sims" in sys.argv

    # 1. Load metadata (small CSVs — always reloaded)
    run_metadata            = real.load_sample_metadata(real.SAMPLE_METADATA_PATH)
    pos_filter              = real.load_filter_metadata(real.FILTER_METADATA_PATH)
    pos_gene, gene_lengths  = real.load_gene_metadata(real.GENE_METADATA_PATH)
    filter_lengths          = real.compute_filter_lengths(pos_filter)

    # 2. Real (subsampled) data — stats, threshold_stats, aggregates,
    #    crossover, and mle_mode_t_idx all live in one cache, since they're
    #    all produced together in one --rebuild-real pass.
    if not rebuild_real and os.path.exists(REAL_CACHE_PATH):
        print(f"Loading real-data cache from {REAL_CACHE_PATH} "
              f"(use --rebuild-real to force re-read) ...")
        with open(REAL_CACHE_PATH, "rb") as _f:
            _cache = pickle.load(_f)
        stats                 = _cache["stats"]
        threshold_stats       = _cache["threshold_stats"]
        mean_stats             = _cache["mean_stats"]
        mean_thr_stats         = _cache["mean_thr_stats"]
        mean_filter_stats      = _cache["mean_filter_stats"]
        mean_filter_thr_stats  = _cache["mean_filter_thr_stats"]
        mean_gene_stats        = _cache["mean_gene_stats"]
        mean_gene_thr_stats    = _cache["mean_gene_thr_stats"]
        all_genes               = _cache["all_genes"]
        crossover                = _cache["crossover"]
        mle_mode_t_idx          = _cache["mle_mode_t_idx"]
    else:
        if rebuild_real:
            print("--rebuild-real: re-reading real condition files ...")
        stats            = {}; threshold_stats  = {}
        filter_stats     = {}; filter_thr_stats = {}
        gene_stats       = {}; gene_thr_stats   = {}
        mean_stats       = {}; mean_thr_stats   = {}
        mean_filter_stats = {}; mean_filter_thr_stats = {}
        mean_gene_stats  = {}; mean_gene_thr_stats = {}
        crossover        = {}
        all_genes        = []

        for condition, filepath in real.CONDITIONS.items():
            s, ts, fs, fts, gs, gts, _nrs, _nrts, _gnrs, _gnrts = real.read_one_condition(
                condition, filepath, pos_filter, pos_gene
            )
            stats[condition]            = s
            threshold_stats[condition]  = ts

            ms, mts = real.aggregate_stats({condition: s}, {condition: ts})
            mean_stats.update(ms)
            mean_thr_stats.update(mts)

            mfs, mfts = real.aggregate_filter_stats({condition: fs}, {condition: fts}, filter_lengths)
            mean_filter_stats.update(mfs)
            mean_filter_thr_stats.update(mfts)

            mgs, mgts, all_genes = real.aggregate_gene_stats(
                {condition: gs}, {condition: gts}, gene_lengths
            )
            mean_gene_stats.update(mgs)
            mean_gene_thr_stats.update(mgts)

            crossover.update(real.find_pp_sa_crossover(
                {condition: ms[condition]}, {condition: mts[condition]}
            ))

        real.write_crossover_csv(crossover, filename=_out("pp_sa_crossover_callrate.csv"))
        real.write_crossover_csv_tidy(crossover, filename=_out("pp_sa_crossover_callrate_tidy.csv"))
        mle_mode_t_idx = real._mode_crossover_threshold(crossover["pr_scaling=MLE"])

        print(f"Saving real-data cache to {REAL_CACHE_PATH} ...")
        with open(REAL_CACHE_PATH, "wb") as _f:
            pickle.dump({
                "stats":                 stats,
                "threshold_stats":       threshold_stats,
                "mean_stats":            mean_stats,
                "mean_thr_stats":        mean_thr_stats,
                "mean_filter_stats":     mean_filter_stats,
                "mean_filter_thr_stats": mean_filter_thr_stats,
                "mean_gene_stats":       mean_gene_stats,
                "mean_gene_thr_stats":   mean_gene_thr_stats,
                "all_genes":             all_genes,
                "crossover":             crossover,
                "mle_mode_t_idx":        mle_mode_t_idx,
            }, _f)

    # 3. Replicate data (sub20_only/) — separate cache, independent files.
    if not rebuild_reps and os.path.exists(REPS_CACHE_PATH):
        print(f"Loading replicate cache from {REPS_CACHE_PATH} "
              f"(use --rebuild-reps to force re-read) ...")
        with open(REPS_CACHE_PATH, "rb") as _f:
            _rep_cache = pickle.load(_f)
        gene_rep_stats     = _rep_cache["gene_rep_stats"]
        gene_rep_thr_stats = _rep_cache["gene_rep_thr_stats"]
    else:
        if rebuild_reps:
            print("--rebuild-reps: re-reading replicate files ...")
        gene_rep_stats, gene_rep_thr_stats = real.read_replicate_data(pos_gene)

        print(f"Saving replicate cache to {REPS_CACHE_PATH} ...")
        with open(REPS_CACHE_PATH, "wb") as _f:
            pickle.dump({
                "gene_rep_stats":     gene_rep_stats,
                "gene_rep_thr_stats": gene_rep_thr_stats,
            }, _f)

    # 4. Simulated data — separate cache, independent files. Crossover isn't
    #    computed here since sim_sa_pp_1v5 doesn't use it (unlike the real
    #    pipeline's mle_mode_t_idx dependency).
    if not rebuild_sims and os.path.exists(SIMS_CACHE_PATH):
        print(f"Loading simulated-data cache from {SIMS_CACHE_PATH} "
              f"(use --rebuild-sims to force re-read) ...")
        with open(SIMS_CACHE_PATH, "rb") as _f:
            _sim_cache = pickle.load(_f)
        sim_stats           = _sim_cache["stats"]
        sim_threshold_stats = _sim_cache["threshold_stats"]
        sim_mean_stats       = _sim_cache["mean_stats"]
        sim_mean_thr_stats   = _sim_cache["mean_thr_stats"]
    else:
        if rebuild_sims:
            print("--rebuild-sims: re-reading simulated condition files ...")
        sim_stats           = {}
        sim_threshold_stats = {}

        for condition, filepath in sim.CONDITIONS.items():
            s, ts = sim.read_one_condition(condition, filepath)
            sim_stats[condition]           = s
            sim_threshold_stats[condition] = ts

        sim_mean_stats, sim_mean_thr_stats = sim.aggregate_stats(sim_stats, sim_threshold_stats)

        print(f"Saving simulated-data cache to {SIMS_CACHE_PATH} ...")
        with open(SIMS_CACHE_PATH, "wb") as _f:
            pickle.dump({
                "stats":           sim_stats,
                "threshold_stats": sim_threshold_stats,
                "mean_stats":      sim_mean_stats,
                "mean_thr_stats":  sim_mean_thr_stats,
            }, _f)

    # 5. Plots
    sa_pp_1v5(mean_stats, mean_thr_stats)
    sim_sa_pp_1v5(sim_mean_stats, sim_mean_thr_stats)
    sample_metadata_1v5(stats, threshold_stats, run_metadata)
    quality_flag_1v5(mean_filter_stats, mean_filter_thr_stats)
    gene_region_1v5(mean_gene_stats, mean_gene_thr_stats, all_genes)
    replicate_acc_diff(gene_rep_stats, gene_rep_thr_stats, all_genes, mle_mode_t_idx, run_metadata)
