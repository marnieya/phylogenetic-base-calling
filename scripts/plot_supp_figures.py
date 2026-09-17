"""
Supplemental ROC plots: "spaghetti" figures overlaying every Run's (or
sample_id's, for simulated data) own SA point and PP curve directly, with no
cross-run averaging -- this turned out to be more informative than a single
mean curve per condition (which is what an earlier version of this script,
and every figure in plot_main_figures.py, plots instead), since it shows the
actual spread across runs rather than collapsing it away. Panels use
independent axes (own x/y range, no shared scale across panels) and no
axis breaks, unlike plot_main_figures.py's figures.

Reuses roc_real_data.py / roc_sim_data.py (data pipelines) and
plot_common.py (generic plotting utilities) directly, exactly like
plot_main_figures.py does -- and also imports plot_main_figures.py itself,
to reuse its _make_row_grid/_add_row_labels layout helpers and its pickle
caches (REAL_CACHE_PATH, SIMS_CACHE_PATH) rather than re-reading the call
files or duplicating layout code.

Caching: this script builds no cache of its own. It expects
plot_main_figures.py to have already been run at least once and simply
loads its two pickles directly -- mainly their raw, per-Run/per-sample_id
`stats`/`threshold_stats` accumulators (for the *_by_run figures), plus the
real cache's `mean_stats`/`mean_thr_stats` (aggregate_stats()'s cross-run
means, for pp_sa_crossover_accuracy.csv below -- the *_by_run figures don't
use these). Run `python plot_main_figures.py` (add --rebuild-real and/or
--rebuild-sims the first time, or whenever the underlying call files
change) before running this script.

Also produces pp_sa_crossover_accuracy.csv (real data): repurposes
plot_main_figures.py's crossover-CSV machinery
(roc_real_data.write_crossover_csv) for a different metric --
find_pp_sa_crossover_accuracy's least-stringent-threshold-meeting-SA's-
accuracy definition instead of find_pp_sa_crossover's most-stringent-
threshold-meeting-SA's-call-rate one (see roc_real_data.py's docstrings and
CLAUDE.md's "Crossover analysis" section for how the two differ). The
plain, call-rate-based version is pp_sa_crossover_callrate.csv, produced by
plot_main_figures.py.

Also produces a "coverage plot" figure (position vs. avg_depth, plus a
coverage-by-quality-flag row) from coverage_only/'s aggregated coverage
tables, via coverage_data.py -- a separate, plot-agnostic data pipeline
paralleling roc_real_data.py/roc_sim_data.py but for
coverage_only/aggregate_coverage.sh's output rather than the call files.
See coverage_data.py and coverage_only/CLAUDE.md for how those tables are
built/read; see the "Coverage plots" section below for how it's drawn.
This part of the script builds its own small pickle cache
(COVERAGE_CACHE_PATH), separate from plot_main_figures.py's, since reading
the two coverage tables it needs (real, sim -- the replicate table is no
longer read here, see "Coverage plots" below) from scratch takes a while
(they're large: ~2.4M / ~5.9M rows).

Produces plot files written to OUTPUT_DIR ("supp_figures/"):
sa_pp_all_by_run.pdf, one sim_sa_pp_err_<error_rate>_by_run.pdf per
error_rate the simulated data has (three as of this writing -- beta0.001,
beta0.004, beta0.007 -- but discovered from the data, not hardcoded),
pp_sa_crossover_accuracy.csv, and coverage_all.pdf (simulated/real coverage
stacked into one 3-row figure -- the third row splitting real coverage by
quality-flag category (pass/caution/mask); skipped, with a message, if
coverage_only/'s aggregated tables aren't present). More supplemental
figures (drawing on other, currently-dormant plot_roc.py /
plot_roc_simulated.py functions) are expected to be added here later.
"""

import os
import random
import sys

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D

import roc_real_data as real
import roc_sim_data as sim
import plot_main_figures as mainfig
import coverage_data
from plot_common import N_THRESHOLDS, _PP_MARKER, GENOME_LENGTH

OUTPUT_DIR = "supp_figures"
COVERAGE_CACHE_PATH = "plot_supp_coverage_cache.pkl"  # --rebuild-coverage to force re-read

# Max columns for these grid figures -- matches replicate_acc_diff's column
# count (both go through mainfig._make_row_grid with the same n_cols=3),
# so this project's whole-genome/all-conditions supplemental figures share
# the same panel size and column width as replicate_acc_diff.pdf. Rows grow
# as needed to fit however many subsamples/depths are actually present --
# unlike the fixed 2-column _1v5 figures.
_ALL_GRID_N_COLS = 3

# "full" is real.SUBSAMPLE_ORDER's full-depth (unsubsampled) entry, not a
# true subsample -- out of scope for this analysis entirely, not just the
# _1v5 figures (which never included it in the first place, since they only
# ever select sub1/sub5 explicitly). sa_pp_all_by_run is the only place in
# the project that would otherwise sweep in every entry of SUBSAMPLE_ORDER,
# so this exclusion is applied locally here rather than touching the shared
# SUBSAMPLE_ORDER constant itself (still an accurate list of every
# subsample value the data can contain).
_REAL_SUBSAMPLES = [s for s in real.SUBSAMPLE_ORDER if s != "full"]

# Every run's PP line is thin (lw=0.8) and semi-transparent (alpha=0.3) so up
# to ~81 overlapping runs per condition per panel stay legible as a
# density/spread cloud rather than a solid smear; the diamond markers on
# those curves share the line's color+alpha but are kept small
# (_BY_RUN_DIAMOND_MS) so they mark each threshold without overpowering the
# line itself -- no gradient coloring here, unlike plot_main_figures.py's
# sa_pp family, since with this many overlapping lines a per-threshold color
# gradient wouldn't be legible anyway, and that encoding is already
# established by the main figures. SA gets its own, higher alpha
# (_BY_RUN_SA_ALPHA) and a larger size (_BY_RUN_SC) than the PP lines/
# diamonds -- SA is one point per run rather than a whole curve, so it needs
# to read clearly against the PP spaghetti rather than blend into it.
_BY_RUN_LW          = 0.8
_BY_RUN_SC          = 40
_BY_RUN_ALPHA        = 0.3
_BY_RUN_SA_ALPHA     = 0.5
_BY_RUN_DIAMOND_MS   = 4


def _out(filename):
    """Path for a plot output file, under OUTPUT_DIR (created on first use)."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    return os.path.join(OUTPUT_DIR, filename)


def _subsample_has_data_by_run(per_run, subsample):
    """Whether any run has non-NaN SA or PP data for this subsample, across
    any condition in `per_run` ({condition: {run: (mean_s, mean_ts)}}, each
    (mean_s, mean_ts) from real._aggregate_for_run_subset called on a
    singleton run_subset) -- mean_s/mean_ts always carry a key for every
    entry in real.SUBSAMPLE_ORDER regardless of whether any run actually has
    that subsample (aggregate_stats's _safe_mean([]) fills a NaN sentinel
    instead), so a plain key check isn't enough to tell whether a subsample
    is really present."""
    for by_run in per_run.values():
        for mean_s, mean_ts in by_run.values():
            sa_acc, _ = mean_s[subsample].get("SA", (float("nan"), float("nan")))
            if sa_acc == sa_acc:
                return True
            for t_idx in range(N_THRESHOLDS):
                acc, _ = mean_ts[subsample].get(t_idx, (float("nan"), float("nan")))
                if acc == acc:
                    return True
    return False


def _draw_sa_pp_panel_by_run(ax, sa_points, curves_by_condition, sa_color, condition_colors):
    """Draw one independent-axis 'by run' sa_pp panel onto `ax`: every run's
    own SA point and PP curve overlaid directly (no cross-run averaging) --
    one thin, semi-transparent (_BY_RUN_ALPHA) line per run per condition
    (with a small diamond marker (_BY_RUN_DIAMOND_MS) at each threshold,
    same color and alpha as its line -- for visibility only, not a
    gradient: with this many overlapping runs a per-threshold color
    encoding wouldn't be legible anyway, and it's already established by
    plot_main_figures.py's sa_pp family), one point per run for SA, larger
    and more opaque (_BY_RUN_SC, _BY_RUN_SA_ALPHA) than the PP lines/
    diamonds so it still reads clearly as one point per run against the PP
    spaghetti rather than blending into it. Legend proxies (one zero-length
    SA scatter, one zero-length PP line+diamond per condition, all
    full-opacity) are added last so their handles/labels are independent of
    how many runs got plotted. Axis limits come from every point/curve
    value actually plotted on this panel (5% pad).

    sa_points: list of (acc, cr), one per run.
    curves_by_condition: {condition: [(curve_acc, curve_cr), ...]} -- one
    (curve_acc, curve_cr) pair per run, each N_THRESHOLDS long.
    """
    panel_y, panel_x = [], []
    for sa_acc, sa_cr in sa_points:
        if sa_acc == sa_acc and sa_cr == sa_cr:
            ax.scatter(sa_cr, sa_acc, marker="o", color=sa_color, s=_BY_RUN_SC,
                       alpha=_BY_RUN_SA_ALPHA, linewidths=0, zorder=5)
            panel_y.append(sa_acc)
            panel_x.append(sa_cr)

    for condition, run_curves in curves_by_condition.items():
        color = condition_colors[condition]
        for curve_acc, curve_cr in run_curves:
            ax.plot(curve_cr, curve_acc, color=color, linewidth=_BY_RUN_LW,
                     alpha=_BY_RUN_ALPHA, marker=_PP_MARKER, markersize=_BY_RUN_DIAMOND_MS,
                     zorder=4)
            panel_y.extend(a for a in curve_acc if a == a)
            panel_x.extend(c for c in curve_cr if c == c)

    # Legend proxies -- full opacity, independent of the real (semi-transparent)
    # draws above, same zero-length-artist technique the sa_pp family's own
    # legend proxies use. Kept at the project's standard legend marker sizes
    # (not _BY_RUN_SC/_BY_RUN_DIAMOND_MS) since a legend swatch is read in
    # isolation, not against a cloud of overlapping points.
    ax.scatter([], [], marker="o", color=sa_color, s=_BY_RUN_SC * 2, label="SA")
    for condition, color in condition_colors.items():
        ax.plot([], [], color=color, linewidth=1.5, marker=_PP_MARKER, markersize=6,
                label=f"PP ({condition})")

    if panel_y:
        lo, hi = min(panel_y), max(panel_y)
        pad = (hi - lo) * 0.05 or 0.01
        ax.set_ylim(lo - pad, hi + pad)
    if panel_x:
        lo, hi = min(panel_x), max(panel_x)
        pad = (hi - lo) * 0.05 or 0.01
        ax.set_xlim(lo - pad, hi + pad)
    ax.ticklabel_format(useOffset=False)


def _finish_sa_pp_grid(fig, last_ax, legend_y, lfs, filename):
    """Shared legend + save step: bottom-anchored (mainfig._make_row_grid's
    plain convention), matching replicate_acc_diff's width exactly. Every
    panel produces identically-labeled proxy handles (SA, PP per
    condition), so handles are pulled from whichever panel was drawn last."""
    handles, labels = last_ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), fontsize=lfs, loc="upper center",
               bbox_to_anchor=(0.5, legend_y), ncol=len(by_label))
    plt.savefig(filename, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# sa_pp_all_by_run (real data)
# ---------------------------------------------------------------------------

def sa_pp_all_by_run(stats, threshold_stats):
    """
    SA + PP curves across every subsample real data has (_REAL_SUBSAMPLES --
    real.SUBSAMPLE_ORDER minus "full", which isn't a true subsample and is
    out of scope for this analysis -- filtered to those actually present:
    sub0.8, sub1, sub3, sub5, sub10, sub20, in that order -- real.SUBSAMPLE_ORDER
    is deliberately just these 6, not the wider 11-depth sweep the coverage
    pipeline uses, since the real call files were only ever base-called at
    these depths), with every Run's own SA point and PP curve plotted
    directly rather than one cross-run mean per condition (see
    _draw_sa_pp_panel_by_run). Panels wrap at
    _ALL_GRID_N_COLS (3) columns, rows growing as needed; panel sizing
    comes from mainfig._make_row_grid, so panels render at exactly
    (mainfig._ROW_PANEL_WIDTH_IN, mainfig._ROW_PANEL_HEIGHT_IN), matching
    every other figure in the project. No y-axis break, no shared range
    across panels -- each panel is a single, independent axis (own x/y
    range) -- this is a supplemental figure meant to show all the
    underlying data, so forcing a shared range/break convention across 6
    subsamples spanning very different accuracy levels would obscure more
    than it reveals.
    Uses the *raw*, per-Run stats/threshold_stats accumulators (already
    saved in plot_main_figures.py's real-data cache, just not previously
    loaded by this script's __main__) via
    real._aggregate_for_run_subset(..., {run}) per run -- a singleton
    run_subset makes its existing cross-run averaging degenerate to that
    one run's own (acc, cr), with all the same NaN/missing-data handling
    already built in.
    Output: sa_pp_all_by_run.pdf
    """
    conditions = {k: v for k, v in real.CONDITIONS.items() if k != "pr_scaling=0.001"}
    condition_colors = {c: real.ALL_CONDITION_COLORS[c] for c in conditions}
    first_cond = next(iter(conditions))
    sa_color = cm.tab10.colors[1]
    lfs = 8

    # One _aggregate_for_run_subset({run}) call per (condition, run) --
    # already loops every subsample internally, so this covers every
    # panel's data for that run in one call rather than recomputing per
    # (run, subsample) pair.
    per_run = {
        condition: {
            run: real._aggregate_for_run_subset(stats[condition], threshold_stats[condition], {run})
            for run in stats.get(condition, {})
        }
        for condition in conditions
    }
    if not any(per_run[c] for c in conditions):
        print("sa_pp_all_by_run: no data available, skipping.")
        return

    subsamples = [s for s in _REAL_SUBSAMPLES if _subsample_has_data_by_run(per_run, s)]
    if not subsamples:
        print("sa_pp_all_by_run: no data available, skipping.")
        return

    n_cols = min(_ALL_GRID_N_COLS, len(subsamples))
    n_rows = (len(subsamples) + n_cols - 1) // n_cols
    fig, axes, legend_y = mainfig._make_row_grid(n_rows, n_cols)

    last_ax = None
    for i, subsample in enumerate(subsamples):
        row, col = divmod(i, n_cols)
        ax = axes[row, col]

        sa_points = [
            mean_s[subsample].get("SA", (float("nan"), float("nan")))
            for mean_s, _ in per_run[first_cond].values()
        ]
        curves_by_condition = {}
        for condition in conditions:
            run_curves = []
            for _, mean_ts in per_run[condition].values():
                curve_acc, curve_cr = [], []
                for t_idx in range(N_THRESHOLDS):
                    acc, cr = mean_ts[subsample].get(t_idx, (float("nan"), float("nan")))
                    curve_acc.append(acc)
                    curve_cr.append(cr)
                run_curves.append((curve_acc, curve_cr))
            curves_by_condition[condition] = run_curves

        _draw_sa_pp_panel_by_run(ax, sa_points, curves_by_condition, sa_color, condition_colors)
        ax.set_title(real._fmt_subsample(subsample), fontsize=lfs + 1)
        ax.set_xlabel("Mean call rate", fontsize=lfs)
        ax.tick_params(labelsize=lfs - 1)
        if col == 0:
            ax.set_ylabel("Mean accuracy", fontsize=lfs)
        last_ax = ax

    # Hide any unfilled trailing panels (e.g. 5 subsamples in a 3-col grid
    # leaves the last row's 3rd panel empty).
    for i in range(len(subsamples), n_rows * n_cols):
        row, col = divmod(i, n_cols)
        axes[row, col].set_visible(False)

    filename = _out("sa_pp_all_by_run.pdf")
    _finish_sa_pp_grid(fig, last_ax, legend_y, lfs, filename)
    print(f"Written {filename}")


# ---------------------------------------------------------------------------
# sim_sa_pp_err_<error_rate>_by_run (simulated data, one file per error_rate)
# ---------------------------------------------------------------------------

def _draw_sim_sa_pp_for_error_rate_by_run(stats, threshold_stats, conditions, condition_colors,
                                           error_rate, filename):
    """Shared body for one sim_sa_pp_err_<error_rate>_by_run.pdf file: every
    sample's own SA point and PP curve overlaid directly (no cross-run
    averaging), across every mean_depth present for this one error_rate
    (columns wrap at _ALL_GRID_N_COLS, rows growing as needed) -- see
    sa_pp_all_by_run / _draw_sa_pp_panel_by_run for the shared per-run
    drawing convention. Unlike the real pipeline, no
    _aggregate_for_run_subset-style helper is needed to pull one run's own
    values: roc_sim_data's raw stats/threshold_stats are already keyed
    directly by run, so sim.per_run_stats/per_run_thresholds just do the
    safe division directly.
    A single row label above the top row names the error_rate
    (mainfig._add_row_labels, reusing the same mechanism sim_sa_pp_1v5's
    block labels and gene_region_1v5's row headers use), since the figure
    would otherwise carry no on-page indication of which error_rate it is
    once separated from its filename (e.g. once embedded elsewhere)."""
    sa_color = cm.tab10.colors[1]
    lfs = 8
    first_cond = next(iter(conditions))

    depths = sorted(
        {d for c in conditions for run in stats.get(c, {})
         for d in stats[c].get(run, {}).get(error_rate, {})},
        key=sim._depth_sort_key,
    )
    if not depths:
        print(f"{filename}: no data available for error_rate={error_rate}, skipping.")
        return

    n_cols = min(_ALL_GRID_N_COLS, len(depths))
    n_rows = (len(depths) + n_cols - 1) // n_cols
    fig, axes, legend_y = mainfig._make_row_grid(n_rows, n_cols)

    last_ax = None
    for i, depth in enumerate(depths):
        row, col = divmod(i, n_cols)
        ax = axes[row, col]

        sa_points = [
            sim.per_run_stats(stats[first_cond], run, error_rate, depth)["SA"]
            for run in stats.get(first_cond, {})
        ]
        curves_by_condition = {}
        for condition in conditions:
            run_curves = []
            for run in stats.get(condition, {}):
                per_t = sim.per_run_thresholds(threshold_stats[condition], run, error_rate, depth)
                curve_acc = [per_t[t_idx][0] for t_idx in range(N_THRESHOLDS)]
                curve_cr  = [per_t[t_idx][1] for t_idx in range(N_THRESHOLDS)]
                run_curves.append((curve_acc, curve_cr))
            curves_by_condition[condition] = run_curves

        _draw_sa_pp_panel_by_run(ax, sa_points, curves_by_condition, sa_color, condition_colors)
        ax.set_title(sim._fmt_depth(depth), fontsize=lfs + 1)
        ax.set_xlabel("Mean call rate", fontsize=lfs)
        ax.tick_params(labelsize=lfs - 1)
        if col == 0:
            ax.set_ylabel("Mean accuracy", fontsize=lfs)
        last_ax = ax

    for i in range(len(depths), n_rows * n_cols):
        row, col = divmod(i, n_cols)
        axes[row, col].set_visible(False)

    mainfig._add_row_labels(fig, axes, lfs, [sim._fmt_error_rate(error_rate)])
    _finish_sa_pp_grid(fig, last_ax, legend_y, lfs, filename)
    print(f"Written {filename}")


def sim_sa_pp_by_error_rate_by_run(stats, threshold_stats):
    """
    One sim_sa_pp_err_<error_rate>_by_run.pdf file per error_rate the
    simulated data has -- discovered from the data, not hardcoded (three as
    of this writing: beta0.001, beta0.004, beta0.007) -- rather than one
    figure stacking every error_rate's block, which got unwieldy once there
    turned out to be three error rates. Each file covers every mean_depth
    present for that one error_rate (roc_sim_data._depth_sort_key, not
    assumed to match the real pipeline's subsamples), wrapped at 3 columns,
    with every sample_id's own curve plotted directly rather than one
    cross-run mean per condition. See _draw_sim_sa_pp_for_error_rate_by_run.
    """
    conditions = {k: v for k, v in sim.CONDITIONS.items() if k != "pr_scaling=0.001"}
    all_condition_colors = {
        cond: col for cond, col in zip(sim.CONDITIONS.keys(), cm.tab10.colors[:len(sim.CONDITIONS)])
    }
    condition_colors = {c: all_condition_colors[c] for c in conditions}

    error_rates = sorted({
        er for c in conditions for run in stats.get(c, {}) for er in stats[c][run]
    })
    if not error_rates:
        print("sim_sa_pp_by_error_rate_by_run: no data available, skipping.")
        return

    for error_rate in error_rates:
        filename = _out(f"sim_sa_pp_err_{error_rate}_by_run.pdf")
        _draw_sim_sa_pp_for_error_rate_by_run(stats, threshold_stats, conditions,
                                               condition_colors, error_rate, filename)


# ---------------------------------------------------------------------------
# pp_sa_crossover_accuracy.csv (real data) -- accuracy-based counterpart to
# plot_main_figures.py's pp_sa_crossover_callrate.csv
# ---------------------------------------------------------------------------

def pp_sa_crossover_accuracy(mean_stats, mean_thr_stats):
    """
    Repurposes plot_main_figures.py's crossover-CSV machinery
    (roc_real_data.write_crossover_csv) for a different metric:
    roc_real_data.find_pp_sa_crossover_accuracy instead of
    find_pp_sa_crossover -- least-stringent threshold whose *accuracy*
    meets or exceeds SA's, rather than most-stringent threshold whose
    *call rate* meets or exceeds SA's. See that function's docstring, and
    CLAUDE.md's "Crossover analysis" section, for the two definitions and
    why they scan in opposite directions. Both files carry the same four
    value columns -- pp_call_rate, sa_call_rate, pp_accuracy, sa_accuracy --
    at each file's own crossover threshold; they differ only in which metric
    drove the threshold search (and thus which threshold each row reports).
    Also writes pp_sa_crossover_accuracy_tidy.csv via
    real.write_crossover_csv_tidy -- a condensed, presentation-ready
    version: pr_scaling=MLE only, "full" excluded, threshold as an exact
    "1-Ce-N" delta string, call rates rounded to 4 decimals and accuracies
    to 6 (accuracy values here all sit within ~3e-5 of 1.0, so 4 decimals
    would round every row to an identical-looking 1.0000).
    Output: pp_sa_crossover_accuracy.csv, pp_sa_crossover_accuracy_tidy.csv
    """
    crossover = real.find_pp_sa_crossover_accuracy(mean_stats, mean_thr_stats)
    real.write_crossover_csv(crossover, filename=_out("pp_sa_crossover_accuracy.csv"))
    real.write_crossover_csv_tidy(crossover, filename=_out("pp_sa_crossover_accuracy_tidy.csv"))


# ---------------------------------------------------------------------------
# Coverage plots (position vs. avg_depth, from coverage_only/'s aggregated
# tables -- see coverage_data.py and coverage_only/CLAUDE.md)
# ---------------------------------------------------------------------------

_COVERAGE_FIG_WIDTH_IN  = 10.0
_COVERAGE_FIG_HEIGHT_IN = 3.0
_COVERAGE_LW          = 1.0

# Shared x-axis for every coverage plot -- built once since it's identical
# (position 1..GENOME_LENGTH) across all three data sources.
_COVERAGE_X = list(range(1, GENOME_LENGTH + 1))

# One color per depth column, shared between coverage_real and coverage_sim
# (see coverage_data.ALL_DEPTH_COLS) so e.g. avg_depth_5x is the same color
# in both plots -- built once, keyed by column name rather than by each
# plot's own local index.
_COVERAGE_DEPTH_CMAP = cm.tab20 if len(coverage_data.ALL_DEPTH_COLS) > 10 else cm.tab10
_COVERAGE_DEPTH_COLORS = {
    c: _COVERAGE_DEPTH_CMAP.colors[i % len(_COVERAGE_DEPTH_CMAP.colors)]
    for i, c in enumerate(coverage_data.ALL_DEPTH_COLS)
}


def _draw_coverage_real(ax, depth_cols, mean_by_depth):
    """Draws coverage_real's content onto `ax`: one line per depth (color =
    depth), y = mean avg_depth across Runs at each position.
    "avg_depth_full" (the complete, unsubsampled depth) is excluded, same
    as sa_pp_all_by_run and every other _REAL_SUBSAMPLES-based figure --
    not a true subsample, out of scope for this analysis; its magnitude
    (thousands of reads) is also orders of magnitude above every subsampled
    depth's, which would otherwise squash them flat on a shared linear
    y-axis. Used by coverage_all as its "Real subsampled data" row (no
    standalone real-only coverage figure is produced)."""
    lfs = 8
    plot_cols = sorted(
        (c for c in depth_cols if c != "avg_depth_full"),
        key=coverage_data.depth_col_sort_key,
    )

    for c in plot_cols:
        ax.plot(_COVERAGE_X, mean_by_depth[c], color=_COVERAGE_DEPTH_COLORS[c],
                 linewidth=_COVERAGE_LW, label=coverage_data.fmt_depth_col(c))
    ax.set_xlabel("Position", fontsize=lfs)
    ax.set_ylabel("Mean avg_depth (across Runs)", fontsize=lfs)
    ax.tick_params(labelsize=lfs - 1)
    ax.legend(fontsize=lfs, title="Depth", title_fontsize=lfs, ncol=2,
              loc="upper left", bbox_to_anchor=(1.01, 1.0))


def _draw_coverage_sim(ax, depth_cols, error_rates, mean_by_error_depth):
    """Draws coverage_sim's content onto `ax`: color = depth, linestyle =
    error_rate, y = mean avg_depth across sample_ids at each position. Two
    independent legends (one per encoding, same convention the rest of
    this project uses for dual color/marker encodings -- e.g.
    replicate_acc_diff's CenterName color + replicate marker) rather than
    one combined legend with an entry per (depth, error_rate) combination.
    Used by coverage_all as its "Simulated data" row (no standalone
    sim-only coverage figure is produced)."""
    lfs = 8
    plot_cols = sorted(depth_cols, key=coverage_data.depth_col_sort_key)

    linestyle_cycle = ["-", "--", ":", "-."]
    linestyles = {er: linestyle_cycle[i % len(linestyle_cycle)] for i, er in enumerate(error_rates)}

    for er in error_rates:
        for c in plot_cols:
            ax.plot(_COVERAGE_X, mean_by_error_depth[er][c], color=_COVERAGE_DEPTH_COLORS[c],
                     linestyle=linestyles[er], linewidth=_COVERAGE_LW)
    ax.set_xlabel("Position", fontsize=lfs)
    ax.set_ylabel("Mean avg_depth (across sample_ids)", fontsize=lfs)
    ax.tick_params(labelsize=lfs - 1)

    color_handles = [Line2D([0], [0], color=_COVERAGE_DEPTH_COLORS[c], linewidth=_COVERAGE_LW,
                              label=coverage_data.fmt_depth_col(c)) for c in plot_cols]
    style_handles = [Line2D([0], [0], color="black", linewidth=_COVERAGE_LW,
                              linestyle=linestyles[er], label=er) for er in error_rates]
    color_legend = ax.legend(handles=color_handles, fontsize=lfs, title="Depth",
                              title_fontsize=lfs, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.add_artist(color_legend)
    ax.legend(handles=style_handles, fontsize=lfs, title="Error rate", title_fontsize=lfs,
              loc="lower left", bbox_to_anchor=(1.01, 0.0))


def _subsample_to_depth_col(subsample):
    """'sub0.8' -> 'avg_depth_0.8x' -- the coverage file's column-naming
    convention for that same depth (coverage_data.ALL_DEPTH_COLS). Used by
    the "coverage by quality flag" row to read/plot only the subsample
    depths real.SUBSAMPLE_ORDER actually uses (the overlap between
    coverage_only's wider depth sweep and the depths base-calling was
    actually run at -- same _REAL_SUBSAMPLES list every other real-data
    figure in this project uses), not coverage_data.ALL_DEPTH_COLS' full
    set."""
    return f"avg_depth_{subsample[3:]}x"


# Each panel covers one pair of subsample depths -- a deliberate content
# grouping (roughly "low/medium/high" of the 6 real depths), not a generic
# every-Nth slice of _REAL_SUBSAMPLES, so it's spelled out explicitly rather
# than sliced programmatically (stays correct/obvious even if
# _REAL_SUBSAMPLES's own order or membership ever changes).
_COVERAGE_DEPTH_GROUPS = [
    ("sub0.8", "sub1"),
    ("sub3", "sub5"),
    ("sub10", "sub20"),
]

# Quality-flag color: same tab10-by-FILTER_VALUES-order assignment
# plot_main_figures._draw_quality_flag_figure uses for quality_flag_1v5.pdf
# (flag_colors), recomputed locally (that dict is private to
# plot_main_figures) so a given FILTER value renders in the same color in
# both figures.
_COVERAGE_FILTER_COLORS = {v: cm.tab10.colors[i] for i, v in enumerate(real.FILTER_VALUES)}

_COVERAGE_VIOLIN_ALPHA         = 0.3
_COVERAGE_VIOLIN_WIDTH         = 0.22
_COVERAGE_VIOLIN_POINT_SC      = 14
_COVERAGE_VIOLIN_POINT_ALPHA   = 0.7
_COVERAGE_VIOLIN_X_JITTER      = 0.06
_COVERAGE_VIOLIN_GROUP_GAP     = 1.3   # x spacing between a panel's two depths
_COVERAGE_VIOLIN_FILTER_OFFSET = 0.28  # x spacing between the 3 filters at one depth


def _draw_coverage_by_filter(axes, per_run_filter_depth, depth_groups=_COVERAGE_DEPTH_GROUPS):
    """Draws the "coverage by quality flag" row onto `axes` (one ax per
    entry of `depth_groups`, in that order): each panel covers one pair of
    subsample depths (_COVERAGE_DEPTH_GROUPS -- roughly low/medium/high of
    the 6 real depths), and at each depth, one violin+jittered-scatter pair
    per real.FILTER_VALUES entry (pass/caution/mask), offset side by side
    (_COVERAGE_VIOLIN_FILTER_OFFSET) so the three quality-flag categories
    are directly comparable at that same target depth. Color = quality flag
    (_COVERAGE_FILTER_COLORS) for both the violin body (semi-transparent,
    _COVERAGE_VIOLIN_ALPHA) and its points -- unlike an earlier version of
    this row, which colored points by CenterName and had one panel per
    filter value instead of per depth-group.

    Each point is one Run's own mean avg_depth at that subsample depth,
    averaged across every position carrying that FILTER value
    (per_run_filter_depth, from coverage_data.read_real_coverage_by_filter)
    -- i.e. one summary value per (Run, filter_value, depth), not a
    per-position value, so a violin has as many points as there are Runs
    with data for that (filter_value, depth) combination. random.seed(42)
    before jittering, same convention plot_main_figures.replicate_acc_diff
    uses for its own x-jitter, so re-running reproduces the same layout.

    Each panel keeps its own independent y-axis (no sharey= across panels,
    unlike this row's first version) -- the three depth-groups can sit at
    very different absolute depths, so a shared scale would flatten the
    higher-depth panels' own pass/caution/mask spread; direct comparison
    across quality-flag categories still holds *within* each panel, which
    is where it matters (comparing pass vs. mask coverage at the same
    target depth), not across depth-groups."""
    lfs = 8
    random.seed(42)
    n_filters = len(real.FILTER_VALUES)
    # Symmetric offsets around each depth's own x position, e.g. n=3 ->
    # [-offset, 0, +offset], one per FILTER_VALUES entry in order.
    filter_offsets = [
        (j - (n_filters - 1) / 2) * _COVERAGE_VIOLIN_FILTER_OFFSET for j in range(n_filters)
    ]

    for ax, depths in zip(axes, depth_groups):
        depth_cols = [_subsample_to_depth_col(d) for d in depths]
        base_x = [i * _COVERAGE_VIOLIN_GROUP_GAP for i in range(len(depths))]

        for j, filt in enumerate(real.FILTER_VALUES):
            color = _COVERAGE_FILTER_COLORS[filt]
            violin_data = []
            violin_positions = []
            for i, depth_col in enumerate(depth_cols):
                values = [
                    filt_dict[filt][depth_col]
                    for filt_dict in per_run_filter_depth.values()
                    if filt in filt_dict and depth_col in filt_dict[filt]
                ]
                if not values:
                    continue
                x0 = base_x[i] + filter_offsets[j]
                violin_data.append(values)
                violin_positions.append(x0)
                xs = [x0 + random.uniform(-_COVERAGE_VIOLIN_X_JITTER, _COVERAGE_VIOLIN_X_JITTER)
                      for _ in values]
                ax.scatter(xs, values, color=color, s=_COVERAGE_VIOLIN_POINT_SC,
                           alpha=_COVERAGE_VIOLIN_POINT_ALPHA, zorder=3)

            if violin_data:
                parts = ax.violinplot(
                    violin_data, positions=violin_positions,
                    widths=_COVERAGE_VIOLIN_WIDTH, showmeans=False, showmedians=False,
                    showextrema=False,
                )
                for body in parts["bodies"]:
                    body.set_facecolor(color)
                    body.set_edgecolor(color)
                    body.set_alpha(_COVERAGE_VIOLIN_ALPHA)

        ax.set_xticks(base_x)
        ax.set_xticklabels([f"{d[3:]}x" for d in depths])
        ax.set_xlabel("Subsample depth", fontsize=lfs)
        ax.set_title(" / ".join(f"{d[3:]}x" for d in depths), fontsize=lfs + 1)
        ax.tick_params(labelsize=lfs - 1)

    axes[0].set_ylabel("Mean avg_depth (per sample)", fontsize=lfs)

    handles = [Line2D([0], [0], marker="o", color="w",
                       markerfacecolor=_COVERAGE_FILTER_COLORS[v], markersize=lfs, label=v)
               for v in real.FILTER_VALUES]
    axes[-1].legend(handles=handles, fontsize=lfs, title="Quality flag", title_fontsize=lfs,
                     loc="upper left", bbox_to_anchor=(1.01, 1.0))


_COVERAGE_ALL_LABELS = [
    "Simulated data",
    "Real subsampled data",
]
_COVERAGE_BY_FILTER_LABEL = "Real subsampled data, by quality-flag category"

# Row-height weights passed to add_gridspec -- the filter row is taller than
# the other two (the user asked for it explicitly, once it grew 3 side-by-
# side violin+scatter groups per panel instead of a single one).
_COVERAGE_ALL_ROW_HEIGHT_RATIOS = [1, 1, 1.8]


def coverage_all(real_depth_cols, real_mean_by_depth,
                  sim_depth_cols, sim_error_rates, sim_mean_by_error_depth,
                  per_run_filter_depth):
    """
    Simulated / real coverage, plus a "coverage by quality flag" row,
    stacked into one 3-row figure -- the only coverage-plot output this
    script produces (earlier versions also wrote the first two rows as
    their own standalone coverage_sim.pdf/coverage_real.pdf, and included a
    third "Replicates ..." row via coverage_reps.pdf/_draw_coverage_reps;
    both dropped -- the standalone files as redundant once this combined
    figure covered the same content, the replicates row at the user's
    request to keep just the simulated/real rows above the filter row
    below). The first two rows are each drawn by their own single-ax helper
    (_draw_coverage_sim/_draw_coverage_real) rather than duplicating any
    drawing logic -- same per-row width as those helpers' original
    standalone figures (_COVERAGE_FIG_WIDTH_IN), with a section label above
    each row (_COVERAGE_ALL_LABELS) naming which data source it is, via a
    plain ax.set_title(loc="left") since each of those rows is a single ax
    spanning the full row width. Row heights follow
    _COVERAGE_ALL_ROW_HEIGHT_RATIOS, not a flat 1:1:1 split -- the filter
    row is taller, at the user's request, to give its 3x2-group violin
    panels more room.

    The third row is instead three side-by-side panels (one per
    _COVERAGE_DEPTH_GROUPS entry), drawn by _draw_coverage_by_filter --
    since ax.set_title(loc="left") isn't available for a row with no single
    spanning ax, its own row label (_COVERAGE_BY_FILTER_LABEL) is placed on
    a fourth, invisible, full-width ax stacked above the three panels (a
    nested subgridspec, height_ratios=[label strip, panels]) rather than
    the pixel-measured fig.text() placement sim_sa_pp_1v5's block labels /
    _add_row_labels use elsewhere in this project -- this reserved-strip
    approach needs no post-draw pixel measurement, at the cost of only
    working for a short one-line label (a taller label could overflow its
    fixed-height strip, which pixel measurement would size correctly).
    Output: coverage_all.pdf
    """
    lfs = 8
    height_ratios = _COVERAGE_ALL_ROW_HEIGHT_RATIOS
    n_rows = len(height_ratios)
    # constrained_layout, not tight_layout() -- tight_layout() doesn't
    # correctly size a mix of plain gridspec cells (rows 0-1) and a nested
    # subgridspec (row 2's label-strip-over-panels split below), and was
    # observed to collapse rows on top of each other; constrained_layout is
    # matplotlib's layout engine built to handle nested gridspecs correctly.
    fig = plt.figure(
        figsize=(_COVERAGE_FIG_WIDTH_IN, _COVERAGE_FIG_HEIGHT_IN * sum(height_ratios)),
        constrained_layout=True,
    )
    gs = fig.add_gridspec(n_rows, 3, height_ratios=height_ratios)

    ax_sim  = fig.add_subplot(gs[0, :])
    ax_real = fig.add_subplot(gs[1, :])
    axes = (ax_sim, ax_real)

    _draw_coverage_sim(ax_sim, sim_depth_cols, sim_error_rates, sim_mean_by_error_depth)
    _draw_coverage_real(ax_real, real_depth_cols, real_mean_by_depth)

    for ax, label in zip(axes, _COVERAGE_ALL_LABELS):
        ax.set_title(label, fontsize=lfs + 2, loc="left")

    gs3 = gs[2, :].subgridspec(2, 3, height_ratios=[0.08, 1], hspace=0.05)
    ax_label = fig.add_subplot(gs3[0, :])
    ax_label.axis("off")
    ax_label.set_title(_COVERAGE_BY_FILTER_LABEL, fontsize=lfs + 2, loc="left")
    ax_f1 = fig.add_subplot(gs3[1, 0])
    ax_f2 = fig.add_subplot(gs3[1, 1])
    ax_f3 = fig.add_subplot(gs3[1, 2])
    _draw_coverage_by_filter((ax_f1, ax_f2, ax_f3), per_run_filter_depth)

    filename = _out("coverage_all.pdf")
    plt.savefig(filename, bbox_inches="tight")
    plt.close(fig)
    print(f"Written {filename}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    missing = [p for p in (mainfig.REAL_CACHE_PATH, mainfig.SIMS_CACHE_PATH)
               if not os.path.exists(p)]
    if missing:
        print("plot_supp_figures.py builds no cache of its own -- it expects "
              "plot_main_figures.py to have already been run (missing: "
              f"{', '.join(missing)}).")
        print("Run `python plot_main_figures.py` first (with --rebuild-real/"
              "--rebuild-sims the first time, or after the call files "
              "change), then re-run this script.")
        sys.exit(1)

    import pickle

    with open(mainfig.REAL_CACHE_PATH, "rb") as _f:
        _real_cache = pickle.load(_f)
    stats           = _real_cache["stats"]
    threshold_stats = _real_cache["threshold_stats"]
    mean_stats      = _real_cache["mean_stats"]
    mean_thr_stats  = _real_cache["mean_thr_stats"]

    with open(mainfig.SIMS_CACHE_PATH, "rb") as _f:
        _sims_cache = pickle.load(_f)
    sim_stats           = _sims_cache["stats"]
    sim_threshold_stats = _sims_cache["threshold_stats"]

    sa_pp_all_by_run(stats, threshold_stats)
    sim_sa_pp_by_error_rate_by_run(sim_stats, sim_threshold_stats)
    pp_sa_crossover_accuracy(mean_stats, mean_thr_stats)

    # Coverage plots -- separate, optional: skipped (not a hard error) if
    # coverage_only/'s aggregated tables aren't there yet.
    rebuild_coverage = "--rebuild-coverage" in sys.argv
    coverage_ready = os.path.exists(COVERAGE_CACHE_PATH) or all(
        os.path.exists(p) for p in (coverage_data.REAL_COVERAGE_PATH,
                                     coverage_data.SIM_COVERAGE_PATH)
    )
    if not coverage_ready:
        print("Skipping coverage plots -- coverage_only/'s aggregated tables not found "
              "(run coverage_only/aggregate_coverage.sh first).")
    else:
        # The flat depth list read_real_coverage_by_filter needs -- every
        # depth _COVERAGE_DEPTH_GROUPS' panels use, derived from that
        # grouping (not a separate list) so it can't drift out of sync with
        # what _draw_coverage_by_filter actually plots.
        filter_depths = [d for group in _COVERAGE_DEPTH_GROUPS for d in group]

        if not rebuild_coverage and os.path.exists(COVERAGE_CACHE_PATH):
            print(f"Loading coverage cache from {COVERAGE_CACHE_PATH} "
                  f"(use --rebuild-coverage to force re-read) ...")
            with open(COVERAGE_CACHE_PATH, "rb") as _f:
                _cov_cache = pickle.load(_f)
            real_depth_cols         = _cov_cache["real_depth_cols"]
            real_mean_by_depth      = _cov_cache["real_mean_by_depth"]
            sim_cov_depth_cols      = _cov_cache["sim_depth_cols"]
            sim_error_rates         = _cov_cache["sim_error_rates"]
            sim_mean_by_error_depth = _cov_cache["sim_mean_by_error_depth"]
            per_run_filter_depth    = _cov_cache["per_run_filter_depth"]
        else:
            if rebuild_coverage:
                print("--rebuild-coverage: re-reading coverage tables ...")
            real_depth_cols, real_mean_by_depth = coverage_data.read_real_coverage()
            sim_cov_depth_cols, sim_error_rates, sim_mean_by_error_depth = \
                coverage_data.read_sim_coverage()

            pos_filter = real.load_filter_metadata(real.FILTER_METADATA_PATH)
            filter_depth_cols = [_subsample_to_depth_col(d) for d in filter_depths]
            per_run_filter_depth = coverage_data.read_real_coverage_by_filter(
                pos_filter, filter_depth_cols
            )

            print(f"Saving coverage cache to {COVERAGE_CACHE_PATH} ...")
            with open(COVERAGE_CACHE_PATH, "wb") as _f:
                pickle.dump({
                    "real_depth_cols":         real_depth_cols,
                    "real_mean_by_depth":      real_mean_by_depth,
                    "sim_depth_cols":          sim_cov_depth_cols,
                    "sim_error_rates":         sim_error_rates,
                    "sim_mean_by_error_depth": sim_mean_by_error_depth,
                    "per_run_filter_depth":    per_run_filter_depth,
                }, _f)

        coverage_all(real_depth_cols, real_mean_by_depth,
                     sim_cov_depth_cols, sim_error_rates, sim_mean_by_error_depth,
                     per_run_filter_depth)
