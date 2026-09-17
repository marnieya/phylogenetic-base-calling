"""
Coverage-data pipeline: reads the three aggregated coverage tables produced
by coverage_only/aggregate_coverage.sh (see coverage_only/CLAUDE.md for how
those are built and their exact column layout) and computes the per-position
summaries plot_supp_figures.py's coverage plots need. Plot-agnostic, same
separation as roc_real_data.py / roc_sim_data.py.

All three input files are large (real: ~2.4M rows / 81 Runs, sim: ~5.9M rows
/ ~80 sample_ids x 3 error_rates, reps: ~23.9M rows / 800 (Run, replicate)
samples), so reading is manual-split streaming rather than csv.reader/
pandas, and the real/sim readers accumulate running (sum, count) per
position rather than holding every Run's/sample's own value -- only the
reps reader needs to keep every individual (Run, replicate) trace, since
those are plotted unaveraged.

Assumes (not re-validated row by row) that every Run/sample contributes
exactly one row per position 1..GENOME_LENGTH, matching how
aggregate_coverage.sh builds these files from per-position coverage_plot.txt
inputs.

read_real_coverage_by_filter is the exception to the "running (sum, count)
per position" streaming pattern above: it collapses across positions within
each (Run, FILTER value) group instead of across Runs within each position,
so its accumulators are bounded by Runs x filter values x depth columns
(tiny) rather than GENOME_LENGTH, and it needs no array('d')/array('L').
"""

from array import array

from plot_common import GENOME_LENGTH

COVERAGE_DIR = "coverage_only"
REAL_COVERAGE_PATH = f"{COVERAGE_DIR}/all_coverage_plots_real.txt"
SIM_COVERAGE_PATH  = f"{COVERAGE_DIR}/all_coverage_plots_sim.txt"
REPS_COVERAGE_PATH = f"{COVERAGE_DIR}/all_coverage_plots_reps.txt"

REP_SUBSAMPLE = "sub20"  # matches roc_real_data.REPLICATE_SUBSAMPLE

# Canonical depth-column order -- real's is the superset (sim's depths are
# a strict subset, sharing the same "avg_depth_<D>x" naming: sim never has
# 15x/10x/20x/full). plot_supp_figures.py assigns one color per entry here
# and looks colors up by name in both coverage_real and coverage_sim, so a
# shared depth (e.g. avg_depth_5x) gets the same color in both plots --
# assigning colors independently per plot by local enumeration index would
# silently drift apart once the two datasets have different column counts,
# even for depths they share.
ALL_DEPTH_COLS = [
    "avg_depth_20x", "avg_depth_15x", "avg_depth_10x", "avg_depth_5x",
    "avg_depth_4x", "avg_depth_3x", "avg_depth_2x", "avg_depth_1x",
    "avg_depth_0.8x", "avg_depth_0.5x", "avg_depth_0.2x", "avg_depth_full",
]


def depth_col_sort_key(col):
    """Sort key for an 'avg_depth_<D>x' / 'avg_depth_full' column name --
    numeric depths ascending, 'full' last (mirrors
    roc_sim_data._depth_sort_key's 'full'-goes-last convention, though the
    source coverage files themselves already list depths loosest-first;
    this is for when column order needs to be re-derived, e.g. from a
    dict's keys)."""
    if col == "avg_depth_full":
        return float("inf")
    try:
        return float(col[len("avg_depth_"):-1])  # strip 'avg_depth_' and trailing 'x'
    except ValueError:
        return float("inf")


def fmt_depth_col(col):
    """'avg_depth_20x' -> '20x', 'avg_depth_full' -> 'full'."""
    if col == "avg_depth_full":
        return "full"
    return col[len("avg_depth_"):]


def read_real_coverage(path=REAL_COVERAGE_PATH):
    """
    Returns (depth_cols, mean_by_depth):
      depth_cols     : the header's depth column names, in file order.
      mean_by_depth  : {depth_col: [mean avg_depth across Runs at pos]},
                       one list per depth_col, length GENOME_LENGTH
                       (index 0 = pos 1). NaN at any position no Run
                       contributed a value for.
    """
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        depth_cols = header[2:]
        n_depths = len(depth_cols)
        sums   = [array("d", [0.0]) * GENOME_LENGTH for _ in range(n_depths)]
        counts = [array("L", [0]) * GENOME_LENGTH for _ in range(n_depths)]

        for line in f:
            fields = line.rstrip("\n").split("\t")
            idx = int(fields[1]) - 1
            if not (0 <= idx < GENOME_LENGTH):
                continue
            for i in range(n_depths):
                sums[i][idx] += float(fields[2 + i])
                counts[i][idx] += 1

    mean_by_depth = {
        depth_cols[i]: [
            (sums[i][p] / counts[i][p]) if counts[i][p] else float("nan")
            for p in range(GENOME_LENGTH)
        ]
        for i in range(n_depths)
    }
    return depth_cols, mean_by_depth


def read_sim_coverage(path=SIM_COVERAGE_PATH):
    """
    Returns (depth_cols, error_rates, mean_by_error_depth):
      depth_cols          : the header's depth column names, in file order.
      error_rates          : error_rate values, discovered from the data
                              (sorted), not hardcoded.
      mean_by_error_depth  : {error_rate: {depth_col: [mean avg_depth
                              across sample_ids at pos]}}, same shape/NaN
                              convention as read_real_coverage's
                              mean_by_depth, one level deeper.
    """
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        depth_cols = header[3:]
        n_depths = len(depth_cols)
        sums   = {}  # error_rate -> [array('d', ...) per depth]
        counts = {}  # error_rate -> [array('L', ...) per depth]

        for line in f:
            fields = line.rstrip("\n").split("\t")
            error_rate = fields[1]
            idx = int(fields[2]) - 1
            if not (0 <= idx < GENOME_LENGTH):
                continue
            if error_rate not in sums:
                sums[error_rate]   = [array("d", [0.0]) * GENOME_LENGTH for _ in range(n_depths)]
                counts[error_rate] = [array("L", [0]) * GENOME_LENGTH for _ in range(n_depths)]
            for i in range(n_depths):
                sums[error_rate][i][idx] += float(fields[3 + i])
                counts[error_rate][i][idx] += 1

    error_rates = sorted(sums.keys())
    mean_by_error_depth = {
        error_rate: {
            depth_cols[i]: [
                (sums[error_rate][i][p] / counts[error_rate][i][p])
                if counts[error_rate][i][p] else float("nan")
                for p in range(GENOME_LENGTH)
            ]
            for i in range(n_depths)
        }
        for error_rate in error_rates
    }
    return depth_cols, error_rates, mean_by_error_depth


def read_real_coverage_by_filter(pos_filter, depth_cols, path=REAL_COVERAGE_PATH):
    """
    Per-Run, per-FILTER-category mean avg_depth at each of `depth_cols`, for
    plot_supp_figures.py's "coverage by quality flag" violin row (the bottom
    row of coverage_all.pdf) -- unlike read_real_coverage's cross-Run,
    per-position means, this keeps one value per (Run, filter_value,
    depth_col): the mean avg_depth, *for that Run alone*, across every
    position carrying that FILTER value. `pos_filter` is
    roc_real_data.load_filter_metadata's {pos: FILTER value} dict; positions
    absent from it are implicitly "pass", same convention as the rest of the
    real pipeline (roc_real_data.py, compute_filter_lengths). `depth_cols`
    restricts accumulation to just the requested columns (e.g. the subset of
    ALL_DEPTH_COLS overlapping real.SUBSAMPLE_ORDER's depths) -- unlike
    read_real_coverage, which always reads every depth column the file has.

    Every row contributes to exactly one (filter_value, depth_col) count per
    Run (one row = one position), so a single per-(Run, filter_value) count
    covers all of that Run's requested depth columns -- unlike
    read_real_coverage's per-position accumulators, this never needs
    GENOME_LENGTH-sized arrays, since the output is bounded by
    Runs x filter values x depth_cols (tiny), not by position.

    Returns {Run: {filter_value: {depth_col: mean_avg_depth}}}.
    """
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        file_depth_cols = header[2:]
        col_idx = [file_depth_cols.index(c) for c in depth_cols]
        n = len(depth_cols)

        sums   = {}  # Run -> filter_value -> [sum per depth_cols entry]
        counts = {}  # Run -> filter_value -> position count contributing to that sum

        for line in f:
            fields = line.rstrip("\n").split("\t")
            run = fields[0]
            pos = int(fields[1])
            filt = pos_filter.get(pos, "pass")

            run_sums = sums.setdefault(run, {})
            s = run_sums.setdefault(filt, [0.0] * n)
            for i, ci in enumerate(col_idx):
                s[i] += float(fields[2 + ci])

            run_counts = counts.setdefault(run, {})
            run_counts[filt] = run_counts.get(filt, 0) + 1

    mean_by_run_filter = {}
    for run, run_sums in sums.items():
        run_counts = counts[run]
        mean_by_run_filter[run] = {
            filt: {depth_cols[i]: s[i] / run_counts[filt] for i in range(n)}
            for filt, s in run_sums.items()
        }
    return mean_by_run_filter


def sample_to_run(sample):
    """'SRR26069013_sub20_rep3' -> 'SRR26069013' (strips the trailing
    '_sub<REP_SUBSAMPLE>_rep<N>')."""
    marker = f"_{REP_SUBSAMPLE}_rep"
    idx = sample.rfind(marker)
    return sample[:idx] if idx != -1 else sample


def read_reps_coverage(path=REPS_COVERAGE_PATH):
    """
    Returns {sample: avg_depth_array}, one entry per (Run, replicate)
    sample (e.g. "SRR26069013_sub20_rep3") -- NOT averaged across
    replicates or Runs, since these are meant to be plotted individually.
    avg_depth_array is an array('d') of length GENOME_LENGTH (index 0 =
    pos 1); use sample_to_run() to recover the Run for CenterName lookup.
    """
    reps_by_sample = {}
    with open(path) as f:
        next(f)  # header
        for line in f:
            sample, pos, avg_depth = line.rstrip("\n").split("\t")
            idx = int(pos) - 1
            if not (0 <= idx < GENOME_LENGTH):
                continue
            if sample not in reps_by_sample:
                reps_by_sample[sample] = array("d", [0.0]) * GENOME_LENGTH
            reps_by_sample[sample][idx] = float(avg_depth)
    return reps_by_sample
