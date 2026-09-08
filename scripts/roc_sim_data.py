"""
Simulated-data pipeline: call-file reading, accumulation, and aggregation for
the phyloBC ROC analysis. Plot-agnostic — this module knows about the
simulated call-file format and sample-name convention, but nothing about how
any particular figure draws them. Used by plot_main_figures.py, and also by
plot_supp_figures.py (which draws other simulated-data plots from this same
pipeline).

Extracted from plot_roc_simulated.py — see that script for the original
derivation notes carried forward here unchanged.
"""

import csv
import matplotlib.cm as cm
from compare_nucleotides import compare_nucleotides
from plot_common import (
    METHODS,
    GENOME_LENGTH,
    THRESHOLDS,
    THRESHOLDS_LOG10,
    N_THRESHOLDS,
    _new_method_acc,
    _new_threshold_acc,
    _add_score,
    _safe_mean,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CONDITIONS = {
    "pr_scaling=0.0":   "all_sim_calls_pr_scaling0.0_precise.txt",
    "pr_scaling=0.001": "all_sim_calls_pr_scaling0.001_precise.txt",
    "pr_scaling=MLE":   "all_sim_calls_pr_scalingMLE_precise.txt",
}


# ---------------------------------------------------------------------------
# Sample name parsing
# ---------------------------------------------------------------------------

def parse_sample(sample):
    """
    Parse simulated sample name into (sample_id, error_rate, mean_depth).

    Example: "EPI_ISL_18106665_err_beta0.001.simulated_5x_cov"
      parts  :  [0]  [1]   [2]   [3]         [4]           [5] [6]
      → sample_id  = "EPI_ISL_18106665"
      → error_rate = "beta0.001"   (strip ".simulated" from parts[4])
      → mean_depth = "sub5"        (strip "x" from parts[-2], prepend "sub")
    """
    parts      = sample.split("_")
    sample_id  = "_".join(parts[0:3])
    error_rate = parts[4].replace(".simulated", "")
    mean_depth = "sub" + parts[-2].rstrip("x")
    return sample_id, error_rate, mean_depth


# ---------------------------------------------------------------------------
# Main read loop
# ---------------------------------------------------------------------------

def read_one_condition(condition, filepath):
    """
    Read one simulated call file.

    Returns
    -------
    stats[run][error_rate][mean_depth][method]    → [acc_sum, cr_sum, n]
    thr_stats[run][error_rate][mean_depth][t_idx] → [acc_sum, cr_sum, n_passing]
    """
    s  = {}
    ts = {}

    print(f"Reading {filepath} ...")

    with open(filepath, newline="") as f:
        reader = csv.reader(f)
        next(reader)

        skipped = 0
        for row in reader:
            if len(row) != 11 or any(field == "" for field in row):  # TODO: remove once confirmed clean
                skipped += 1
                continue
            try:  # TODO: remove once confirmed clean
                sample = row[0]
                gt     = row[3]
                sa     = row[4]
                pr     = row[5]
                ll     = row[6]
                pp     = row[7]
                pp_val = float(row[10])
            except ValueError:
                skipped += 1
                continue

            run, error_rate, mean_depth = parse_sample(sample)

            if run not in s:
                s[run]  = {}
                ts[run] = {}
            if error_rate not in s[run]:
                s[run][error_rate]  = {}
                ts[run][error_rate] = {}
            if mean_depth not in s[run][error_rate]:
                s[run][error_rate][mean_depth]  = {m: _new_method_acc() for m in METHODS}
                ts[run][error_rate][mean_depth] = _new_threshold_acc()

            method_scores = {
                "SA": compare_nucleotides(sa, gt),
                "PR": compare_nucleotides(pr, gt),
                "LL": compare_nucleotides(ll, gt),
                "PP": compare_nucleotides(pp, gt),
            }

            for method in METHODS:
                acc, cr = method_scores[method]
                _add_score(s[run][error_rate][mean_depth][method], acc, cr)

            pp_acc, pp_cr = method_scores["PP"]
            for t_idx in range(N_THRESHOLDS):
                if pp_val > THRESHOLDS_LOG10[t_idx]:
                    _add_score(ts[run][error_rate][mean_depth][t_idx], pp_acc, pp_cr)

    if skipped:
        print(f"  WARNING: skipped {skipped} malformed rows in {filepath}")

    return s, ts


# ---------------------------------------------------------------------------
# Per-run access (no cross-run averaging)
# ---------------------------------------------------------------------------

def per_run_stats(stats_cond, run, error_rate, mean_depth, denom=GENOME_LENGTH):
    """
    Compute per-method (acc, cr) for exactly one run/sample_id -- the
    single-run analogue of roc_real_data._agg_methods (which averages over a
    run_subset instead). Added for plot_supp_figures.py's *_by_run figures,
    which overlay every run's own curve rather than the cross-run mean.
    Unlike the real pipeline, no _aggregate_for_run_subset-style helper is
    needed here: roc_sim_data's raw stats/threshold_stats are already keyed
    directly by run, so this only needs the same safe-division/n>0 gate
    aggregate_stats() applies before averaging.

    Returns {method: (acc, cr)} -- NaN for every method if this run has no
    data at all for (error_rate, mean_depth), NaN for just one method if
    that method never got any calls (n == 0) despite the combination
    existing.
    """
    run_dict = stats_cond.get(run, {}).get(error_rate, {}).get(mean_depth)
    result = {}
    for method in METHODS:
        if run_dict is None:
            result[method] = (float("nan"), float("nan"))
            continue
        acc_sum, cr_sum, n = run_dict[method]
        if n > 0:
            result[method] = (acc_sum / cr_sum if cr_sum > 0 else float("nan"), cr_sum / denom)
        else:
            result[method] = (float("nan"), float("nan"))
    return result


def per_run_thresholds(threshold_stats_cond, run, error_rate, mean_depth, denom=GENOME_LENGTH):
    """Threshold analogue of per_run_stats -- the single-run version of what
    aggregate_stats() averages over every run. Returns {t_idx: (acc, cr)}."""
    run_list = threshold_stats_cond.get(run, {}).get(error_rate, {}).get(mean_depth)
    result = {}
    for t_idx in range(N_THRESHOLDS):
        if run_list is None:
            result[t_idx] = (float("nan"), float("nan"))
            continue
        acc_sum, cr_sum, n_passing = run_list[t_idx]
        cr = n_passing / denom
        acc = acc_sum / n_passing if n_passing > 0 else float("nan")
        result[t_idx] = (acc, cr)
    return result


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate_stats(stats, threshold_stats):
    """
    Average per-run (sample_id) stats for each (condition, error_rate, mean_depth).

    Parameters
    ----------
    stats[condition][run][error_rate][mean_depth][method]    → [acc_sum, cr_sum, n]
    threshold_stats[condition][run][error_rate][mean_depth][t_idx] → [acc_sum, cr_sum, n_passing]

    Returns
    -------
    mean_stats[condition][error_rate][mean_depth][method]    → (mean_acc, mean_cr)
    mean_thr_stats[condition][error_rate][mean_depth][t_idx] → (mean_acc, mean_cr)
    """
    mean_s  = {}
    mean_ts = {}

    for condition, run_dict in stats.items():
        mean_s[condition]  = {}
        mean_ts[condition] = {}
        thr_dict = threshold_stats[condition]

        # Discover all (error_rate, mean_depth) pairs present in the data
        pairs = {
            (error_rate, mean_depth)
            for run, err_d in run_dict.items()
            for error_rate, dep_d in err_d.items()
            for mean_depth in dep_d
        }

        for error_rate, mean_depth in pairs:
            mean_s[condition].setdefault(error_rate, {})[mean_depth] = {}
            mean_ts[condition].setdefault(error_rate, {})[mean_depth] = {}

            runs = [
                r for r in run_dict
                if error_rate in run_dict[r] and mean_depth in run_dict[r][error_rate]
            ]

            for method in METHODS:
                run_accs, run_crs = [], []
                for run in runs:
                    acc_sum, cr_sum, n = run_dict[run][error_rate][mean_depth][method]
                    if n > 0:
                        run_accs.append(acc_sum / cr_sum if cr_sum > 0 else float("nan"))
                        run_crs.append(cr_sum / GENOME_LENGTH)
                mean_s[condition][error_rate][mean_depth][method] = (
                    _safe_mean(run_accs), _safe_mean(run_crs)
                )

            for t_idx in range(N_THRESHOLDS):
                run_accs, run_crs = [], []
                for run in runs:
                    if (run not in thr_dict
                            or error_rate not in thr_dict[run]
                            or mean_depth not in thr_dict[run][error_rate]):
                        continue
                    acc_sum, cr_sum, n_passing = thr_dict[run][error_rate][mean_depth][t_idx]
                    run_crs.append(n_passing / GENOME_LENGTH)
                    if n_passing > 0:
                        run_accs.append(acc_sum / n_passing)
                mean_ts[condition][error_rate][mean_depth][t_idx] = (
                    _safe_mean(run_accs), _safe_mean(run_crs)
                )

    return mean_s, mean_ts


# ---------------------------------------------------------------------------
# Crossover analysis
# ---------------------------------------------------------------------------

def find_pp_sa_crossover(mean_stats, mean_thr_stats):
    """
    For each (condition, error_rate, mean_depth), find the least stringent threshold
    (highest index, lowest value) at which PP accuracy > SA accuracy.

    Returns
    -------
    crossover[condition][error_rate][mean_depth] = (threshold_value, threshold_idx)
        or (None, None) if PP never exceeds SA.
    """
    crossover = {}
    for condition, err_dict in mean_stats.items():
        crossover[condition] = {}
        for error_rate, dep_dict in err_dict.items():
            crossover[condition][error_rate] = {}
            for mean_depth, method_dict in dep_dict.items():
                sa_acc, _ = method_dict.get("SA", (float("nan"), float("nan")))
                result = (None, None)
                for t_idx in range(N_THRESHOLDS - 1, -1, -1):
                    pp_acc, _ = mean_thr_stats[condition][error_rate][mean_depth].get(
                        t_idx, (float("nan"), float("nan"))
                    )
                    if pp_acc > sa_acc:
                        result = (THRESHOLDS[t_idx], t_idx)
                        break
                crossover[condition][error_rate][mean_depth] = result
    return crossover


def write_crossover_csv(crossover, filename="sim_pp_sa_crossover.csv"):
    rows = []
    for condition, err_dict in crossover.items():
        for error_rate, dep_dict in err_dict.items():
            for mean_depth, (thr_val, t_idx) in dep_dict.items():
                rows.append({
                    "condition":     condition,
                    "error_rate":    error_rate,
                    "mean_depth":    mean_depth,
                    "threshold":     thr_val,
                    "threshold_idx": t_idx,
                })
    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["condition", "error_rate", "mean_depth",
                           "threshold", "threshold_idx"]
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Written {filename}")


# ---------------------------------------------------------------------------
# Formatting / color helpers
# ---------------------------------------------------------------------------

def _depth_sort_key(d):
    if d == "full":
        return float("inf")
    try:
        return float(d[3:])  # strip leading "sub"
    except ValueError:
        return float("inf")


def _fmt_depth(d):
    """Format mean_depth label for panel titles: 'sub5' → 'mean depth = 5x'."""
    if d == "full":
        return "full"
    return f"mean depth = {d[3:]}x"


def _fmt_error_rate(e):
    """Format error_rate label for row/block titles: 'beta0.001' → 'error
    rate = 0.001'. Added for sim_sa_pp_all (plot_supp_figures.py), which
    discovers error_rate values from the data rather than hardcoding the two
    values sim_sa_pp_1v5 knows about (_ERROR_RATE_LABELS)."""
    if e.startswith("beta"):
        return f"error rate = {e[len('beta'):]}"
    return e


def _cond_colors():
    """Return {condition: RGB color} using tab10, indexed by CONDITIONS order."""
    cmap = cm.tab10
    return {c: cmap.colors[i % len(cmap.colors)] for i, c in enumerate(CONDITIONS)}
