"""
Real (subsampled sequencing) data pipeline: call-file reading, accumulation,
and aggregation for the phyloBC ROC analysis. Plot-agnostic — this module
knows about the call-file format, sample-name convention, and metadata CSVs
for the *real* data, but nothing about how any particular figure draws them.
Used by plot_main_figures.py, and also by plot_supp_figures.py (which draws
other real-data plots from this same pipeline).

Extracted from plot_roc.py — see that script and its project CLAUDE.md for
the original derivation notes (sample-name parsing, PP_val log10 format,
accuracy/call-rate definitions, crossover algorithm history, etc.) carried
forward here unchanged.
"""

import csv
import os
import statistics
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

import matplotlib.cm as cm

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CONDITIONS = {
    "pr_scaling=0.0":   "all_calls_pr_scaling0.0_precise.txt",
    "pr_scaling=0.001": "all_calls_pr_scaling0.001_precise.txt",
    "pr_scaling=MLE":   "all_calls_pr_scalingMLE_precise.txt",
}

# Subset of conditions used in plots 2 and 3 (pr_scaling=0.0 excluded)
CONDITIONS_23 = {k: v for k, v in CONDITIONS.items() if k != "pr_scaling=0.0"}

# Condition subset for plot1c (MLE only)
CONDITIONS_1C = {"pr_scaling=MLE": CONDITIONS["pr_scaling=MLE"]}

# Colors assigned by position in the full CONDITIONS dict, so a given
# condition's color stays consistent across every plot/script that draws
# from this module — not just within one figure.
ALL_CONDITION_COLORS = {
    cond: col for cond, col in zip(CONDITIONS.keys(), cm.tab10.colors[:len(CONDITIONS)])
}

SUBSAMPLE_ORDER = ["sub0.8", "sub1", "sub3", "sub5", "sub10", "sub20", "full"]
# Deliberately 6 depths, not the full 11 aggregate_coverage.sh's REAL_DEPTHS
# sweeps (sub0.2/sub0.5/sub2/sub4/sub15 also excluded) -- a scope decision
# made back when the real-data call files themselves were compiled, not a
# plotting-side filter: the 11-depth series got unwieldy and repetitive
# (several depths show very similar patterns) and didn't line up well
# against the simulated pipeline's own depth set, so base-calling was only
# ever run for these 6 real depths. The call files have no rows at all for
# the other 5 -- confirmed by briefly widening this list to all 11, which
# produced entirely empty crossover/aggregate entries for exactly those 5,
# for every condition. The coverage pipeline (aggregate_coverage.sh /
# coverage_data.py) intentionally keeps its own, wider REAL_DEPTHS
# independent of this list -- coverage was measured at all 11 depths even
# though base-calling wasn't, and showing that extra spread in the
# supplemental coverage plots is useful on its own, with no need to match
# the ROC-analysis depth set exactly.

FILTER_VALUES = ["PASS", "mask", "caution"]

STRUCTURAL_GENES = {"Envelope", "Membrane", "Nucleocapsid", "Spike"}
GENE_REGION_TYPES = ["structural protein genes", "nonstructural genes"]

SAMPLE_METADATA_PATH = "sample_sequencing_metadata.csv"  # columns: Run, CenterName, primer_scheme
FILTER_METADATA_PATH = "position_masks_metadata.csv"     # columns: pos, FILTER
GENE_METADATA_PATH   = "position_regions_metadata.csv"   # columns: pos, gene_region, gene_region_length

# Replicate data (sub20 only). Same 11-column format as the main CONDITIONS files
# and living in their own subdirectory, but prefixed "all_prob_calls_" rather than
# "all_calls_". A pr_scaling=0.001 replicate file also exists there, but is
# intentionally not read/used for now. Replicates are embedded in the `sample`
# field (e.g. "SRR30176165_sub20_rep3"). See read_replicate_data().
REPLICATE_DIR = "sub20_only"
REPLICATE_CONDITIONS = {
    k: v.replace("all_calls_", "all_prob_calls_")
    for k, v in CONDITIONS.items() if k in ("pr_scaling=0.0", "pr_scaling=MLE")
}
REPLICATE_SUBSAMPLE = "sub20"


# ---------------------------------------------------------------------------
# Metadata loading
# ---------------------------------------------------------------------------

def load_sample_metadata(path):
    """
    Returns {Run: {"CenterName": str, "primer_scheme": str, ...}}
    """
    metadata = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            metadata[row["Run"]] = dict(row)
    return metadata


def load_filter_metadata(path):
    """
    Returns {pos (int): filter_value (str)}.
    Positions absent from this dict are implicitly "PASS".
    """
    pos_filter = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pos_filter[int(row["pos"])] = row["FILTER"]
    return pos_filter


def load_gene_metadata(path):
    """
    Returns
    -------
    pos_gene : {pos (int): gene_region (str)}
    gene_lengths : {gene_region (str): gene_region_length (int)}
        Constant denominators for call_rate in gene-region-stratified plots.
    """
    pos_gene = {}
    gene_lengths = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pos_gene[int(row["pos"])] = row["gene_region"]
            gene_lengths[row["gene_region"]] = int(row["gene_region_length"])
    return pos_gene, gene_lengths


def compute_filter_lengths(pos_filter):
    """
    Derive constant position counts per FILTER value from the metadata dict.
    Used as denominators for call_rate in filter-stratified plots.

    Returns {filter_value: n_positions}
    """
    n_mask    = sum(1 for v in pos_filter.values() if v == "mask")
    n_caution = sum(1 for v in pos_filter.values() if v == "caution")
    return {
        "PASS":    GENOME_LENGTH - n_mask - n_caution,
        "mask":    n_mask,
        "caution": n_caution,
    }


# ---------------------------------------------------------------------------
# Main read loop
# ---------------------------------------------------------------------------

def read_one_condition(condition, filepath, pos_filter, pos_gene):
    """
    Read one call file, accumulating summary statistics for that condition.

    For condition == "pr_scaling=MLE", also populates nreads_stats and
    nreads_thr_stats (both None for all other conditions).

    Returns
    -------
    s, ts, fs, fts, gs, gts  — per-condition accumulators (keyed by Run):
        s   : [Run][subsample][method]       -> [acc_sum, cr_sum, n]
        ts  : [Run][subsample][t_idx]        -> [acc_sum, cr_sum, n]
        fs  : [Run][subsample][filt][method] -> [acc_sum, cr_sum, n]
        fts : [Run][subsample][filt][t_idx]  -> [acc_sum, cr_sum, n]
        gs  : [Run][subsample][gene][method] -> [acc_sum, cr_sum, n]
        gts : [Run][subsample][gene][t_idx]  -> [acc_sum, cr_sum, n]
    nreads_stats          : [subsample][n_reads_cat][method][Run] -> [acc_sum, cr_sum, n]
                            (None unless condition == "pr_scaling=MLE")
    nreads_thr_stats      : [subsample][n_reads_cat][t_idx][Run]  -> [acc_sum, cr_sum, n]
                            (None unless condition == "pr_scaling=MLE")
    gene_nreads_stats     : [subsample][gene_region][n_reads_cat][method][Run] -> [acc_sum, cr_sum, n]
                            (None unless condition == "pr_scaling=MLE")
    gene_nreads_thr_stats : [subsample][gene_region][n_reads_cat][t_idx][Run]  -> [acc_sum, cr_sum, n]
                            (None unless condition == "pr_scaling=MLE")
    """
    s   = {}; ts  = {}
    fs  = {}; fts = {}
    gs  = {}; gts = {}

    is_mle = (condition == "pr_scaling=MLE")
    nreads_stats          = {} if is_mle else None
    nreads_thr_stats      = {} if is_mle else None
    gene_nreads_stats     = {} if is_mle else None
    gene_nreads_thr_stats = {} if is_mle else None

    print(f"Reading {filepath} ...")

    with open(filepath, newline="") as f:
        reader = csv.reader(f)
        next(reader)  # skip header

        skipped = 0
        for row in reader:
            if len(row) != 11 or any(f == "" for f in row):  # TODO: remove once file corruption is resolved
                skipped += 1
                continue
            try:  # TODO: remove once file corruption is resolved
                sample   = row[0]
                pos      = int(float(row[1]))
                n_reads  = int(float(row[2]))
                gt       = row[3]
                sa       = row[4]
                pr       = row[5]
                ll       = row[6]
                pp       = row[7]
                # pr_val = float(row[8])
                # ll_val = float(row[9])
                pp_val   = float(row[10])
            except ValueError:
                skipped += 1
                continue

            parts = sample.split("_", 1)
            if len(parts) == 2:
                run, subsample = parts
            else:
                run, subsample = parts[0], "full"

            filt_val = pos_filter.get(pos, "PASS")
            gene_val = pos_gene.get(pos)

            if run not in s:
                s[run]   = {}; ts[run]  = {}
                fs[run]  = {}; fts[run] = {}
                gs[run]  = {}; gts[run] = {}
            if subsample not in s[run]:
                s[run][subsample]   = {m: _new_method_acc() for m in METHODS}
                ts[run][subsample]  = _new_threshold_acc()
                fs[run][subsample]  = {}
                fts[run][subsample] = {}
                gs[run][subsample]  = {}
                gts[run][subsample] = {}

            method_scores = {
                "SA": compare_nucleotides(sa, gt),
                "PR": compare_nucleotides(pr, gt),
                "LL": compare_nucleotides(ll, gt),
                "PP": compare_nucleotides(pp, gt),
            }

            for method in METHODS:
                acc, cr = method_scores[method]
                _add_score(s[run][subsample][method], acc, cr)

            pp_acc, pp_cr = method_scores["PP"]
            for t_idx in range(N_THRESHOLDS):
                if pp_val > THRESHOLDS_LOG10[t_idx]:
                    # thresholds descending: all subsequent also satisfied
                    _add_score(ts[run][subsample][t_idx], pp_acc, pp_cr)

            if filt_val not in fs[run][subsample]:
                fs[run][subsample][filt_val]  = {m: _new_method_acc() for m in METHODS}
                fts[run][subsample][filt_val] = _new_threshold_acc()
            for method in METHODS:
                acc, cr = method_scores[method]
                _add_score(fs[run][subsample][filt_val][method], acc, cr)
            for t_idx in range(N_THRESHOLDS):
                if pp_val > THRESHOLDS_LOG10[t_idx]:
                    _add_score(fts[run][subsample][filt_val][t_idx], pp_acc, pp_cr)

            if gene_val is not None:
                if gene_val not in gs[run][subsample]:
                    gs[run][subsample][gene_val]  = {m: _new_method_acc() for m in METHODS}
                    gts[run][subsample][gene_val] = _new_threshold_acc()
                for method in METHODS:
                    acc, cr = method_scores[method]
                    _add_score(gs[run][subsample][gene_val][method], acc, cr)
                for t_idx in range(N_THRESHOLDS):
                    if pp_val > THRESHOLDS_LOG10[t_idx]:
                        _add_score(gts[run][subsample][gene_val][t_idx], pp_acc, pp_cr)

            # n_reads accumulators (MLE only, sub0.8–sub20 only, n_reads>0, no per-run averaging)
            if is_mle and subsample != "full" and n_reads > 0:
                nr_cat = f"n_{n_reads}"
                if subsample not in nreads_stats:
                    nreads_stats[subsample]          = {}
                    nreads_thr_stats[subsample]      = {}
                    gene_nreads_stats[subsample]     = {}
                    gene_nreads_thr_stats[subsample] = {}
                if nr_cat not in nreads_stats[subsample]:
                    nreads_stats[subsample][nr_cat]     = {m: {} for m in METHODS}
                    nreads_thr_stats[subsample][nr_cat] = [{} for _ in range(N_THRESHOLDS)]
                for method in METHODS:
                    acc, cr = method_scores[method]
                    if run not in nreads_stats[subsample][nr_cat][method]:
                        nreads_stats[subsample][nr_cat][method][run] = _new_method_acc()
                    _add_score(nreads_stats[subsample][nr_cat][method][run], acc, cr)
                for t_idx in range(N_THRESHOLDS):
                    if pp_val > THRESHOLDS_LOG10[t_idx]:
                        if run not in nreads_thr_stats[subsample][nr_cat][t_idx]:
                            nreads_thr_stats[subsample][nr_cat][t_idx][run] = _new_method_acc()
                        _add_score(nreads_thr_stats[subsample][nr_cat][t_idx][run], pp_acc, pp_cr)

                # gene × n_reads accumulators
                if gene_val is not None:
                    if gene_val not in gene_nreads_stats[subsample]:
                        gene_nreads_stats[subsample][gene_val]     = {}
                        gene_nreads_thr_stats[subsample][gene_val] = {}
                    if nr_cat not in gene_nreads_stats[subsample][gene_val]:
                        gene_nreads_stats[subsample][gene_val][nr_cat]     = {m: {} for m in METHODS}
                        gene_nreads_thr_stats[subsample][gene_val][nr_cat] = [{} for _ in range(N_THRESHOLDS)]
                    for method in METHODS:
                        acc, cr = method_scores[method]
                        if run not in gene_nreads_stats[subsample][gene_val][nr_cat][method]:
                            gene_nreads_stats[subsample][gene_val][nr_cat][method][run] = _new_method_acc()
                        _add_score(gene_nreads_stats[subsample][gene_val][nr_cat][method][run], acc, cr)
                    for t_idx in range(N_THRESHOLDS):
                        if pp_val > THRESHOLDS_LOG10[t_idx]:
                            if run not in gene_nreads_thr_stats[subsample][gene_val][nr_cat][t_idx]:
                                gene_nreads_thr_stats[subsample][gene_val][nr_cat][t_idx][run] = _new_method_acc()
                            _add_score(gene_nreads_thr_stats[subsample][gene_val][nr_cat][t_idx][run], pp_acc, pp_cr)

    if skipped:
        print(f"  WARNING: skipped {skipped} malformed rows in {filepath}")

    return s, ts, fs, fts, gs, gts, nreads_stats, nreads_thr_stats, gene_nreads_stats, gene_nreads_thr_stats


# ---------------------------------------------------------------------------
# Per-position diff read (for gene_acc_diff_kde)
# ---------------------------------------------------------------------------

def read_mle_per_position_diffs(filepath, pos_gene, mle_mode_t_idx):
    """
    Second pass over the MLE file. For each (pos, run) pair, accumulates per-subsample
    accuracy and call-rate differences (PP_thresh − SA). Only positions with a known gene
    region are kept.

    Returns
    -------
    pos_run_diffs : {pos: {"gene": str, "runs": {run: [(acc_diff, cr_diff), ...]}}}
        One (acc_diff, cr_diff) tuple per subsample observation.
    """
    threshold = THRESHOLDS_LOG10[mle_mode_t_idx]
    pos_run_diffs = {}

    print(f"Reading {filepath} for per-position diffs (t_idx={mle_mode_t_idx}) ...")

    with open(filepath, newline="") as f:
        reader = csv.reader(f)
        next(reader)

        skipped = 0
        for row in reader:
            if len(row) != 11 or any(field == "" for field in row):  # TODO: remove once file corruption is resolved
                skipped += 1
                continue
            try:  # TODO: remove once file corruption is resolved
                sample = row[0]
                pos    = int(float(row[1]))
                gt     = row[3]
                sa     = row[4]
                pp     = row[7]
                pp_val = float(row[10])
            except ValueError:
                skipped += 1
                continue

            gene_val = pos_gene.get(pos)
            if gene_val is None:
                continue

            run = sample.split("_", 1)[0]

            sa_acc, sa_cr = compare_nucleotides(sa, gt)
            pp_acc, pp_cr = compare_nucleotides(pp, gt)
            pp_thresh_acc = pp_acc if pp_val > threshold else 0.0
            pp_thresh_cr  = pp_cr  if pp_val > threshold else 0.0
            acc_diff = pp_thresh_acc - sa_acc
            cr_diff  = pp_thresh_cr  - sa_cr

            if pos not in pos_run_diffs:
                pos_run_diffs[pos] = {"gene": gene_val, "runs": {}}
            pos_run_diffs[pos]["runs"].setdefault(run, []).append((acc_diff, cr_diff))

    if skipped:
        print(f"  WARNING: skipped {skipped} malformed rows in {filepath}")

    return pos_run_diffs


# ---------------------------------------------------------------------------
# Replicate data read (sub20-only replicate call files)
# ---------------------------------------------------------------------------

def read_replicate_data(pos_gene):
    """
    Reads the sub20-only replicate call files (one per condition in
    REPLICATE_CONDITIONS, living in REPLICATE_DIR, same 11-column format and
    filename as the corresponding main condition file). Replicates are embedded
    in the `sample` field, e.g. "SRR30176165_sub20_rep3" -> run="SRR30176165",
    subsample="sub20", replicate="rep3" (split on "_": first part = run, last
    part = replicate). Rows are restricted to subsample == REPLICATE_SUBSAMPLE
    and n_reads > 0.

    Returns
    -------
    gene_rep_stats     : [condition][gene][n_reads_cat]["SA"|"PP"][(run, replicate)] -> [acc_sum, cr_sum, n]
    gene_rep_thr_stats : [condition][gene][n_reads_cat][t_idx][(run, replicate)]      -> [acc_sum, cr_sum, n]
        Conditions with no file present in REPLICATE_DIR are simply absent from
        both dicts.
    """
    gene_rep_stats     = {}
    gene_rep_thr_stats = {}

    for condition, basename in REPLICATE_CONDITIONS.items():
        filepath = os.path.join(REPLICATE_DIR, basename)
        if not os.path.exists(filepath):
            print(f"read_replicate_data: {filepath} not found, skipping {condition}.")
            continue

        gs  = gene_rep_stats[condition]     = {}
        gts = gene_rep_thr_stats[condition] = {}
        print(f"Reading {filepath} for replicate data ...")

        with open(filepath, newline="") as f:
            reader = csv.reader(f)
            next(reader)  # skip header

            skipped = 0
            for row in reader:
                if len(row) != 11 or any(field == "" for field in row):  # TODO: remove once file corruption is resolved
                    skipped += 1
                    continue
                try:  # TODO: remove once file corruption is resolved
                    sample  = row[0]
                    pos     = int(float(row[1]))
                    n_reads = int(float(row[2]))
                    gt      = row[3]
                    sa      = row[4]
                    pp      = row[7]
                    pp_val  = float(row[10])
                except ValueError:
                    skipped += 1
                    continue

                parts = sample.split("_")
                if len(parts) < 3:
                    skipped += 1
                    continue
                run, subsample, replicate = parts[0], parts[1], parts[-1]
                if subsample != REPLICATE_SUBSAMPLE or n_reads <= 0:
                    continue

                gene_val = pos_gene.get(pos)
                if gene_val is None:
                    continue

                nr_cat  = f"n_{n_reads}"
                rep_key = (run, replicate)

                sa_acc, sa_cr = compare_nucleotides(sa, gt)
                pp_acc, pp_cr = compare_nucleotides(pp, gt)

                if gene_val not in gs:
                    gs[gene_val]  = {}
                    gts[gene_val] = {}
                if nr_cat not in gs[gene_val]:
                    gs[gene_val][nr_cat]  = {"SA": {}, "PP": {}}
                    gts[gene_val][nr_cat] = [{} for _ in range(N_THRESHOLDS)]

                bucket = gs[gene_val][nr_cat]
                if rep_key not in bucket["SA"]:
                    bucket["SA"][rep_key] = _new_method_acc()
                _add_score(bucket["SA"][rep_key], sa_acc, sa_cr)
                if rep_key not in bucket["PP"]:
                    bucket["PP"][rep_key] = _new_method_acc()
                _add_score(bucket["PP"][rep_key], pp_acc, pp_cr)

                thr_list = gts[gene_val][nr_cat]
                for t_idx in range(N_THRESHOLDS):
                    if pp_val > THRESHOLDS_LOG10[t_idx]:
                        if rep_key not in thr_list[t_idx]:
                            thr_list[t_idx][rep_key] = _new_method_acc()
                        _add_score(thr_list[t_idx][rep_key], pp_acc, pp_cr)

        if skipped:
            print(f"  WARNING: skipped {skipped} malformed rows in {filepath}")

    return gene_rep_stats, gene_rep_thr_stats


# ---------------------------------------------------------------------------
# Aggregation helpers
# ---------------------------------------------------------------------------

def _agg_methods(stats_cond, run_subset, subsample, denom):
    """
    Compute per-method (mean_acc, mean_cr) averaged over run_subset for one subsample.
    call_rate = cr_sum / denom  (constant denominator).
    Returns {method: (mean_acc, mean_cr)}.
    """
    result = {}
    for method in METHODS:
        run_accs, run_crs = [], []
        for run in run_subset:
            if run not in stats_cond or subsample not in stats_cond[run]:
                continue
            acc_sum, cr_sum, n = stats_cond[run][subsample][method]
            if n > 0:
                run_accs.append(acc_sum / cr_sum if cr_sum > 0 else float("nan"))
                run_crs.append(cr_sum / denom)
        result[method] = (_safe_mean(run_accs), _safe_mean(run_crs))
    return result


def _agg_thresholds(threshold_stats_cond, run_subset, subsample, denom):
    """
    Compute per-threshold (mean_acc, mean_cr) averaged over run_subset for one subsample.
    call_rate = n_passing / denom  (constant denominator).
    Returns {t_idx: (mean_acc, mean_cr)}.
    """
    result = {}
    for t_idx in range(N_THRESHOLDS):
        run_accs, run_crs = [], []
        for run in run_subset:
            if run not in threshold_stats_cond or subsample not in threshold_stats_cond[run]:
                continue
            acc_sum, cr_sum, n_passing = threshold_stats_cond[run][subsample][t_idx]
            run_crs.append(n_passing / denom)
            if n_passing > 0:
                run_accs.append(acc_sum / n_passing)
        result[t_idx] = (_safe_mean(run_accs), _safe_mean(run_crs))
    return result


def _aggregate_for_run_subset(stats_cond, threshold_stats_cond, run_subset, denom=GENOME_LENGTH):
    """
    Aggregate base and threshold stats for a given run subset and single condition.

    Returns
    -------
    mean_s  : {subsample: {method: (mean_acc, mean_cr)}}
    mean_ts : {subsample: {t_idx:  (mean_acc, mean_cr)}}
    """
    mean_s  = {}
    mean_ts = {}
    for subsample in SUBSAMPLE_ORDER:
        mean_s[subsample]  = _agg_methods(stats_cond, run_subset, subsample, denom)
        mean_ts[subsample] = _agg_thresholds(threshold_stats_cond, run_subset, subsample, denom)
    return mean_s, mean_ts


def aggregate_stats(stats, threshold_stats):
    """
    [condition][subsample][method]  -> (mean_acc, mean_cr)
    [condition][subsample][t_idx]   -> (mean_acc, mean_cr)
    """
    mean_stats     = {}
    mean_thr_stats = {}
    for condition in stats:
        all_runs = set(stats[condition].keys())
        mean_stats[condition], mean_thr_stats[condition] = _aggregate_for_run_subset(
            stats[condition], threshold_stats[condition], all_runs
        )
    return mean_stats, mean_thr_stats


def aggregate_filter_stats(filter_stats, filter_thr_stats, filter_lengths):
    """
    [condition][subsample][filter_value][method]  -> (mean_acc, mean_cr)
    [condition][subsample][filter_value][t_idx]   -> (mean_acc, mean_cr)
    call_rate denominator = filter_lengths[filter_value].
    """
    mean_fs  = {c: {s: {} for s in SUBSAMPLE_ORDER} for c in filter_stats}
    mean_fts = {c: {s: {} for s in SUBSAMPLE_ORDER} for c in filter_stats}

    for condition in filter_stats:
        for subsample in SUBSAMPLE_ORDER:
            for filt_val in FILTER_VALUES:
                denom = filter_lengths.get(filt_val, 1)

                for method in METHODS:
                    run_accs, run_crs = [], []
                    for run, sub_dict in filter_stats[condition].items():
                        if subsample not in sub_dict or filt_val not in sub_dict[subsample]:
                            continue
                        acc_sum, cr_sum, n = sub_dict[subsample][filt_val][method]
                        if n > 0:
                            run_accs.append(acc_sum / cr_sum if cr_sum > 0 else float("nan"))
                            run_crs.append(cr_sum / denom)
                    mean_fs[condition][subsample][filt_val] = mean_fs[condition][subsample].get(filt_val, {})
                    mean_fs[condition][subsample][filt_val][method] = (
                        _safe_mean(run_accs), _safe_mean(run_crs)
                    )

                for t_idx in range(N_THRESHOLDS):
                    run_accs, run_crs = [], []
                    for run, sub_dict in filter_thr_stats[condition].items():
                        if subsample not in sub_dict or filt_val not in sub_dict[subsample]:
                            continue
                        acc_sum, cr_sum, n_passing = sub_dict[subsample][filt_val][t_idx]
                        run_crs.append(n_passing / denom)
                        if n_passing > 0:
                            run_accs.append(acc_sum / n_passing)
                    mean_fts[condition][subsample][filt_val] = mean_fts[condition][subsample].get(filt_val, {})
                    mean_fts[condition][subsample][filt_val][t_idx] = (
                        _safe_mean(run_accs), _safe_mean(run_crs)
                    )

    return mean_fs, mean_fts


def aggregate_gene_stats(gene_stats, gene_thr_stats, gene_lengths):
    """
    [condition][subsample][gene_region][method]  -> (mean_acc, mean_cr)
    [condition][subsample][gene_region][t_idx]   -> (mean_acc, mean_cr)
    call_rate denominator = gene_lengths[gene_region].
    Also returns sorted list of all gene region labels found in the data.
    """
    all_genes = sorted({
        gene
        for condition in gene_stats
        for run_dict in gene_stats[condition].values()
        for sub_dict in run_dict.values()
        for gene in sub_dict.keys()
    })

    mean_gs  = {c: {s: {} for s in SUBSAMPLE_ORDER} for c in gene_stats}
    mean_gts = {c: {s: {} for s in SUBSAMPLE_ORDER} for c in gene_stats}

    for condition in gene_stats:
        for subsample in SUBSAMPLE_ORDER:
            for gene in all_genes:
                denom = gene_lengths.get(gene, 1)

                for method in METHODS:
                    run_accs, run_crs = [], []
                    for run, sub_dict in gene_stats[condition].items():
                        if subsample not in sub_dict or gene not in sub_dict[subsample]:
                            continue
                        acc_sum, cr_sum, n = sub_dict[subsample][gene][method]
                        if n > 0:
                            run_accs.append(acc_sum / cr_sum if cr_sum > 0 else float("nan"))
                            run_crs.append(cr_sum / denom)
                    mean_gs[condition][subsample][gene] = mean_gs[condition][subsample].get(gene, {})
                    mean_gs[condition][subsample][gene][method] = (
                        _safe_mean(run_accs), _safe_mean(run_crs)
                    )

                for t_idx in range(N_THRESHOLDS):
                    run_accs, run_crs = [], []
                    for run, sub_dict in gene_thr_stats[condition].items():
                        if subsample not in sub_dict or gene not in sub_dict[subsample]:
                            continue
                        acc_sum, cr_sum, n_passing = sub_dict[subsample][gene][t_idx]
                        run_crs.append(n_passing / denom)
                        if n_passing > 0:
                            run_accs.append(acc_sum / n_passing)
                    mean_gts[condition][subsample][gene] = mean_gts[condition][subsample].get(gene, {})
                    mean_gts[condition][subsample][gene][t_idx] = (
                        _safe_mean(run_accs), _safe_mean(run_crs)
                    )

    return mean_gs, mean_gts, all_genes


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def _fmt_subsample(s):
    """Format a subsample key for panel titles: 'sub1' → 'mean depth = 1x'.
    'full' is passed through as-is (s[3:] would otherwise slice it to 'l') --
    this never came up while every caller only passed sub1/sub5, but
    sa_pp_all (plot_supp_figures.py) plots every subsample in
    SUBSAMPLE_ORDER, 'full' included."""
    if s == "full":
        return "full"
    return f"mean depth = {s[3:]}x"


# ---------------------------------------------------------------------------
# Crossover analysis
# ---------------------------------------------------------------------------

def find_pp_sa_crossover(mean_stats, mean_thr_stats):
    """
    For each (condition, subsample) in the supplied mean_stats, find the most
    stringent PP threshold whose call rate still meets or exceeds SA's call
    rate (not accuracy).

    Scans threshold indices from N_THRESHOLDS-1 (least stringent, threshold
    0.0, call rate 1.0) toward 0 (most stringent, ~1-1e-15) -- i.e. toward
    higher thresholds, since call rate falls as the threshold gets stricter.
    The scan keeps going as long as call rate is still >= SA's call rate; as
    soon as a threshold's call rate dips below SA's, the scan stops and steps
    one threshold back (looser) to the last index that was still >= SA's call
    rate, guaranteeing the result is never below SA's call rate (it can land
    exactly on it, but need not). If call rate never dips below SA's even at
    the strictest threshold (t_idx=0), that index is returned directly.

    Returns crossover[condition][subsample] = (threshold_value, threshold_idx,
    pp_call_rate, sa_call_rate), or (None, None, None, None) if no call rate
    data is available at any threshold. pp_call_rate is PP's call rate at
    threshold_idx -- the value actually compared against sa_call_rate to find
    the crossover, kept here so it doesn't have to be looked back up from
    mean_thr_stats later. sa_call_rate is SA's own call rate for that
    (condition, subsample), included alongside it for reference even though
    it doesn't vary by threshold.
    """
    crossover = {}
    for condition in mean_stats:
        crossover[condition] = {}
        for subsample in SUBSAMPLE_ORDER:
            _, sa_cr = mean_stats[condition][subsample].get(
                "SA", (float("nan"), float("nan"))
            )
            last_valid_idx = None
            last_valid_cr = None
            for t_idx in range(N_THRESHOLDS - 1, -1, -1):
                _, pp_cr = mean_thr_stats[condition][subsample].get(
                    t_idx, (float("nan"), float("nan"))
                )
                if pp_cr != pp_cr:  # NaN: no data at this threshold, skip it
                    continue
                if pp_cr < sa_cr:
                    break
                last_valid_idx = t_idx
                last_valid_cr = pp_cr
            if last_valid_idx is None:
                crossover[condition][subsample] = (None, None, None, None)
            else:
                crossover[condition][subsample] = (
                    THRESHOLDS[last_valid_idx], last_valid_idx, last_valid_cr, sa_cr
                )
    return crossover


def find_pp_sa_crossover_accuracy(mean_stats, mean_thr_stats):
    """
    For each (condition, subsample) in the supplied mean_stats, find the least
    stringent PP threshold whose accuracy still meets or exceeds SA's accuracy
    (not call rate) -- the accuracy-based counterpart to find_pp_sa_crossover's
    call-rate-based definition above (see CLAUDE.md's "Crossover analysis"
    section for how the two differ and why both exist; this one feeds
    plot_supp_figures.py's pp_sa_crossover_accuracy.csv rather than any plot).

    Scans threshold indices from 0 (most stringent, ~1-1e-15) toward
    N_THRESHOLDS-1 (least stringent, threshold 0.0) -- i.e. toward looser
    thresholds, since accuracy falls as the threshold gets looser (the
    opposite trend from call rate, hence the reversed scan direction from
    find_pp_sa_crossover). The scan keeps going as long as accuracy is still
    >= SA's accuracy; as soon as a threshold's accuracy dips below SA's, the
    scan stops and steps one threshold back (stricter) to the last index that
    was still >= SA's accuracy, guaranteeing the result is never below SA's
    accuracy (it can land exactly on it, but need not -- a tie still counts
    as "meets", so the scan keeps loosening through a plateau of equal
    accuracy values rather than stopping at the first one it meets). If
    accuracy never dips below SA's even at the loosest threshold
    (t_idx=N_THRESHOLDS-1), that index is returned directly.

    Returns crossover[condition][subsample] = (threshold_value, threshold_idx,
    pp_accuracy, sa_accuracy), or (None, None, None, None) if no accuracy
    data is available at any threshold. pp_accuracy is PP's accuracy at
    threshold_idx -- the value actually compared against sa_accuracy to find
    the crossover, kept here so it doesn't have to be looked back up from
    mean_thr_stats later. sa_accuracy is SA's own accuracy for that
    (condition, subsample), included alongside it for reference even though
    it doesn't vary by threshold.
    """
    crossover = {}
    for condition in mean_stats:
        crossover[condition] = {}
        for subsample in SUBSAMPLE_ORDER:
            sa_acc, _ = mean_stats[condition][subsample].get(
                "SA", (float("nan"), float("nan"))
            )
            last_valid_idx = None
            last_valid_acc = None
            for t_idx in range(N_THRESHOLDS):
                pp_acc, _ = mean_thr_stats[condition][subsample].get(
                    t_idx, (float("nan"), float("nan"))
                )
                if pp_acc != pp_acc:  # NaN: no data at this threshold, skip it
                    continue
                if pp_acc < sa_acc:
                    break
                last_valid_idx = t_idx
                last_valid_acc = pp_acc
            if last_valid_idx is None:
                crossover[condition][subsample] = (None, None, None, None)
            else:
                crossover[condition][subsample] = (
                    THRESHOLDS[last_valid_idx], last_valid_idx, last_valid_acc, sa_acc
                )
    return crossover


def write_crossover_csv(crossover, filename="pp_sa_crossover.csv", value_label="value"):
    """Write the full crossover dict to CSV (all conditions, all subsamples
    with data). Each row also includes the PP/SA values that were actually
    compared to find the crossover -- named "pp_<value_label>"/
    "sa_<value_label>" so the same writer serves both find_pp_sa_crossover
    (value_label="call_rate") and find_pp_sa_crossover_accuracy
    (value_label="accuracy") without the column names being generic/
    ambiguous in either file. (condition, subsample) entries with no
    crossover found (threshold_idx is None -- no data at all for that
    combination, as opposed to a threshold search that came up empty) are
    skipped rather than written as a blank row, matching
    write_crossover_csv_tidy's own behavior -- there's nothing meaningful
    to report for a combination the underlying call files never had."""
    pp_col = f"pp_{value_label}"
    sa_col = f"sa_{value_label}"
    rows = []
    for condition, sub_dict in crossover.items():
        for subsample, (threshold, t_idx, pp_value, sa_value) in sub_dict.items():
            if t_idx is None:
                continue
            rows.append({
                "condition":     condition,
                "subsample":     subsample,
                "threshold":     threshold,
                "threshold_idx": t_idx,
                pp_col:          pp_value,
                sa_col:          sa_value,
            })
    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["condition", "subsample", "threshold", "threshold_idx",
                           pp_col, sa_col]
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Written {filename}")


def _crossover_delta_str(t_idx):
    """Format threshold_idx as an exact "1-Ce-N" delta string (delta = 1 -
    threshold, C in {1,5}), reconstructed from the THRESHOLDS construction
    formula (plot_common.py: 1 - c*10**-n for n=15..1, two steps per decade,
    plus a t=0.0/idx=N_THRESHOLDS-1 sentinel) rather than subtracting the
    printed decimal threshold from 1 -- that subtraction cancels almost all
    of a float64's precision for thresholds this close to 1, since idx=0 is
    ~1-1e-15. The sentinel (t=0.0, "always passes") has no meaningful delta;
    returned as "0" rather than a bogus "1-5e-1"."""
    t_idx = int(t_idx)
    if t_idx == N_THRESHOLDS - 1:
        return "0"
    n = 15 - t_idx // 2
    c = 1 if t_idx % 2 == 0 else 5
    return f"1-{c}e-{n}"


def write_crossover_csv_tidy(crossover, filename, value_label, decimals,
                              condition="pr_scaling=MLE", exclude_subsamples=("full",)):
    """
    Condensed, presentation-ready counterpart to write_crossover_csv's full
    table -- one condition only (default pr_scaling=MLE), exclude_subsamples
    dropped (default just "full", which isn't a true subsample and is out of
    scope for this analysis), threshold expressed as an exact "1-Ce-N" delta
    string (_crossover_delta_str, reconstructed from threshold_idx rather
    than the printed decimal) instead of the raw decimal, and
    pp_<value_label>/sa_<value_label> rounded to `decimals` places.
    `condition`/`threshold_idx` columns are dropped entirely -- condition is
    constant across rows once filtered, and threshold_idx is an internal
    detail a reader of a results table doesn't need (the delta string next
    to it already conveys the threshold itself). Row order follows
    `crossover[condition]`'s own key order, which is SUBSAMPLE_ORDER's order
    (both crossover-finding functions build it by iterating SUBSAMPLE_ORDER).
    Subsamples with no crossover found (threshold_idx is None) are skipped.
    """
    pp_col = f"pp_{value_label}"
    sa_col = f"sa_{value_label}"
    rows = []
    for subsample, (_threshold, t_idx, pp_value, sa_value) in crossover.get(condition, {}).items():
        if subsample in exclude_subsamples or t_idx is None:
            continue
        rows.append({
            "subsample": subsample,
            "threshold": _crossover_delta_str(t_idx),
            pp_col:      round(pp_value, decimals),
            sa_col:      round(sa_value, decimals),
        })
    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["subsample", "threshold", pp_col, sa_col])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Written {filename}")


def _mode_crossover_threshold(condition_crossover):
    """
    Given crossover[subsample] = (threshold_value, threshold_idx, pp_value,
    sa_value), return the mode threshold_idx across subsamples. On a tie,
    returns the least stringent (highest index) among the tied values.
    Returns None if no valid entries. Indexes into each tuple by position
    (entry[1]) rather than destructuring all four fields, since this only
    ever needs threshold_idx.
    """
    indices = [entry[1] for entry in condition_crossover.values() if entry[1] is not None]
    if not indices:
        return None
    try:
        return statistics.mode(indices)
    except statistics.StatisticsError:
        # Tie: return the highest (least stringent) index among the tied modes
        counts = {}
        for i in indices:
            counts[i] = counts.get(i, 0) + 1
        max_count = max(counts.values())
        return max(k for k, v in counts.items() if v == max_count)
