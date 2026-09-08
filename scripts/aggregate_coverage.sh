#!/bin/bash
# Aggregates coverage_plot.txt files into combined tables. See CLAUDE.md in
# this directory for output formats and the rationale behind each choice
# here. Run from (or point at) the directory holding the per-Run/sample
# *_files/ subdirectories: ./aggregate_coverage.sh [real|sim|reps|all] [dir]
set -euo pipefail

mode="${1:-all}"
reads_dir="${2:-.}"
reads_dir="${reads_dir%/}/"

out_real="${reads_dir}all_coverage_plots_real.txt"
out_sim="${reads_dir}all_coverage_plots_sim.txt"
out_reps="${reads_dir}all_coverage_plots_reps.txt"

REAL_DEPTHS=(20 15 10 5 4 3 2 1 0.8 0.5 0.2 full)
SIM_DEPTHS=(5 4 3 2 1 0.8 0.5 0.2)
SIM_ERROR_RATES=(beta0.001 beta0.004 beta0.007)
REP_DEPTH="20"
REP_IDS=($(seq 1 10))

_real_depth_file() {
    local dir="$1" sra="$2" d="$3"
    if [[ "${d}" == "full" ]]; then
        printf '%s/%s.coverage_plot.txt' "${dir}" "${sra}"
    else
        printf '%s/%s_sub%s.coverage_plot.txt' "${dir}" "${sra}" "${d}"
    fi
}

_real_depth_label() {
    local d="$1"
    if [[ "${d}" == "full" ]]; then
        printf 'avg_depth_full'
    else
        printf 'avg_depth_%sx' "${d}"
    fi
}

_paste_depth_row() {
    local id="$1"; shift
    local files=("$@")
    local tmpdir; tmpdir="$(mktemp -d)"

    cut -f2,3 "${files[0]}" > "${tmpdir}/0"
    local tmp_files=("${tmpdir}/0")
    local i=1
    for f in "${files[@]:1}"; do
        cut -f3 "${f}" > "${tmpdir}/${i}"
        tmp_files+=("${tmpdir}/${i}")
        i=$((i + 1))
    done

    paste -d $'\t' "${tmp_files[@]}" | sed "s/^/${id}\t/"
    rm -rf "${tmpdir}"
}

aggregate_real() {
    {
        printf 'Run\tpos'
        for d in "${REAL_DEPTHS[@]}"; do printf '\t%s' "$(_real_depth_label "${d}")"; done
        printf '\n'
    } > "${out_real}"

    local n_skipped=0
    for dir in "${reads_dir}"SRR*_files; do
        [[ -d "${dir}" ]] || continue
        local sra="${dir%_files}"; sra="${sra##*/}"

        local files=() missing=0
        for d in "${REAL_DEPTHS[@]}"; do
            local f; f="$(_real_depth_file "${dir}" "${sra}" "${d}")"
            if [[ ! -f "${f}" ]]; then
                echo "aggregate_real: missing ${f}, skipping ${sra}" >&2
                missing=1
                break
            fi
            files+=("${f}")
        done
        if (( missing )); then
            n_skipped=$((n_skipped + 1))
            continue
        fi

        _paste_depth_row "${sra}" "${files[@]}" >> "${out_real}"
    done
    echo "Written ${out_real} (${n_skipped} Run(s) skipped for missing files)"
}

aggregate_sim() {
    {
        printf 'sample_id\terror_rate\tpos'
        for d in "${SIM_DEPTHS[@]}"; do printf '\tavg_depth_%sx' "${d}"; done
        printf '\n'
    } > "${out_sim}"

    local n_skipped=0
    for dir in "${reads_dir}"EPI_ISL*_files; do
        [[ -d "${dir}" ]] || continue
        local sample_id="${dir%_files}"; sample_id="${sample_id##*/}"

        for err in "${SIM_ERROR_RATES[@]}"; do
            local id="${sample_id}_err_${err}"
            local files=() missing=0
            for d in "${SIM_DEPTHS[@]}"; do
                local f="${dir}/${id}.simulated_${d}x_cov.coverage_plot.txt"
                if [[ ! -f "${f}" ]]; then
                    echo "aggregate_sim: missing ${f}, skipping ${id}" >&2
                    missing=1
                    break
                fi
                files+=("${f}")
            done
            if (( missing )); then
                n_skipped=$((n_skipped + 1))
                continue
            fi

            _paste_depth_row "${sample_id}	${err}" "${files[@]}" >> "${out_sim}"
        done
    done
    echo "Written ${out_sim} (${n_skipped} (sample_id, error_rate) combination(s) skipped for missing files)"
}

aggregate_reps() {
    printf 'sample\tpos\tavg_depth\n' > "${out_reps}"

    local n_skipped=0
    for dir in "${reads_dir}"SRR*_files; do
        [[ -d "${dir}" ]] || continue
        local sra="${dir%_files}"; sra="${sra##*/}"

        for i in "${REP_IDS[@]}"; do
            local f="${dir}/${sra}_sub${REP_DEPTH}_rep${i}.coverage_plot.txt"
            if [[ ! -f "${f}" ]]; then
                echo "aggregate_reps: missing ${f}, skipping ${sra} rep${i}" >&2
                n_skipped=$((n_skipped + 1))
                continue
            fi
            cut -f2,3 "${f}" | sed "s/^/${sra}_sub${REP_DEPTH}_rep${i}\t/" >> "${out_reps}"
        done
    done
    echo "Written ${out_reps} (${n_skipped} replicate file(s) skipped for missing files)"
}

case "${mode}" in
    real) aggregate_real ;;
    sim)  aggregate_sim ;;
    reps) aggregate_reps ;;
    all)  aggregate_real; aggregate_sim; aggregate_reps ;;
    *)    echo "Usage: $0 [real|sim|reps|all] [dir]" >&2; exit 1 ;;
esac
