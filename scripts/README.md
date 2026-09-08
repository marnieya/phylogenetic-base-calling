# README


## Pipeline

```
step0_run_all_steps.sh		      # outlines how the shell scripts fit together
├─ step1_sra_to_fastq.sh              # download/prepare real SRA reads
├─ step1_gisaid_to_fastq.sh           # simulate reads from GISAID genomes
├─ step2_fastq_to_ubam.sh             # fastq -> UBAM (can be modified/uncommented to handle real, simulated, and sub20 replicates)
├─ get_subsample_n_reads.py           # compute read counts for each subsampling depth
├─ step3_call_bases_from_pileup.sh
│   └─ call_bases.py                  # per-position base calling from pileups
├─ python plot_main_figures.py        # main ROC-like curves/coverage-analysis figures
├─ python plot_supp_figures.py        # supplemental figures
└─ ./aggregate_coverage.sh {real|reps|sim} <dir>   # has to be run three times, once per data source
```

### Plotting module dependencies

```
plot_main_figures.py  ─┬─> roc_real_data.py ─┬─> compare_nucleotides.py
plot_supp_figures.py  ─┤                     └─> plot_common.py
                        ├─> roc_sim_data.py  ───> compare_nucleotides.py, plot_common.py
                        ├─> coverage_data.py ───> plot_common.py
                        └─> plot_main_figures.py   (plot_supp_figures.py only)
```
