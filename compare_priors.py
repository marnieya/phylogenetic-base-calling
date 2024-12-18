#!/home/marniella/miniconda3/envs/py310/bin/python

import os
os.environ["OMP_NUM_THREADS"] = '4'
os.environ["OPENBLAS_NUM_THREADS"] = '4'
os.environ["MKL_NUM_THREADS"] = '4'

import numpy as np
import pandas as pd
import re
import sys
import math

home_dir = "/space/s1/marniella/phylogenetic-base-calling/sra_raw_read_results/"
sample_names = [str.split(f.name,"_")[0] for f in os.scandir(home_dir) if (f.is_dir() and ("_files" in f.name))]

# look up the call set encoded by each letter
iupac_lookup = {
    "A":("A"), "C":("C"), "G":("G"), "T":("T"),
    "R":("A","G"), "Y":("C","T"), "S":("G","C"),
    "W":("A","T"), "K":("G","T"), "M":("A","C"),
    "B":("C", "G", "T"), "D":("A", "G", "T"),
    "H":("A", "C", "T"), "V":("A", "C", "G")
}

# distinguish between nucleotides and ambiguous encoding
nuc = ("A","C","G","T")
amb_nuc = ("R", "Y", "S", "W", "K", "M", "B", "D", "H", "V")

def compare_nucleotides(base1, base2):
    # base1 method call, base2 ground truth
    # return accurate/called/positions
    if (base2=="N") | (base2=="-") | (base1=="-"):
        return (0,0,0)
    elif (base1=="N"):
        return (0,0,1)
    call_set = iupac_lookup[base1]
    gt_set = iupac_lookup[base2]
    overlap = len(set(call_set) & set(gt_set))
    call_rate = 1/len(call_set)
    return((overlap*call_rate, 1, 1))

msa_paths = [(f.name + "/" + str.split(f.name,"_")[0] + "_pr_compare_msa.fasta") for f in os.scandir(home_dir) if (f.is_dir() and ("_files" in f.name))]

sample_prior_accuracy = {}
sample_ref_accuracy = {}
for sample,msa_path in zip(sample_names, msa_paths):
    # if sample != "SRR26069449": # check one of the least accurate tronko priors
    #     continue
    ref_total_accurate = 0
    ref_total_positions = 0 # fraction accurate for Broad pipeline
    pr_total_accurate = 0
    pr_total_positions = 0 # fraction accurate for Tronko 

    # compare pr, gt, msa_ref, and ref
    sample_seqs = ['']*5
    line_num = 0
    with open(home_dir + msa_path, 'r') as f:
        for line in f: 
            if line_num%2 == 1:
                sample_seqs[math.floor(line_num/2)] = line.strip('\n')
            line_num = line_num+1

    for i,(pr_b,gt_b,ref_b) in enumerate(zip(list(sample_seqs[0]), list(sample_seqs[1]), list(sample_seqs[2]))):
        # for each sample, we want to know how many bases the Wuhan-Hu genome matches
        # and how many bases the tronko prior matches
        pr_scores = compare_nucleotides(pr_b, gt_b)
        if pr_scores[1] > 0:
            pr_total_positions = pr_total_positions + pr_scores[2]
            pr_total_accurate = pr_total_accurate + pr_scores[0]
        ref_scores = compare_nucleotides(ref_b, gt_b)
        if ref_scores[1] > 0:
            ref_total_positions = ref_total_positions + ref_scores[2]
            ref_total_accurate = ref_total_accurate + ref_scores[0]

    sample_prior_accuracy[sample] = [pr_total_accurate, pr_total_positions]
    sample_ref_accuracy[sample] = [ref_total_accurate, ref_total_positions]

accuracies_df = pd.DataFrame.from_dict(sample_ref_accuracy,orient="index")
accuracies_df.columns = ["accurate","total"]
accuracies_df.to_csv(home_dir + "ref_prior_accuracies.csv")

# ref_genome_paths = [(f.name + "/" + str.split(f.name,"_")[0] + "_assemblies_msa_refonly.fasta") for f in os.scandir(home_dir) if (f.is_dir() and ("_files" in f.name))]
# ground_truth_paths = [(f.name + "/" + str.split(f.name,"_")[0] + "_sa.fasta") for f in os.scandir(home_dir) if (f.is_dir() and ("_files" in f.name))]
# prior_paths = [(f.name + "/tronko_PP/" + str.split(f.name,"_")[0] + ".PP.txt") for f in os.scandir(home_dir) if (f.is_dir() and ("_files" in f.name))]
# prior_genome_dict = {}
# header = ['A','C','G','T']
# for sample,prior_path,gt_path,ref_path in zip(sample_names, prior_paths, ground_truth_paths,ref_genome_paths):
#     ref_total_accurate = 0
#     ref_total_positions = 0 # fraction accurate for Broad pipeline
#     pr_total_accurate = 0
#     pr_total_positions = 0 # fraction accurate for Tronko 
#     with open(home_dir + prior_path, 'r') as f:
#         next(f)
#         for line in f: 
#             pos,pp_A,pp_C,pp_G,pp_T = str.split(line.strip('\n'))
#             prior_genome_dict[int(pos)] = header[np.nanargmax(np.array([pp_A,pp_C,pp_G,pp_T]))]
#     with open(home_dir + re.sub("_sa","_pr",gt_path), 'a') as f:
#             f.write(''.join([b for i,b in prior_genome_dict.items]) + '\n')

