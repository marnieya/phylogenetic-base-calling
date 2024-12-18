#!/home/marniella/miniconda3/envs/py310/bin/python

import os
os.environ["OMP_NUM_THREADS"] = '4' 
os.environ["OPENBLAS_NUM_THREADS"] = '4'
os.environ["MKL_NUM_THREADS"] = '4' 

import numpy as np
import pandas as pd
import re
import sys

calls_folder_path = sys.argv[1] # "/space/s1/marniella/phylogenetic-base-calling/sra_raw_read_results/"
call_mode = sys.argv[2] # "normal"

calls_file_path_list = [[(f.name + "/" + str.split(f.name,"_")[0] + "_prob_calls_pr_scaling" + str(pr_scaling) + ".txt") for pr_scaling in [0.1,0.01,0.001,0.0001] if (f.is_dir() and ("_files" in f.name))] for f in os.scandir(calls_folder_path)]
calls_file_path_list = [calls_folder_path + c for c in [c for cl in calls_file_path_list for c in cl]]

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

with open(calls_file_path_list[0]) as f:
    first_line = f.readline()
call_idxs_dict = {header_name:i for i,header_name in enumerate(str.split(first_line.strip("\n"),","))}

# exclude positions where even at full coverage we would not get it right
# from manual review, these are usually indels 
indel_regions = {}
for calls_file_path in calls_file_path_list:
    # only using prior scaling factor 0.001, which performs indististinguishably from 0.001
    if "0.001" not in calls_file_path:
        continue
    sample_indel_regions = set()
    with open(calls_file_path) as f:
        next(f)
        for line in f: 
            line_list = str.split(line.strip(), ",")
            depth = str.split(line_list[call_idxs_dict['sample']],"_")[1]
            sample = str.split(line_list[call_idxs_dict['sample']],"_")[0]
            if depth != "full":
                continue
            else:
                sa_call = compare_nucleotides(line_list[call_idxs_dict['SA']],
                                            line_list[call_idxs_dict['GT']])
                pp_call = compare_nucleotides(line_list[call_idxs_dict['PP']],
                                            line_list[call_idxs_dict['GT']])
                if (sa_call[0] == 1) & (pp_call[0] == 0):
                    sample_indel_regions = {int(line_list[call_idxs_dict['pos']])} | sample_indel_regions
        if len(sample_indel_regions) > 0:
            indel_regions[sample] = sample_indel_regions

# now manual review, outside of indel positions, where do we get it wrong?
for calls_file_path in calls_file_path_list:
    # only using prior scaling factor 0.001, which performs indististinguishably from 0.001
    if "0.001" not in calls_file_path:
        continue
    with open(calls_file_path) as f:
        next(f)
        for line in f: 
            line_list = str.split(line.strip(), ",")
            depth = str.split(line_list[call_idxs_dict['sample']],"_")[1]
            sample = str.split(line_list[call_idxs_dict['sample']],"_")[0]
            pos = int(line_list[call_idxs_dict['pos']])
            if ((sample in indel_regions) and (pos in indel_regions[sample])) or (depth == "full"):
                continue
            else:
                sa_call = compare_nucleotides(line_list[call_idxs_dict['SA']],
                                            line_list[call_idxs_dict['GT']])
                pp_call = compare_nucleotides(line_list[call_idxs_dict['PP']],
                                            line_list[call_idxs_dict['GT']])
                if (sa_call[0] == 1) & (pp_call[0] == 0):
                    print(','.join(line_list))


problem_sites = pd.read_table("/space/s1/marniella/phylogenetic-base-calling/SARS_CoV_2_problem_sites.txt",header=0,sep="\t")
caution_sites = problem_sites[problem_sites.FILTER == 'caution']['POS'].tolist()
mask_sites = problem_sites[problem_sites.FILTER == 'mask']['POS'].tolist() 

center_samples = pd.read_table("/space/s1/marniella/phylogenetic-base-calling/sra_raw_read_results/SRA_Illumina_Runs_230814-230919_random_sample.csv",header=0,sep=",")
outlier_center_samples = center_samples[center_samples.CenterName == 'LAOPHWGS']['Run'].tolist() 

methods = ['SA','PR','LL','PP']

# dictionaries of aggregated scores by sample or by position
# i.e. it will always be either sample-subsample (over all pos), or pos-subsample (over all samples)
samp_sub_scores_dict = {}
pos_sub_scores_dict = {}

# tally up accurate/called/positions
method_totals = {} 
for calls_file_path in calls_file_path_list:
    with open(calls_file_path) as f:
        next(f)
        for line in f: 
            if depth == "full": # look at subsamples only, if full depth, skip
                break
            else:
                # parse line and get ground truth
                line_list = str.split(line.strip("\n"),",")
                gt = line_list[call_idxs_dict['GT']]
                sa,pr,ll,pp = (line_list[call_idxs_dict[m]] for m in methods)
                cov_met = ( int(line_list[call_idxs_dict['mincov3_met']]) == 1 )

                sample_subsample = line_list[call_idxs_dict['sample']]
                if sample_subsample.startswith("EPI"):
                    [sample,subsample] = [re.sub("_cov.*","",s.strip("_")) for s in sample_subsample.split("err")]
                else:
                    [sample,subsample] = sample_subsample.split("_")
                sample_center = sample in outlier_center_samples
                subsample = re.sub("_sub","",subsample)
                sample_subsample = (sample + "_" + subsample) + "_" + str(sample_center) + "_" + str(cov_met)
                
                pos = line_list[call_idxs_dict['pos']]
                if int(pos) in indel_regions[sample]: # if indel, skip
                    break
                pos_mask = "" 
                if int(pos) in mask_sites:
                    pos_mask = "mask" 
                elif int(pos) in caution_sites:
                    pos_mask = "caution"
                else:
                    pos_mask = "all_other_sites"
                pos_subsample = pos + "_" + subsample + "_" + pos_mask + "_" + str(cov_met)
                
                new_scores = [compare_nucleotides(gt,method_call) for method_call in [sa,pr,ll,pp]]

                if call_mode == "roc":
                    # we add roc keys for only the thresholds this position's calculation passed
                    roc_key = [re.sub("pp_thresh_pass_","",k) for k,i in call_idxs_dict.items() if k.startswith("pp_thresh_pass_") and int(line_list[i])==1]
                    for rk in roc_key:
                        # sample/subsample/position_mask/mincov3_met/posterior_threshold
                        sample_subsample_rk = sample_subsample + "_" + rk
                        # position/subsample/sample_center/mincov3_met/posterior_threshold
                        pos_subsample_rk = pos_subsample + "_" + rk
                        if sample_subsample_rk not in samp_sub_scores_dict:
                            # initialize the sample dictionaries containing aggregated scores
                            # aggregate by sample/subsample/position_mask/mincov3_met (add over positions)
                            samp_sub_scores_dict[sample_subsample_rk] = new_scores
                            # aggregate by position/subsample/sample_center/mincov3_met (add over samples)  
                            pos_sub_scores_dict[pos_subsample_rk] = new_scores
                        else:
                            # # already initialized
                            # samp_sub_scores_dict = scores_dict[sample][0]
                            # pos_sub_scores_dict = scores_dict[sample][1]
                            # summing over positions
                            if sample_subsample_rk in samp_sub_scores_dict:
                                samp_sub_scores_dict[sample_subsample_rk] = np.add(samp_sub_scores_dict[sample_subsample_rk], new_scores)
                            else:
                                samp_sub_scores_dict[sample_subsample_rk] = new_scores
                            # or summing over samples
                            if pos_subsample_rk in pos_sub_scores_dict:
                                pos_sub_scores_dict[pos_subsample_rk] = np.add(pos_sub_scores_dict[pos_subsample_rk], new_scores)
                            else:
                                pos_sub_scores_dict[pos_subsample_rk] = new_scores
                else:
                    if sample_subsample not in samp_sub_scores_dict:
                        # initialize the sample dictionaries containing aggregated scores
                        # aggregate by sample/subsample/position_mask/mincov3_met (add over positions)
                        samp_sub_scores_dict[sample_subsample] = new_scores
                        # aggregate by position/subsample/sample_center/mincov3_met (add over samples)  
                        pos_sub_scores_dict[pos_subsample] = new_scores
                    else:
                        # # already initialized
                        # samp_sub_scores_dict = scores_dict[sample][0]
                        # pos_sub_scores_dict = scores_dict[sample][1]
                        # summing over positions
                        if sample_subsample in samp_sub_scores_dict:
                            samp_sub_scores_dict[sample_subsample] = np.add(samp_sub_scores_dict[sample_subsample], new_scores)
                        else:
                            samp_sub_scores_dict[sample_subsample] = new_scores
                        # or summing over samples
                        if pos_subsample in pos_sub_scores_dict:
                            pos_sub_scores_dict[pos_subsample] = np.add(pos_sub_scores_dict[pos_subsample], new_scores)
                        else:
                            pos_sub_scores_dict[pos_subsample] = new_scores

with open(calls_folder_path + "all_calls_summarized_by_sample.txt",'a') as f:
    for samp_sub,scores_list in samp_sub_scores_dict.items():
        for method_scores,method_name in zip(scores_list,methods):
            f.write(','.join([samp_sub, method_name] + [str(s) for s in method_scores]) + '\n')
with open(calls_folder_path + "all_calls_summarized_by_pos.txt",'a') as f:
    for pos_subsample,scores_list in pos_sub_scores_dict.items():
        for method_scores,method_name in zip(scores_list,methods):
            f.write(','.join([pos_subsample, method_name] + [str(s) for s in method_scores]) + '\n')


