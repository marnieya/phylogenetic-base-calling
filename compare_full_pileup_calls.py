#!/home/marniella/miniconda3/envs/py310/bin/python

import os
os.environ["OMP_NUM_THREADS"] = '4' 
os.environ["OPENBLAS_NUM_THREADS"] = '4'
os.environ["MKL_NUM_THREADS"] = '4' 

import numpy as np
import pandas as pd
import re
import sys
import itertools

calls_folder_path = sys.argv[1] # "/space/s1/marniella/phylogenetic-base-calling/sra_raw_read_results/"
call_mode = sys.argv[2] # "normal"

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

def flatten(xss):
    return [x for xs in xss for x in xs]

gt_complete = list((pd.read_table('/space/s1/marniella/phylogenetic-base-calling/sra_complete.txt', header = None))[0])

if (call_mode == "roc") | (call_mode == "roc_overgenes"):
    calls_file_path_list = [[(f.name + "/" + str.split(f.name,"_files")[0] + "_prob_calls_pr_scaling" + str(pr_scaling) + "_ROC.txt") for pr_scaling in ["0.0"] if 
                             (f.is_dir() and ("_files" in f.name) and ( (str.split(f.name,"_")[0] in gt_complete) | ( (f.name).startswith("EPI") ) ) )] for f in os.scandir(calls_folder_path)]
else:
    calls_file_path_list = [[(f.name + "/" + str.split(f.name,"_files")[0] + "_prob_calls_pr_scaling" + str(pr_scaling) + ".txt") for pr_scaling in [0.0001, 0.001, 0.01, 0.1, 0.0] 
                             if (f.is_dir() and ("_files" in f.name) and ( (str.split(f.name,"_")[0] in gt_complete) | ( (f.name).startswith("EPI") ) ) )] for f in os.scandir(calls_folder_path)]

calls_file_path_list = [calls_folder_path + c for c in flatten(calls_file_path_list)]
if "gisaid" in calls_folder_path:
    calls_file_path_list = [re.sub("_prob_calls_pr_scaling", "_err_beta_prob_calls_pr_scaling", c) for c in calls_file_path_list]
    calls_file_path_list = flatten([[re.sub("beta", "beta" + str(e),c) for c in calls_file_path_list] for e in [0.001, 0.004, 0.007]])

with open(calls_file_path_list[0]) as f:
    first_line = f.readline()
call_idxs_dict = {header_name:i for i,header_name in enumerate(str.split(first_line.strip("\n"),","))}

# files to be read
problem_sites = pd.read_table("/space/s1/marniella/phylogenetic-base-calling/SARS_CoV_2_problem_sites.txt",header=0,sep="\t")
caution_sites = problem_sites[problem_sites.FILTER == 'caution']['POS'].tolist()
mask_sites = problem_sites[problem_sites.FILTER == 'mask']['POS'].tolist() 

methods = ['SA','PR','LL','PP']

# dictionaries of aggregated scores by sample or by position
# i.e. it will always be either sample-subsample (over all pos), or pos-subsample (over all samples)
samp_sub_scores_dict = {}
pos_sub_scores_dict = {}

# tally up accurate/called/positions
method_totals = {} 
sample_started = False

for calls_file_path in calls_file_path_list:
    with open(calls_file_path) as f:
        next(f)
        pr_scaling = float(re.sub("_ROC","",re.sub(".txt","",str.split(calls_file_path,"pr_scaling")[1])))
        for line in f: 
            # process sample name
            line_list = str.split(line.strip(), ",")
            if "gisaid" in calls_folder_path:
                sample = '_'.join(str.split(line_list[call_idxs_dict['sample']],"_")[0:3]) 
                depth = str.split(line_list[call_idxs_dict['sample']],"_")[5] 
                sim_err = str.split(line_list[call_idxs_dict['sample']],"_")[4] 
            else:
                sample = str.split(line_list[call_idxs_dict['sample']],"_")[0]
                depth = str.split(line_list[call_idxs_dict['sample']],"_")[1]
            
            pos = int(line_list[call_idxs_dict['pos']])
            gt = line_list[call_idxs_dict['GT']]
            sa,pr,ll,pp = (line_list[call_idxs_dict[m]] for m in methods)
            
            pileup_reads = int(line_list[call_idxs_dict['n_reads']]) # if summarizing without thresholding, this is binary cov=3 Y/N
            # if summarizing in "ROC" fashion, this is number of reads

            # for SRA: build a map of full coverage skipped positions 
            if sample.startswith("SRR"):
                # the flow of this block relies on full calls being the first alphabetically (i.e. SRR* vs SRR*_sub)
                if depth == "full":
                    if sample_started == False: 
                        skipped_regions = [0]*30000
                        sample_started = True
                    else: # continue adding to it
                        sa_call = compare_nucleotides(line_list[call_idxs_dict['SA']],line_list[call_idxs_dict['GT']])
                        pp_call = compare_nucleotides(line_list[call_idxs_dict['PP']],line_list[call_idxs_dict['GT']])
                        if (sa_call[0] == 1) & (pp_call[0] == 0):
                            skipped_regions[pos-1] = 1
                else:
                    sample_started = False # done building 
                    if skipped_regions[pos-1] == 1: # skip positions that full coverage cannot get right
                        continue 
            
            # parse calls this line, get ground truth
            sample_subsample = line_list[call_idxs_dict['sample']]
            if sample_subsample.startswith("EPI"):
                [sample,subsample] = [re.sub("_cov.*","",s.strip("_")) for s in sample_subsample.split("err")]
            else:
                [sample,subsample] = sample_subsample.split("_")
            subsample = re.sub("_sub","",subsample)
            sample_subsample = (sample + "_" + subsample) + "_" + str(pileup_reads) + "_" + str(pr_scaling)

            pos_mask = "" 
            if int(pos) in mask_sites:
                pos_mask = "mask" 
            elif int(pos) in caution_sites:
                pos_mask = "caution"
            else:
                pos_mask = "all_other_sites"
            pos_subsample = str(pos) + ";" + subsample + ";" + pos_mask + ";" + str(pileup_reads) + ";" + str(pr_scaling)
            
            new_scores = [compare_nucleotides(method_call,gt) for method_call in [sa,pr,ll,pp]]

            if call_mode == "roc":
                # we add roc keys for only the thresholds this position's calculation passed
                pp_thresholds = [re.sub("pp_thresh_pass_","",k) for k,i in call_idxs_dict.items() if k.startswith("pp_thresh_pass_")]
                pp_thresholds_passed = [int(line_list[i]) for k,i in call_idxs_dict.items() if k.startswith("pp_thresh_pass_")]
                for rk,agg_bool in zip(pp_thresholds, pp_thresholds_passed):
                    # structure of the keys to the dict with the score aggregation: sample/subsample/pr_scaling/posterior_threshold
                    # this is the roc key but we only want to aggregate if this particular site passed
                    sample_subsample_rk = sample_subsample + "_" + rk
                    if sample_subsample_rk not in samp_sub_scores_dict:
                        if agg_bool:
                            samp_sub_scores_dict[sample_subsample_rk] = new_scores
                        else: 
                            samp_sub_scores_dict[sample_subsample_rk] = (0,0,1)
                    else: 
                        # if this particular threshold was not passed, we still want to keep a tally of the sites visited
                        # so we will add it as (0,0,1) to indicate a position was visited
                        if agg_bool:
                            samp_sub_scores_dict[sample_subsample_rk] = np.add(samp_sub_scores_dict[sample_subsample_rk], new_scores)
                        else:
                            samp_sub_scores_dict[sample_subsample_rk] = np.add(samp_sub_scores_dict[sample_subsample_rk], (0,0,1))
            elif call_mode == "roc_overgenes":
                pp_thresholds = [re.sub("pp_thresh_pass_","",k) for k,i in call_idxs_dict.items() if k.startswith("pp_thresh_pass_")]
                pp_thresholds_passed = [int(line_list[i]) for k,i in call_idxs_dict.items() if k.startswith("pp_thresh_pass_")]
                for rk,agg_bool in zip(pp_thresholds, pp_thresholds_passed):
                    if pos > 13468 and pos <= 16236:
                        rk = rk + "_rdrp"
                    elif pos > 21563 and pos <= 25384:
                        rk = rk + "_spike"
                    else:
                        continue
                    sample_subsample_rk = sample_subsample + "_" + rk
                    if sample_subsample_rk not in samp_sub_scores_dict:
                        if agg_bool:
                            samp_sub_scores_dict[sample_subsample_rk] = new_scores
                        else: 
                            samp_sub_scores_dict[sample_subsample_rk] = (0,0,1)
                    else: 
                        if agg_bool:
                            samp_sub_scores_dict[sample_subsample_rk] = np.add(samp_sub_scores_dict[sample_subsample_rk], new_scores)
                        else:
                            samp_sub_scores_dict[sample_subsample_rk] = np.add(samp_sub_scores_dict[sample_subsample_rk], (0,0,1))
            else:
                if sample_subsample in samp_sub_scores_dict:
                    samp_sub_scores_dict[sample_subsample] = np.add(samp_sub_scores_dict[sample_subsample], new_scores)
                else:
                    samp_sub_scores_dict[sample_subsample] = new_scores

                if pos_subsample in pos_sub_scores_dict:
                    pos_sub_scores_dict[pos_subsample] = np.add(pos_sub_scores_dict[pos_subsample], new_scores)
                else:
                    pos_sub_scores_dict[pos_subsample] = new_scores

if call_mode == "roc":
    samp_sub_scores_df = pd.DataFrame.from_dict({ k:itertools.chain(*itertools.chain(*[[([tr[0], tr[1]]) for tr in v if not np.isscalar(tr)]])) for k,v in samp_sub_scores_dict.items() }, orient = "index")
    samp_sub_scores_df.columns = itertools.chain(*[['accurate_total_' + m, 'called_total_' + m] for m in methods])
    if sample.startswith("SRR"):
        samp_sub_scores_df.to_csv(calls_folder_path + "sra_gt_complete_calls_pr_scaling0.0_summarized_by_sample_totals_ROC.csv", index_label='id' )
    else:
        samp_sub_scores_df.to_csv(calls_folder_path + "gisaid_sim_calls_pr_scaling0.0_summarized_by_sample_totals_ROC.csv", index_label='id' )

elif call_mode == "roc_overgenes":
    samp_sub_scores_df = pd.DataFrame.from_dict({ k:itertools.chain(*itertools.chain(*[[([tr[0], tr[1]]) for tr in v if not np.isscalar(tr)]])) for k,v in samp_sub_scores_dict.items() }, orient = "index")
    samp_sub_scores_df.columns = itertools.chain(*[['accurate_total_' + m, 'called_total_' + m] for m in methods])
    if sample.startswith("SRR"):
        samp_sub_scores_df.to_csv(calls_folder_path + "sra_gt_complete_calls_summarized_by_sample_totals_ROC_bygene.csv", index_label='id' )
    else:
        samp_sub_scores_df.to_csv(calls_folder_path + "gisaid_sim_calls_summarized_by_sample_totals_ROC_bygene.csv", index_label='id' )

elif call_mode == "normal":
    samp_sub_scores_df = pd.DataFrame.from_dict({ k:itertools.chain(*itertools.chain(*[[([tr[0], tr[1]]) for tr in v if not np.isscalar(tr)]])) for k,v in samp_sub_scores_dict.items() }, orient = "index")
    samp_sub_scores_df.columns = itertools.chain(*[['accurate_total_' + m, 'called_total_' + m] for m in methods])
    if sample.startswith("SRR"):
        samp_sub_scores_df.to_csv(calls_folder_path + "sra_gt_complete_calls_summarized_by_sample_totals.csv", index_label='id' )
    else:
        samp_sub_scores_df.to_csv(calls_folder_path + "gisaid_sim_calls_summarized_by_sample_totals.csv", index_label='id' )

    pos_sub_scores_df = pd.DataFrame.from_dict({ k:itertools.chain(*itertools.chain(*[[([tr[0], tr[1]]) for tr in v if not np.isscalar(tr)]])) for k,v in pos_sub_scores_dict.items() }, orient = "index")
    pos_sub_scores_df.columns = itertools.chain(*[['accurate_total_' + m, 'called_total_' + m] for m in methods])
    if sample.startswith("SRR"):
        pos_sub_scores_df.to_csv(calls_folder_path + "sra_gt_complete_calls_summarized_by_pos_totals.csv", index_label='id' )
    else:
        pos_sub_scores_df.to_csv(calls_folder_path + "gisaid_sim_complete_calls_summarized_by_pos_totals.csv", index_label='id' )
    

# with open(calls_folder_path + "all_prob_calls_summarized_by_sample" + call_mode + ".txt",'a') as f:
#     for samp_sub,scores_list in samp_sub_scores_dict.items():
#         for method_scores,method_name in zip(scores_list,methods):
#             f.write(','.join([samp_sub, method_name] + [str(s) for s in method_scores]) + '\n')
            
# if call_mode == "normal":
#     with open(calls_folder_path + "all_prob_calls_summarized_by_pos.txt",'a') as f:
#         for pos_subsample,scores_list in pos_sub_scores_dict.items():
#             for method_scores,method_name in zip(scores_list,methods):
#                 f.write(','.join([pos_subsample, method_name] + [str(s) for s in method_scores]) + '\n')

