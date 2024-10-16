#!/home/marniella/miniconda3/envs/py310/bin/python

import os
os.environ["OMP_NUM_THREADS"] = '4'
os.environ["OPENBLAS_NUM_THREADS"] = '4'
os.environ["MKL_NUM_THREADS"] = '4'

import numpy as np
import pandas as pd
import re
import sys

msat_path = sys.argv[1] #"/space/s1/marniella/phylogenetic-base-calling/sra_raw_read_results_min_cov_0/SRR26070223_files/SRR26070223_calls_msat.txt"

problem_sites = pd.read_table("/space/s1/marniella/phylogenetic-base-calling/SARS_CoV_2_problem_sites.txt",header=0,sep="\t")
caution_sites = problem_sites[problem_sites.FILTER == 'caution']['POS'].tolist()
mask_sites = problem_sites[problem_sites.FILTER == 'mask']['POS'].tolist() 

msat_df = pd.read_csv(msat_path, skiprows = 1, header = None)
msat_df = msat_df[0].str.split('', expand=True)
msat_df = msat_df.drop(labels = [0,msat_df.shape[1]-1], axis = 1)

msat_header = pd.read_csv(msat_path, nrows = 1, header = None)
msat_df.columns = msat_header.values[0]

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
    if (overlap > 0):
        return (1, call_rate, 1)
    else:
        return (0, call_rate, 1)

call_idxs_dict = {}
ref_idx = 0
sample = ""
for idx_first_line,sample_name in enumerate(msat_header.values[0]):
    # parse sample name
    sample_name_list = sample_name.split("_")
    sample = sample_name_list[0]
    if sample.startswith("EPI"):
        ref_idx = idx_first_line
        continue
    else:
        if sample_name == sample:
            depth = "full"
            calltype = "GT"
        else:
            if len(sample_name_list) == 2:
                calltype = "sa"
            else:
                calltype = '_'.join(sample_name_list[2:])
            depth = sample_name_list[1]
        # shorten call type
        if 'like' in calltype:
            calltype = "ll"
        elif 'prior' in calltype:
            calltype = 'pr'
        elif 'post' in calltype:
            thresh = float(re.sub('_postprob_calls','',(re.sub('pr_thresh','',calltype))))
            calltype = 'pp'
        # construct dictionaries for where each call is
        if depth not in call_idxs_dict:
            if calltype == 'pp':
                call_idxs_dict[depth] = {calltype:{thresh:idx_first_line}}
            else:
                call_idxs_dict[depth] = {calltype:idx_first_line}
        else:
            if (calltype == 'pp') and (calltype in call_idxs_dict[depth]):
                call_idxs_dict[depth]['pp'][thresh] = idx_first_line
            elif (calltype == 'pp') and ('pp' not in call_idxs_dict[depth]):
                call_idxs_dict[depth]['pp'] = {thresh:idx_first_line}
            elif (calltype != 'pp'):
                call_idxs_dict[depth][calltype] = idx_first_line
if 'full' in call_idxs_dict:
    call_idxs_dict['full']['sa'] = call_idxs_dict['full']['GT']

# calls_dict = {}
full_pp = msat_df[msat_header.values[0][call_idxs_dict['full']['pp'][0.001]]]
full_sa = msat_df[msat_header.values[0][call_idxs_dict['full']['sa']]]
wrong_positions = [pos for pos,s in enumerate([compare_nucleotides(p,s) for p,s in zip(full_pp,full_sa)]) if ((s[0]==0) & (s[2]==1))]
poly_a_pos = [((full == 'A') and (posterior == '-')) for posterior,full in zip(full_pp, full_sa) ] and msat_df.iloc[wrong_positions][[sample,msat_header.values[0][call_idxs_dict['full']['pp'][0.001]]]].duplicated(keep=False).index.tolist()
poly_a_df = msat_df.iloc[poly_a_pos]
msat_df = msat_df.drop(poly_a_pos)
poly_a_df.to_csv(re.sub('_calls_msat.txt','_calls_polyA_rem.csv',msat_path))

gt = msat_df[msat_header.values[0][call_idxs_dict['full']['GT']]]
gt_len = np.sum([1 if (b!='-') else 0 for b in gt])
gt_ll = msat_df[msat_header.values[0][call_idxs_dict['full']['ll']]] 
gt_pr = msat_df[msat_header.values[0][call_idxs_dict['full']['pr']]] 

ref_bias_sites = [i for i,call in enumerate(gt_ll) if call == "N"] # sites we cannot evaluate fairly, they are missing in ground truth pileup
indel_sites = [i for i,call in enumerate(gt_pr) if call == "-"] # sites we will not evaluate, as they are indels
problem_sites_mask = [i for n,i in enumerate([i for i,b in enumerate(msat_df[msat_header.values[0][ref_idx]]) if b!="-"]) if (n in mask_sites)]
caution_sites_mask = [i for n,i in enumerate([i for i,b in enumerate(msat_df[msat_header.values[0][ref_idx]]) if b!="-"]) if (n in caution_sites)]
keep_out_sites = set(ref_bias_sites + indel_sites)
keep_sites_mask = [i for i in range(len(gt)) if (i not in keep_out_sites)]

with open(re.sub('_calls_msat.txt','_calls_summarized.txt',msat_path),'a') as f:
    for depth in call_idxs_dict:
        if 'pp' in call_idxs_dict[depth]:
            # calls_dict[depth] = {}
            sa,pr,ll = (msat_df[msat_header.values[0][call_idxs_dict[depth][m]]] for m in ['sa','pr','ll'])
            pr_scores = np.array([compare_nucleotides(method_call,ground_truth) for method_call,ground_truth in zip(pr,gt)])
            ll_scores = np.array([compare_nucleotides(method_call,ground_truth) for method_call,ground_truth in zip(ll,gt)])
            sa_scores = np.array([compare_nucleotides(method_call,ground_truth) for method_call,ground_truth in zip(sa,gt)])
            
            pp1,pp2,pp3,pp4 = (msat_df[msat_header.values[0][call_idxs_dict[depth]['pp'][thresh]]] for thresh in [0.1,0.01,0.001,0.0001])
            pp_scores1 = np.array([compare_nucleotides(method_call,ground_truth) for method_call,ground_truth in zip(pp1,gt)])
            pp_scores2 = np.array([compare_nucleotides(method_call,ground_truth) for method_call,ground_truth in zip(pp2,gt)])
            pp_scores3 = np.array([compare_nucleotides(method_call,ground_truth) for method_call,ground_truth in zip(pp3,gt)])
            pp_scores4 = np.array([compare_nucleotides(method_call,ground_truth) for method_call,ground_truth in zip(pp4,gt)])

            for scores,score_type in zip(
                [pr_scores, ll_scores, sa_scores, pp_scores1, pp_scores2, pp_scores3, pp_scores4],
                ['pr','ll','sa','pp_th0.1', 'pp_th0.01', 'pp_th0.001', 'pp_th0.0001']):
                masked_scores = scores[keep_sites_mask]
                caution_sites_scores = scores[caution_sites_mask]
                problem_sites_scores = scores[problem_sites_mask]
                total_accurate, total_called, total_positions = (sum([s[0] for s in masked_scores]), sum([s[1] for s in masked_scores]), sum([s[2] for s in masked_scores]))
                caution_accurate, caution_called, caution_positions = (sum([s[0] for s in caution_sites_scores]), sum([s[1] for s in caution_sites_scores]), sum([s[2] for s in caution_sites_scores]))
                problem_accurate, problem_called, problem_positions = (sum([s[0] for s in problem_sites_scores]), sum([s[1] for s in problem_sites_scores]), sum([s[2] for s in problem_sites_scores]))
                # calls_dict[depth][score_type] = [(total_accurate, total_called, total_positions),
                #                                 (caution_accurate, caution_called, caution_positions),
                #                                 (problem_accurate, problem_called, problem_positions)]
                f.write(','.join([sample, depth, score_type, 'total', ','.join([str(n) for n in [total_accurate, total_called, total_positions, gt_len]]), '\n']))
                f.write(','.join([sample, depth, score_type, 'caution', ','.join([str(n) for n in [caution_accurate, caution_called, caution_positions, gt_len]]), '\n']))
                f.write(','.join([sample, depth, score_type, 'problem', ','.join([str(n) for n in [problem_accurate, problem_called, problem_positions, gt_len]]), '\n']))

