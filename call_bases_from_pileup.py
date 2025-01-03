#!/home/marniella/miniconda3/envs/py310/bin/python

import os
os.environ["OMP_NUM_THREADS"] = '4' 
os.environ["OPENBLAS_NUM_THREADS"] = '4'
os.environ["MKL_NUM_THREADS"] = '4' 

import numpy as np
import re
import sys

pileup_file = sys.argv[1] 
pileup_results_path = sys.argv[2] 
call_mode = sys.argv[3]

# pileup_file = "SRR26069449_piledup.txt"
# pileup_results_path = "/space/s1/marniella/phylogenetic-base-calling/sra_raw_read_results/SRR26069449_files/"

prior_scaling = [0.0001, 0.001, 0.01, 0.1]

if "SRR" in pileup_file:
    # for SRA
    prefix = re.sub('_piledup.txt','', pileup_file)
    prefix_base = re.sub('_sub.+', '', prefix)
else:
    # for GISAID
    prefix = re.sub('_piledup.txt','', pileup_file)
    prefix_base = re.sub('_[0-9]*\.*[0-9]+x_cov.sim.+', '', prefix) # re.sub('.sim.+', '', prefix)

with open(pileup_results_path + prefix + '_sa.fasta') as f: 
    subsample_assembly = f.readlines()
subsample_assembly = subsample_assembly[1].strip('\n')

# read the ground truth in order to skip indels
with open(pileup_results_path + re.sub("_err.*$","",prefix_base) + '_gt.fasta') as f: 
    ground_truth_assembly = f.readlines()
ground_truth_assembly = ground_truth_assembly[1].strip('\n')

# read the reference in order to get the indexing right for indels
with open(pileup_results_path + prefix_base + '_assemblies_msa_refonly.fasta') as f: 
    assembly_msa_ref = f.readlines()
assembly_msa_ref = assembly_msa_ref[1].strip('\n')

# msa_dict_idx[9] will give us the index in the MSA corresponding to the 10th position of Wuhan-Hu-1
# so then msa_dict_idx.index(9) will give us the position of the 10th Wuhan-Hu-1 base in the MSA, allowing us to connect pileup/prior to assembly and ground truth 
# need to access using msa_dict_idx.index(9)+1 in the pileup_dict/prior_dict because of indexing difference
msa_dict_idx = [p for p,b in enumerate(assembly_msa_ref) if b!= "-"] 

header = ['A','C','G','T','-']
pattern1 = re.compile('-[0-9]+[ACGTNacgtn]+')
pattern2 = re.compile('\\+[0-9]+[ACGTNacgtn*#]+')
pattern3 = re.compile('\\$')
pattern4 = re.compile('\\^.{1}')

def get_other_match(m,e):
    # probs             A       C       G       T       -
    match m:
        case 'R': #     X               X
            return [(1-e)/2, e/2    , (1-e)/2, e/2    ,  0]
        case 'Y': #             X               X
            return [e/2    , (1-e)/2, e/2    , (1-e)/2,  0]
        case 'S': #             X       X   
            return [e/2    , (1-e)/2, (1-e)/2, e/2    ,  0]
        case 'W': #     X                       X
            return [(1-e)/2, e/2    , e/2    , (1-e)/2,  0]
        case 'K': #                     X       X
            return [e/2    , e/2    , (1-e)/2, (1-e)/2,  0]
        case 'M': #     X       X 
            return [(1-e)/2, (1-e)/2, e/2    , e/2    ,  0]
        case 'B': #             X       X       X
            return [e      , (1-e)/3, (1-e)/3, (1-e)/3,  0]
        case 'D': #     X               X       X
            return [(1-e)/3, e      , (1-e)/3, (1-e)/3,  0]
        case 'H': #     X       X               X
            return [(1-e)/3, (1-e)/3, e      , (1-e)/3,  0]
        case 'V': #     X       X       X
            return [(1-e)/3, (1-e)/3, (1-e)/3, e      ,  0]

def get_q(ref, m, e):
    qualities = [e/3]*4 + [0]     
    if m==',' or m=='.': 
        m = ref
    if m.upper() in header: 
        qualities[header.index(m.upper())] = 1-e
        return qualities
    elif m=='*' or m=='N': 
        return [1]*5 
    elif m.upper() in ['R','Y','S','W','K','M','B','D','H','V']:
        return get_other_match(m.upper(),e)
    return [np.nan]*5 

bases_prior_dict = {}
bases_prior_path = pileup_results_path + "/tronko_PP/" + prefix + ".PP.txt"
with open(bases_prior_path) as f:
    next(f)
    for line in f:
        (pos,prior_A,prior_C,prior_G,prior_T) = (float(i) for i in str.split(line.strip("\n"), '\t'))
        bases_prior_dict[int(pos)] = np.array([prior_A,prior_C,prior_G,prior_T,0],dtype='float64')
                                              
# where does the ground truth assembly start? where does tronko prior end?
ground_truth_pos_start = [(i+1) for i,(gt,ref) in enumerate(zip(ground_truth_assembly, assembly_msa_ref)) if (gt != "-") & (ref != "-")][0]
tronko_pos_max = np.max(list(bases_prior_dict.keys()))

pileup_dict = {}
pileup_path = pileup_results_path + pileup_file

# pileups are in terms of pos
with open(pileup_path) as f:
    for line in f: 
        # read and store pileup info; we will need it later
        (pos, ref_base, n_reads, match_str, phred_list) = str.split(line.strip('\n'), '\t')[1:6]
        pos = int(pos)
        if (pos < ground_truth_pos_start) or (pos >= tronko_pos_max):
            continue
        is_insertion = (assembly_msa_ref[pos-1] == '-')
        is_deletion = (ground_truth_assembly[pos-1] == '-')
        if is_deletion or is_insertion:
            continue
        else:
            # translate pos to be wrt Wuhan-Hu ref so we can get tronko prior
            pr_list = bases_prior_dict[pos]
            # read likelihoods & calculate posteriors
            match_str = pattern1.sub('', pattern2.sub('', pattern3.sub('', pattern4.sub('', match_str))))
            q_list = np.array([get_q(ref_base, m, pow(10, -(ord(p) - 33) / 10.0)) for m, p in zip(match_str, phred_list)], dtype='float')
            ll_list = np.prod(q_list, axis=0)
            for c in prior_scaling:
                # prior "scaling" : Pi* = (Pi+c/4)/(c+ P1+P2+P3+P4)
                pr_list = [(p + c/4)/(np.sum(pr_list) + c) for p in pr_list]
                pp_list = pr_list
                pp_list = np.multiply(pr_list, ll_list) # try normal calculations
                if (np.sum(ll_list) == 0) or (np.sum(pp_list) == 0): # we have underflow 
                    log_pr_list = np.append(np.log10(pr_list[:4]),0) # get log priors
                    log_ll_list = np.nansum( # get log likelihoods
                        [ np.log10(np.append(q[:4], np.min(q[:4]) / 10.0)) if (np.sum(q) < 5) else np.array([0]*5) for q in q_list ], axis = 0) 
                    if (np.sum(log_ll_list) == 0.0) and (subsample_assembly[pos-1] == '-'): 
                        log_pp_list = np.array([-1]*4 + [0]) # we do not make a call if there is a gap 
                    else:
                        log_pp_list = np.add(log_pr_list, log_ll_list) # get log posteriors if no gap
                    ll_list = [np.power(10, ll - np.max(log_ll_list)) for ll in log_ll_list] # normalize log likelihoods, exponent
                    pp_list = [np.power(10, pp - np.max(log_pp_list)) for pp in log_pp_list] # normalize log posteriors, exponent
                else: # normal calculations can proceed
                    if (np.sum(ll_list) == 5.0) and (subsample_assembly[pos-1] == '-'): 
                        pp_list = np.array([0]*4 + [1]) # we do not make a call if there is a gap 
                        ll_list = [0]*4 + [1]
                    else: 
                        pp_list = np.multiply(pr_list, ll_list) # get posteriors if no gap
                # probabilities need to sum to 1
                pp_list = np.divide(pp_list,sum(pp_list))
                pr_list = np.divide(pr_list,sum(pr_list))
                if c in pileup_dict:
                    pileup_dict[c][pos] = (pr_list, ll_list, pp_list, match_str, ref_base)
                else:
                    pileup_dict[c] = {pos:(pr_list, ll_list, pp_list, match_str, ref_base)}

if call_mode != "roc":
    for c in prior_scaling:
        with open(pileup_results_path + prefix_base + '_prob_calls_pr_scaling' + str(c) + '.txt', 'a') as f:
            for msa_idx,(s, g, ref) in enumerate(zip(subsample_assembly, ground_truth_assembly, assembly_msa_ref)):
                if ref == "-": # this is an insertion, no prior, we skip
                    continue
                elif g == "-": # this is a deletion, no read info, we skip
                    continue
                elif g == "N": # if no call at full depth, also skip
                    continue
                else:
                    pileup_pos = msa_dict_idx.index(msa_idx) + 1
                    if pileup_pos not in pileup_dict[c]:
                        # this position had no reads
                        # call N in LL, then PR for both PR/PP
                        scaled_pr = [(p + c/4)/(np.sum(pr_list) + c) for p in bases_prior_dict[pileup_pos]]
                        r = header[np.nanargmax(scaled_pr)]
                        (r,l,p) = (r, "N", r)
                        mincov3_met = 0
                    else:
                        # this position is in the pileup, pileup dict contains scaled PR, LL, and PP
                        # sometimes pileup has * (no reads), so use prior in these cases (force LL to be N)
                        pileup_entry = pileup_dict[c][pileup_pos]
                        (r,l,p) = (header[np.nanargmax(e)] if np.nanargmax(e) < 4 else "N" for e in [pileup_entry[0], pileup_entry[1], pileup_entry[2]])
                        if l == "N" and p == "N":
                            l,p = (r,r)
                        mincov3_met = int(len(pileup_entry[3]) >= 3)
                        if prefix == prefix_base: # only for the full-coverage base calls
                            f.write(','.join([prefix + "_full", str(pileup_pos), str(mincov3_met), g, s, r, l, p]) + '\n')
                        else:
                            f.write(','.join([prefix, str(pileup_pos), str(mincov3_met), g, s, r, l, p]) + '\n')
elif call_mode == "roc":
    c = 0.001 # choose one prior scaling factor
    posterior_thesholds = [ # add additional columns to mask posterior calls
        0.999999999999, 0.999999999995, 0.99999999999, 0.99999999995, 0.9999999999, 0.9999999995, 0.999999999, 0.999999995, 0.99999999,0.99999995,
        0.9999999, 0.9999995, 0.999999, 0.999995, 0.99999, 0.9999,0.999,0.99,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5]
    with open(pileup_results_path + prefix_base + '_prob_calls_pr_scaling' + str(c) + '_ROC.txt', 'a') as f:
        # header written separately but here to remember order
        # f.write(','.join(['sample', 'pos','mincov3_met'] + ["pp_thresh_pass_" + str(p) for p in posterior_thesholds] + ['GT','SA','PR','LL','PP']) + '\n')
        for msa_idx,(s, g, ref) in enumerate(zip(subsample_assembly, ground_truth_assembly, assembly_msa_ref)):
            if ref == "-": # this is an insertion, no prior, we skip
                    continue
            elif g == "-": # this is a deletion, no read info, we skip
                continue
            elif g == "N": # if no call at full depth, also skip
                continue
            else:
                pileup_pos = msa_dict_idx.index(msa_idx) + 1
                if pileup_pos not in pileup_dict[c]:
                    # this position had no reads
                    # call N in LL, then PR for both PR/PP
                    scaled_pr = [(p + c/4)/(np.sum(pr_list) + c) for p in bases_prior_dict[pileup_pos]]
                    r = header[np.nanargmax(scaled_pr)]
                    thresholds_passed = [int(np.nanmax(scaled_pr) >= pp_thresh) for pp_thresh in posterior_thesholds]
                    (r,l,p) = (r, "N", r)
                    mincov3_met = 0
                else:
                    #this position had reads, pileup dict contains scaled PR, LL, and PP
                    # sometimes pileup has * (no reads), so use prior in these cases (force to be N)
                    pileup_entry = pileup_dict[c][pileup_pos]
                    (r,l,p) = (header[np.nanargmax(e)] if np.nanargmax(e) < 4 else "N" for e in [pileup_entry[0], pileup_entry[1], pileup_entry[2]])
                    if l == "N" and p == "N":
                        l,p = (r,r)
                    thresholds_passed = [int(np.nanmax(pileup_entry[2]) >= pp_thresh) if (np.nanargmax(pileup_entry[2]) < 4) else (int(np.nanmax(pileup_entry[0]) >= pp_thresh)) for pp_thresh in posterior_thesholds]
                    mincov3_met = int(len(pileup_entry[3]) >= 3)
                if prefix == prefix_base: # only for the full-coverage base calls
                    f.write(','.join([prefix + "_full", str(pileup_pos), str(mincov3_met)] + [str(i) for i in thresholds_passed] + [g, s, r, l, p]) + '\n')
                else:
                    f.write(','.join([prefix, str(pileup_pos), str(mincov3_met)] + [str(i) for i in thresholds_passed] + [g, s, r, l, p]) + '\n')

#### previous way, had some bugs
# pr = [header[np.nanargmax(bases_prior_dict[i+1])] if (i+1 in bases_prior_dict) else "-" for i,sa_base in enumerate(sa)]
# ll = [header[np.nanargmax(pileup_dict[0.1][p+1][1])] if ((p+1 in pileup_dict[0.1]) and (np.sum(pileup_dict[0.1][p+1][1]) != 5.0)) else "N" for p in range(len(sa))]

# # posterior calls with no minimum n reads
# pp = { # c is scaling factor, p is position, b is base call
#     c:[header[np.nanargmax(pileup_dict[c][p+1][2])] if ((p+1 in pileup_dict[c]) & (b!='-')) else b for p,b in enumerate(pr)] 
#     for c in prior_scaling}

# sa = [(msa_idx,b) for msa_idx,b in enumerate(subsample_assembly) if b!="-"] # length of SA genome
# gt = [ground_truth_assembly[msa_idx] for msa_idx,b in enumerate(subsample_assembly) if b!="-"] # length of SA genome

# # posterior calls masked with N if n reads < 3
# pp_nreads_mask = [len(pileup_dict[c][p+1][3])>=3 if ((p+1 in pileup_dict[c]) & (b!='-')) else True for p,b in enumerate(pr)]

# else: # to align calls 
#     with open(pileup_results_path + prefix + '_prob_calls.fasta', 'a') as f:
#         for c in prior_scaling:
#             f.write('>' + prefix + '_pr_scale' + str(c) + '_postprob_calls\n')
#             f.write(''.join(pp[c]) + '\n')
#             f.write('>' + prefix + '_pr_scale' + str(c) + '_minreads3' + '_postprob_calls\n')
#             f.write(''.join([b if keep else 'N' for b,keep in zip(pp[c],pp_nreads_mask)]) + '\n')
#         f.write('>' + prefix + '_priorprob_calls\n')
#         f.write(''.join(pr) + '\n')
#         f.write('>' + prefix + '_likelihood_calls\n')
#         f.write(''.join(ll) + '\n')

#### previous debugging
# with open(pileup_results_path + prefix_base + '_calculations_by_pos_pr_thresh' + str(prior_threshold) + '.txt', 'a') as f:
#     f.write(';'.join(["subsample", "pileup_pos", "tronko_pos", 
#                       "GT_call", "SA_call", "PR_call", "LL_call", "PP_call",
#                       "PR_probs", "LL_probs", "PP_probs"]) + '\n')
#     for (pos,g,s,r,l,p) in zip(positions_idx, gt, sa, pr, ll, pp):
#         if g != 'N':
#             f.write(';'.join([prefix, str(pos), str(msa_idx_ref[pos]), 
#                               g, s, r, l, p,
#                               ",".join([str(p) for p in pileup_dict[pos][0]]), 
#                               ",".join([str(p) for p in pileup_dict[pos][1]]), 
#                               ",".join([str(p) for p in pileup_dict[pos][2]])
#                               ]) + '\n')