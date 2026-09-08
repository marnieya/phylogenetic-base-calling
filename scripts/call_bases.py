#!/home/marniella/miniconda3/bin/python

import os
os.environ["OMP_NUM_THREADS"] = '4' 
os.environ["OPENBLAS_NUM_THREADS"] = '4'
os.environ["MKL_NUM_THREADS"] = '4'

import numpy as np
import re
import sys

pileup_file = sys.argv[1]
pileup_results_path = sys.argv[2]

#pileup_file = "SRR26069290_sub20_piledup.txt"
#pileup_results_path = "/path_here/SRR26069290_files/"

prior_scaling = ["MLE", 0.001, 0.0]

if "SRR" in pileup_file:
    # for SRA
    prefix = re.sub('_piledup.txt','', pileup_file)
    prefix_base = re.sub('_sub.+', '', prefix)
else:
    # for GISAID
    prefix = re.sub('_piledup.txt','', pileup_file)
    prefix_base =  re.sub('.sim.+', '', prefix) # re.sub('_[0-9]*\.*[0-9]+x_cov.sim.+', '', prefix)

if prefix != prefix_base:
    with open(pileup_results_path + prefix + '_gt_ref_msa.fasta') as f: 
        gt_ref_msa = f.readlines()
    ground_truth_assembly = gt_ref_msa[1].strip('\n') # the genome assembled from original reads
    wh1_gisaid_msa_ref = gt_ref_msa[3].strip('\n') # the Wuhan-Hu-1 reference genome
    subsample_assembly = gt_ref_msa[5].strip('\n') # the genome assembled from subsampled reads
else:
    with open(pileup_results_path + prefix + '_pwise_msa.fasta') as f: 
        gt_ref_msa = f.readlines()
    ground_truth_assembly = gt_ref_msa[1].strip('\n')
    subsample_assembly = gt_ref_msa[1].strip('\n')
    wh1_gisaid_msa_ref = gt_ref_msa[3].strip('\n')

# msa_dict_idx[9] will give us the index in the MSA corresponding to the 10th position of Wuhan-Hu-1
msa_dict_idx = [p for p,b in enumerate(wh1_gisaid_msa_ref) if b!= "-"] 

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

const_bases_prior_dict = {}
with open(pileup_results_path + "/tronko_PP/" + re.sub('_rep.*','',prefix) + ".PP.txt") as f:
    next(f)
    for line in f:
        (pos,prior_A,prior_C,prior_G,prior_T) = (float(i) for i in str.split(line.strip("\n"), '\t'))
        const_bases_prior_dict[int(pos)] = np.array([prior_A,prior_C,prior_G,prior_T,0],dtype='float64')

mle_bases_prior_dict = {}
with open(pileup_results_path + "/tronko_PP/corrected." + re.sub('_rep.*','',prefix) + ".PP.txt") as f:
    next(f)
    for line in f:
        (pos,prior_A,prior_C,prior_G,prior_T) = (float(i) for i in str.split(line.strip("\n"), '\t'))
        mle_bases_prior_dict[int(pos)] = np.array([prior_A,prior_C,prior_G,prior_T,0],dtype='float64')

# where does the ground truth assembly start? where does tronko prior end?
ground_truth_pos_start = [(i+1) for i,(gt,ref) in enumerate(zip(ground_truth_assembly, wh1_gisaid_msa_ref)) if (gt != "-") & (ref != "-")][0]

pileup_dict = {}
pileup_path = pileup_results_path + pileup_file

# pileups are in terms of pos
with open(pileup_path) as f:
    for line in f: 
        # read and store pileup info; we will need it later
        (pos, ref_base, n_reads, match_str, phred_list) = str.split(line.strip('\n'), '\t')[1:6]
        pos = int(pos)
        if (pos < ground_truth_pos_start) or (pos >= len(msa_dict_idx)):
            continue
        msa_idx = msa_dict_idx[pos-1]
        is_insertion = (wh1_gisaid_msa_ref[msa_idx] == '-')
        is_deletion = (ground_truth_assembly[msa_idx] == '-')
        if is_deletion or is_insertion:
            continue
        else:
            # read likelihoods & calculate posteriors
            match_str = pattern1.sub('', pattern2.sub('', pattern3.sub('', pattern4.sub('', match_str))))
            q_list = np.array([get_q(ref_base, m, pow(10, -(ord(p) - 33) / 10.0)) for m, p in zip(match_str, phred_list)], dtype='float')
            ll_list_raw = np.prod(q_list, axis=0, dtype=float)
            for c_name in prior_scaling:
                ll_list = ll_list_raw.copy() # fresh linear copy each iteration; ll_list gets overwritten below
                # translate pos to be wrt Wuhan-Hu ref so we can get tronko prior
                if c_name != "MLE":
                    pr_list = const_bases_prior_dict[pos] #unaltered priors, to be scaled
                    c=c_name
                else:
                    pr_list = mle_bases_prior_dict[pos] #priors already corrected, treat scaling as c=0
                    c=0.0
                # prior "scaling" : Pi* = (Pi+c/4)/(c+ P1+P2+P3+P4)
                pr_list = [(p + c/4)/(np.sum(pr_list) + c) for p in pr_list[:4]] + [pr_list[4]]
                pp_list = pr_list
                pp_list = np.multiply(pr_list, ll_list, dtype=float) # try normal calculations
                if (np.sum(ll_list) == 0) or (np.sum(pp_list) == 0): # we have underflow 
                    log_pr_list = np.append(np.log10(pr_list[:4]),0) # get log priors
                    log_ll_list = np.nansum( # get log likelihoods
                        [ np.log10(np.append(q[:4], np.min(q[:4]) / 10.0)) if (np.sum(q) < 5) else np.array([0]*5) for q in q_list ], axis = 0) 
                    if (np.sum(log_ll_list) == 0.0) and (subsample_assembly[pos-1] == '-'): 
                        log_pp_list = np.array([-1]*4 + [0]) # we do not make a call if there is a gap 
                    else:
                        log_pp_list = np.add(log_pr_list, log_ll_list) # get log posteriors if no gap
                    ll_list = log_ll_list # keep raw log likelihoods as-is, no shift/exponentiate (that round trip is what collapsed them to 0)
                    pp_list = [np.power(10, pp - np.max(log_pp_list), dtype=float) for pp in log_pp_list] # normalize log posteriors, exponent
                else: # normal calculations can proceed
                    if (np.sum(ll_list) == 5.0) and (subsample_assembly[pos-1] == '-'):
                        pp_list = np.array([0]*4 + [1]) # we do not make a call if there is a gap
                        ll_list = np.array([-np.inf]*4 + [0])
                    else:
                        pp_list = np.multiply(pr_list, ll_list) # get posteriors if no gap
                        with np.errstate(divide='ignore'):
                            ll_list = np.append(np.log10(ll_list[:4]), -np.inf) # convert to log space to match underflow branch; index 4 is a placeholder that must never win the argmax
                # probabilities need to sum to 1
                pp_list = np.append(np.divide(pp_list[:4],sum(pp_list[:4]), dtype=float), pp_list[4])
                pr_list = np.append(np.divide(pr_list[:4],sum(pr_list[:4]), dtype=float), pr_list[4])
                if c_name in pileup_dict:
                    pileup_dict[c_name][pos] = (pr_list, ll_list, pp_list, match_str, ref_base)
                else:
                    pileup_dict[c_name] = {pos:(pr_list, ll_list, pp_list, match_str, ref_base)}

sa_nogaps = [sa for sa,ref in zip(subsample_assembly, wh1_gisaid_msa_ref) if ref!= "-"]
gt_nogaps = [gt for gt,ref in zip(ground_truth_assembly, wh1_gisaid_msa_ref) if ref!="-"]
ref_nogaps = [ref for ref in wh1_gisaid_msa_ref if ref!="-"]

for c_name in prior_scaling:
    with open(pileup_results_path + prefix_base + '_prob_calls_pr_scaling' + str(c_name) + '_precise.txt', 'a') as f:
        for pos,(s, g, ref) in enumerate(zip(sa_nogaps, gt_nogaps, ref_nogaps)):
            if g == "N": # if no call at full depth, skip
                continue
            else:
                pileup_pos = pos + 1
                if pileup_pos not in pileup_dict[c_name]:
                    # position not in the pileup; therefore we skip
                    continue
                else:
                    # position in the pileup; pileup dict contains scaled PR, LL, and PP
                    pileup_entry = pileup_dict[c_name][pileup_pos]
                    pileup_reads = pileup_reads = len(pileup_entry[3])
                    (r,l,p) = (header[np.nanargmax(e)] if np.nanargmax(e) < 4 else "N" for e in [pileup_entry[0], pileup_entry[1], pileup_entry[2] ])
                    pr = np.log10(np.nanmax(pileup_entry[0])) if np.nanargmax(pileup_entry[0]) < 4 else np.nan
                    ll = np.nanmax(pileup_entry[1]) if np.nanargmax(pileup_entry[1]) < 4 else np.nan # already log-scale, no further log10
                    pp = np.log10(np.nanmax(pileup_entry[2])) if np.nanargmax(pileup_entry[2]) < 4 else np.nan
                    if (l == "N" and p == "N"):
                        # position in the pileup but has * (no reads) 
                        continue
                    else:
                        # compare methods only if a position was both: (1) in the pileup and (2) had reads
                        if prefix == prefix_base:
                            # only for the full datasets
                            f.write(','.join([prefix + "_full", str(pileup_pos), str(pileup_reads), g, s, r, l, p, str(pr), str(ll), str(pp)]) + '\n')
                        else:
                            # for the subsamples 
                            f.write(','.join([prefix, str(pileup_pos), str(pileup_reads), g, s, r, l, p, str(pr), str(ll), str(pp)]) + '\n')
