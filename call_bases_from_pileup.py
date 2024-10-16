#!/Users/marniella/anaconda3/envs/py310/bin/python

import numpy as np
import re
import sys

pileup_file = sys.argv[1] 
pileup_results_path = sys.argv[2] 

prior_threshold = [0.0001, 0.001, 0.01, 0.1]

# for SRA
prefix = re.sub('_piledup.txt','', pileup_file)
prefix_base = re.sub('_sub.+', '', prefix)

# for GISAID
# prefix = re.sub('_piledup.txt','', pileup_file)
# prefix_base = re.sub('_[0-9]*\.*[0-9]+x_cov.sim.+', '', prefix) # re.sub('.sim.+', '', prefix)

with open(pileup_results_path + prefix + '_sa.fasta') as f: 
    subsample_assembly = f.readlines()
subsample_assembly = subsample_assembly[1].strip('\n')

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
                                              
pileup_dict = {}
pileup_path = pileup_results_path + pileup_file

# pileups are in terms of pos
with open(pileup_path) as f:
    for line in f: 
        # read and store pileup info
        (pos, ref_base, n_reads, match_str, phred_list) = str.split(line.strip('\n'), '\t')[1:6]
        pos = int(pos) 
        # if pos does not have a tronko prior, skip
        if pos not in bases_prior_dict:
            continue
        pr_list = bases_prior_dict[pos]
        # read likelihoods & calculate posteriors
        match_str = pattern1.sub('', pattern2.sub('', pattern3.sub('', pattern4.sub('', match_str))))
        q_list = np.array([get_q(ref_base, m, pow(10, -(ord(p) - 33) / 10.0)) for m, p in zip(match_str, phred_list)], dtype='float')
        ll_list = np.prod(q_list, axis=0)
        for thresh in prior_threshold:
            # upper limit on priors Pi* = (Pi+c/4)/(c+ P1+P2+P3+P4)
            c = thresh 
            pr_list = [(p + c/4)/(np.sum(pr_list) + c) for p in pr_list]
            pp_list = pr_list
            pp_list = np.multiply(pr_list, ll_list) # try normal calculations
            if (np.sum(ll_list) == 0) or (np.sum(pp_list) == 0): # we have underflow 
                log_pr_list = np.append(np.log10(pr_list[:4]),0) # get log priors
                log_ll_list = np.nansum( # get log likelihoods
                    [ np.log10(np.append(q[:4], np.min(q[:4]) / 10.0)) if (np.sum(q) < 5) else np.array([0]*5) for q in q_list ], axis = 0) 
                if (np.sum(log_ll_list) == 0.0) and (subsample_assembly[pos-1] == '-'): 
                    log_pp_list = np.array([-1]*4 + [0]) # we do not make a call if there is a gap in the assembly
                else:
                    log_pp_list = np.add(log_pr_list, log_ll_list) # get log posteriors if no gap
                ll_list = [np.power(10, ll - np.max(log_ll_list)) for ll in log_ll_list] # normalize log likelihoods, exponent
                pp_list = [np.power(10, pp - np.max(log_pp_list)) for pp in log_pp_list] # normalize log posteriors, exponent
            else: # normal calculations can proceed
                if (np.sum(ll_list) == 5.0) and (subsample_assembly[pos-1] == '-'): 
                    pp_list = np.array([0]*4 + [1]) # we do not make a call if there is a gap in the assembly
                    ll_list = [0]*4 + [1]
                else: 
                    pp_list = np.multiply(pr_list, ll_list) # get posteriors if no gap
            # probabilities need to sum to 1
            pp_list = np.divide(pp_list,sum(pp_list))
            pr_list = np.divide(pr_list,sum(pr_list))
            if thresh in pileup_dict:
                pileup_dict[thresh][pos] = (pr_list, ll_list, pp_list, match_str, ref_base)
            else:
                pileup_dict[thresh] = {pos:(pr_list, ll_list, pp_list, match_str, ref_base)}

positions_idx = sorted(list(pileup_dict.keys()))

pr = [header[np.nanargmax(bases_prior_dict[p])] if (subsample_assembly[p-1] != "-") else "-" for p in bases_prior_dict]
ll = [header[np.nanargmax(pileup_dict[0.1][p+1][1])] if (p+1 in pileup_dict[0.1]) else "N" for p in range(len(pr))]

pp = {
    thresh:[header[np.nanargmax(pileup_dict[thresh][p+1][2])] if ((p+1 in pileup_dict[thresh]) & (b!='-')) else b for p,b in enumerate(pr)] 
    for thresh in prior_threshold}

# only for the full-coverage base calls
if prefix == prefix_base:
    prefix = prefix + "_full"
    prefix_base = prefix_base + "_full"

with open(pileup_results_path + prefix + '_prob_calls.fasta', 'a') as f:
    for thresh in prior_threshold:
        f.write('>' + prefix + '_pr_thresh' + str(thresh) + '_postprob_calls\n')
        f.write(''.join(pp[thresh]) + '\n')
    f.write('>' + prefix + '_priorprob_calls\n')
    f.write(''.join(pr) + '\n')
    f.write('>' + prefix + '_likelihood_calls\n')
    f.write(''.join(ll) + '\n')

# now instead of the calls_by_pos we need to deal with msat for scoring


