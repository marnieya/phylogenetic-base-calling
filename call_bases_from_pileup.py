#!/Users/marniella/anaconda3/envs/py310/bin/python

import numpy as np
import re
import sys

pileup_file = sys.argv[1] 
pileup_results_path = sys.argv[2] 
prefix = re.sub('_piledup.txt','', pileup_file)

# sra
# pileup_file = "SRR25117579_piledup.txt"
# pileup_results_path = "/Users/marniella/research/nielsen_lab/project-extra-files/phylobc/raw_read_results/"

home_dir = '/Users/marniella/research/nielsen_lab/'

with open(home_dir + 'Wuhan_Hu_reference_MSA.fasta') as f: 
    reference = f.readlines()
reference = reference[1].strip('\n')

with open(pileup_results_path + prefix + '_aligned.fasta') as f: 
    og_reference = f.readlines()
og_reference = og_reference[1].strip('\n')

msa_og_j = [j for j,ref_base in enumerate(reference) if ref_base!='-']

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
bases_prior_path = home_dir + msa_dir + 'OU_PP_indexed.txt'
with open(bases_prior_path) as f:
    next(f)
    for line in f:
        (msa_pos,prior_A,prior_C,prior_G,prior_T) = str.split(line.strip("\n"), '\t')
        bases_prior_dict[int(msa_pos)] = np.array([prior_A,prior_C,prior_G,prior_T,0],dtype='float64')

pileup_dict = {}
pileup_path = pileup_results_path + pileup_file

with open(pileup_path) as f:
    for line in f: 
        (pos, ref_base, n_reads, match_str, phred_list) = str.split(line.strip('\n'), '\t')[1:6]
        pos = int(pos) - 1
        match_str = pattern1.sub('', pattern2.sub('', pattern3.sub('', pattern4.sub('', match_str))))
        pr_list = bases_prior_dict[msa_og_j[pos]]
        q_list = np.array([get_q(ref_base, m, pow(10, -(ord(p) - 33) / 10.0)) for m, p in zip(match_str, phred_list)])
        ll_list = np.prod(q_list, axis=0)
        pp_list = np.array([0.0]*len(ll_list))
        if (0 in ll_list[:4]) : # we might have underflow in likelihood of reads; careful with log of 0
            log_pr_list = np.append(np.log10(pr_list[:4]),0) # get log priors
            log_ll_list = np.nansum( # get log likelihoods
                [ np.log10(np.append(q[:4], np.min(q[:4]) / 10.0)) if (np.sum(q) < 5) else np.array([0]*5) for q in q_list ], axis = 0) 
            if (np.sum(log_ll_list) == 0.0) and (og_reference[msa_og_j[pos]] == '-'): 
                log_pp_list = np.array([-1]*4 + [0])
            else:
                log_pp_list = np.add(log_pr_list, log_ll_list) # get log posteriors
            ll_list = [np.power(10, ll - np.max(log_ll_list)) for ll in log_ll_list] # normalize log likelihoods, exponent
            pp_list = [np.power(10, pp - np.max(log_pp_list)) for pp in log_pp_list] # normalize log posteriors, exponent
        else: # can resume normal calculations
            if (np.sum(ll_list) == 5.0) and (og_reference[msa_og_j[pos]] == '-'): 
                pp_list = np.array([0]*4 + [1])
            else: 
                pp_list = np.multiply(pr_list, ll_list) 
        pileup_dict[pos] = (pr_list, ll_list, pp_list) 
        

positions_idx = sorted(list(pileup_dict.keys()))

called_bases_posterior = [header[np.nanargmax(p)] for p in [np.divide(pileup_dict[p][2],sum(pileup_dict[p][2])) for p in positions_idx]]
called_bases_likelihood = [header[np.nanargmax(p)] for p in [np.divide(pileup_dict[p][1],sum(pileup_dict[p][1])) for p in positions_idx]]

pr = [header[np.nanargmax(pileup_dict[p][0])] for p in positions_idx]
pp = called_bases_posterior
ll = called_bases_likelihood

with open(pileup_results_path + prefix + '_calls.fasta', 'a') as f:
    f.write('>' + prefix + '_postprob_calls\n')
    f.write(''.join(pp) + '\n')
    f.write('>' + prefix + '_priorprob_calls\n')
    f.write(''.join(pr) + '\n')
    f.write('>' + prefix + '_likelihood_calls\n')
    f.write(''.join(ll) + '\n')