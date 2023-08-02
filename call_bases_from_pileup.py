#!/Users/marniella/anaconda3/envs/py310/bin/python

import numpy as np
import re
import sys

pileup_file = sys.argv[1] 
pileup_results_path = sys.argv[2] 

# og
# pileup_file = "10x_reads_50bp_verr_is50_sim_piledup.txt"
# pileup_results_path = "step3_phylobc/sim_results_verr/sim_results_10x_reads_50bp_verr_is50/"

# # sra
# pileup_file = "SRR25117579_piledup.txt"
# pileup_results_path = "step3_phylobc/raw_read_results/"

home_dir = '/Users/marniella/research/nielsen_lab/phylogenetic-base-calling/'
msa_dir = 'ou_pruned/'

with open(home_dir + 'justmn.fasta') as f: 
    reference = f.readlines()
reference = reference[1].strip('\n')

with open(home_dir + 'justou.fasta') as f: 
    og_reference = f.readlines()
og_reference = og_reference[1].strip('\n')

msa_og_j = [j for j,ref_base in enumerate(reference) if ref_base!='-']

header = ['A','C','G','T','-']
pattern1 = re.compile('-[0-9]+[ACGTNacgtn]+')
pattern2 = re.compile('\\+[0-9]+[ACGTNacgtn*#]+')
pattern3 = re.compile('\\$')
pattern4 = re.compile('\\^.{1}')

def get_other_match(m,e):
    # probs          A     C    G    T  -
    match m:
        case 'R': #  X          X
            return [1-e, e/2, 1-e, e/2, 0]
        case 'Y': #        X        X
            return [e/2, 1-e, e/2, 1-e, 0]
        case 'S': #        X    X
            return [e/2, 1-e, 1-e, e/2, 0]
        case 'W': #  X              X
            return [1-e, e/2, e/2, 1-e, 0]
        case 'K': #             X   X
            return [e/2, e/2, 1-e, 1-e, 0]
        case 'M': # X      X 
            return [1-e, 1-e, e/2, e/2, 0]
        case 'B': #        X    X   X
            return [e/2, 1-e, 1-e, 1-e, 0]
        case 'D': # X           X   X
            return [1-e, e/2, 1-e, 1-e, 0]
        case 'H': # X      X        X
            return [1-e, 1-e, e/2, 1-e, 0]
        case 'V': # X      X    X
            return [1-e, 1-e, 1-e, e/2, 0]

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

# TODO: review how i'm doing the logs and the underflow -- can't have a 0 then call prob=1 after exponent
bases_prior_dict = {}
bases_prior_path = home_dir + msa_dir + 'OU_PP_indexed.txt'
with open(bases_prior_path) as f:
    next(f)
    for line in f:
        (msa_pos,prior_A,prior_C,prior_G,prior_T) = str.split(line.strip("\n"), '\t')
        bases_prior_dict[int(msa_pos)] = np.array([prior_A,prior_C,prior_G,prior_T,0],dtype='float64')

pileup_dict = {}
debug_dict = {}
pileup_path = pileup_results_path + pileup_file

with open(pileup_path) as f:
    for line in f: 
        (pos, ref_base, n_reads, match_str, phred_list) = str.split(line.strip('\n'), '\t')[1:6]
        pos = int(pos) - 1
        match_str = pattern1.sub('', pattern2.sub('', pattern3.sub('', pattern4.sub('', match_str))))
        pr_list = bases_prior_dict[msa_og_j[pos]]
        q_list = np.array([get_q(ref_base, m, pow(10, -(ord(p) - 33) / 10.0)) for m, p in zip(match_str, phred_list)])
        ll_list = np.prod(q_list, axis=0)
        debug_dict[pos] = (n_reads, match_str, phred_list, q_list)
        if(len(ll_list) == 4):
            print(pos)
        pp_list = np.array([0.0]*len(ll_list))
        if (0 in ll_list[:4]): # we might have underflow in likelihood of reads
            log_pr_list = np.log10(pr_list)
            log_ll_list = np.nansum([np.append(np.log10(q[:4]), q[4]) for q in q_list], axis = 0) # get log likelihoods
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

def compare_bases(a,b):
    if a=='N' or b=='N':
        return -1
    if (b=='-' and a !='-') or (b!='-' and a =='-'): 
        return -1 
    return a==b

called_bases_posterior = [header[np.nanargmax(p)] for p in [np.divide(pileup_dict[p][2],sum(pileup_dict[p][2])) for p in positions_idx]]
called_bases_likelihood = [header[np.nanargmax(p)] for p in [np.divide(pileup_dict[p][1],sum(pileup_dict[p][1])) for p in positions_idx]]

pr = [header[np.nanargmax(pileup_dict[p][0])] for p in positions_idx]
pp = called_bases_posterior
ll = called_bases_likelihood
og = [og_reference[msa_og_j[pos]] for pos in positions_idx]

def get_all_comparisons(list1, list2):
    return([compare_bases(l1,l2) for l1,l2 in zip(list1, list2)])

pr_og_all_comp = get_all_comparisons(pr, og)
ll_og_all_comp = get_all_comparisons(pp, og)
pp_og_all_comp = get_all_comparisons(ll, og)

def get_agreements(list1, list2, pos_list):
    ag_all = [(compare_bases(l1, l2), pos) for l1, l2, pos in zip(list1, list2, pos_list) if compare_bases(l1, l2) > -1]
    ag = [a for a, pos in ag_all]
    metrics = [np.mean(ag), np.sum(ag), len(ag), np.sum([not a for a in ag])]
    return [str(m) for m in metrics]

pr_og_agreements = get_agreements(pr, og, positions_idx) 
pp_og_agreements = get_agreements(pp, og, positions_idx) 
ll_og_agreements = get_agreements(ll, og, positions_idx) 
print(','.join([pileup_file, pr_og_agreements[0], pp_og_agreements[0], ll_og_agreements[0]]))