#!/usr/bin/env python3

import numpy as np
import re
import sys

pileup_file = sys.argv[1] 
pileup_results_path = sys.argv[2] 

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

def get_q(ref, match, e):
    qualities = [e/3]*4 + [0] 
    if match==',' or match=='.': 
        match = ref
    if match.upper() in header: 
        qualities[header.index(match.upper())] = 1-e
        return qualities
    if match=='*': 
        return [1]*5 
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
        q_list = np.array([get_q(ref_base, m, pow(10, -(ord(p) - 33) / 10.0)) for m, p in zip(match_str, phred_list)])
        ll_list = np.nanprod(q_list, axis = 0)
        msa_pos = msa_og_j[pos]
        pr_list = bases_prior_dict[msa_pos]
        pp_list = np.array([0.0]*len(ll_list))
        if (np.sum(ll_list) == 5.0) and (og_reference[msa_pos] == '-'): 
            pp_list = np.array([0]*4 + [1])
        else: 
            pp_list = np.multiply(pr_list,ll_list)
        pileup_dict[pos] = (pr_list, ll_list, pp_list) # prior, likelihood, posterior

positions_idx = sorted(list(pileup_dict.keys()))

def compare_bases(a,b):
    if a=='N' or b=='N':
        return -1
    if b=='-': 
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