#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import sys

pileup_file = sys.argv[1] # pileup file
pileup_results_path = sys.argv[2] # folder containing pileup file
pileup_prefix = re.sub('_piledup.txt','',pileup_file)

### ============ get prior probabilities ============

home_dir = '/Users/marniella/research/rotations/r3_nielsen/'
bases_prior_df = pd.read_csv(home_dir + 'OU061397_1_PP.txt', sep = '\t', index_col='position')
bases_prior_df.columns = ['A','C','G','T']

### ============ get likelihood from raw reads ============

# read in quality scores / alignment info at all positions
# positions are 1-indexed, don't have gaps
with open(home_dir + pileup_path) as f:
    positions_temp = f.readlines()
    f.close()
positions_temp = [str.split(p.strip('\n'), '\t')[1:6] for p in positions_temp]
positions_df = pd.DataFrame(data=positions_temp, columns=['pos', 'ref_base', 'n_reads','match_str', 'read_base_qualities'])

def clean_up_match_str(ms):
    # we care about matches or mismatches (,.actgnACTGN)
    # * are deleted bases, and $ or indel strings are referring to what comes after this position-- we will not use that
    new_ms = re.sub('-[0-9]+[ACGTNacgtn]+', '', re.sub('\\+[0-9]+[ACGTNacgtn*#]+', '', re.sub('\\$','', re.sub('\\^.{1}', '', ms))))
    return new_ms

positions_df['clean_match_str'] = [clean_up_match_str(s) for s in positions_df['match_str']]
positions_df['pos'] = [int(p)-1 for p in positions_df['pos']]
positions = positions_df.values.tolist()

# build a dict positions_dict of position:(reference base, CLEAN match string, error rate) 
unique_phred = list(set([i for s in [list(set(qs)) for qs in positions_df['read_base_qualities']] for i in s]) )
phred_e_dict = {phred:pow(10, -(ord(phred)-33)/10.0) for phred in unique_phred} # calculate error rate from phred score
positions_dict = {int(p[0]):(p[1], p[5], [phred_e_dict[phred] for phred in list(p[4]) if phred in phred_e_dict]) for p in positions}
positions_idx = sorted(list(positions_dict.keys())) # bases with reads, in order

# read in the gapped MSA
with open(home_dir + 'alignment101_mafft/msa_ref.fasta') as f: 
    reference = f.readlines()
    f.close()
reference = reference[1].strip('\n')

with open(home_dir + 'alignment101_mafft/msa_og.fasta') as f: 
    og_reference = f.readlines()
    f.close()
og_reference = og_reference[1].strip('\n')

# 'pos' refers to the "absolute" position on the Wuhan-Hu reference genome where the read aligned
# 'i' refers to position on the MSA where the ground-truth OU genome does not have a gap (i.e. content of array is MSA pos, index is absolute pos)
# 'j' refers to position on the MSA where the Wuhan-Hu reference genome does not have a gap (i.e. content of array is MSA pos, index is absolute pos)

# all are 0-indexed at this point in the program
msa_og_i = [i for i,ref_base in enumerate(og_reference) if ref_base!='-']
msa_og_j = [j for j,ref_base in enumerate(reference) if ref_base!='-']


### ============ calculate posteriors ============

# algorithm: get pos from positions_dict 
#            get msa_pos from msa_og_j[pos]
#            get prior from bases_prior_df.loc[[msa_pos]]
#            get likelihood from all read info in positions_dict[pos]
#            get unscaled posteriors for each base [A,C,G,T]
#            sum across and scale posteriors

header = ['A','C','G','T','-']
def get_q(ref, match, e):
    # get the read's adjusted quality score from the match string and error rate (e) from the phred score
    qualities = [e/3]*4 + [0] # assume there is a match
    if match==',' or match=='.': # read base is same as reference
        match = ref
    if match.upper() in header: # change quality score based on the position of the base
        qualities[header.index(match.upper())] = 1-e
        return qualities
    if match=='*': # if no match, we weigh everything equally 
        return [1]*5 
    return [None]*5 # if match string not a valid base / does not match ref

def get_likelihood_prod(pos):
    r = positions_dict[pos][0] # get reference base at pos
    # then get the likelihoods for each base (q), using match string (m) & error rate (e)
    ll = [get_q(r,m,e) for m,e in zip(positions_dict[pos][1], positions_dict[pos][2])] 
    # put in df for easy multiply across all reads
    ll_df = pd.DataFrame(columns = header, data = ll).dropna().reset_index(drop=True)
    ll_list = [np.prod(ll_df['A']), np.prod(ll_df['C']), np.prod(ll_df['G']), np.prod(ll_df['T']), np.prod(ll_df['-'])]
    return ll_list

def get_unscaled_posterior(pos, use_prior):
    # map raw read pos --> MSA pos
    msa_pos = msa_og_j[pos]
    og_msa_base = og_reference[msa_pos]
    # get the likelihoods for A,C,G,T,-
    ll = get_likelihood_prod(pos) 
    pr = [float(p) for p in [bases_prior_df.loc[[msa_pos]]['A'], bases_prior_df.loc[[msa_pos]]['C'], 
            bases_prior_df.loc[[msa_pos]]['G'], bases_prior_df.loc[[msa_pos]]['T']]]
    if (np.sum(ll) == 5.0) and (og_msa_base == '-'): # no reads aligned and it's a gap in the MSA
        return [0]*4 + [1]
    else: # if we either had reads, or the MSA is not a gap, we use the prior or the likelihood
        if use_prior:
            # fetch prior from tronko_build
            pr = pr + [0]
            # return prior*ll for each base in that same order
            return [prior*likelihood for prior,likelihood in zip(pr,ll)] 
        else:
            return ll

### ============ compare phylogenetic vs. reads-only base calling ============

# get posterior for each position
unscaled_posterior = [get_unscaled_posterior(pos=p,use_prior=True) for p in positions_idx]
unscaled_likelihood = [get_unscaled_posterior(pos=p,use_prior=False) for p in positions_idx]

# make into dataframes for easy scaling
unscaled_posterior_df = pd.DataFrame(columns = header, data = unscaled_posterior, index = positions_idx)
scaled_posterior_df = unscaled_posterior_df.div(unscaled_posterior_df.sum(axis=1), axis=0)

unscaled_likelihood_df = pd.DataFrame(columns = header, data = unscaled_likelihood, index = positions_idx)
scaled_likelihood_df = unscaled_likelihood_df.div(unscaled_likelihood_df.sum(axis=1), axis=0)

called_bases_posterior = scaled_posterior_df.dropna().idxmax(axis=1)
called_bases_likelihood = scaled_likelihood_df.dropna().idxmax(axis=1)

# how similar are likelihood and posterior base calling?
## first make sure our accuracy isn't affected by N's
def compare_bases(a,b):
    if a=='N' or b=='N':
        return -1
    if b=='-': # assuming b is the og/'target' genome re-constructed from MSA
        return -1 # we will ignore positions that are gaps in the MSA
    return a==b

pp = [called_bases_posterior[p] for p in called_bases_posterior.index.tolist()]
ll = [called_bases_likelihood[p] for p in called_bases_likelihood.index.tolist()]
pp_ll_agreements = [a for a in [compare_bases(p,l) for p,l in zip(pp, ll)] if a > -1]

# how similar are they to the original called bases?
og = [og_reference[msa_og_j[pos]] for pos in called_bases_posterior.index.tolist()]
pp_og_agreements = [a for a in [compare_bases(p,o) for p,o in zip(pp, og)] if a > -1]
ll_og_agreements = [a for a in [compare_bases(l,o) for l,o in zip(ll, og)] if a > -1]

# output stats
# agreement b/ posterior prob (pp) & likelihood-only (ll) calling
print(', '.join(pileup_prefix.split('_') + ['pp_ll_agreements'] + [str(i) for i in (np.mean(pp_ll_agreements), np.sum(pp_ll_agreements), len(pp_ll_agreements))]))
# agreement b/ consensus & phylo pp calling
print(', '.join(pileup_prefix.split('_') + ['pp_og_agreements'] + [str(i) for i in (np.mean(pp_og_agreements), np.sum(pp_og_agreements), len(pp_og_agreements))]))
# agreement b/ consensus & ll calling
print(', '.join(pileup_prefix.split('_') + ['ll_og_agreements'] + [str(i) for i in (np.mean(ll_og_agreements), np.sum(ll_og_agreements), len(ll_og_agreements))]))

