
## slightlly different approach: as soon as we detect different bases, move on
## i.e. only compute seq sim to check for exactly identical sequences

## ABANDONED: I'm gonna write something like this but in C++ instead; 
## the idea is to have more efficient 2d array navigation using pointers
## and do away with python type-less shenanigans, my own mem alloc

## but do the crop to avoid n's and gaps
cropped_seqs = [seq[1000:-1000] for seq in seqs]
nucs = ['a','c','g','t','-']

def compare_bases(base1,base2):
    # returns bool of whether the two bases are the same nucleotide
    # but if it's n and a nucleotide, we don't necessarily wanna say that is distinct
    if (base1.lower() in nucs) & (base2.lower() in nucs):
        return int(base1==base2)
    elif (base1.lower() != 'n') & (base2.lower() == 'n'):
        return 1
    elif (base1.lower() == 'n') & (base2.lower() != 'n'):
        return 1
    else:
        return 0

def get_seq_sim(seq1, seq2):
    # this function goes through the two sequences
    # returns 0 if they are distinct, or 1 if they're identical
    for b1,b2 in zip(list(seq1),list(seq2)):
        #print(b1 + " vs. " + b2 + str(compare_bases(b1,b2)))
        if compare_bases(b1,b2) == 0:
            return 0
    return 1

l = len(seqs)
seq_sim_mat = np.zeros((l,l))
for i,s1 in enumerate(cropped_seqs):
    for j,s2 in enumerate(cropped_seqs):
        if j <= l - i: # we really only care about the lower triangle
            seq_sim_mat[i,j] = get_seq_sim(s1,s2)
