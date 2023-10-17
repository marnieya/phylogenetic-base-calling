#!/Users/marniella/anaconda3/envs/py310/bin/python
import sys
dir = sys.argv[1]
sra = sys.argv[2]

def compare_bases(a,b):
    if (b=='N') or ((a == '-') & (b == '-')): 
        return 'N' # cases of gaps or ground truth = N skipped
    return a==b

pos=1
with open(dir + sra + "_calls_msat.txt","r") as msat_f:
    h = msat_f.readline().strip('\n').split(',') 
    with open(dir + sra + "_calls_metrics.txt", "a") as metrics_f:
        for line in msat_f:
            # each line is a position in the genome
            calls = list(line)
            gt = calls[h.index(sra)]
            if compare_bases('check', gt) == 'N':
                next # skip position if ground truth = N
            for subsample in [(sra + "_sub" + str(s)) for s in [0.2, 0.5, 0.8, 1, 2, 3, 4, 5, 10, 15, 20]]:
                pp = calls[h.index(subsample + "_postprob_calls")]
                ll = calls[h.index(subsample + "_likelihood_calls")]
                if compare_bases(pp, gt) == 'N':
                    next # skip position if it's a gap
                # metrics are TP, FP, FN, and TN
                metrics = ",".join([ str(int(b)) for b in [
                    ((compare_bases(pp, gt)) & (not compare_bases(ll, gt))), 
                    ((not compare_bases(pp, gt)) & (compare_bases(ll, gt))),
                    ((not compare_bases(pp, gt)) & (not compare_bases(ll, gt))), 
                    ((compare_bases(pp, gt)) & (compare_bases(ll, gt))) ] ])
                metrics_f.write(subsample + "," + str(pos) + "," + metrics + "\n")
            pos=pos+1
    metrics_f.close()
msat_f.close()