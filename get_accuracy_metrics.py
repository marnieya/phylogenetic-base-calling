#!/Users/marniella/anaconda3/envs/py310/bin/python
import sys
dir = sys.argv[1]
sra = sys.argv[2]

pos=1
with open(dir + sra + "_calls_msat.txt","r") as msat_f:
    h = msat_f.readline().strip('\n').split(',')
    with open(dir + sra + "_calls_metrics.txt", "a") as metrics_f:
        for line in msat_f:
            # each line is a position in the genome
            all_calls = list(line)
            gt = all_calls[h.index(sra)]
            # we do not include positions where the ground truth is an N (because we are unable to assess accuracy)
            if (gt == 'N'):
                pos=pos+1
                continue 
            # if not, let's break down call accuracy for this position
            for subsample in [(sra + "_sub" + str(s)) for s in [0.2, 0.5, 0.8, 1, 2, 3, 4, 5, 10, 15, 20]]:
                subsample_assembly =  all_calls[h.index(subsample)]
                pp = all_calls[h.index(subsample + "_postprob_calls")]
                ll = all_calls[h.index(subsample + "_likelihood_calls")]
                pr = all_calls[h.index(subsample + "_priorprob_calls")]
                
                # if the assembly was N, we are interested on whether our method called it correctly
                if (subsample_assembly == 'N') or (subsample_assembly == '-'):
                    tp = 0
                    fp = 0
                    fn = 0
                    tn = 0
                else:
                    # if the assembly was not N, errors go into our gain calculation (% errors removed from the typical assembly)
                    tp = (pp == gt) & (subsample_assembly != gt)
                    fp = (pp != gt) & (subsample_assembly == gt)
                    fn = (pp != gt) & (subsample_assembly != gt)
                    tn = (pp == gt) & (subsample_assembly == gt)
                # calls are ground truth call, subsample assembly call, prior call, likelihood call, posterior call
                # metrics are TP, FP, FN, and TN
                calls = ",".join([gt, subsample_assembly, pr, ll, pp])
                metrics = ",".join([tp, fp, fn, tn])
                metrics_f.write(subsample + "," + str(pos) + "," + calls + "," + metrics + "\n")
            pos=pos+1
    metrics_f.close()
msat_f.close()