#!/bin/bash

# this script assumes runs have been downloaded and subsampled (or simulated), and consensus assemblies are done (after step0-2.sh)
# all of this will be in a directory titled "[id]_files" for each subsampled sequencing run, or simulated set of raw reads
# this identifier is the first input, the full path to the directory is the second input

id=$1
dir=$2  # e.g. for SRA, "sra_raw_read_results/"${id}"_files/"
	# for GISAID, "gisaid_sim_read_results/"${id}"_files/"

home_path="/space/s1/marniella/phylogenetic-base-calling/"
ref_path=${home_path}SARS_CoV_2_Wuhan_ref_genome.fasta
phylobc_path=${home_path}"scripts/call_bases.py"

# we use the pileups to call bases in the subsamples 
# we use subsample assemblies & ground truth from MSA done in previous step for indexing 
if [[ ${id} == SRR* ]]; then
	echo "sample,pos,n_reads,GT,SA,PR,LL,PP,PR_val,LL_val,PP_val" > ${dir}${id}_prob_calls_pr_scaling0.0_precise.txt
	echo "sample,pos,n_reads,GT,SA,PR,LL,PP,PR_val,LL_val,PP_val" > ${dir}${id}_prob_calls_pr_scaling0.001_precise.txt
	echo "sample,pos,n_reads,GT,SA,PR,LL,PP,PR_val,LL_val,PP_val" > ${dir}${id}_prob_calls_pr_scalingMLE_precise.txt
	${phylobc_path} ${id}_piledup.txt ${dir} 
 	for rep in {1..10}; do #d in 0.8 1 3 5 10 20
		${phylobc_path} ${id}_sub20_rep${rep}_piledup.txt ${dir} #${id}_sub${d}_piledup.txt ${dir}
	done
else
	for e in "beta0.001" "beta0.004" "beta0.007"; do
		echo "sample,pos,n_reads,GT,SA,PR,LL,PP,PR_val,LL_val,PP_val" > ${dir}${id}_err_${e}_prob_calls_pr_scaling0.001_precise.txt
		echo "sample,pos,n_reads,GT,SA,PR,LL,PP,PR_val,LL_val,PP_val" > ${dir}${id}_err_${e}_prob_calls_pr_scaling0.0_precise.txt
		echo "sample,pos,n_reads,GT,SA,PR,LL,PP,PR_val,LL_val,PP_val" > ${dir}${id}_err_${e}_prob_calls_pr_scalingMLE_precise.txt
		for d in 0.2 0.5 0.8 1 2 3 4 5; do
			${phylobc_path} ${id}_err_${e}.simulated_${d}x_cov_piledup.txt ${dir} 
		done	
	done
fi

### we have (1): MSA of assemblies with the (possibly gapped) reference genome from the global database
### this is necessary for matching up indexes of priors, to the pileup indexes for a particular assembly

### previous approach: add calls to (1), then MSA and transpose, then compare
# cat ${dir}${id}*prob_calls.fasta > ${dir}${id}_calls_concat.fasta
# mafft --thread 4 --maxiterate 1001 --preservecase --globalpair --keeplength --add ${dir}${id}_calls_concat.fasta ${dir}${id}_assemblies_msa.fasta > ${dir}${id}_calls_msa.fasta

# transpose_path=${home_path}"transpose.awk"
# seqtk seq ${dir}${id}_calls_msa.fasta | sed '/^>/s/.*//' | awk -f ${transpose_path} > ${dir}${id}_calls_msat.txt
# sed -n -e '/^>/p' ${dir}${id}_calls_msa.fasta | sed 's/>//g' | sed -z 's/\n/,/g' | sed 's/,$//g' > ${dir}${id}_tmp.txt
# echo "" >> ${dir}${id}_tmp.txt
# cat ${dir}${id}_calls_msat.txt >> ${dir}${id}_tmp.txt
# mv ${dir}${id}_tmp.txt ${dir}${id}_calls_msat.txt

### new approach: ground truth, reference, and assemblies are already matched up
# base calling script can still refer to the reference and its particular assembly
# and spit out the corresponding pr,ll,pp, plus ref and assembly, plus ground truth (for accuracy evaluation only)

# so now the base calling script does more -- creating a separate file for each subsample's calls, and appending to the overall accuracy metrics file


