#!/bin/bash

# this script assumes runs have been downloaded, subsampled, and consensus calls are done (after from_sra_to_pileups.sh)
# all of this will be in a directory titled "SRRXXXXX_files" for each run downloaded from SRA
# this identifier in the form of SRRXXXXX is the first input

sra=$1

dir="/Users/marniella/research/nielsen_lab/project-extra-files/phylobc/step3_phylobc/raw_read_results/"
ref_path="/Users/marniella/research/nielsen_lab/phylogenetic-base-calling/SARS_CoV_2_Wuhan_ref_genome.fasta"

dir=${dir}"ex5_after814/"${sra}"_files/"

# this part of the script uses the downsampled pileups for phylogenetic base calling
# then compares to the consensus sequence obtained from the full set of reads ("ground truth")

phylobc_path="/Users/marniella/research/nielsen_lab/phylogenetic-base-calling/call_bases_from_pileup.py"

{ time ${phylobc_path} ${sra}_piledup.txt ${dir} ; } 2>> phylobc_times.txt
for d in 0.2 0.5 0.8 1 2 3 4 5 10 15 20; do
	{ time ${phylobc_path} ${sra}_sub${d}_piledup.txt ${dir} ; } 2>> phylobc_times.txt
done

# then we concatenate the consensus fasta with the phylobc calls and run mafft
cat ${dir}${sra}*.fasta >> ${dir}${sra}_calls.fasta
mafft --maxiterate 1000 --preservecase --globalpair ${dir}${sra}_calls.fasta | seqtk seq > ${dir}${sra}_calls_msa.fasta

# rm ${dir}${sra}*.fasta

# from msa we transpose to make it easier to compute accuracy stats
sed '/^>/s/.*//' ${dir}${sra}_calls_msa.fasta | awk -f transpose.awk > ${dir}${sra}_calls_msat.txt
