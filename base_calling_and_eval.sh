#!/bin/bash

# first run get_pileup_from_sra.sh, wherein: 
# given an SRA + URL, download the SRA file --> fastq --> align to reference genome --> consensus fasta
#						      --> downsample fastq --> align to reference --> get pileups

sra=$1
url=$2

mm_path="/Users/marniella/anaconda3/bin/minimap2"
st_path="/Users/marniella/anaconda3/bin/samtools"
fqd_path="/Users/marniella/Downloads/sratoolkit/bin/fasterq-dump"
bcf_path="/usr/local/Cellar/bcftools/1.17/bin/bcftools"

dir="/Users/marniella/research/nielsen_lab/project-extra-files/phylobc/step3_phylobc/raw_read_results/"
ref_path="/Users/marniella/research/nielsen_lab/phylogenetic-base-calling/SARS_CoV_2_Wuhan_ref_genome.fasta"

dir=${dir}${sra}"_files/"

# second, this part uses the downsampled pileups for phylogenetic base calling
# then compares to the consensus sequence obtained from the full set of reads

echo "**************** phylogenetic base calling *********************"
# from previous part of the script, we have a consensus fasta for raw reads
# now we use the pileup from that same sra's downsampled fastq/bam to run phylobc
phylobc_path="/Users/marniella/research/nielsen_lab/phylogenetic-base-calling/call_bases_from_pileup.py"
{ time ${phylobc_path} ${sra}_sub100_piledup.txt ${dir} ; } 2>> phylobc_times.txt

echo "******************** evaluating accuracies *********************"
# then we concatenate the consensus fasta with the phylobc calls and run mafft
cat ${dir}${sra}_consensus.fasta >> ${dir}${sra}_calls.fasta
mafft --maxiterate 1000 --preservecase --globalpair ${dir}${sra}_calls.fasta | seqtk seq > ${dir}${sra}_calls_msa.fasta

rm ${dir}${sra}_consensus.fasta

# from msa we transpose to make it easier to compute accuracy stats
sed '/^>/s/.*//' ${dir}${sra}_calls_msa.fasta | awk -f transpose.awk > ${dir}${sra}_calls_msat.txt
