#!/bin/bash

# given an SRA + URL, download the SRA file --> fastq --> align to reference genome --> get pileup

sra=$1
url=$2

mm_path="/Users/marniella/anaconda3/bin/minimap2"
st_path="/Users/marniella/anaconda3/bin/samtools"
fqd_path="/Users/marniella/Downloads/sratoolkit/bin/fasterq-dump"
bcf_path="/usr/local/Cellar/bcftools/1.17/bin/bcftools"

dir="/Users/marniella/research/nielsen_lab/raw_read_results/"
ref_path="/Users/marniella/research/nielsen_lab/phylogenetic-base-calling/SARS_CoV_2_Wuhan_ref_genome.fasta"

dir=${dir}${sra}"/"
mkdir ${dir}

wget --directory-prefix=${dir} ${url}
${fqd_path} --split-files --skip-technical --outdir ${dir}fastq ${sra}
${mm_path} -t 8 -a -x sr ${ref_path} ${dir}fastq/${sra}_1.fastq ${dir}fastq/${sra}_2.fastq -o ${dir}${sra}.sam
${st_path} view -bS ${dir}${sra}.sam > ${dir}${sra}.bam
${st_path} sort ${dir}${sra}.bam -o ${dir}${sra}.sorted.bam
${st_path} mpileup -f ${ref_path} ${dir}${sra}.sorted.bam > ${dir}${sra}_piledup.txt

# after getting the pileup, use bcftools to call variants and create a consensus sequence
${bcf_path} mpileup -Ou -f ${ref_path} ${dir}${sra}.sorted.bam | ${bcf_path} call -mv -Oz -o ${dir}${sra}_calls.vcf.gz
${bcf_path} index ${dir}${sra}_calls.vcf.gz
cat ${ref_path} | ${bcf_path} consensus ${dir}${sra}_calls.vcf.gz | seqkit seq -w 0 > ${dir}${sra}_consensus.fasta

# from previous part of the script, we have a consensus fasta for raw reads
# now we use the pileup from that same sra to run phylobc

phylobc_path="/Users/marniella/research/nielsen_lab/phylogenetic-base-calling/call_bases_from_pileup.py"
time ${phylobc_path} ${sra}_piledup.txt ${dir}

# then we concatenate the consensus fasta with the phylobc calls and run mafft
cat ${dir}${sra}_consensus.fasta >> ${dir}${sra}_calls.fasta
mafft --maxiterate 1000 --preservecase --globalpair ${dir}${sra}_calls.fasta > ${dir}${sra}_calls_msa.fasta

# from msa we transpose to make it easier to compute accuracy stats
sed '/^>/s/.*//' ${dir}${sra}_calls_msa.fasta | awk -f transpose.awk > ${dir}${sra}_calls_msat.txt
