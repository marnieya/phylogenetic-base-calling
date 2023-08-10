#!/bin/bash

# given an SRA + URL, download the SRA file --> fastq --> align to reference genome --> get pileup

sra=$1
url=$2

mm_path="/Users/marniella/anaconda3/bin/minimap2"
st_path="/Users/marniella/anaconda3/bin/samtools"
fqd_path="/Users/marniella/Downloads/sratoolkit/bin/fasterq-dump"
bcf_path="/usr/local/Cellar/bcftools/1.17/bin/bcftools"

dir="/Users/marniella/research/nielsen_lab/project-extra-files/phylobc/step3_phylobc/raw_read_results/"
ref_path="/Users/marniella/research/nielsen_lab/phylogenetic-base-calling/SARS_CoV_2_Wuhan_ref_genome.fasta"

dir=${dir}${sra}"_files/"
mkdir ${dir}

wget --directory-prefix=${dir} ${url}
${fqd_path} --split-files --skip-technical --outdir ${dir}fastq ${sra}

echo "********************** downsampling fastq **********************"
# downsample fastq
seqtk sample -s100 ${dir}fastq/${sra}_1.fastq 100 > ${dir}fastq/${sra}_sub100_1.fq
seqtk sample -s100 ${dir}fastq/${sra}_2.fastq 100 > ${dir}fastq/${sra}_sub100_2.fq

echo "******************** aligning to reference *********************"
${mm_path} -t 8 -a -x sr ${ref_path} ${dir}fastq/${sra}_1.fastq ${dir}fastq/${sra}_2.fastq -o ${dir}${sra}.sam
${mm_path} -t 8 -a -x sr ${ref_path} ${dir}fastq/${sra}_sub100_1.fq ${dir}fastq/${sra}_sub100_2.fq -o ${dir}${sra}_sub100.sam

rm -r ${dir}fastq 

${st_path} view -bS ${dir}${sra}.sam > ${dir}${sra}.bam
${st_path} view -bS ${dir}${sra}_sub100.sam > ${dir}${sra}_sub100.bam

rm ${dir}${sra}.sam
rm ${dir}${sra}_sub100.sam

${st_path} sort ${dir}${sra}.bam -o ${dir}${sra}.sorted.bam
${st_path} sort ${dir}${sra}_sub100.bam -o ${dir}${sra}_sub100.sorted.bam

rm ${dir}${sra}.bam
rm ${dir}${sra}_sub100.bam

echo "********************** generating pileup ***********************"
${st_path} mpileup -f ${ref_path} ${dir}${sra}_sub100.sorted.bam > ${dir}${sra}_sub100_piledup.txt

# after getting the pileup, use bcftools to call variants and create a consensus sequence
${bcf_path} mpileup -Ou -f ${ref_path} ${dir}${sra}.sorted.bam | ${bcf_path} call -mv -Oz -o ${dir}${sra}_calls.vcf.gz
${bcf_path} index ${dir}${sra}_calls.vcf.gz
cat ${ref_path} | ${bcf_path} consensus ${dir}${sra}_calls.vcf.gz | seqtk seq > ${dir}${sra}_consensus.fasta

rm ${dir}${sra}.sorted.bam
rm ${dir}${sra}_sub100.sorted.bam
rm ${dir}${sra}_calls.vcf.gz
