#!/bin/bash

# given number of base pairs (bp), intended sequencing depth (depth), and sequencing error per site (e)
# use Lenore's script to simulate reads from OU061397.1 fasta & generate pileup file
# this version has a constant insert size
# 'verr' stands for 'variable error', example of errors per site is error_per_site_29884sites.txt
## for that file, error at each site was generated with R: sample(rbeta(29884,0.05,10))

bp=50
depth=$1
e='error_per_site_29884sites.txt'

n_reads="$(( (${depth}*30000)/${bp} ))"
prefix="${depth}x_reads_${bp}bp_verr_is50"
res_dir="sim_results_${prefix}"

mkdir ${res_dir}

mm_path="/Users/marniella/anaconda3/bin/minimap2"
st_path="/Users/marniella/anaconda3/bin/samtools"

echo "#!/bin/bash" > ${prefix}_commands.sh
echo "./simulate_paired_end_reads_fastq_error.pl SARS_CoV_2_genome.fasta ${e} ${prefix} ${bp} 50 ${n_reads}" >> ${prefix}_commands.sh
echo "${mm_path} -t 8 -a -x sr SARS_CoV_2_Wuhan_ref_genome.fasta ${e}.simulated_${prefix}_read1.fq ${e}.simulated_${prefix}_read2.fq -o ${prefix}_sim.sam" >> ${prefix}_commands.sh
echo "${st_path} view -bS ${prefix}_sim.sam > ${prefix}_sim.bam" >> ${prefix}_commands.sh
echo "${st_path} sort ${prefix}_sim.bam -o ${prefix}_sim.sorted.bam" >> ${prefix}_commands.sh
echo "${st_path} mpileup -f SARS_CoV_2_Wuhan_ref_genome.fasta ${prefix}_sim.sorted.bam > ${prefix}_sim_piledup.txt" >> ${prefix}_commands.sh

chmod +x ${prefix}_commands.sh
./${prefix}_commands.sh
mv ${prefix}_* ${res_dir} 
