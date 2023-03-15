#!/bin/bash

# given number of base pairs (bp), error rate (e), and intended sequencing depth (depth)
# use Lenore's script to simulate reads from OU061397.1 fasta & generate pileup file

bp=50
e=0.005
i=50
depth=$1

n_reads="$(( (${depth}*30000)/${bp} ))"
prefix="${depth}x_reads_${bp}bp_e${e}_is${i}"
res_dir="sim_results_${prefix}"

mkdir ${res_dir}

mm_path="/Users/marniella/anaconda3/bin/minimap2"
st_path="/Users/marniella/anaconda3/bin/samtools"

echo "#!/bin/bash" > ${prefix}_commands.sh
echo "./simulate_paired_end_reads_fastq.pl SARS_CoV_2_genome.fasta ${prefix} ${bp} ${e} ${i} ${n_reads}" >> ${prefix}_commands.sh
echo "${mm_path} -t 8 -a -x sr SARS_CoV_2_Wuhan_ref_genome.fasta ${prefix}.simulated_${bp}_read1.fq ${prefix}.simulated_${bp}_read2.fq -o ${prefix}_sim.sam" >> ${prefix}_commands.sh
echo "${st_path} view -bS ${prefix}_sim.sam > ${prefix}_sim.bam" >> ${prefix}_commands.sh
echo "${st_path} sort ${prefix}_sim.bam -o ${prefix}_sim.sorted.bam" >> ${prefix}_commands.sh
echo "${st_path} mpileup -f SARS_CoV_2_Wuhan_ref_genome.fasta ${prefix}_sim.sorted.bam > ${prefix}_sim_piledup.txt" >> ${prefix}_commands.sh

chmod +x ${prefix}_commands.sh
./${prefix}_commands.sh
mv ${prefix}_* ${res_dir} 
