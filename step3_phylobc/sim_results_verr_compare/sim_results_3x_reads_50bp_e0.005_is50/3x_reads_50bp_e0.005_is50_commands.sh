#!/bin/bash
./simulate_paired_end_reads_fastq.pl SARS_CoV_2_genome.fasta 3x_reads_50bp_e0.005_is50 50 0.005 50 1800
/Users/marniella/anaconda3/bin/minimap2 -t 8 -a -x sr SARS_CoV_2_Wuhan_ref_genome.fasta 3x_reads_50bp_e0.005_is50.simulated_50_read1.fq 3x_reads_50bp_e0.005_is50.simulated_50_read2.fq -o 3x_reads_50bp_e0.005_is50_sim.sam
/Users/marniella/anaconda3/bin/samtools view -bS 3x_reads_50bp_e0.005_is50_sim.sam > 3x_reads_50bp_e0.005_is50_sim.bam
/Users/marniella/anaconda3/bin/samtools sort 3x_reads_50bp_e0.005_is50_sim.bam -o 3x_reads_50bp_e0.005_is50_sim.sorted.bam
/Users/marniella/anaconda3/bin/samtools mpileup -f SARS_CoV_2_Wuhan_ref_genome.fasta 3x_reads_50bp_e0.005_is50_sim.sorted.bam > 3x_reads_50bp_e0.005_is50_sim_piledup.txt
