#!/bin/bash
./simulate_paired_end_reads_fastq.pl SARS_CoV_2_genome.fasta 5x_reads_50bp_e0.005_is50 50 0.005 50 3000
/Users/marniella/anaconda3/bin/minimap2 -t 8 -a -x sr SARS_CoV_2_Wuhan_ref_genome.fasta 5x_reads_50bp_e0.005_is50.simulated_50_read1.fq 5x_reads_50bp_e0.005_is50.simulated_50_read2.fq -o 5x_reads_50bp_e0.005_is50_sim.sam
/Users/marniella/anaconda3/bin/samtools view -bS 5x_reads_50bp_e0.005_is50_sim.sam > 5x_reads_50bp_e0.005_is50_sim.bam
/Users/marniella/anaconda3/bin/samtools sort 5x_reads_50bp_e0.005_is50_sim.bam -o 5x_reads_50bp_e0.005_is50_sim.sorted.bam
/Users/marniella/anaconda3/bin/samtools mpileup -f SARS_CoV_2_Wuhan_ref_genome.fasta 5x_reads_50bp_e0.005_is50_sim.sorted.bam > 5x_reads_50bp_e0.005_is50_sim_piledup.txt
