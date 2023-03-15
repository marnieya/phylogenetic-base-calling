#!/bin/bash
./simulate_paired_end_reads_fastq_error.pl SARS_CoV_2_genome.fasta error_per_site_29884sites.txt 4x_reads_50bp_verr_is50 50 50 2400
/Users/marniella/anaconda3/bin/minimap2 -t 8 -a -x sr SARS_CoV_2_Wuhan_ref_genome.fasta error_per_site_29884sites.txt.simulated_4x_reads_50bp_verr_is50_read1.fq error_per_site_29884sites.txt.simulated_4x_reads_50bp_verr_is50_read2.fq -o 4x_reads_50bp_verr_is50_sim.sam
/Users/marniella/anaconda3/bin/samtools view -bS 4x_reads_50bp_verr_is50_sim.sam > 4x_reads_50bp_verr_is50_sim.bam
/Users/marniella/anaconda3/bin/samtools sort 4x_reads_50bp_verr_is50_sim.bam -o 4x_reads_50bp_verr_is50_sim.sorted.bam
/Users/marniella/anaconda3/bin/samtools mpileup -f SARS_CoV_2_Wuhan_ref_genome.fasta 4x_reads_50bp_verr_is50_sim.sorted.bam > 4x_reads_50bp_verr_is50_sim_piledup.txt
