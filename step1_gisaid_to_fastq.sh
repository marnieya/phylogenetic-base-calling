#!/bin/bash

# given a genome to simulate from (id)
# generate a random distribution of error rates per site, and use Lenore's script to simulate reads at varying depths
# this version has a constant insert size (50 bp) and a constant read length based on avgLength from SRA metadata (215 bp)

input_file=$1

home_path="/space/s1/marniella/phylogenetic-base-calling/"
simulate_reads="/space/s1/marniella/lenore-scripts/simulate_paired_end_reads_fastq_error.pl"

while IFS="," read -r id strain remainder; do
	dir="${home_path}gisaid_sim_read_results/${id}_files/"
	mkdir ${dir} 

	grep -A1 "${id}" ${home_path}"gisaid_sim_read_results/GISAID_Genomes_230814-230919.fasta" | seqtk seq | tr -s ' ' | tr ' ' '_' > ${dir}${id}.fasta
	sed -i 's/hCoV-19\/.*EPI/EPI/g' ${dir}${id}.fasta
	sed -i 's/|2023-.*-.*//g' ${dir}${id}.fasta
	genome_length=$(cat ${dir}${id}.fasta | awk 'NR==2 {print length}')

	# generate error rates files for the specific number of sites of this genome
	low_err="beta0.001"
	med_err="beta0.004"
	high_err="beta0.007"

	Rscript <(echo "invisible(lapply(rbeta(${genome_length}, 10, 9990) + 1e-04, write, \"${dir}${id}_err_${low_err}\", append=TRUE, ncolumns=1))")
	Rscript <(echo "invisible(lapply(rbeta(${genome_length}, 40, 9960) + 1e-04, write, \"${dir}${id}_err_${med_err}\", append=TRUE, ncolumns=1))")
	Rscript <(echo "invisible(lapply(rbeta(${genome_length}, 700, 99300) + 1e-04, write, \"${dir}${id}_err_${high_err}\", append=TRUE, ncolumns=1))")
	
	fastq_dir="${dir}fastq/"
	mkdir ${fastq_dir}

	for d in 5 4 3 2 1 0.8 0.5 0.2; do
		n_reads=$(bc <<< "scale=0; (${d}*${genome_length})/430")
		cmnd1="${simulate_reads} ${dir}${id}.fasta ${dir}${id}_err_${low_err} ${fastq_dir}${id}_err_${low_err}_${d}x_cov 215 50 ${n_reads}"
		cmnd2="${simulate_reads} ${dir}${id}.fasta ${dir}${id}_err_${med_err} ${fastq_dir}${id}_err_${med_err}_${d}x_cov 215 50 ${n_reads}"
		cmnd3="${simulate_reads} ${dir}${id}.fasta ${dir}${id}_err_${high_err} ${fastq_dir}${id}_err_${high_err}_${d}x_cov 215 50 ${n_reads}"
		${cmnd1} & ${cmnd2} & ${cmnd3}
		# mv ${dir}${id}*.fq ${fastq_dir}
		sed -i 's/_[1-2]$//g' ${fastq_dir}*.fq
	done
done < ${input_file}
