#!/bin/bash

home_path="/space/s1/marniella/phylogenetic-base-calling/"
sra_runs_path=${home_path}"sra_raw_read_results/SRA_Illumina_Runs_230814-230919_random_sample_inputs.csv"
gisaid_genomes_path=${home_path}"gisaid_sim_read_results/GISAID_Illumina_Genomes_230814-230919_random_sample_inputs.csv"
ref_path=${home_path}"SARS_CoV_2_Wuhan_ref_genome.fasta"

fastq_to_ubam="/space/s1/marniella/viral-pipelines/pipes/WDL/workflows/fastq_to_ubam.wdl"
assemble_refbased="/space/s1/marniella/viral-pipelines/pipes/WDL/workflows/assemble_refbased.wdl"

# note, to keep miniwdl from burning down the server, use a modest N for assemblies (e.g. N=5) and:
export OMP_NUM_THREADS=1; export MINIWDL__SCHEDULER__TASK_CONCURRENCY=10

# download or simulate fastq
./step1_sra_to_fastq.sh ${sra_runs_path} 
./step1_gisaid_to_fastq.sh ${gisaid_genomes_path}

# convert these to ubam
for sra_dir in ${home_path}sra_raw_read_results/*_files; do
	sra_path=$(echo ${sra_dir} | sed 's/_files//g')
	sra=${sra_path##*/}
	./step2_fastq_to_ubam.sh ${sra_dir} ${sra} "_" "fastq" & 
done

N=20
for epi_dir in ${home_path}gisaid_sim_read_results/*_files; do
	(
	epi_path=$(echo ${epi_dir} | sed 's/_files//g')
	epi=${epi_path##*/}
	for d in 5 4 3 2 1 0.8 0.5 0.2; do
		for e in "beta0.007" "beta0.004" "beta0.001"; do
			./step2_fastq_to_ubam.sh ${epi_dir}/ "${epi}_err_${e}_${d}x_cov.simulated_215" "_read" "fq" 
		done
	done
	) &
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		wait -n
	fi
done 

# for all generated full depth ubam, run assembly pipeline
N=25
while IFS="," read -r sra ubam_path bed_file remainder; do
        (
        cmnd="miniwdl run ${assemble_refbased} reads_unmapped_bams=${ubam_path} reference_fasta=${ref_path} min_coverage=20 skip_mark_dupes=true trim_coords_bed=${bed_file}"
        ${cmnd}
        ) &
       if [[ $(jobs -r -p | wc -l) -ge $N ]]; then wait -n; fi
done < full_ubam_paths+primers.csv

rm -r ${home_path}*_fastq_to_ubam

# for all assemblies, move relevant output to SRA-specific folder
# use outputs to calculate downsample reads, and save those to same folder
# then re-run the ubam conversion and assembly for each downsample of SRA reads 
# (the simulated reads should already be done)

for assembly_dir in ${home_path}*assemble_refbased; do
	fasta_path=$(ls ${assembly_dir}/call-call_consensus/out/refined_assembly_fasta/*.fasta)
	id_path=$(echo ${fasta_path} | sed 's/.fasta//g')
        id=${id_path##*/}
	id_dir=${home_path}gisaid_sim_read_results/$(echo ${id} | sed 's/_err_beta.*//g')_files/
	if [[ ${id} == SRR* ]]
	then
		id_dir=${home_path}sra_raw_read_results/${id}_files/
	fi
	cp ${fasta_path} ${id_dir}
	cp $(ls ${assembly_dir}/call-plot_ref_coverage/out/coverage_tsv/*.txt) ${id_dir} 
	cp ${assembly_dir}/outputs.json ${id_dir}${id}_assembly_outputs.json
	bam_path=$(ls ${assembly_dir}/call-ivar_trim-0/out/aligned_trimmed_bam/*.bam)
	samtools mpileup -f ${ref_path} ${bam_path} > ${id_dir}${id}_piledup.txt
	if [[ ${id} == SRR* ]]
	then
		for d in 20 15 10 5 4 3 2 1 0.8 0.5 0.2; do
			n_reads=`python get_subsample_n_reads.py ${d} ${id_dir}/${id}_assembly_outputs.json`
			seqtk sample -s1000 ${id_dir}fastq/${id}_1.fastq ${n_reads} > ${id_dir}fastq/${id}_sub${d}_1.fq
			seqtk sample -s1000 ${id_dir}fastq/${id}_2.fastq ${n_reads} > ${id_dir}fastq/${id}_sub${d}_2.fq
			./step2_fastq_to_ubam.sh ${id_dir} ${id}_sub${d} "_" "fq"
		done &
	fi
done

# run subsample assemblies with paths and their corresponding primer .bed files:
while IFS="," read -r sra ubam_path bed_file remainder; do
	(
	cmnd="miniwdl run ${assemble_refbased} reads_unmapped_bams=${ubam_path} reference_fasta=${ref_path} min_coverage=3 skip_mark_dupes=true trim_coords_bed=${bed_file}"
	${cmnd}
	) &
       if [[ $(jobs -r -p | wc -l) -ge $N ]]; then wait -n; fi
done < subsample_ubam_paths+primers.csv  


# now, all assemblies (subsamples only), move outputs etc. 
for assembly_dir in ${home_path}*assemble_refbased; do
        fasta_path=$(ls ${assembly_dir}/out/assembly_fasta/*.fasta)
        sra_path=$(echo ${fasta_path} | sed 's/.fasta//g')
        sra_w_sub=${sra_path##*/}
	sra=$(echo ${sra_w_sub} | sed 's/_sub.*//g')
	sra_dir=${home_path}sra_raw_read_results/${sra}_files/
        cp ${fasta_path} ${sra_dir}
        cp $(ls ${assembly_dir}/out/align_to_ref_merged_coverage_tsv/*.txt) ${sra_dir}
	cp ${assembly_dir}/outputs.json ${sra_dir}/${sra_w_sub}_assembly_outputs.json 
	bam_path=$(ls ${assembly_dir}/out/align_to_ref_merged_aligned_trimmed_only_bam/*.bam)
        cp ${bam_path} ${sra_dir}
	samtools mpileup -f ${ref_path} ${bam_path} > ${sra_dir}${sra_w_sub}_piledup.txt 
done

# at this step, all assemblies are done; place resulting assembled genomes using tronko
# correct tronko priors with error_estimate.c to obtain MLE scaling
# lastly, go through each folder, and call bases for each SRA and its subsets
# (base calling script automatically uses unscaled, constant scaled, and MLE scaled priors)

# lastly, go through each folder, and call bases for each SRA and its subsets
N=5
for sra_dir in ${home_path}sra_raw_read_results/*_files; do
	(
	sra_path=$(echo ${sra_dir} | sed 's/_files//g')
	sra=${sra_path##*/}
	# formatting files containing the priors
	for pp in ${sra_dir}/tronko_PP/*PP.txt; do
		echo "position        PP_A    PP_C    PP_G    PP_T" >> ${sra_dir}/tronko_PP/pp_temp.txt
		nl -w4 -s $'\t' ${pp} >> ${sra_dir}/tronko_PP/pp_temp.txt
	      	mv ${sra_dir}/tronko_PP/pp_temp.txt ${pp}
       	done
	# aligning subsampled assemblies with the global MSA reference genome
	cat ${sra_dir}/${sra}.fasta > ${sra_dir}/${sra}_pwise.fasta
	cat ${home_path}Wuhan_Hu_reference_MSA.fasta >> ${sra_dir}/${sra}_pwise.fasta
	mafft --thread 4 --maxiterate 1001 --preservecase --globalpair ${sra_dir}/${sra}_pwise.fasta | seqtk seq > ${sra_dir}/${sra}_pwise_msa.fasta
       	for sub in 0.2 0.5 0.8 1 2 3 4 5 10 15 20; do
	(
		cat ${sra_dir}/${sra}.fasta > ${sra_dir}/${sra}_sub${sub}_gt_ref.fasta
		cat ${home_path}Wuhan_Hu_reference_MSA.fasta >> ${sra_dir}/${sra}_sub${sub}_gt_ref.fasta
		cat ${sra_dir}/${sra}_sub${sub}.fasta >> ${sra_dir}/${sra}_sub${sub}_gt_ref.fasta
		mafft --thread 4 --maxiterate 1001 --preservecase --globalpair ${sra_dir}/${sra}_sub${sub}_gt_ref.fasta | seqtk seq > ${sra_dir}/${sra}_sub${sub}_gt_ref_msa.fasta
	) &
       	# allow to execute up to $N jobs in parallel
        if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
                # now there are $N jobs already running, so wait here for any job
                # to be finished so there is a place to start next one.
                wait -n
        fi
	./step3_call_bases_from_pileup.sh ${sra} ${sra_dir}/ 
done

N=5
for epi_dir in ${home_path}gisaid_sim_read_results/*_files; do 
	epi_path=$( echo ${epi_dir} | sed 's/_files//g' )
	epi=${epi_path##*/}
	for sim_genome in ${epi_dir}/*_215.fasta; do 
		sim_genome_id=$( echo ${sim_genome##*/} | sed 's/.fasta//g' )
		cat ${epi_dir}/${epi}.fasta | seqtk seq > ${epi_dir}/${sim_genome_id}_gt_ref.fasta
		cat Wuhan_Hu_reference_MSA.fasta | seqtk seq >> ${epi_dir}/${sim_genome_id}_gt_ref.fasta
		cat ${epi_dir}/${sim_genome_id}.fasta | seqtk seq >> ${epi_dir}/${sim_genome_id}_gt_ref.fasta
	done

	for epi_msa in ${epi_dir}/*_gt_ref.fasta; do 
		(
		epi_msa_id=$( echo ${epi_msa##*/} | sed 's/_gt_ref.fasta//g' )
		mafft --thread 4 --maxiterate 1001 --preservecase --globalpair ${epi_msa} | seqtk seq > ${epi_dir}/${epi_msa_id}_gt_ref_msa.fasta 
		) & 
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then 
			wait -n 
		fi
	done

	( 
	./step3_call_bases_from_pileup.sh ${epi} ${epi_dir}/ 
	) & 
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		wait -n
	fi
done

# for the sub20 replicates analysis, start with the full-depth fastq's
# and instead of different depths, subsample with different seeds
N=20
for sra_dir in ${home_path}sra_raw_read_results/*_files; do
        (
        sra_path=$(echo ${sra_dir} | sed 's/_files//g')
        sra=${sra_path##*/}
	for rep in {1..10}; do
		n_reads=`python get_subsample_n_reads.py 20 ${sra_dir}/${sra}_assembly_outputs.json`
		seqtk sample -s${rep} ${sra_dir}/fastq/${sra}_1.fastq ${n_reads} > ${home_path}sra_raw_read_results_sub20/${sra}_files/fastq/${sra}_sub20_rep${rep}_1.fq
		seqtk sample -s${rep} ${sra_dir}/fastq/${sra}_2.fastq ${n_reads} > ${home_path}sra_raw_read_results_sub20/${sra}_files/fastq/${sra}_sub20_rep${rep}_2.fq
	done
) & if [[ $(jobs -r -p | wc -l) -ge $N ]]; then wait -n; fi 
done

# all replicates to ubam
for sra_dir in ${home_path}sra_raw_read_results_sub20/*_files; do
	(
	sra_path=$(echo ${sra_dir} | sed 's/_files//g')
        sra=${sra_path##*/}
	for rep in {1..10}; do
		./step2_fastq_to_ubam.sh ${sra_dir}/ ${sra} "_sub20_rep${rep}_" "fq"
	done	
	) & if [[ $(jobs -r -p | wc -l) -ge $N ]]; then wait -n; fi
done

# now run the replicate assemblies like any subsample, with corresponding primers
while IFS="," read -r sra ubam_path bed_file remainder; do
        (
        cmnd="miniwdl run ${assemble_refbased} reads_unmapped_bams=${ubam_path} reference_fasta=${ref_path} min_coverage=3 skip_mark_dupes=true trim_coords_bed=${bed_file}"
        ${cmnd}
        ) &
       if [[ $(jobs -r -p | wc -l) -ge $N ]]; then wait -n; fi
done < sub20_ubam_paths+primers.csv

# move outputs, call bases
for assembly_dir in ${home_path}*assemble_refbased; do
        fasta_path=$(ls ${assembly_dir}/out/assembly_fasta/*.fasta)
        sra_path=$(echo ${fasta_path} | sed 's/.fasta//g')
        sra_w_sub=${sra_path##*/}
        sra=$(echo ${sra_w_sub} | sed 's/_sub.*//g')
        sra_dir=${home_path}sra_raw_read_results_sub20/${sra}_files/
        cp ${fasta_path} ${sra_dir}
        cp $(ls ${assembly_dir}/out/align_to_ref_merged_coverage_tsv/*.txt) ${sra_dir}
        cp ${assembly_dir}/outputs.json ${sra_dir}/${sra_w_sub}_assembly_outputs.json
        bam_path=$(ls ${assembly_dir}/out/align_to_ref_merged_aligned_trimmed_only_bam/*.bam)
        cp ${bam_path} ${sra_dir}
        samtools mpileup -f ${ref_path} ${bam_path} > ${sra_dir}${sra_w_sub}_piledup.txt
done

for sra_dir in ${home_path}sra_raw_read_results_sub20/*_files; do
        (
        sra_path=$(echo ${sra_dir} | sed 's/_files//g')
        sra=${sra_path##*/}
        for rep in {1..10}; do
                cat ${home_path}sra_raw_read_results_bed_dupest/${sra}_files/${sra}.fasta > ${sra_dir}/${sra}_sub20_rep${rep}_gt_ref.fasta
                cat ${home_path}Wuhan_Hu_reference_MSA.fasta >> ${sra_dir}/${sra}_sub20_rep${rep}_gt_ref.fasta
                cat ${sra_dir}/${sra}_sub20_rep${rep}.fasta >> ${sra_dir}/${sra}_sub20_rep${rep}_gt_ref.fasta
                mafft --thread 4 --maxiterate 1001 --preservecase --globalpair ${sra_dir}/${sra}_sub20_rep${rep}_gt_ref.fasta | seqtk seq > ${sra_dir}/${sra}_sub20_rep${rep}_gt_ref_msa.fasta	
	done
	./step3_call_bases_from_pileup.sh ${sra} ${sra_dir}/
	) & if [[ $(jobs -r -p | wc -l) -ge $N ]]; then wait -n; fi
done

# summarize the accuracies and plot
python plot_main_figures.py # --rebuild-real / --rebuild-sims / --rebuild-reps to update the cache 
python plot_supp_figures.py
./aggregate_coverage.sh real sra_raw_read_results_bed_dupest
./aggregate_coverage.sh reps sra_raw_read_results_sub20
./aggregate_coverage.sh sim  gisaid_sim_read_results

