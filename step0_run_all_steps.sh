#!/bin/bash

home_path="/space/s1/marniella/phylogenetic-base-calling/"
sra_runs_path=${home_path}"sra_raw_read_results/SRA_Illumina_Runs_230814-230919_random_sample_inputs.csv"
gisaid_genomes_path=${home_path}"gisaid_sim_read_results/GISAID_Illumina_Genomes_230814-230919_random_sample_inputs.csv"
ref_path=${home_path}"SARS_CoV_2_Wuhan_ref_genome.fasta"

fastq_to_ubam="/space/s1/marniella/viral-pipelines/pipes/WDL/workflows/fastq_to_ubam.wdl"
assemble_refbased="/space/s1/marniella/viral-pipelines/pipes/WDL/workflows/assemble_refbased.wdl"

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

# for all generated ubam, run assembly pipeline
N=25
for ubam_dir in ${home_path}*_fastq_to_ubam; do
        (
	ubam_path=$(ls ${ubam_dir}/out/unmapped_bam/*.bam)
	id_path=$(echo ${ubam_path} | sed 's/.bam//g')
        id=${id_path##*/}
	align="bwa"
	if [[ ${id} == SRR* ]]
	then
		align="minimap2"
	fi
	cmnd="miniwdl run ${assemble_refbased} reads_unmapped_bams=${ubam_path} reference_fasta=${ref_path} aligner=${align} min_coverage=0"
	${cmnd} 
	) & 
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		wait -n
	fi
done

rm -r ${home_path}*_fastq_to_ubam

# for all assemblies, move relevant output to SRA-specific folder
# use outputs to calculate downsample reads, and save those to same folder
# then re-run the ubam conversion and assembly for each downsample of SRA reads 
# (the simulated reads should already be done)

for assembly_dir in ${home_path}*assemble_refbased; do
	fasta_path=$(ls ${assembly_dir}/out/assembly_fasta/*.fasta)
	id_path=$(echo ${fasta_path} | sed 's/.fasta//g')
        id=${id_path##*/}
	id_dir=${home_path}gisaid_sim_read_results/$(echo ${id} | sed 's/_err_beta.*//g')_files/
	if [[ ${id} == SRR* ]]
	then
		id_dir=${home_path}sra_raw_read_results/${id}_files/
	fi
	cp ${fasta_path} ${id_dir}
	cp $(ls ${assembly_dir}/out/align_to_ref_merged_coverage_tsv/*.txt) ${id_dir} 
	cp ${assembly_dir}/outputs.json ${id_dir}${id}_assembly_outputs.json
	bam_path=$(ls ${assembly_dir}/out/align_to_ref_merged_aligned_trimmed_only_bam/*.bam)
	samtools mpileup -f ${ref_path} ${bam_path} > ${id_dir}${id}_piledup.txt
	if [[ ${id} == SRR* ]]
	then
		for d in 20 15 10 5 4 3 2 1 0.8 0.5 0.2; do
			n_reads=`python ./get_subsample_n_reads.py ${d} ${assembly_dir}/outputs.json`
			# random_seed=`python ./get_subsample_random_seed.py ${id}`
			seqtk sample -s1000 ${id_dir}fastq/${id}_1.fastq ${n_reads} > ${id_dir}fastq/${id}_sub${d}_1.fq
			seqtk sample -s1000 ${id_dir}fastq/${id}_2.fastq ${n_reads} > ${id_dir}fastq/${id}_sub${d}_2.fq
			./step2_fastq_to_ubam.sh ${id_dir} ${id}_sub${d} "_" "fq"
		done &
	fi
done

rm -r ${home_path}*assemble_refbased

# now, for all generated ubam (subsamples only), do same as above
# parallelized 
N=25
for ubam_dir in ${home_path}*_fastq_to_ubam; do
        (
        # parallel code here
        ubam_path=$(ls ${ubam_dir}/out/unmapped_bam/*.bam)
        cmnd="miniwdl run ${assemble_refbased} reads_unmapped_bams=${ubam_path} reference_fasta=${ref_path} min_coverage=0"
        ${cmnd}
        ) &
        # allow to execute up to $N jobs in parallel
        if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
                # now there are $N jobs already running, so wait here for any job
                # to be finished so there is a place to start next one.
                wait -n
        fi
done
wait

rm -r ${home_path}*fastq_to_ubam

# now, all assemblies (subsamples only), do same as above
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

rm -r ${home_path}*assemble_refbased
rm _LAST

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
       	for sub in 0.2 0.5 0.8 1 2 3 4 5 10 15 20; do
	(
		cat ${sra_dir}/${sra} > ${sra_dir}/${sra}_sub${sub}_gt_ref.fasta
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
	./step3_call_bases_from_pileup.sh ${sra} ${sra_dir}/ "normal"
	./step3_call_bases_from_pileup.sh ${sra} ${sra_dir}/ "roc"
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
	./step3_call_bases_from_pileup.sh ${epi} ${epi_dir}/ "normal" 
	./step3_call_bases_from_pileup.sh ${epi} ${epi_dir}/ "roc" 
	) & 
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		wait -n
	fi
done

# merge the outputs 
./step4_coverage_plots.sh ${home_path}"sra_raw_read_results/"

# summarize the accuracies
./compare_full_pileup_calls.py "${home_path}sra_raw_read_results/" "normal"
./compare_full_pileup_calls.py "${home_path}sra_raw_read_results/" "roc"
./compare_full_pileup_calls.py "${home_path}sra_raw_read_results/" "roc_overgenes"

#for thresh in 0.1 0.01 0.001 0.0001; do
#	all_calls="${home_path}sra_raw_read_results_min_cov_0/all_prob_calls_pr_scaling${thresh}.txt"
#	echo "id,method,accurate,called,positions" > "${home_path}sra_raw_read_results_min_cov_0/all_prob_calls_summarized_by_sample_pr_scaling${thresh}.txt"
#	echo "id,method,accurate,called,positions" > "${home_path}sra_raw_read_results_min_cov_0/all_prob_calls_summarized_by_pos_pr_scaling${thresh}.txt"
#	./compare_called_bases.py "${all_calls}" "normal"
#done

# generate ROC curves
roc_header="sample,pos,n_reads,pp_thresh_pass_0.999999999999,pp_thresh_pass_0.999999999995,pp_thresh_pass_0.99999999999,pp_thresh_pass_0.99999999995,pp_thresh_pass_0.9999999999,pp_thresh_pass_0.9999999995,pp_thresh_pass_0.999999999,pp_thresh_pass_0.999999995,pp_thresh_pass_0.99999999,pp_thresh_pass_0.99999995,pp_thresh_pass_0.9999999,pp_thresh_pass_0.9999995,pp_thresh_pass_0.999999,pp_thresh_pass_0.999995,pp_thresh_pass_0.99999,pp_thresh_pass_0.9999,pp_thresh_pass_0.999,pp_thresh_pass_0.99,pp_thresh_pass_0.9,pp_thresh_pass_0.85,pp_thresh_pass_0.8,pp_thresh_pass_0.75,pp_thresh_pass_0.7,pp_thresh_pass_0.65,pp_thresh_pass_0.6,pp_thresh_pass_0.55,pp_thresh_pass_0.5,GT,SA,PR,LL,PP"

for sra_dir in ${home_path}sra_raw_read_results_min_cov_0/*_files; do 
	( 
	sra_path=$(echo ${sra_dir} | sed 's/_files//g')
	sra=${sra_path##*/}
	rm ${sra_dir}/*_ROC.txt
	for pileupfile in ${sra_dir}/*_piledup.txt; do 
		p=${pileupfile##*/}
		./call_bases_from_pileup.py "${p}" "${sra_dir}/" "roc"
	done
	) & if [[ $(jobs -r -p | wc -l) -ge $N ]]; then wait -n; fi 
done

echo ${roc_header} > ${home_path}sra_raw_read_results_min_cov_0/all_prob_calls_pr_scaling0.001_ROC.txt
for sra_dir in ${home_path}sra_raw_read_results_min_cov_0/*_files; do
	sra_path=$(echo ${sra_dir} | sed 's/_files//g')
	sra=${sra_path##*/}
	cat ${sra_dir}/${sra}*_ROC.txt >> ${home_path}sra_raw_read_results_min_cov_0/all_prob_calls_pr_scaling0.001_ROC.txt
done

# TODO: add a header before running this
./compare_called_bases.py "${home_path}sra_raw_read_results_min_cov_0/all_prob_calls_pr_scaling0.001_ROC.txt" "roc"

# Rscript step6_plots.R

# done??


