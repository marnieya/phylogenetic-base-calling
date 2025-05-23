#!/bin/bash

home_path="/space/s1/marniella/phylogenetic-base-calling/"

# for SRA
# reads_dir=${home_path}"sra_raw_read_results_min_cov_0/"
# echo -e 'Run\tpos\tavg_depth_20x\tavg_depth_15x\tavg_depth_10x\tavg_depth_5x\tavg_depth_4x\tavg_depth_3x\tavg_depth_2x\tavg_depth_1x\tavg_depth_0.8x\tavg_depth_0.5x\tavg_depth_0.2x' > ${reads_dir}all_coverage_plots.txt

reads_dir=$1
echo -e 'Run\tpos\tavg_depth_5x\tavg_depth_4x\tavg_depth_3x\tavg_depth_2x\tavg_depth_1x\tavg_depth_0.8x\tavg_depth_0.5x\tavg_depth_0.2x' > ${reads_dir}all_coverage_plots.txt

for dir in ${reads_dir}*_files; do
       sra_path=$(echo ${dir} | sed 's/_files//g')
       sra=${sra_path##*/}
       if [[ ${sra} == SRR* ]]; then
	       paste -d"\t" <(cut -f2,3 ${dir}/${sra}_sub20.coverage_plot.txt) <(cut -f3 ${dir}/${sra}_sub15.coverage_plot.txt) \
		       <(cut -f3 ${dir}/${sra}_sub10.coverage_plot.txt) \
		       <(cut -f3 ${dir}/${sra}_sub5.coverage_plot.txt) \
		       <(cut -f3 ${dir}/${sra}_sub4.coverage_plot.txt) \
		       <(cut -f3 ${dir}/${sra}_sub3.coverage_plot.txt) \
		       <(cut -f3 ${dir}/${sra}_sub2.coverage_plot.txt) \
		       <(cut -f3 ${dir}/${sra}_sub1.coverage_plot.txt) \
		       <(cut -f3 ${dir}/${sra}_sub0.8.coverage_plot.txt) \
		       <(cut -f3 ${dir}/${sra}_sub0.5.coverage_plot.txt) \
		       <(cut -f3 ${dir}/${sra}_sub0.2.coverage_plot.txt) > ${reads_dir}${sra}_tmp.txt
	else
		for err in "beta0.001" "beta0.004" "beta0.007"; do
			id="${sra}_err_${err}"
			paste -d"\t" <(cut -f2,3 ${dir}/${id}.simulated_5x_cov.coverage_plot.txt) <(cut -f3 ${dir}/${id}.simulated_4x_cov.coverage_plot.txt) \
				<(cut -f3 ${dir}/${id}.simulated_3x_cov.coverage_plot.txt) \
				<(cut -f3 ${dir}/${id}.simulated_2x_cov.coverage_plot.txt) \
				<(cut -f3 ${dir}/${id}.simulated_1x_cov.coverage_plot.txt) \
				<(cut -f3 ${dir}/${id}.simulated_0.8x_cov.coverage_plot.txt) \
				<(cut -f3 ${dir}/${id}.simulated_0.5x_cov.coverage_plot.txt) \
				<(cut -f3 ${dir}/${id}.simulated_0.2x_cov.coverage_plot.txt) > ${reads_dir}${id}_tmp.txt 
		done
       fi
done

awk -i inplace '{print FILENAME"\t"$0}' ${reads_dir}*_tmp.txt
cat ${reads_dir}*_tmp.txt | sed 's/_tmp.txt//g' >> ${reads_dir}all_coverage_plots.txt
rm ${reads_dir}*_tmp.txt
