#!/bin/bash

# given a file with inputs (sra accession, url)
# run commands to download fastq

input_file=$1

fqd_path="/space/s1/marniella/sratoolkit/bin/fasterq-dump"
home_path="/space/s1/marniella/phylogenetic-base-calling/sra_raw_read_results/"

while IFS="," read -r sra url remainder; do
	dir=${home_path}${sra}"_files/"
	mkdir ${dir}
	wget --directory-prefix=${dir} ${url}
	${fqd_path} --split-files --skip-technical --outdir ${dir}fastq ${sra} &
done < ${input_file}
