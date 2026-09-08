#!/bin/bash

# given a path where we can find a "fastq" folder, 
# as well as the SRA "prefix" (i.e. SRRXXXXX or SRRXXXXX_subYY) 
# convert to ubam 

dir=$1
sra=$2
sep=$3
ext=$4

echo $0

miniwdl run \
  /space/s1/marniella/viral-pipelines/pipes/WDL/workflows/fastq_to_ubam.wdl \
  FastqToUBAM.fastq_1="${dir}fastq/${sra}${sep}1.${ext}" \
  FastqToUBAM.fastq_2="${dir}fastq/${sra}${sep}2.${ext}" \
  FastqToUBAM.library_name=${sra} \
  FastqToUBAM.platform_name="ILLUMINA" \
  FastqToUBAM.sample_name=${sra} \
  FastqToUBAM.additional_picard_options="QUALITY_FORMAT=Standard"
