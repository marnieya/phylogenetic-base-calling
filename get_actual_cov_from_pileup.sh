#!/bin/bash

# this script calculates actual mean coverage from a pileup file

sim_folder=$1
d=$2
bp=50
e=$3 # can be 'verr' or 'e0.005'

sim_nreads_sum=$( cut -f 4 ${sim_folder}/sim_results_${d}x_reads_${bp}bp_${e}_is50/${d}x_reads_${bp}bp_${e}_is50_sim_piledup.txt | awk '{s+=$1}END{print s}' )
sim_nreads_avg=$( echo "scale=2 ; ${sim_nreads_sum} / 29884" | bc ); 
echo ${d} ${e} ${sim_nreads_avg}
