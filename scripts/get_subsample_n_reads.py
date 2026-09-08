#!/home/marniella/miniconda3/bin/python

import sys
import json

read_depth = sys.argv[1]
outputs_path = sys.argv[2]

with open(outputs_path, 'r') as f:
    outputs = json.load(f)

assembly_len = outputs['assemble_refbased.assembly_length']
avg_read_len = float(outputs['assemble_refbased.align_to_ref_merged_bases_aligned']) / float(outputs['assemble_refbased.align_to_ref_merged_reads_aligned'])

n_reads = ( float(read_depth) * int(assembly_len) ) / avg_read_len
print(int(n_reads))
