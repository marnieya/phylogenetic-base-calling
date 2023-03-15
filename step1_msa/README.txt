// explanation of the files in this directory
// author: monica 


step-by-step

1. random_seqs.fasta		500 sequences collected by Lenore from GISAID
2. SARS-CoV-2_500genomes.fasta	linearized fasta
3. SARS-CoV-2_501genomes.fasta	same as above but + our target genome
4. SARS-CoV-2_316genomes.fasta	excluding samples that have 50 or more N's 
5. SARS-CoV-2_317genomes.fasta	same as above but + our target genome
6. alignment317.fasta		alignment with MAFFT
7. msa317.fasta			same as above but linearized + all caps
-----
8. SARS-CoV-2_282genomes.fasta	excluding samples that have 10 or more N's
9. SARS-CoV-2_234genomes.fasta	same as above but + our target genome + Wuhan-Hu reference
10. alignment284.fasta		alignment with MAFFT
11. msa284.fasta		same as above but linearized + all caps + spaces removed + ambiguous bases replaced
12. msa283.phylmAln		same as above but converted using Lenore's perl script

files below the dashed line used in the next step for tronko
