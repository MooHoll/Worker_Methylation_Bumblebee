### Running command line Blast on ALICE2

#!/bin/bash

#PBS -q devel
#PBS -N blast
#PBS -l walltime=00:05:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load blast+/2.5.0

blastx \
-query upreg_repro_genes.fasta \
-db "nr" \
-entrez_query "Drosophila melanogaster [organism]" \
-remote \
-max_target_seqs 1 \
-outfmt 5 \
-out upreg_repro_genes.xml \
-evalue 1e-3

## After this put the fasta and xml file into Blast2GO and put through all steps to get
## a prorogation report for GO enrichment analysis. 