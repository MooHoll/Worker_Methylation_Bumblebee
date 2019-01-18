#!/bin/bash

#PBS -N bismark_alignment
#PBS -l walltime=18:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load bismark/0.16.3
module load bowtie2/2.2.9
module load samtools/1.3.2

# Define file paths
REF_FA=/scratch/monoallelic/hm257/liverpool_tests/genome

# Align all fasta files to the reference
for file in $(ls *1.fq)
do
    base=$(basename $file "1.fq")
    bismark \
    ${REF_FA} \
    -1 ${base}1.fq \
    -2 ${base}2.fq
done
                                                                                                                                                                                                                                                                                                                                                                      


     