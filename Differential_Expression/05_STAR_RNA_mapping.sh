#!/bin/bash

#PBS -N alignments_RNA
#PBS -l walltime=00:30:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=2:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load star/2.5.2b
                                                                                      
# create the directory where the output files are to be written  
OUTPUT=alignments                                                                                                                                 
if [ ! -d "$OUTPUT" ]; then                                                                                 
    mkdir -p ${OUTPUT} 
fi

######################################################################
REF=/scratch/monoallelic/hm257/repro_transcription/genome

# Run STAR alignments 
for file in $(ls *1.fq.gz)
do
    base=$(basename $file "1.fq.gz")
    STAR \
    --runThreadN 16 \
    --genomeDir ${REF} \
    --readFilesCommand gunzip -c \
    --readFilesIn ${base}1.fq.gz ${base}2.fq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${base}
done


