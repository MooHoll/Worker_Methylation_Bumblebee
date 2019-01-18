#!/bin/bash

#PBS -N STAR_genome_index
#PBS -l walltime=00:10:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load star/2.5.2b

# create the directory where the output files are to be written (vital for this)   
OUTPUT=genome_dir                                                                                                                                 
if [ ! -d "$OUTPUT" ]; then                                                                                 
    mkdir -p ${OUTPUT} 
fi

# Define file paths
REF_FA=/scratch/monoallelic/hm257/repro_transcription/genome/GCF_000214255.1_Bter_1.0_genomic.fa
REF_GFF=/scratch/monoallelic/hm257/repro_transcription/genome/ref_Bter_1.0_top_level.gff3
                                                                                      

# Run STAR to index new reference genomes
    STAR \
    --runMode genomeGenerate \
    --genomeFastaFiles ${REF_FA} \
    --sjdbGTFfile ${REF_GFF} \
    --sjdbGTFtagExonParentTranscript Parent \
    --runThreadN 8 \
    --genomeDir ./${OUTPUT}      
