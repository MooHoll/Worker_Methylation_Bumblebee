#!/bin/bash

#PBS -N fastqc
#PBS -l walltime=05:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Change directory to the one the job was submitted in
cd $PBS_O_WORKDIR 

# Load required modules
module load fastqc/0.11.5                                                                

# create the directory where the output files are to be written 
OUTPUT=fastqc                                                                                                                                     
if [ ! -d "$OUTPUT" ]; then                                                                                 
    mkdir -p ${OUTPUT} 
fi

# Create a list of the files to be called                                                                                                            
FILES=*.fq.gz

for file in ${FILES}
do
	fastqc -o ${OUTPUT} ${file}
done


### NOTE: prior to this script all files were downloaded from BGI servers and concatenated (cat command) into
# one forward and one reverse read file, then zipped (gzip command).