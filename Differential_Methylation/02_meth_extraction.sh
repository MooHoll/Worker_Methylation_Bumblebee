#!/bin/bash

#PBS -N m_extraction
#PBS -l walltime=16:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run in current working directory
cd $PBS_O_WORKDIR

# Load modules
module load bismark/0.18.1
module load samtools/1.3.2

# Run full methylation extraction including read trimming

bismark_methylation_extractor -p --no_overlap --comprehensive --bedgraph --ignore 3 --ignore_r2 7 --ignore_3prime_r2 2 --report j1nr.deduplicated.$
bismark_methylation_extractor -p --no_overlap --comprehensive --bedgraph --ignore 5 --ignore_r2 8 --ignore_3prime_r2 5 --report j1r.deduplicated.b$
bismark_methylation_extractor -p --no_overlap --comprehensive --bedgraph --ignore 8 --ignore_r2 8 --report j5nr.deduplicated.bam
bismark_methylation_extractor -p --no_overlap --comprehensive --bedgraph --ignore 5 --ignore_r2 8 --ignore_3prime_r2 25 --report j5r.deduplicated.$
bismark_methylation_extractor -p --no_overlap --comprehensive --bedgraph --ignore 3 --ignore_r2 8 --ignore_3prime_r2 2 --report j8nr.deduplicated.$
bismark_methylation_extractor -p --no_overlap --comprehensive --bedgraph --ignore 5 --ignore_r2 8 --ignore_3prime_r2 5 --report j8r.deduplicated.b$

## Need to update this script with actual methylation extraction and deduplication