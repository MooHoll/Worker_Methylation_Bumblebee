# Running HTSeq for counting reads per genes for transcription data (needed for differential gene expression analysis)

# First need to make a virtual python environment on ALICE2 (HPC)
# Making it inside the directory called 'ptyhon_environ'
module load python/gcc/27
virtualenv --system-site-packages --prompt=$USER-python_environ ./

# Then activate the environment before use 
source /scratch/monoallelic/hm257/repro_transcription/python_enviro/bin/activate

# Need PySam to read a .bam file
pip install PySam

# Don't need to install HTSeq as already available (11:50am)
# Tested code with one file:
python //
-m HTSeq.scripts.count //
-i gene //
-f bam //
-o annotated_j1_17.bam //
/scratch/monoallelic/hm257/repro_transcription/bams/trim_j1_17_Aligned.sortedByCoord.out.bam //
/scratch/monoallelic/hm257/repro_transcription/genome/ref_Bter_1.0_top_level.gff3 //
> j1_17.counts

# It works so put the following in a shell script and called from the activated python environment
source /scratch/monoallelic/hm257/repro_transcription/python_enviro/bin/activate

pip install PySam

REF_GFF=/scratch/monoallelic/hm257/repro_transcription/genome/ref_Bter_1.0_top_level.gff3

for file in $(ls /scratch/monoallelic/hm257/repro_transcription/bams/*.bam)
do
    base=$(basename ${file} ".bam")
    python \
    -m HTSeq.scripts.count \
    -i gene \
    -f bam \
    -o ${base}_annotated.bam \
    ${file} \
    ${REF_GFF} \
    > ${base}.counts
done

# Then put this in a .PBS script:

#!/bin/bash

#PBS -N read_counts_for_genes
#PBS -l walltime=06:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=2:ppn=8

# Run script in the working directory it was submitted in
cd $PBS_O_WORKDIR

# Load software needed
module load python/gcc/27

bash htseq.sh