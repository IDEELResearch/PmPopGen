#!/bin/bash
##############################################################

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=100G
#SBATCH -t 24:00:00

module load samtools

cd /work/users/z/p/zpopkinh/Pm_rerun/

cd Picard_output

for i in *dedupped.bam
do samtools coverage ${i} > ${i%_dedupped.bam}_coverage.txt
done