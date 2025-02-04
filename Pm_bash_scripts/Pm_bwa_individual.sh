#!/bin/bash
##############################################################

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH -t 1:00:00

module load bwa-mem2

module load samtools

cd /work/users/z/p/zpopkinh/Pm_rerun/

cd bbsplit_output

for i in *Pm.fq;
do sbatch -p general -N 1 -n 12 --mem=100g -t 24:00:00 --wrap="bwa-mem2 mem -M -t 12 /proj/ideel/resources/genomes/Pmalariae/PlasmoDB-67_PmalariaeUG01_Genome.fasta  ${i} | samtools view -bS - > ${i%\w[10]-\w[10]_Pm.fq}.bam";
done
