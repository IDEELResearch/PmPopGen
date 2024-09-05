#!/bin/bash
##############################################################

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=1G
#SBATCH -t 1:00:00

module load samtools

module load gatk

module load picard

cd /work/users/z/p/zpopkinh/Pm_rerun/Picard_output/

for i in *dedupped.bam
do
sbatch -p general -N 1 -n 12 --mem=100G -t 24:00:00 --wrap="gatk HaplotypeCaller -R /proj/ideel/resources/genomes/Pmalariae/PlasmoDB-67_PmalariaeUG01_Genome.fasta -I ${i} -O ${i%_dedupped.bam}.g.vcf.gz -ERC GVCF"
done