#!/bin/bash
##############################################################

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH -t 4:00:00

module load samtools

#module load gatk

module load picard

cd /work/users/z/p/zpopkinh/Pm_rerun/bwa_output/

for i in *.bam
do
sbatch -p general -N 1 -n 12 --mem=100G -t 24:00:00 --wrap="picard AddOrReplaceReadGroups I=${i} O=${i%.fq.bam}_sorted.bam SORT_ORDER=coordinate RGID=${i%.fq.bam}.RGID RGLB=${i%.fq.bam}.RGID RGPL=illumina RGPU=${i%.fq.bam}.RGPU RGSM=${i%.fq.bam}.RGSM;picard MarkDuplicates I=${i%.fq.bam}_sorted.bam O=${i%.fq.bam}_dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${FILE}_dedupped.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000"
done