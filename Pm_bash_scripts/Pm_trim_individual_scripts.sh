#!/bin/bash
##############################################################

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=1G
#SBATCH -t 1:00:00

module load trim_galore

cd /work/users/z/p/zpopkinh/Pm_full_HC/

mkdir NovaSeq

cd NovaSeq/

mkdir trim_galore_output

mkdir bbsplit_output

for i in /proj/ideel/julianog/HTSF/230914_UNC41-A00434_0694_AHJKN7DSX7/*S*_L004_R*.fastq.gz;
do sbatch -p general -N 1 -n 4 --mem=100g -t 4:00:00 --wrap="trim_galore ${i%[12]_001.fastq.gz}1_001.fastq.gz ${i%[12]_001.fastq.gz}2_001.fastq.gz --illumina --paired --fastqc -o /work/users/z/p/zpopkinh/Pm_full_HC/NovaSeq/trim_galore_output/";
done
