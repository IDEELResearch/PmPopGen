#!/bin/bash
##############################################################

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=100G
#SBATCH -t 8:00:00

module load trim_galore

module load bwa-mem2

module load samtools

cd /work/users/z/p/zpopkinh/Pm_full_HC/

#cp -r ref/ NovaSeq/trim_galore_output/

#mkdir NovaSeq

cd NovaSeq/

#mkdir trim_galore_output

#mkdir bbsplit_output

#shopt -s globstar

#for i in /proj/ideel/julianog/HTSF/230914_UNC41-A00434_0694_AHJKN7DSX7/*S*_L004_R1.fastq.gz;
#do sbatch -p general -N 1 -n 4 --mem=100g -t 24:00:00 --wrap="trim_galore ${i%[12]_001.fastq.gz}1_001.fastq.gz ${i%[12]_001.fastq.gz}2_001.fastq.gz --illumina --paired --fastqc -o /work/users/z/p/zpopkinh/Pm_full_HC/NovaSeq/trim_galore_output/";
#done

#/nas/longleaf/apps/bbmap/38.96/bbmap/bbsplit.sh build=1 ref_Pm=/proj/ideel/resources/genomes/Pmalariae/PmalariaeUG01.fasta ref_Hs=/proj/ideel/resources/genomes/Hsapiens/hg38.fa ref_Pf=/proj/ideel/resources/genomes/Pfalciparum/genomes/Pf3D7.fasta

cd trim_galore_output/
for i in /work/users/z/p/zpopkinh/Pm_full_HC/NovaSeq/trim_galore_output/*S*_L004_R1*.fq.gz; 
do #mkdir ${i%_L004*.fq.gz}
cd ${i%_L004*.fq.gz}
#cp -r /work/users/z/p/zpopkinh/Pm_full_HC/ref/ .
sbatch -p general -N 1 -n 4 --mem=100g -t 24:00:00 --wrap="/nas/longleaf/apps/bbmap/38.96/bbmap/bbsplit.sh build=1 in=${i%[12]_001_val_[12].fq.gz}1_001_val_1.fq.gz in2=${i%[12]_001_val_[12].fq.gz}2_001_val_2.fq.gz basename=${i%_L004*.fq.gz}_%.fq ziplevel=2 refstats=${i%_L004*.fq.gz}_refstats.txt"; 
done; 

#mv *_\w\w.fq bbsplit_output/

#for i in bbsplit_output/*Pm.fq;
#do bwa-mem2 mem -M -t 12 /proj/ideel/resources/genomes/Pmalariae/PmalariaeUG01.fasta ${i} | samtools view -bS - > ${i%\w[10]-\w[10]_Pm.fq}.bam;
#done
