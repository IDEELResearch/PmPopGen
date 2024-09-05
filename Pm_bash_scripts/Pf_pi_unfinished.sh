#!/bin/bash
##############################################################

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=100G
#SBATCH -t 4:00:00

cd /work/users/z/p/zpopkinh/Pm_rerun/Variants/

while IFS=$'\t' read -r f1 f2 f3

do  CHROM=("$f1") START=("$f2") STOP=("$f3") LENGTH=$(($STOP-$START))
vcftools --gzvcf Pf_ortholog_samples_orthos_only.vcf.gz --chr $CHROM --from-bp $START --to-bp $STOP --window-pi $LENGTH --window-pi-step 1 --out Pf_pi/$CHROM-$START-$STOP
done < Pf_masked_orthologs_unfinished.bed

#append these into one file for each species then use R to extract only the lines that actually matter by comparing it to the ortholog intervals

cd Pf_pi/
cat *.windowed.pi >> all_Pf_pi.txt