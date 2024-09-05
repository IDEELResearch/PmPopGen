#!/bin/bash
##############################################################

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=100G
#SBATCH -t 24:00:00
cd /work/users/z/p/zpopkinh/Pm_rerun/Pf_VCFs/ortholog_samples/

for i in ../*.vcf.gz

do bcftools view -S Pf_ortholog_samples.txt -O z -o ./${i%.vcf.gz}_ortholog_samples.vcf.gz ${i}
done

bcftools concat *ortholog_samples.vcf.gz -Oz -o Pf_ortholog_samples.vcf.gz

bcftools index Pf_ortholog_samples.vcf.gz