###############################################################
################# Pm_VariantFiltration ########################
###############################################################
#Description: applies hard quality filtering and missingness thresholds to VCF file,
#limits to SNPs, and excludes hypervariable regions


#! /bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -t 8:00:00

module add gatk

module add samtools

module add vcftools

cd /work/users/z/p/zpopkinh/Pm_rerun/Variants/

#apply quality filter threshoolds
	gatk VariantFiltration \
	-R /proj/ideel/resources/genomes/Pmalariae/PlasmoDB-67_PmalariaeUG01_Genome.fasta/ \
	-V Pm_HC_raw.vcf.gz \
	-O Pm_HC_hard_filtered.vcf.gz \
	#each filter expression must be updated with final threshold
	--filter-name "lowQD" \
	--filter-expression "QD<2.5" \
	--filter-name "highFS" \
	--filter-expression "FS>10.0" \
	--filter-name "lowMQ" \
	--filter-expression "MQ<50.0" \
	--filter-name "lowMQRankSum" \
	--filter-expression "MQRankSum<-2.5" \
	--filter-name "lowReadPosRankSum" \
	--filter-expression "ReadPosRankSum<-2.5"

#limit variants to SNPs only
bcftools view -m2 -M2 -v snps Pm_HC_hard_filtered.vcf.gz -O z -o Pm_HC_hard_filtered_biallelic_snps_only.vcf.gz

#output table of SNPs
gatk VariantsToTable -V Pm_HC_hard_filtered_biallelic_snps_only.vcf.gz -F QD -F FS -F MQ -F MQRankSum -F ReadPosRankSum -O filtered_biallelics.table

#mask tandem repeats
bcftools view -T ^Pmalariae_trf_sort.bed -Oz -o Pm_TRs_masked.vcf.gz Pm_HC_hard_filtered_biallelic_snps_only.vcf.gz

#mask hypervariable Plasmodium Interspersed Repeat (pir) regions
bcftools view -T ^Pm_PIRs.bed -Oz -o Pm_PIRs_masked.vcf.gz Pm_TRs_masked.vcf.gz

#excludes variant sites with missing data in >20% of samples
vcftools --gzvcf Pm_PIRs_masked.vcf.gz --max-missing 0.8 --recode --recode-INFO-all --stdout | gzip -c > Pm_HC_missingness_filtered_first.vcf.gz
