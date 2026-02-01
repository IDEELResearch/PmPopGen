###############################################################
############### Pm_HC_gentype_gVCFs ###########################
###############################################################
#Description: uses GATK GenotypeGVCFs to genotype individual gVCF 
#files across entire sample pool, yielding vcf file showing variant sites

#!/bin/bash
##############################################################

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=200G
#SBATCH -t 24:00:00

module load samtools

module load gatk

module load picard

cd /work/users/z/p/zpopkinh/Pm_rerun/

mkdir Variants

cd Picard_output

#generate DBImport database of samples, reference genome intervals, and output directory
gatk GenomicsDBImport --sample-name-map ../Pm.sample_map --genomicsdb-workspace-path ../Variants/gVCF -L /proj/ideel/resources/genomes/Pmalariae/PlasmoDB-67_PmalariaeUG01_Genome.fasta.bed

#produce VCF file containing all variants in Pm genome across all samples
gatk GenotypeGVCFs -R /proj/ideel/resources/genomes/Pmalariae/PlasmoDB-67_PmalariaeUG01_Genome.fasta -V gendb://../Variants/gVCF -O ../Variants/Pm_HC_raw.vcf.gz

#output table of variant position, statistics, and quality metrics
gatk VariantsToTable -V Pm_HC_raw.vcf.gz -F QD -F FS -F MQ -F MQRankSum -F ReadPosRankSum -O raw.table
