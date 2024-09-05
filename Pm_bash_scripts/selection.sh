#!/bin/bash
##############################################################

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=100G
#SBATCH -t 8:00:00

cd /work/users/z/p/zpopkinh/Pm_rerun/Variants/

#bcftools norm -f /proj/ideel/resources/genomes/Pmalariae/PlasmoDB-67_PmalariaeUG01_Genome.fasta -m +any -Oz -o Pm_monoclonal_normalized.vcf.gz Pm_monoclonals_missingness_only.vcf.gz

while IFS=$'\t' read -r f1 f2 f3 f4
do CHROM=("$f1") START=("$f2") STOP=("$f3") NAME=("$f4")
echo $CHROM $START $STOP $NAME
bcftools view -r $CHROM:$START-$STOP -Oz -o sample_VCFs/Pm/$NAME.vcf.gz Pm_monoclonals_missingness_only.vcf.gz
done < target_genes.bed

for i in sample_VCFs/Pm/*.vcf.gz
do
mkdir sample_VCFs/Pm/${i%.vcf.gz}
bcftools +split -Oz -o ${i%.vcf.gz} ${i}

done

#bedtools getfasta -split -fo target_genes.fa -fi /proj/ideel/resources/genomes/Pmalariae/PlasmoDB-67_PmalariaeUG01_Genome.fasta -bed target_genes.bed 

#cat target_genes.fa | awk '{
#        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}
#        print $0 > filename
#}'

#for i in *.fa

#do samtools faidx ${i}

#done


while IFS=$'\t' read -r f1 f2 f3 f4
do CHROM=("$f1") START=("$f2") STOP=("$f3") NAME=("$f4")

for i in sample_VCFs/Pm/$NAME/*.vcf.gz

do
bcftools index ${i}
samtools faidx /proj/ideel/resources/genomes/Pmalariae/PlasmoDB-67_PmalariaeUG01_Genome.fasta $CHROM:$START-$STOP | bcftools consensus -H I -o ${i%.vcf.gz}_$NAME.fa ${i}
done

done < target_genes.bed

mv sample_VCFs/Pm/*/*.fa fastas/

cd fastas/

for i in AMA1 CRT CSP DHFR-TS Kelch13 LSA1 MDR1 MDR2 MRP1 MRP2 MSP1 P25 P48-45 PPPK-DHPS TRAP

do

awk '/^>/{print ">" substr(FILENAME,1,length(FILENAME)-3); next} 1' *${i}.fa > ${i}.fa

done

#for i in AMA1 CRT CSP DHFR-TS Kelch13 LSA1 MDR1 MDR2 MRP1 MRP2 MSP1 P25 P48-45 PPPK-DHPS TRAP

#do

#awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' ${i}.fa > ${i}_noempties.fa

#done

cd Pm/work/users/z/p/zpopkinh/Pm_rerun/Variants/

#bcftools norm -f /proj/ideel/resources/genomes/Pmalariae/PlasmoDB-67_PmalariaeUG01_Genome.fasta -m +any -Oz -o Pm_monoclonal_normalized.vcf.gz Pm_monoclonals_missingness_only.vcf.gz

while IFS=$'\t' read -r f1 f2 f3 f4
do CHROM=("$f1") START=("$f2") STOP=("$f3") NAME=("$f4")
echo $CHROM $START $STOP $NAME
bcftools view -r $CHROM:$START-$STOP -Oz -o sample_VCFs/Pf/$NAME.vcf.gz Pf_ortholog_samples.vcf.gz
done < Pf_target_genes.bed

for i in sample_VCFs/Pf/*.vcf.gz
do
mkdir sample_VCFs/Pf/${i%.vcf.gz}
bcftools +split -Oz -o ${i%.vcf.gz} ${i}

done

#bedtools getfasta -split -fo Pf_target_genes.fa -fi /proj/ideel/resources/genomes/Pfalciparum/genomes/PlasmoDB-47_Pfalciparum3D7_Genome.fasta -bed Pf_target_genes.bed 

#cat Pf_target_genes.fa | awk '{
#        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}
#        print $0 > filename
#}'

#for i in *.fa

#do samtools faidx ${i}

#done

while IFS=$'\t' read -r f1 f2 f3 f4
do CHROM=("$f1") START=("$f2") STOP=("$f3") NAME=("$f4")

for i in sample_VCFs/Pf/$NAME/*.vcf.gz

do
bcftools index ${i}
samtools faidx /proj/ideel/resources/genomes/Pfalciparum/genomes/PlasmoDB-47_Pfalciparum3D7_Genome.fasta | bcftools consensus -H I -o ${i%.vcf.gz}_$NAME.fa ${i}
done

done < Pf_target_genes.bed