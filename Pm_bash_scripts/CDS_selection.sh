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
bcftools view -r $CHROM:$START-$STOP -Oz -o sample_VCFs/Pm/CDS/$NAME.vcf.gz Pm_monoclonals_missingness_only.vcf.gz
done < target_CDS.bed

for i in sample_VCFs/Pm/CDS/*.vcf.gz
do
mkdir sample_VCFs/Pm/CDS/${i%.vcf.gz}
bcftools +split -Oz -o ${i%.vcf.gz} ${i}

done


while IFS=$'\t' read -r f1 f2 f3 f4
do CHROM=("$f1") START=("$f2") STOP=("$f3") NAME=("$f4")

for i in sample_VCFs/Pm/CDS/$NAME/*.vcf.gz

do
bcftools index ${i}
samtools faidx /proj/ideel/resources/genomes/Pmalariae/PlasmoDB-67_PmalariaeUG01_Genome.fasta $CHROM:$START-$STOP | bcftools consensus -H R -o ${i%.vcf.gz}_$NAME.fa ${i}
done

done < target_CDS.bed

mkdir sample_VCFs/Pm/CDS/CRT_concat

mv sample_VCFs/Pm/CDS/CRT*/*.fa sample_VCFs/Pm/CDS/CRT_concat

cd sample_VCFs/Pm/CDS/CRT_concat

for i in *CRT1.fa; do tail -q -n +2 ${i%1.fa}*.fa > ${i%1.fa}_concat.fa

done
for i in *concat.fa;
do 
sed -i "1s/^/>${i%_concat.fa}\n/" ${i} 

done

cd ..

mkdir CSP_concat

cd CSP_concat/

mv ../CSP*/*.fa .

for i in *CSP1.fa; do tail -q -n +2 ${i%1.fa}*.fa > ${i%1.fa}_concat.fa

done
for i in *concat.fa;
do 
sed -i "1s/^/>${i%_concat.fa}\n/" ${i} 

done

cd ..

mkdir MRP2_concat

cd MRP2_concat/

mv ../MRP2*/*.fa .

for i in *MRP21.fa; do tail -q -n +2 ${i%1.fa}*.fa > ${i%1.fa}_concat.fa

done
for i in *concat.fa;
do 
sed -i "1s/^/>${i%_concat.fa}\n/" ${i} 

done

cd ..

mkdir MSP1_concat

cd MSP1_concat/

mv ../MSP1*/*.fa .

for i in *MSP11.fa; do tail -q -n +2 ${i%1.fa}*.fa > ${i%1.fa}_concat.fa

done
for i in *concat.fa;
do 
sed -i "1s/^/>${i%_concat.fa}\n/" ${i} 

done

cd ..

mkdir PPPK-DHPS_concat

cd PPPK-DHPS_concat/

mv ../PPPK-DHPS*/*.fa .

for i in *PPPK-DHPS1.fa; do tail -q -n +2 ${i%1.fa}*.fa > ${i%1.fa}_concat.fa

done
for i in *concat.fa;
do 
sed -i "1s/^/>${i%_concat.fa}\n/" ${i} 

done

cd /work/users/z/p/zpopkinh/Pm_rerun/Variants/

mkdir fastas/CDS/Pm/

mv sample_VCFs/Pm/CDS/*/*.fa fastas/CDS/Pm/

cd fastas/CDS/Pm/

for i in AMA1 CRT_concat CSP_concat DHFR-TS Kelch13 LSA1 MDR1 MDR2 MRP1 MRP2_concat MSP1_concat P25 P48-45 PPPK-DHPS_concat TRAP

do

awk '/^>/{print ">" substr(FILENAME,1,length(FILENAME)-3); next} 1' *${i}.fa > ${i}.fa

done


cd /work/users/z/p/zpopkinh/Pm_rerun/Variants/

mkdir sample_VCFs/Pf/CDS/

while IFS=$'\t' read -r f1 f2 f3 f4
do CHROM=("$f1") START=("$f2") STOP=("$f3") NAME=("$f4")
echo $CHROM $START $STOP $NAME
bcftools view -r $CHROM:$START-$STOP -Oz -o sample_VCFs/Pf/CDS/$NAME.vcf.gz Pf_ortholog_samples.vcf.gz
done < Pf_target_CDS.bed

for i in sample_VCFs/Pf/CDS/*.vcf.gz
do
mkdir sample_VCFs/Pf/CDS/${i%.vcf.gz}
bcftools +split -Oz -o ${i%.vcf.gz} ${i}

done

while IFS=$'\t' read -r f1 f2 f3 f4
do CHROM=("$f1") START=("$f2") STOP=("$f3") NAME=("$f4")

for i in sample_VCFs/Pf/CDS/$NAME/*.vcf.gz

do
bcftools index ${i}
samtools faidx PlasmoDB-47_Pfalciparum3D7_Genome.fasta $CHROM:$START-$STOP | bcftools consensus -H R -o ${i%.vcf.gz}_$NAME.fa ${i}
done

done < Pf_target_CDS.bed

mkdir sample_VCFs/Pf/CDS/CRT_concat

mv sample_VCFs/Pf/CDS/CRT*/*.fa sample_VCFs/Pf/CDS/CRT_concat

cd sample_VCFs/Pf/CDS/CRT_concat

for i in *CRT1.fa; do tail -q -n +2 ${i%1.fa}*.fa > ${i%1.fa}_concat.fa

done
for i in *concat.fa;
do 
sed -i "1s/^/>${i%_concat.fa}\n/" ${i} 

done

cd ..

mkdir PPPK-DHPS_concat

cd PPPK-DHPS_concat/

mv ../PPPK-DHPS*/*.fa .

for i in *PPPK-DHPS1.fa; do tail -q -n +2 ${i%1.fa}*.fa > ${i%1.fa}_concat.fa

done
for i in *concat.fa;
do 
sed -i "1s/^/>${i%_concat.fa}\n/" ${i} 

done

cd /work/users/z/p/zpopkinh/Pm_rerun/Variants/

mkdir fastas/CDS/Pf/

mv sample_VCFs/Pf/CDS/*/*.fa fastas/CDS/Pf/

cd fastas/CDS/Pf/

for i in AMA1 CRT_concat CSP DHFR-TS Kelch13 LSA1 MDR1 MDR2 MRP1 MRP2 MSP1 P25 P48-45 PPPK-DHPS_concat TRAP

do

awk '/^>/{print ">" substr(FILENAME,1,length(FILENAME)-3); next} 1' *${i}.fa > ${i}.fa

done

cd ..

for i in **/CRT_concat.fa **/CSP_concat.fa **/MRP2_concat.fa **/MSP1_concat.fa **/PPPK-DHPS_concat.fa; do mv ${i} ${i%_concat.fa}.fa; done

for i in AMA1 CRT CSP DHFR-TS Kelch13 LSA1 MDR1 MDR2 MRP1 MRP2 MSP1 P25 P48-45 PPPK-DHPS TRAP

do cat Pm/${i}.fa Pf/${i}.fa > ${i}.fa

done

for i in AMA1 CRT CSP DHFR-TS Kelch13 LSA1 MDR1 MDR2 MRP1 MRP2 MSP1 P25 P48-45 PPPK-DHPS TRAP

do /nas/longleaf/rhel8/apps/mafft/7.490/bin/mafft  --auto --reorder ${i}.fa > ${i}_CDS_alignment.fa
done