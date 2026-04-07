#### SPLIT VCF BY CHR ####

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
do bcftools view -r PmUG01_${i}_v1 -O z -o nsl/${i}.vcf.gz Pm_monoclonals_missingness_only.vcf.gz
done

#### RUN nSL ####

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
do

bcftools +missing2ref -O z -o ${i}_fixed.vcf.gz ${i}.vcf.gz

selscan --nsl --vcf ${i}_fixed.vcf.gz --out ${i}.res

done


#### NORMALIZE nSL ####
for i in *.res.nsl.out

do

norm --ihs --files ${i}

done