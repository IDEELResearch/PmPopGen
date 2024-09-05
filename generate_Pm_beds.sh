#!/bin/bash
##############################################################

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=100G
#SBATCH -t 24:00:00

cd /work/users/z/p/zpopkinh/Pm_rerun/Picard_output/

### generate bed files for each sample by species call that show the intervals of the genome with over 5X covergae
for i in *dedupped.bam
do bedtools genomecov -ibam ${i} -bg | awk '$4>=5' > /work/users/z/p/zpopkinh/Pm_rerun/Variants/${i%_dedupped.bam}_cov5.bed
done

cd /work/users/z/p/zpopkinh/Pm_rerun/Variants/

mkdir Pm_cov5_beds

mv *cov5.bed Pm_cov5_beds/

### generate a bed file reflecting intervals where {prop} of samples had at least {cov_filter} coverage (see config.yaml)

cd Pm_cov5_beds

bedtools multiinter -i *.bed | awk '$4>=0.6' > Pm_60percentcov5_unmerged.bed

bedtools merge -i Pm_60percentcov5_unmerged.bed -d 10 > Pm_60percentcov5_merged.bed

### generate bed files for each sample by species call that show the intervals of the genome with over 5X covergae Pf on Oscar
#for i in $(ls /users/zpopkinh/data/shared/plasmodium/falciparum/pfpubdata/bams/*sorted.bam)
#do bedtools genomecov -ibam ${i} -bg | awk '$4>=5' > ./${i%.mdups.sorted.bam}_cov5.bed
#done
