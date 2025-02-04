cd /mnt/c/Users/zacha/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/1D_demographics

for i in bottlegrowth_1d growth snm_1d three_epoch two_epoch

do dadi-cli StatDM --fs ../dadi_1pop_wsaf.fs --model ${i} --demo-popt ${i}.demog.params.InferDM.bestfits --grids 80 90 100 --bootstrapping-dir bootstraps/ --output ${i}.godambe.ci --nomisid > ${i}_matrices.txt

done