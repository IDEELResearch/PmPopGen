# PmPopGen

This repository contains code used to process and analyze 81 <i>Plasmodium malariae</i> genomes generated via hybrid capture sequencing. The original data files are available via SRA (BioProject ID PRJNA1157442).

Shell scripts are used to process sequence read files and are optimized for SLURM sbatch submission. While it is possible to adapt them to run on a Unix desktop, it will take much longer and require a large quantity of memory and storage space.

Initial processing using these shell scripts is expected to take up to a week on a high-performance computing cluster.

In order to run the subsequent analysis, you will need a working installation of R and RStudio. Code was optimized to run in R 4.2.2 and RStudio 2022.07.2. Individual R packages are detailed within each respective script.

Shell scripts rely on functioning installations of trim_galore, BBMap, bwa-mem2, GATK4, bcftools, vcftools, Tandem Repeats Finder, bedtools, python3, dadi-cli, donni, ADMIXTURE, PLINK v1.9, and RAxML Next Generation.
