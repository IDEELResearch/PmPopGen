# PmPopGen

This repository contains code used to process and analyze 81 <i>Plasmodium malariae</i> genomes generated via hybrid capture sequencing. The original data files are available via SRA (BioProject ID PRJNA1157442).

Shell scripts are used to process sequence read files and are optimized for SLURM sbatch submission. While it is possible to adapt them to run on a Unix desktop, it will take much longer and require a large quantity of memory and storage space.

Initial processing using these shell scripts is expected to take up to a week on a high-performance computing cluster.

In order to run the subsequent analysis, you will need a working installation of R and RStudio. Code was optimized to run in R 4.2.2 and RStudio 2022.07.2. Individual R packages are detailed within each respective script.

Shell scripts rely on functioning installations of trim_galore, BBMap, bwa-mem2, GATK4, bcftools, vcftools, Tandem Repeats Finder, bedtools, python3, dadi-cli, donni, ADMIXTURE, PLINK v1.9, and RAxML Next Generation.

The following scripts are intended to be run in the order given:

1. Pm_trim_individual_scripts.sh: trims sequencing adapters from fastq files
2. Pm_bbsplit.sh: competitively aligns reads across multiple reference genomes
3. Pm_bwa_individual.sh: select reads best aligned to the Pm reference genome after competitive alignment
4. Pm_Picard.sh: adds readgroup information to Pm-aligned bam files
5. Pm_HC_generate_gVCF.sh: generate gVCF files of variants across Pm genome for each sample
6. Pm_HC_genotype_gVCFs.sh: genotypes individual gVCF files across entire sample pool, yielding vcf file showing variant sites
7. Pm_filtering_determination.R: titrate quality filtering thresholds for raw variant output from unfiltered vcf file
8. Pm_VariantFiltration.sh: applies hard quality filtering and missingness thresholds to VCF file, limits to SNPs, and excludes hypervariable regions
9. COI_Pm_coiaf.R: estimates complexity of infection for each sample using COIAF
10. Pf_sample_picker.R: selects, filters, and compiles metadata for P. falciparum samples from Pf7 database

These scripts must be completed before running any others (excluding scripts assessing sequencing metrics such as coverage, depth, and degree of enrichment) because the other scripts exclusively use monoclonal samples.

11. download_Pf7_vcfs.sh: downloads publicly-available P. falciparum VCF files from Pf7 database
12. Pf_ortholog_samples.sh: subsets VCF to matched samples for nucleotide diversity calculation and comparison
13. generate_Pm_beds.sh: generates bed files of P. malariae genome containing and excluding specific genomic intervals
14. ortholog_masker.R: generates bed files containing only the 1-to-1 orthologs between Pf and Pm genomes
15. Pm_pi.sh: calculates nucleotide diversity of in P. malariae orthologous genes
16. Pf_pi.sh: calculates nucleotide diversity of in P. malariae orthologous genes
17. Pm_pi.R: plots nucloetide diversity values across Pm and Pf orthologues
18. LD_decay.R: calculates and graphs the decay of linkage disequilibrium across genomic intervals in both P. malariae and P. falciparum
19. Pm_hmmibdr.R: uses hidden Markov model to identify genomic segments that are identical by descent
20. Pm_PCA.R: performs and graphs principal components analysis on monoclonal P. malariae isolates
21. Pm_DAPC.R: performs and plots discriminant analysis of principal components to identify samples clustered by genetic similarity
22. ADMIXTURE.R: Uses ADMIXTURE to calculate estimated number of population clusters within sample pool.
23. Pm_FST.R: Calculates weir-Fst between countries of origin
24. MIT_API_processing.sh: extracts mitochondrial and apicoplast sequences for phylogenetic analyses
25. Pm_phylogeny.R: generates maximum likelihood phylogenetic trees of Pm samples using RaxML
26. Pm_nSL.sh: calculates NSL by chromosome
27. Pm_selection.R: calculates, extracts, and visualizes Tajima's D and extracts and visualizes NSL among different genomic regions and intervals
28. selection.sh: generates bed file showing genomic windows containing specific Pm genes of interest
29. selection_Pf.sh: generates bed file showing genomic windows containing specific Pf genes of interest
30. CDS_selection.sh: extracts sequences of complementarity-determining regions of specific genes of interest for selection analysis

Other scripts do not need to be run in any particular order, other than needing to complete alignment and deduplicating before calculating coverage and depth.
