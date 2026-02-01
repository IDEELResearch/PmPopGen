###############################################################
######################### Pm_FST ##############################
###############################################################
#Description: Calculates weir-Fst between countries of origin

setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun")

Pm_wsaf_filtered <- vcfR::read.vcfR("Pm_monoclonals_wsaf_filtered.vcf.gz") |> vcfR::vcfR2tidy()

#generate sample list by country
sample_list <- Pm_wsaf_filtered$gt$Indiv |> unique()

sample_list <- sample_list |> as.data.frame()

names(sample_list) <- "Sample"

sample_list <- sample_list |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                       stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                       stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                       stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                       .default = "Tanzania"))

sample_list |> dplyr::group_by(Country) |> dplyr::group_split()

country_list <- sample_list |> dplyr::group_by(Country) |> dplyr::group_split()

Cameroon_list <- country_list[[1]] |> as.data.frame()

DRC_list <- country_list[[2]] |> as.data.frame()

Nigeria_list <- country_list[[3]] |> as.data.frame()

Tanzania_list <- country_list[[4]] |> as.data.frame()

Cameroon_list$Sample |> as.list() |> data.table::fwrite("Cameroon_samples.txt", sep = "\n", col.names = F)

DRC_list$Sample |> as.list() |> data.table::fwrite("DRC_samples.txt", sep = "\n", col.names = F)

Nigeria_list$Sample |> as.list() |> data.table::fwrite("Nigeria_samples.txt", sep = "\n", col.names = F)

Tanzania_list$Sample |> as.list() |> data.table::fwrite("Tanzania_samples.txt", sep = "\n", col.names = F)

#calculate Weir Fst between each combination of countries
system("vcftools --gzvcf Pm_monoclonals_wsaf_filtered.vcf.gz --weir-fst-pop Cameroon_samples.txt --weir-fst-pop DRC_samples.txt --out Cam_DRC")

system("vcftools --gzvcf Pm_monoclonals_wsaf_filtered.vcf.gz --weir-fst-pop Cameroon_samples.txt --weir-fst-pop Nigeria_samples.txt --out Cam_Nigeria")       

system("vcftools --gzvcf Pm_monoclonals_wsaf_filtered.vcf.gz --weir-fst-pop Cameroon_samples.txt --weir-fst-pop Tanzania_samples.txt --out Cam_Tanzania")       

system("vcftools --gzvcf Pm_monoclonals_wsaf_filtered.vcf.gz --weir-fst-pop DRC_samples.txt --weir-fst-pop Nigeria_samples.txt --out DRC_Nigeria")       

system("vcftools --gzvcf Pm_monoclonals_wsaf_filtered.vcf.gz --weir-fst-pop DRC_samples.txt --weir-fst-pop Tanzania_samples.txt --out DRC_Tanzania")       


system("vcftools --gzvcf Pm_monoclonals_wsaf_filtered.vcf.gz --weir-fst-pop Nigeria_samples.txt --weir-fst-pop Tanzania_samples.txt --out Nigeria_Tanzania")
