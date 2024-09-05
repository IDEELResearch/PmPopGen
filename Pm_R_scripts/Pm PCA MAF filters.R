setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

#MOST RECENT - JUST FILTERING OUT INTERMEDIATE WSAF VALUES

#bcftools view -Oz -o Pm_monoclonals_wsaf_filtered.vcf.gz -i 'FORMAT/WSAF =0 | FORMAT/WSAF =1' Pm_monoclonals_wsaf.vcf.gz

#bcftools +prune -m 0.6 -e 'F_MISSING>=0.02' Pm_monoclonals_wsaf_filtered.vcf.gz -Oz -o Pm_monoclonals_wsaf_filtered.pruned.vcf.gz

#plink --pca var-wts --vcf Pm_monoclonals_wsaf_filtered.vcf.gz --const-fid --allow-extra-chr --out plink_WSAF_filtered --make-rel

#awk '{print $NR}' plink_WSAF_filtered.rel > plink_WSAF_filtered.rel.diag

Pm_wsaf_pruned_eigenvalues <- data.table::fread("plink_WSAF_filtered_pruned.eigenval")
Pm_wsaf_pruned_eigenvectors <- data.table::fread("plink_WSAF_filtered_pruned.eigenvec")

rel_diag <- data.table::fread("plink_WSAF_filtered_pruned.rel.diag")

Pm_wsaf_pruned_eigenvalues <- Pm_wsaf_pruned_eigenvalues |> dplyr::mutate(prop_var = V1/sum(rel_diag), prop_var2 = V1/sum(V1))

Pm_wsaf_table_pruned <- data.frame(sample.id = Pm_wsaf_pruned_eigenvectors$V2, 
                            EV1 = Pm_wsaf_pruned_eigenvectors$V3,   
                            EV2 = Pm_wsaf_pruned_eigenvectors$V4, 
                            stringsAsFactors = FALSE)

Pm_wsaf_table_pruned <- Pm_wsaf_table_pruned |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                           stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                           stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                           .default = "Tanzania"))


plot(Pm_wsaf_table_pruned$EV2, Pm_wsaf_table_pruned$EV1, xlab="eigenvector 2", ylab="eigenvector 1",)

library(ggplot2)
library(viridis)

Pm_wsaf_table_pruned |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC1 (",signif(Pm_wsaf_pruned_eigenvalues$prop_var[1], 4),"%)")) +
  ylab(paste0("PC2 (",signif(Pm_wsaf_pruned_eigenvalues$prop_var[2], 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK WSAF Filtered.png", width = 15, height = 12, units = "in", dpi = 600)

Pm_wsaf_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC1 (",signif(Pm_wsaf_eigenvalues[1], 4),"%)")) +
  ylab(paste0("PC2 (",signif(Pm_wsaf_eigenvalues[2], 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  xlim(-0.2, 0.05) +
  ylim(-0.25, 0.25) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

#####
Pm_wsaf_eigenvalues <- data.table::fread("plink_WSAF_filtered.eigenval")
Pm_wsaf_eigenvectors <- data.table::fread("plink_WSAF_filtered.eigenvec")

rel_diag <- data.table::fread("plink_WSAF_filtered.rel.diag")

Pm_wsaf_eigenvalues <- Pm_wsaf_eigenvalues |> dplyr::mutate(prop_var = V1/sum(rel_diag), prop_var2 = V1/sum(V1))

Pm_wsaf_table <- data.frame(sample.id = Pm_wsaf_eigenvectors$V2, 
                          EV1 = Pm_wsaf_eigenvectors$V3,   
                          EV2 = Pm_wsaf_eigenvectors$V4, 
                          stringsAsFactors = FALSE)

Pm_wsaf_table <- Pm_wsaf_table |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                       stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                       stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                       .default = "Tanzania"))


plot(Pm_wsaf_table$EV2, Pm_wsaf_table$EV1, xlab="eigenvector 2", ylab="eigenvector 1",)

library(ggplot2)
library(viridis)

Pm_wsaf_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC1 (",signif(Pm_wsaf_eigenvalues$prop_var[1], 4),"%)")) +
  ylab(paste0("PC2 (",signif(Pm_wsaf_eigenvalues$prop_var[2], 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK WSAF Filtered.png", width = 15, height = 12, units = "in", dpi = 600)

Pm_wsaf_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC1 (",signif(Pm_wsaf_eigenvalues[1], 4),"%)")) +
  ylab(paste0("PC2 (",signif(Pm_wsaf_eigenvalues[2], 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  xlim(-0.2, 0.05) +
  ylim(-0.25, 0.25) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK WSAF Filtered Zoom.png", width = 15, height = 12, units = "in", dpi = 600)

Pm_wsaf_table <- Pm_wsaf_table |> dplyr::mutate(Capture = dplyr::case_when(stringr::str_detect(sample.id, "HC[:digit:]") ~ stringr::str_extract(sample.id, "HC[:digit:]"),
                                                .default = "Test"))

Pm_wsaf_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Capture)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Capture") +
  xlab(paste0("PC1 (",signif(Pm_wsaf_eigenvalues[1], 4),"%)")) +
  ylab(paste0("PC2 (",signif(Pm_wsaf_eigenvalues[2], 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK WSAF Filtered Capture.png", width = 15, height = 12, units = "in", dpi = 600)

#ADEGENET METHOD

Pm_wsaf <- vcfR::read.vcfR("Pm_monoclonals_wsaf_filtered.vcf.gz") |> vcfR::vcfR2genlight()

gl_PCA <- adegenet::glPca(Pm_wsaf)

#redoing PCA with progressively more stringent MAF filters

#bcftools view -e 'MAF[0]<0.05' -r PmUG01_01_v1,PmUG01_02_v1,PmUG01_03_v1,PmUG01_04_v1,PmUG01_05_v1,PmUG01_06_v1,PmUG01_07_v1,PmUG01_08_v1,PmUG01_09_v1,PmUG01_10_v1,PmUG01_11_v1,PmUG01_12_v1,PmUG01_13_v1,PmUG01_14_v1 -Oz -o Pm_monoclonal_exclude_MAF5.vcf.gz Pm_monoclonals_missingness_only.vcf.gz

# $plink --pca --vcf Pm_monoclonal_exclude_MAF5.vcf.gz --const-fid --allow-extra-chr --out plink_MAF5

MAF5_eigenvalues <- data.table::fread("plink_MAF5.eigenval")
MAF5_eigenvectors <- data.table::fread("plink_MAF5.eigenvec")

MAF5_table <- data.frame(sample.id = MAF5_eigenvectors$V2, 
                         EV1 = MAF5_eigenvectors$V3,   
                         EV2 = MAF5_eigenvectors$V4, 
                         stringsAsFactors = FALSE)

MAF5_table <- MAF5_table |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                                   stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                                   stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                   .default = "Tanzania"))


plot(MAF5_table$EV2, MAF5_table$EV1, xlab="eigenvector 2", ylab="eigenvector 1",)

MAF5_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab("PC1 (29.75%)") +
  ylab("PC2 (5.90%)") +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK MAF5.png", width = 15, height = 12, units = "in", dpi = 600)

#bcftools view -e 'MAF[0]<0.1' -Oz -o Pm_monoclonal_exclude_MAF10.vcf.gz Pm_monoclonal_exclude_MAF5.vcf.gz

#$plink --pca --vcf Pm_monoclonal_exclude_MAF10.vcf.gz --const-fid --allow-extra-chr --out plink_MAF10

#REDID THIS ONE INCORPORATING WSAF FILTERING AS FOLLOWS

#bcftools view Pm_monoclonal_exclude_MAF10.vcf.gz | vcfdo wsaf | bcftools view -i 'FORMAT/WSAF >0.1' -Oz -o Pm_monoclonal_exclude_MAF10_WSAF10.vcf.gz

#$plink --pca --vcf Pm_monoclonal_exclude_MAF10_WSAF10.vcf.gz --const-fid --allow-extra-chr --out plink_MAF10_WSAF10

MAF10_eigenvalues <- data.table::fread("plink_MAF10_WSAF10.eigenval")
MAF10_eigenvectors <- data.table::fread("plink_MAF10_WSAF10.eigenvec")

MAF10_table <- data.frame(sample.id = MAF10_eigenvectors$V2, 
                         EV1 = MAF10_eigenvectors$V3,   
                         EV2 = MAF10_eigenvectors$V4, 
                         stringsAsFactors = FALSE)

MAF10_table <- MAF10_table |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                     stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                     stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                     .default = "Tanzania"))


plot(MAF10_table$EV2, MAF10_table$EV1, xlab="eigenvector 2", ylab="eigenvector 1",)

MAF10_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC1 (",signif(MAF10_eigenvalues[1], 4),"%)")) +
  ylab(paste0("PC2 (",signif(MAF10_eigenvalues[2], 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK MAF10 WSAF10.png", width = 15, height = 12, units = "in", dpi = 600)

#bcftools view -e 'MAF[0]<0.15' -Oz -o Pm_monoclonal_exclude_MAF15.vcf.gz Pm_monoclonal_exclude_MAF10.vcf.gz

#$plink --pca --vcf Pm_monoclonal_exclude_MAF15.vcf.gz --const-fid --allow-extra-chr --out plink_MAF15

MAF15_eigenvalues <- data.table::fread("plink_MAF15.eigenval")
MAF15_eigenvectors <- data.table::fread("plink_MAF15.eigenvec")

MAF15_table <- data.frame(sample.id = MAF15_eigenvectors$V2, 
                          EV1 = MAF15_eigenvectors$V3,   
                          EV2 = MAF15_eigenvectors$V4, 
                          stringsAsFactors = FALSE)

MAF15_table <- MAF15_table |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                       stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                       stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                       .default = "Tanzania"))


plot(MAF15_table$EV2, MAF15_table$EV1, xlab="eigenvector 2", ylab="eigenvector 1",)

MAF15_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC1 (",signif(MAF15_eigenvalues[1], 4),"%)")) +
  ylab(paste0("PC2 (",signif(MAF15_eigenvalues[2], 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK MAF15.png", width = 15, height = 12, units = "in", dpi = 600)

#bcftools view -e 'MAF[0]<0.2' -Oz -o Pm_monoclonal_exclude_MAF20.vcf.gz Pm_monoclonal_exclude_MAF15.vcf.gz

#$plink --pca --vcf Pm_monoclonal_exclude_MAF20.vcf.gz --const-fid --allow-extra-chr --out plink_MAF20

MAF20_eigenvalues <- data.table::fread("plink_MAF20.eigenval")
MAF20_eigenvectors <- data.table::fread("plink_MAF20.eigenvec")

MAF20_table <- data.frame(sample.id = MAF20_eigenvectors$V2, 
                          EV1 = MAF20_eigenvectors$V3,   
                          EV2 = MAF20_eigenvectors$V4, 
                          stringsAsFactors = FALSE)

MAF20_table <- MAF20_table |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                       stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                       stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                       .default = "Tanzania"))


plot(MAF20_table$EV2, MAF20_table$EV1, xlab="eigenvector 2", ylab="eigenvector 1",)

MAF20_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC1 (",signif(MAF20_eigenvalues[1], 4),"%)")) +
  ylab(paste0("PC2 (",signif(MAF20_eigenvalues[2], 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK MAF20.png", width = 15, height = 12, units = "in", dpi = 600)

#bcftools view -e 'MAF[0]<0.25' -Oz -o Pm_monoclonal_exclude_MAF25.vcf.gz Pm_monoclonal_exclude_MAF20.vcf.gz

#$plink --pca --vcf Pm_monoclonal_exclude_MAF25.vcf.gz --const-fid --allow-extra-chr --out plink_MAF25

MAF25_eigenvalues <- data.table::fread("plink_MAF25.eigenval")
MAF25_eigenvectors <- data.table::fread("plink_MAF25.eigenvec")

MAF25_table <- data.frame(sample.id = MAF25_eigenvectors$V2, 
                          EV1 = MAF25_eigenvectors$V3,   
                          EV2 = MAF25_eigenvectors$V4, 
                          stringsAsFactors = FALSE)

MAF25_table <- MAF25_table |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                       stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                       stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                       .default = "Tanzania"))


plot(MAF25_table$EV2, MAF25_table$EV1, xlab="eigenvector 2", ylab="eigenvector 1",)

MAF25_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC1 (",signif(MAF25_eigenvalues[1], 4),"%)")) +
  ylab(paste0("PC2 (",signif(MAF25_eigenvalues[2], 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK MAF25.png", width = 15, height = 12, units = "in", dpi = 600)

#bcftools view -e 'MAF[0]<0.3' -Oz -o Pm_monoclonal_exclude_MAF30.vcf.gz Pm_monoclonal_exclude_MAF25.vcf.gz

#$plink --pca --vcf Pm_monoclonal_exclude_MAF30.vcf.gz --const-fid --allow-extra-chr --out plink_MAF30

MAF30_eigenvalues <- data.table::fread("plink_MAF30.eigenval")
MAF30_eigenvectors <- data.table::fread("plink_MAF30.eigenvec")

MAF30_table <- data.frame(sample.id = MAF30_eigenvectors$V2, 
                          EV1 = MAF30_eigenvectors$V3,   
                          EV2 = MAF30_eigenvectors$V4, 
                          stringsAsFactors = FALSE)

MAF30_table <- MAF30_table |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                       stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                       stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                       .default = "Tanzania"))


plot(MAF30_table$EV2, MAF30_table$EV1, xlab="eigenvector 2", ylab="eigenvector 1",)

MAF30_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC1 (",signif(MAF30_eigenvalues[1], 4),"%)")) +
  ylab(paste0("PC2 (",signif(MAF30_eigenvalues[2], 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK MAF30.png", width = 15, height = 12, units = "in", dpi = 600)

#bcftools view -e 'MAF[0]<0.35' -Oz -o Pm_monoclonal_exclude_MAF35.vcf.gz Pm_monoclonal_exclude_MAF30.vcf.gz

#$plink --pca --vcf Pm_monoclonal_exclude_MAF35.vcf.gz --const-fid --allow-extra-chr --out plink_MAF35

MAF35_eigenvalues <- data.table::fread("plink_MAF35.eigenval")
MAF35_eigenvectors <- data.table::fread("plink_MAF35.eigenvec")

MAF35_table <- data.frame(sample.id = MAF35_eigenvectors$V2, 
                          EV1 = MAF35_eigenvectors$V3,   
                          EV2 = MAF35_eigenvectors$V4, 
                          stringsAsFactors = FALSE)

MAF35_table <- MAF35_table |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                       stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                       stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                       .default = "Tanzania"))


plot(MAF35_table$EV2, MAF35_table$EV1, xlab="eigenvector 2", ylab="eigenvector 1",)

MAF35_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC1 (",signif(MAF35_eigenvalues[1], 4),"%)")) +
  ylab(paste0("PC2 (",signif(MAF35_eigenvalues[2], 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK MAF35.png", width = 15, height = 12, units = "in", dpi = 600)

#bcftools view -e 'MAF[0]<0.4' -Oz -o Pm_monoclonal_exclude_MAF40.vcf.gz Pm_monoclonal_exclude_MAF35.vcf.gz

#$plink --pca --vcf Pm_monoclonal_exclude_MAF40.vcf.gz --const-fid --allow-extra-chr --out plink_MAF40

MAF40_eigenvalues <- data.table::fread("plink_MAF40.eigenval")
MAF40_eigenvectors <- data.table::fread("plink_MAF40.eigenvec")

MAF40_table <- data.frame(sample.id = MAF40_eigenvectors$V2, 
                          EV1 = MAF40_eigenvectors$V3,   
                          EV2 = MAF40_eigenvectors$V4, 
                          stringsAsFactors = FALSE)

MAF40_table <- MAF40_table |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                       stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                       stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                       .default = "Tanzania"))


plot(MAF40_table$EV2, MAF40_table$EV1, xlab="eigenvector 2", ylab="eigenvector 1",)

MAF40_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC1 (",signif(MAF40_eigenvalues[1], 4),"%)")) +
  ylab(paste0("PC2 (",signif(MAF40_eigenvalues[2], 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK MAF40.png", width = 15, height = 12, units = "in", dpi = 600)

#bcftools view -e 'MAF[0]<0.45' -Oz -o Pm_monoclonal_exclude_MAF45.vcf.gz Pm_monoclonal_exclude_MAF40.vcf.gz

#$plink --pca --vcf Pm_monoclonal_exclude_MAF45.vcf.gz --const-fid --allow-extra-chr --out plink_MAF45

MAF45_eigenvalues <- data.table::fread("plink_MAF45.eigenval")
MAF45_eigenvectors <- data.table::fread("plink_MAF45.eigenvec")

MAF45_table <- data.frame(sample.id = MAF45_eigenvectors$V2, 
                          EV1 = MAF45_eigenvectors$V3,   
                          EV2 = MAF45_eigenvectors$V4, 
                          stringsAsFactors = FALSE)

MAF45_table <- MAF45_table |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                       stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                       stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                       .default = "Tanzania"))


plot(MAF45_table$EV2, MAF45_table$EV1, xlab="eigenvector 2", ylab="eigenvector 1",)

MAF45_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC1 (",signif(MAF45_eigenvalues[1], 4),"%)")) +
  ylab(paste0("PC2 (",signif(MAF45_eigenvalues[2], 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK MAF45.png", width = 15, height = 12, units = "in", dpi = 600)

#bcftools view -e 'MAF[0]<0.49' -Oz -o Pm_monoclonal_exclude_MAF50.vcf.gz Pm_monoclonal_exclude_MAF45.vcf.gz

#$plink --pca --vcf Pm_monoclonal_exclude_MAF50.vcf.gz --const-fid --allow-extra-chr --out plink_MAF50

MAF49_eigenvalues <- data.table::fread("plink_MAF50.eigenval")
MAF49_eigenvectors <- data.table::fread("plink_MAF50.eigenvec")

MAF49_table <- data.frame(sample.id = MAF49_eigenvectors$V2, 
                          EV1 = MAF49_eigenvectors$V3,   
                          EV2 = MAF49_eigenvectors$V4, 
                          stringsAsFactors = FALSE)

MAF49_table <- MAF49_table |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                       stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                       stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                       .default = "Tanzania"))


plot(MAF49_table$EV2, MAF49_table$EV1, xlab="eigenvector 2", ylab="eigenvector 1",)

MAF49_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC1 (",signif(MAF49_eigenvalues[1], 4),"%)")) +
  ylab(paste0("PC2 (",signif(MAF49_eigenvalues[2], 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK MAF49.png", width = 15, height = 12, units = "in", dpi = 600)