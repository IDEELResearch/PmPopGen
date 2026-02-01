###############################################################
###################### Pm_PCA #############################
###############################################################
#Description: performs and graphs principal components analysis on monoclonal P. malariae isolates

setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

#subsets to site by missingess
system("bcftools view -e 'F_MISSING>=0.02' -Oz -o Pm_PCA.vcf.gz Pm_IBD_pruned.vcf.gz")

#uses plink to calculate PCA across all filtered sites
system("plink --pca var-wts --vcf Pm_PCA.vcf.gz --const-fid --allow-extra-chr --out plink_IBD_pruned --make-rel")

system("awk '{print $NR}' plink_IBD_pruned.rel > plink_IBD_pruned.rel.diag")

#extracts eigenvalues for Scree plotting
Pm_eigenvalues <- data.table::fread("plink_IBD_pruned.eigenval")

#extracts eigenvectors for PCA plotting
Pm_eigenvectors <- data.table::fread("plink_IBD_pruned.eigenvec")

rel_diag <- data.table::fread("plink_IBD_pruned.rel.diag")

Pm_eigenvalues <- Pm_eigenvalues |> dplyr::mutate(prop_var = V1/sum(rel_diag), prop_var2 = V1/sum(V1))

Pm_table <- data.frame(sample.id = Pm_eigenvectors$V2,
                       EV1 = Pm_eigenvectors$V3,
                       EV2 = Pm_eigenvectors$V4,
                       EV3 = Pm_eigenvectors$V5,
                       EV4 = Pm_eigenvectors$V6,
                       EV5 = Pm_eigenvectors$V7,
                       EV6 = Pm_eigenvectors$V8,
                       stringsAsFactors = FALSE)

Pm_table <- Pm_table |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                 stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                 stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                 .default = "Tanzania"))

library(ggplot2)
library(viridis)

#plot PCA for PC1 and PC2 of Pm samples
Pm_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC1 (",signif(Pm_eigenvalues$prop_var[1]*100, 4),"%)")) +
  ylab(paste0("PC2 (",signif(Pm_eigenvalues$prop_var[2]*100, 4),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA.png", width = 15, height = 12, units = "in", dpi = 600)

#generate scree plot
scree_df <- Pm_eigenvalues |> dplyr::rename(Eigenvalue = V1) |> tibble::rownames_to_column(var = "PC")

scree_df$PC <- as.numeric(scree_df$PC)

scree_plot <- scree_df |> ggplot() + geom_line(aes(x = PC, y = prop_var)) + geom_point(aes(x = PC, y = prop_var)) + labs(x = "Principal Component", y = "Variance Explained") + theme_classic()

ggsave("scree_plot.png", scree_plot, dpi = 600)

#plot PC3 and PC4
PC34 <- Pm_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV3, y = EV4, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC3 (",signif(Pm_eigenvalues$prop_var[3]*100, 2),"%)")) +
  ylab(paste0("PC4 (",signif(Pm_eigenvalues$prop_var[4]*100, 2),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

#Plot PC5 and PC6
PC56 <- Pm_table |> ggplot() + theme_bw() +
  geom_point(aes(x = EV3, y = EV4, color = factor(Country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab(paste0("PC5 (",signif(Pm_eigenvalues$prop_var[5]*100, 2),"%)")) +
  ylab(paste0("PC6 (",signif(Pm_eigenvalues$prop_var[6]*100, 2),"%)")) +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

library(patchwork)

extra_PCs <- PC34 + PC56

ggsave("Extra PCs.png", extra_PCs, width = 30, height = 12, units = "in", dpi = 600)

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


