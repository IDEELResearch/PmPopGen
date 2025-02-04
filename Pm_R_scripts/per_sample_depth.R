setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

Pm_VCF <- vcfR::read.vcfR("Pm_HC_missingness_filtered_first.vcf")

Pm_maf <- Pm_VCF |> vcfR::maf(element = 2)

Pm_maf_filtered <- Pm_maf |> as.data.frame() |> subset(Frequency >= 0.01)

Pm_maf_filtered <- Pm_maf_filtered |> tibble::rownames_to_column(var = "CHROM_POS")

Pm_maf_filtered <- Pm_maf_filtered |> dplyr::mutate(CHROM = dplyr::case_when(stringr::str_detect(CHROM_POS, "archived") ~ stringr::str_extract(CHROM_POS, "PmUG01_00_v1_archived_contig_[:digit:]+"),
                                                                             .default = stringr::str_extract(CHROM_POS, "PmUG01_[:alnum:]+_v1")))

Pm_maf_filtered <- Pm_maf_filtered |> dplyr::mutate(POS = stringr::str_extract(CHROM_POS, "[:digit:]+$"))

Pm_maf_filtered$POS <- Pm_maf_filtered$POS |> as.numeric()

Pm_VCF_tidy <- Pm_VCF |> vcfR::vcfR2tidy()

Pm_filtered_SNPs <- Pm_VCF_tidy$fix |> subset(FILTER == "PASS")

Pm_filtered_SNPs <- dplyr::semi_join(Pm_filtered_SNPs, Pm_maf_filtered)

Pm_filtered_SNPs_stats <- dplyr::semi_join(Pm_VCF_tidy$gt, Pm_filtered_SNPs)

Pm_per_sample_DP <- Pm_filtered_SNPs_stats |> dplyr::group_by(Indiv) |> dplyr::summarize(mean_DP = mean(gt_DP, na.rm = T), median_DP = median(gt_DP, na.rm = T))

Pm_per_sample_genotypes <-Pm_filtered_SNPs_stats |> dplyr::group_by(Indiv) |> dplyr::count(gt_GT)

Pm_GT_NA_per_sample <- Pm_per_sample_genotypes |> subset(is.na(gt_GT) == TRUE)

Pm_GT_SNPs_per_sample <- Pm_per_sample_genotypes |> subset(gt_GT != "0/0" & gt_GT != "0|0" & is.na(gt_GT) == FALSE) |> dplyr::summarize(alt_GT = sum(n))

#Pm_DP <- Pm_VCF |> vcfR::extract.gt(element='DP', as.numeric=TRUE)

#Pm_DP <- Pm_DP |> reshape2::melt(varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)

library(ggplot2)

DP_box <- Pm_filtered_SNPs_stats |> dplyr::filter(is.na(gt_DP) == FALSE) |> ggplot() + geom_boxplot(aes(x = reorder(Indiv, -gt_DP, FUN = median), y = gt_DP)) +scale_y_continuous(breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 300)) + theme_classic() + theme(axis.text.x = element_blank(), axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(), text = element_text(size = 14), title = element_text(size = 24)) + xlab("Sample") + ylab("Depth") + scale_y_break(c(10,20), scales = 0.5)

ggsave("depth_boxplot.png", dpi = 600, width = 10, height = 10, units = "in")

#DP_violin <- Pm_filtered_SNPs_stats |> ggplot() + geom_violin(aes(x = reorder(Indiv, -gt_DP), y = gt_DP)) + theme_classic() + theme(axis.text.x = element_blank()) + xlab("Sample") + ylab("Depth")

#ggsave("depth_violin.png", dpi = 600, width = 10, height = 10, units = "in")

DP_box_2 <- DP_box + ylim(0, 50)

ggsave("depth_boxplot_zoom.png", dpi = 600, width = 10, height = 10, units = "in")

#DP_violin_2 <- DP_violin + ylim(0,50)

#ggsave("depth_violin_zoom.png", dpi = 600, width = 10, height = 10, units = "in")

DP_box_3 <- DP_box + ylim(0, 20)

ggsave("depth_boxplot_zoom20.png", dpi = 600, width = 10, height = 10, units = "in")

#DP_violin_3 <- DP_violin + ylim(0,20)

#ggsave("depth_violin_zoom20.png", dpi = 600, width = 10, height = 10, units = "in")

DP_box_4 <- DP_box + ylim(0, 10)

ggsave("depth_boxplot_zoom10.png", dpi = 600, width = 10, height = 10, units = "in")

#DP_violin_4 <- DP_violin + ylim(0, 10)

#ggsave("depth_violin_zoom10.png", dpi = 600, width = 10, height = 10, units = "in")

Pm_filtered_SNPs_stats_no_zeros <- Pm_filtered_SNPs_stats |> dplyr::filter(gt_DP > 0)

DP_box_5 <- Pm_filtered_SNPs_stats_no_zeros  |> ggplot() + geom_boxplot(aes(x = reorder(Indiv, -gt_DP, FUN = median), y = gt_DP)) + theme_classic() + theme(axis.text.x = element_blank()) + xlab("Sample") + ylab("Depth")

ggsave("depth_boxplot_no_zeros.png", dpi = 600, width = 10, height = 10, units = "in")

DP_box_6 <- DP_box_5 + ylim(0, 50)

ggsave("depth_boxplot_zoom_no_zeros.png", dpi = 600, width = 10, height = 10, units = "in")

DP_box_7 <- DP_box_5 + ylim(0, 20)

ggsave("depth_boxplot_zoom20_no_zeros.png", dpi = 600, width = 10, height = 10, units = "in")

DP_box_8 <- DP_box_5 + ylim(0,10)

ggsave("depth_boxplot_zoom10_no_zeros.png", dpi = 600, width = 10, height = 10, units = "in")

individ_dist <- Pm_filtered_SNPs_stats |> dplyr::group_by(Indiv) |> dplyr::count(gt_DP) |> dplyr::mutate(Prop = n/sum(n))

depth_barplot <- individ_dist |> ggplot() + geom_bar(aes(x = reorder(Indiv, -gt_DP, FUN = median), y = n, fill = gt_DP), position = "fill", stat = "identity") + theme_classic() + scale_fill_distiller(palette = "Dark2", name = "DP") + labs(x = "Sample", y = "Proportion") + theme(axis.text.x = element_blank(), text = element_text(size = 14), title = element_text(size = 24))

ggsave("depth_barplot.png", dpi = 600, width = 10, height = 10, units = "in")

depth_barplot_under100 <- individ_dist |> dplyr::filter(gt_DP < 100) |> ggplot() + geom_bar(aes(x = reorder(Indiv, -gt_DP, FUN = median), y = n, fill = gt_DP), position = "fill", stat = "identity") + theme_classic() + scale_fill_distiller(palette = "Dark2", name = "DP") + labs(x = "Sample", y = "Proportion") + theme(axis.text.x = element_blank())

ggsave("depth_barplot_under100.png", dpi = 600, width = 10, height = 10, units = "in")

depth_barplot_under50 <- individ_dist |> dplyr::filter(gt_DP < 50) |> ggplot() + geom_bar(aes(x = reorder(Indiv, -gt_DP, FUN = median), y = n, fill = gt_DP), position = "fill", stat = "identity") + theme_classic() + scale_fill_distiller(palette = "Dark2", name = "DP") + labs(x = "Sample", y = "Proportion") + theme(axis.text.x = element_blank())

ggsave("depth_barplot_under50.png", dpi = 600, width = 10, height = 10, units = "in")

depth_barplot_under10 <- individ_dist |> dplyr::filter(gt_DP < 10) |> ggplot() + geom_bar(aes(x = reorder(Indiv, -gt_DP, FUN = median), y = n, fill = gt_DP), position = "fill", stat = "identity") + theme_classic() + scale_fill_distiller(palette = "Dark2", name = "DP", breaks = c(0:10)) + labs(x = "Sample", y = "Proportion") + theme(axis.text.x = element_blank(),text = element_text(size = 14), title = element_text(size = 24))

ggsave("depth_barplot_under10.png", dpi = 600, width = 10, height = 10, units = "in")

library(patchwork)

Pm_depth_fig <- DP_box + (depth_barplot / depth_barplot_under10) + plot_annotation(tag_levels = "A") + theme(plot.tag = element_text(size = 24))

ggsave("Pm_depth_fig.png", dpi = 600, width = 20, height = 12, units = "in")

#Pf_VCF <- vcfR::read.vcfR("Pf_depth.vcf.gz") |> vcfR::vcfR2tidy()

#Pf_no_NA <- Pf_VCF$gt |> dplyr::filter(is.na(gt_DP) == FALSE)

#Pf file is too big so we're using an alternative strategy

system("/bin/bcftools query -f ['%CHROM\t%POS0\t%SAMPLE\t%DP\n'] Pf_ortholog_samples_MAF1.vcf.gz > Pf_DP.txt")

Pf_DP <- data.table::fread("Pf_DP.txt", header = F)

colnames(Pf_DP) <- c("CHROM", "POS", "SAMPLE", "DP")

Pf_DP$DP <- as.numeric(Pf_DP$DP)

Pf_DP <- Pf_DP |> subset(is.na(DP) == FALSE)

#test <- Pf_DP |> head()

#test_wider <- test |> tidyr::pivot_wider(names_from = c(CHROM, POS), values_from = DP)

#Pf_DP_wider <- Pf_DP |> tidyr::pivot_wider(names_from = c(CHROM, POS), values_from = DP)

#Pf_DP_longer <- Pf_DP_wider |> tidyr::pivot_longer()

#Pf_DP_wider |> saveRDS("Pf_DP_wider")

#Pf_DP_wider <- readRDS("Pf_DP_wider")

#Pf_DP_wider |> data.table::fwrite("Pf_DP_wider.txt")

library(ggplot2)

DP_box_Pf <- Pf_DP |> ggplot() + geom_boxplot(aes(x = reorder(SAMPLE, -DP, FUN = median), y = DP)) + theme_classic() + theme(axis.text.x = element_blank()) + xlab("Sample") + ylab("Depth")

ggsave("depth_boxplot_Pf.png", dpi = 600, width = 10, height = 10, units = "in")

#DP_violin <- Pm_filtered_SNPs_stats |> ggplot() + geom_violin(aes(x = reorder(Indiv, -gt_DP), y = gt_DP)) + theme_classic() + theme(axis.text.x = element_blank()) + xlab("Sample") + ylab("Depth")

#ggsave("depth_violin.png", dpi = 600, width = 10, height = 10, units = "in")

DP_box_2_Pf <- DP_box_Pf + ylim(0, 50)

ggsave("depth_boxplot_zoom_Pf.png", dpi = 600, width = 10, height = 10, units = "in")

#DP_violin_2 <- DP_violin + ylim(0,50)

#ggsave("depth_violin_zoom.png", dpi = 600, width = 10, height = 10, units = "in")

DP_box_3_Pf <- DP_box_Pf + ylim(0, 20)

ggsave("depth_boxplot_zoom20_Pf.png", dpi = 600, width = 10, height = 10, units = "in")

#DP_violin_3 <- DP_violin + ylim(0,20)

#ggsave("depth_violin_zoom20.png", dpi = 600, width = 10, height = 10, units = "in")

DP_box_4_Pf <- DP_box_Pf + ylim(0, 10)

ggsave("depth_boxplot_zoom10_Pf.png", dpi = 600, width = 10, height = 10, units = "in")

#DP_violin_4 <- DP_violin + ylim(0, 10)

#ggsave("depth_violin_zoom10.png", dpi = 600, width = 10, height = 10, units = "in")

Pf_no_zeros <- Pf_VCF$gt |> dplyr::filter(gt_DP > 0)

DP_box_5_Pf <- Pf_no_zeros  |> ggplot() + geom_boxplot(aes(x = reorder(Indiv, -gt_DP, FUN = median), y = gt_DP)) + theme_classic() + theme(axis.text.x = element_blank()) + xlab("Sample") + ylab("Depth")

ggsave("depth_boxplot_no_zeros_Pf.png", dpi = 600, width = 10, height = 10, units = "in")

DP_box_6_Pf <- DP_box_5_Pf + ylim(0, 50)

ggsave("depth_boxplot_zoom_no_zeros_Pf.png", dpi = 600, width = 10, height = 10, units = "in")

DP_box_7_Pf <- DP_box_5_Pf + ylim(0, 20)

ggsave("depth_boxplot_zoom20_no_zeros_Pf.png", dpi = 600, width = 10, height = 10, units = "in")

DP_box_8_Pf <- DP_box_5_Pf + ylim(0,10)

ggsave("depth_boxplot_zoom10_no_zeros_Pf.png", dpi = 600, width = 10, height = 10, units = "in")

individ_dist_Pf <- Pf_VCF$gt |> dplyr::group_by(Indiv) |> dplyr::count(gt_DP) |> dplyr::mutate(Prop = n/sum(n))

depth_barplot_Pf <- individ_dist_Pf |> ggplot() + geom_bar(aes(x = Indiv, y = n, fill = gt_DP), position = "fill", stat = "identity") + theme_classic() + scale_fill_distiller(palette = "Dark2", name = "DP") + labs(x = "Sample", y = "Proportion") + theme(axis.text.x = element_blank())

ggsave("depth_barplot_Pf.png", dpi = 600, width = 10, height = 10, units = "in")

depth_barplot_under100_Pf <- individ_dist_Pf |> dplyr::filter(gt_DP < 100) |> ggplot() + geom_bar(aes(x = Indiv, y = n, fill = gt_DP), position = "fill", stat = "identity") + theme_classic() + scale_fill_distiller(palette = "Dark2", name = "DP") + labs(x = "Sample", y = "Proportion") + theme(axis.text.x = element_blank())

ggsave("depth_barplot_under100_Pf.png", dpi = 600, width = 10, height = 10, units = "in")

depth_barplot_under50_Pf <- individ_dist_Pf |> dplyr::filter(gt_DP < 50) |> ggplot() + geom_bar(aes(x = Indiv, y = n, fill = gt_DP), position = "fill", stat = "identity") + theme_classic() + scale_fill_distiller(palette = "Dark2", name = "DP") + labs(x = "Sample", y = "Proportion") + theme(axis.text.x = element_blank())

ggsave("depth_barplot_under50_Pf.png", dpi = 600, width = 10, height = 10, units = "in")

depth_barplot_under10_Pf <- individ_dist_Pf |> dplyr::filter(gt_DP < 10) |> ggplot() + geom_bar(aes(x = Indiv, y = n, fill = gt_DP), position = "fill", stat = "identity") + theme_classic() + scale_fill_distiller(palette = "Dark2", name = "DP", breaks = c(0:10)) + labs(x = "Sample", y = "Proportion") + theme(axis.text.x = element_blank())

ggsave("depth_barplot_under10_Pf.png", dpi = 600, width = 10, height = 10, units = "in")
