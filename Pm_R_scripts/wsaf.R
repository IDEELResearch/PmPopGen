setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

#bcftools view Pm_monoclonals_missingness_only.vcf.gz -r PmUG01_01_v1,PmUG01_02_v1,PmUG01_03_v1,PmUG01_04_v1,PmUG01_05_v1,PmUG01_06_v1,PmUG01_07_v1,PmUG01_08_v1,PmUG01_09_v1,PmUG01_10_v1,PmUG01_11_v1,PmUG01_12_v1,PmUG01_13_v1,PmUG01_14_v1 | vcfdo wsaf --ignore-missing | bcftools view -Oz -o Pm_monoclonals_wsaf.vcf.gz

Pm <- vcfR::read.vcfR("Pm_monoclonals_wsaf.vcf.gz") |> vcfR::vcfR2tidy()

library(ggplot2)

wsaf_dist <- Pm$gt |> subset(gt_WSAF != -1) |> ggplot(aes(x = gt_WSAF)) + theme_classic() + geom_density(color = "black", alpha = 0.8)

wsaf_dist_zoom <- wsaf_dist + xlim(0.05, 0.95)

library(patchwork)

wsaf_dist_plots <- wsaf_dist + wsaf_dist_zoom & theme(labs(x = "WSAF"))

ggsave("wsaf_dist.png", wsaf_dist_plots)

wsaf_hist <- data.table::fread("Pm_monclonals_wsaf_hist.txt")

colnames(wsaf_hist) <- c("Sample", "Bin", "Frequency")

wsaf_hist2 <- wsaf_hist |> dplyr::group_by(Sample) |> dplyr::mutate(Sites = sum(Frequency), Prop = Frequency/Sites)

wsaf_hist2 |> ggplot(aes(x = Bin, y = Prop)) + theme_classic() + geom_col()

#essentially I want a per-site metric for how many samples are at each threshold i.e. at each site, how many samples do we lose if we cut it out

#but maybe also just easier to try some scaling WSAF cutoffs and see how many sites we lose?

wsaf_hist_summary <- tidytable::summarize(.by = Bin, Sample_Count = tidytable::n_distinct())

#wsaf_summary <- wsaf_hist |> dplyr::group_by(Bin) |> dplyr::summarize(Mean_Freq = mean(Frequency))

#wsaf_summary |> ggplot(aes(x = Bin, y= Mean_Freq)) + theme_classic() + geom_col()

wsaf_summary <- data.frame(Pm$fix$CHROM, Pm$fix$POS, Pm$gt$Indiv, Pm$gt$gt_WSAF)

colnames(wsaf_summary) <- c("CHROM", "POS", "Sample", "WSAF")

wsaf_summary$WSAF <- wsaf_summary$WSAF |> dplyr::na_if(-1)

wsaf_bins <- wsaf_summary |> dplyr::mutate(WSAF_bin = cut(WSAF, breaks = c (0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1), include.lowest = TRUE))

wsaf_bins_summary <- wsaf_bins |> tidytable::summarize(.by = c(CHROM, POS, WSAF_bin), count = tidytable::n_distinct(Sample)) |> subset(WSAF_bin != "NA") |> dplyr::group_by(WSAF_bin) |> dplyr::summarize(total_count = sum(count), per_site_count = sum(count)/dplyr::n_distinct(CHROM, POS))

wsaf_bins_summary <- wsaf_bins_summary |> dplyr::mutate(normalized_count = total_count/sum(total_count))

barplot(wsaf_bins_summary$normalized_count)

wsaf_bins_summary |> dplyr::group_by(WSAF_bin) |> dplyr::summarize(total_count = sum(count))

wsaf_summary2 <- wsaf_summary |> tidytable::summarize(.by = c(CHROM, POS, WSAF), count = tidytable::n_distinct(Sample)) |> subset(WSAF != "NA")

wsaf_summary3 <- wsaf_summary2 |> dplyr::mutate(CHROM_POS = paste0(CHROM, "_", POS)) |> dplyr::select(CHROM_POS, WSAF, count)

wsaf_summary3 |> dplyr::mutate(row = dplyr::row_number()) |> tidyr::pivot_wider(names_from = c(CHROM_POS, row), values_from = c(WSAF, count))

wsaf_plot2 <- wsaf_summary2 |> ggplot(aes(x = WSAF, y = count)) + theme_classic() + geom_point()

wsaf_summary3 <- wsaf_summary2 |> dplyr::group_by(WSAF) |> dplyr::summarize(avg_samples = mean(count))

wsaf_plot2 <- wsaf_summary3 |> ggplot(aes(x = WSAF, y = avg_samples)) + theme_classic() + geom_point()

wsaf_plot <- wsaf_bins_summary |> ggplot(aes(x = WSAF_bin, y = count)) + theme_classic() + geom_col()

wsaf_scatter <- wsaf_bins_summary |> ggplot(aes(x = WSAF_bin, y = count)) + theme_classic() + geom_point()
