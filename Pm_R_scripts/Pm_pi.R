setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun")

Pm_ortholog_pi <- data.table::fread("all_Pm_pi.txt")
Pm_ortholog_pi$BIN_START <- as.integer(Pm_ortholog_pi$BIN_START)
Pm_ortholog_pi$BIN_END <- as.integer(Pm_ortholog_pi$BIN_END)

masked_orthologs <- readxl::read_xlsx("Pf-Pm_masked_orthologs.xlsx")

Pm_masked_orthos <- masked_orthologs |> dplyr::select(Group_ID, Pm_CHROM, Pm_START, Pm_END, Pm_LENGTH, Pm_ortho) |> dplyr::rename(CHROM = Pm_CHROM, BIN_START = Pm_START, BIN_END = Pm_END)

Pm_pi_orthos_only <- dplyr::left_join(Pm_masked_orthos,Pm_ortholog_pi)

Pf_ortholog_pi <- data.table::fread("all_Pf_pi.txt")

Pf_ortholog_pi$BIN_START <- as.integer(Pf_ortholog_pi$BIN_START)
Pf_ortholog_pi$BIN_END <- as.integer(Pf_ortholog_pi$BIN_END)

Pf_masked_orthos <- masked_orthologs |> dplyr::select(Group_ID, Pf_CHROM, Pf_START, Pf_END, Pf_LENGTH, Pf_ortho) |> dplyr::rename(CHROM = Pf_CHROM, BIN_START = Pf_START, BIN_END = Pf_END)

Pf_pi_orthos_only <- dplyr::left_join(Pf_masked_orthos, Pf_ortholog_pi)

Pm_ortho_pi <- Pm_pi_orthos_only |> dplyr::rename(Pm_CHROM = CHROM, Pm_START = BIN_START, Pm_END = BIN_END, Pm_SNPs = N_VARIANTS, Pm_PI = PI)

Pf_ortho_pi <- Pf_pi_orthos_only |> dplyr::rename(Pf_CHROM = CHROM, Pf_START = BIN_START, Pf_END = BIN_END, Pf_SNPs = N_VARIANTS, Pf_PI = PI)

combined_ortho_pi <- dplyr::left_join(Pm_ortho_pi, Pf_ortho_pi)

combined_ortho_pi |> writexl::write_xlsx("Pm_Pf_Ortholog_pi.xlsx")

missing_orthos <- combined_ortho_pi |> subset(is.na(Pm_PI) == TRUE | is.na(Pf_PI) == TRUE)

pi_df <- combined_ortho_pi |> dplyr::select(Group_ID, Pm_PI, Pf_PI) |> reshape2::melt(id = "Group_ID")

pi_df$value <- as.numeric(pi_df$value)

pi_df <- pi_df |> dplyr::mutate(log_pi = log(value))

pi_plot <- pi_df |> ggplot () +
  geom_boxplot(aes(x = variable, y = value))

pi_violin <- pi_df |> ggplot() +
  geom_violin(aes(x = variable, y = value))

pi_log <- pi_df |> ggplot() +
  geom_boxplot(aes(x = variable, y = value, fill = variable)) + theme_classic() +
  scale_y_log10() + scale_x_discrete(labels = c(expression(italic("P. malariae"), italic("P. falciparum")))) +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 20)) + labs(y = expression(paste("Nucleotide Diversity (", pi, ")")))

saveRDS(pi_log, file = "pi_plot.rds")

ggsave("Pf_Pm_ortholog_pi.png", pi_log, width = 10, height = 10, units = "in", dpi = 600)

pi_log_violin <- pi_df |> ggplot() +
  geom_violin(aes(x = variable, y = log_pi, color = variable), draw_quantiles = c(0.25, 0.5, 0.75)) + geom_jitter(aes(x = variable, y = log_pi, color = variable))

species_t <- t.test(pi_df$value ~ pi_df$variable)

#nuc_div <- data.table::fread("out.sites.pi", header = T)

#library(qqman)

#nuc_div$SNP <- row.names(nuc_div)

#manhattan(nuc_div, chr = "CHROM", bp = "POS", p = "PI", ylim = c(min(nuc_div$PI), max(nuc_div$PI)), logp = FALSE, ylab = expression(pi), ) 

#nuc_div$CHROM <- nuc_div$CHROM |> stringr::str_replace(pattern = "15", replacement = "API")

#library(viridis)
#palette(turbo(15))

#nuc_div$CHROM <- factor(nuc_div$CHROM, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "API"))

#nuc_div |> 
#  ggplot() +
#  geom_violin(aes(x = CHROM, y = PI, color = CHROM, fill = CHROM)) +
#  #geom_jitter(aes(x = Chromosome, y = coverage, color = Chromosome), alpha = 0.3, size = 0.5, height=0.1) +
#  theme_linedraw() +
#  theme(legend.position = "none") + 
  #ylim(0,10) +
#  labs(y = expression(pi), x = "Chromosome")

#ggsave("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/pi_violin.png", device = "png", dpi #= 600, width = 12, height = 10, units = "in")

#nuc_div |> 
#  ggplot() +
#  geom_boxplot(aes(x = CHROM, y = PI, color = CHROM)) +
#  #geom_jitter(aes(x = Chromosome, y = coverage, color = Chromosome), alpha = 0.3, size = 0.5, height=0.1) +
#  theme_linedraw() +
#  theme(legend.position = "none") + 
#  theme(axis.text.y = element_text(size = 30),
#        axis.title.y = element_text(size = 30),
#        axis.title.x = element_text(size = 30),
#        plot.title = element_text(size = 30),
#        axis.text.x = element_text(size = 30)) +
  #ylim(0,10) +
#  labs(y = expression(pi), x = "Chromosome") +
#  ggtitle(expression(paste("Per-Site Nucleotide Diversity ", "(", pi, ,")", " within Tanzanian ", italic("P. malariae"), " Isolates")))

#ggsave("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/pi_boxplot.png", device = "png", dpi = 600, width = 15, height = 8, units = "in")

#summary(nuc_div$PI)

#pi_chrom_aov <- aov(nuc_div$PI ~ nuc_div$CHROM)

#summary(pi_chrom_aov)

#TukeyHSD(pi_chrom_aov)

#pi_stats_by_chrom <- nuc_div |> dplyr::group_by(CHROM) |> dplyr::summarise(Average = mean(PI)) #may add more later

#max(pi_stats_by_chrom$Average)
#min(pi_stats_by_chrom$Average)
