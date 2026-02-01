###############################################################
################### Pm_pi ##########################
###############################################################
#Description: plots nucloetide diversity values across Pm and Pf
#orthologus genes

setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun")

Pm_ortholog_pi <- data.table::fread("all_Pm_pi.txt")
Pm_ortholog_pi$BIN_START <- as.integer(Pm_ortholog_pi$BIN_START)
Pm_ortholog_pi$BIN_END <- as.integer(Pm_ortholog_pi$BIN_END)

masked_orthologs <- readxl::read_xlsx("Pf-Pm_masked_orthologs.xlsx")

colnames(masked_orthologs) <- c("Pm_CHROM", "Pm_START", "Pm_END", "Pm_LENGTH", "Pm_STRAND", "Group_ID", "Pm_ortho", "Pf_CHROM", "Pf_START", "Pf_END", "Pf_LENGTH", "Pf_STRAND", "Pf_ortho")

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

