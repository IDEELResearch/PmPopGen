setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun")

#/home/zpopkinh/plink --vcf Pm_HC_missingness_filtered_first.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.01 --geno 0.1 --mind 0.5 --r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --out Pm_LD

#/home/zpopkinh/plink --vcf Pf_ortholog_samples.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.01 --geno 0.1 --mind 0.5 --r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --out Pf_LD

Pm_LD <- data.table::fread("Pm_LD.ld.gz")

Pm_LD <- Pm_LD |> subset(!stringr::str_detect(CHR_A, "archived")) |> subset(!stringr::str_detect(CHR_B, "archived")) |> subset(!stringr::str_detect(CHR_A, "MIT")) |> subset(!stringr::str_detect(CHR_B, "MIT")) |> subset(!stringr::str_detect(CHR_A, "API")) |> subset(!stringr::str_detect(CHR_B, "API"))

Pm_LD <- Pm_LD |> dplyr::mutate(Distance = BP_B - BP_A)

#Pm_unique_distances <- unique(Pm_LD$Distance)

Pm_decay <- Pm_LD |> dplyr::group_by(Distance) |> dplyr::summarise(Mean_R2 = mean(R2))

library(ggplot2)

Pm_decay |> ggplot() + geom_smooth(aes(x = Distance, y = Mean_R2))

Pf_LD <- data.table::fread("Pf_LD.ld.gz")

Pf_LD <- Pf_LD |> dplyr::mutate(Distance = BP_B - BP_A)

Pf_decay <- Pf_LD |> dplyr::group_by(Distance) |> dplyr::summarise(Mean_R2 = mean(R2))

Pf_decay |> ggplot() + geom_smooth(aes(x = Distance, y = Mean_R2))

Pm_decay$Species <- "Pm"

Pf_decay$Species <- "Pf"

combined_decay <- rbind(Pm_decay, Pf_decay)

combined_decay_plot <- combined_decay |> ggplot() + geom_smooth(aes(x = Distance, y = Mean_R2, color = Species)) + scale_color_brewer(palette = "Dark2") + theme_classic() + ggtitle("LD Decay by Species") + ylab(Mean~R^2)

species_t <- t.test(combined_decay$Mean_R2 ~ combined_decay$Species)

inset_decay <- combined_decay_plot + theme(legend.position = "none") + scale_x_continuous(breaks = c(2, 4, 6, 8, 10), limits = c(0, 10))

library(patchwork)

LD_with_inset <- combined_decay_plot + inset_element(inset_decay, left = 0.5, bottom = 0.5, right = 1, top = 1) & theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20), plot.title = element_blank(), legend.title = element_text(size = 24), legend.text = element_text(size = 20))

ggsave("LD Decay by Species.png", LD_with_inset, dpi = 600)

pi_plot <- readRDS("pi_plot.rds")

Fig2 <- pi_plot + LD_with_inset + plot_annotation(tag_levels = list(c("A", "B"))) & theme(plot.tag = element_text(size = 24, family = "bold"))

ggsave("Fig2.png", Fig2, dpi = 600, height = 10, width = 15, units = "in")

##############0.05 MAF threshold

Pm_LD_0.05 <- data.table::fread("Pm_LD_0.05.ld.gz")

Pm_LD_0.05 <- Pm_LD_0.05 |> subset(!stringr::str_detect(CHR_A, "archived")) |> subset(!stringr::str_detect(CHR_B, "archived")) |> subset(!stringr::str_detect(CHR_A, "MIT")) |> subset(!stringr::str_detect(CHR_B, "MIT")) |> subset(!stringr::str_detect(CHR_A, "API")) |> subset(!stringr::str_detect(CHR_B, "API"))

Pm_LD_0.05 <- Pm_LD_0.05 |> dplyr::mutate(Distance = BP_B - BP_A)

#Pm_unique_distances <- unique(Pm_LD$Distance)

Pm_decay_0.05 <- Pm_LD_0.05 |> dplyr::group_by(Distance) |> dplyr::summarise(Mean_R2 = mean(R2))

library(ggplot2)

Pm_decay_0.05 |> ggplot() + geom_smooth(aes(x = Distance, y = Mean_R2))

Pf_LD_0.05 <- data.table::fread("Pf_LD_0.05.ld.gz")

Pf_LD_0.05 <- Pf_LD_0.05 |> dplyr::mutate(Distance = BP_B - BP_A)

Pf_decay_0.05 <- Pf_LD_0.05 |> dplyr::group_by(Distance) |> dplyr::summarise(Mean_R2 = mean(R2))

Pf_decay_0.05 |> ggplot() + geom_smooth(aes(x = Distance, y = Mean_R2))

Pm_decay_0.05$Species <- "Pm"

Pf_decay_0.05$Species <- "Pf"

combined_decay_0.05 <- rbind(Pm_decay_0.05, Pf_decay_0.05)

combined_decay_plot_0.05 <- combined_decay_0.05 |> ggplot() + geom_smooth(aes(x = Distance, y = Mean_R2, color = Species)) + scale_color_brewer(palette = "Dark2") + theme_classic() + ggtitle("LD Decay by Species") + ylab(Mean~R^2)

species_t_0.05 <- t.test(combined_decay_0.05$Mean_R2 ~ combined_decay_0.05$Species)

inset_decay_0.05 <- combined_decay_plot_0.05 + theme(legend.position = "none") + scale_x_continuous(breaks = c(2, 4, 6, 8, 10), limits = c(0, 10))

library(patchwork)

LD_with_inset_0.05 <- combined_decay_plot_0.05 + inset_element(inset_decay_0.05, left = 0.5, bottom = 0.5, right = 1, top = 1) & theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20), plot.title = element_blank(), legend.title = element_text(size = 24), legend.text = element_text(size = 20))

ggsave("LD Decay by Species 0.05.png", LD_with_inset_0.05, dpi = 600)

pi_plot <- readRDS("pi_plot.rds")

Fig2_0.05 <- pi_plot + LD_with_inset_0.05 + plot_annotation(tag_levels = list(c("A", "B"))) & theme(plot.tag = element_text(size = 24, family = "bold"))

ggsave("Fig2_0.05.png", Fig2_0.05, dpi = 600, height = 10, width = 15, units = "in")
