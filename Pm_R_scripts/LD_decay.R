setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun")

system("plink --vcf Pm_HC_missingness_filtered_first.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.01 --geno 0.1 --mind 0.5 --r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --out Pm_LD")

system("plink --vcf Pf_ortholog_samples.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.01 --geno 0.1 --mind 0.5 --r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --out Pf_LD")

Pm_LD <- data.table::fread("Pm_LD.ld.gz")

Pm_LD <- Pm_LD |> subset(!stringr::str_detect(CHR_A, "archived")) |> subset(!stringr::str_detect(CHR_B, "archived")) |> subset(!stringr::str_detect(CHR_A, "MIT")) |> subset(!stringr::str_detect(CHR_B, "MIT")) |> subset(!stringr::str_detect(CHR_A, "API")) |> subset(!stringr::str_detect(CHR_B, "API"))

Pm_LD <- Pm_LD |> dplyr::mutate(Distance = BP_B - BP_A)

#Pm_unique_distances <- unique(Pm_LD$Distance)

Pm_decay <- Pm_LD |> dplyr::group_by(Distance) |> dplyr::summarise(Mean_R2 = mean(R2), Median_R2 = median(R2))

library(ggplot2)

Pm_decay |> ggplot() + geom_smooth(aes(x = Distance, y = Mean_R2))

Pm_decay |> ggplot() + geom_line(aes(x = Distance, y = Mean_R2))

Pm_decay |> ggplot() + geom_col(aes(x = Distance, y = Mean_R2))

Pm_decay |> ggplot() + geom_step(aes(x = Distance, y = Mean_R2))

Pm_decay |> ggplot() + geom_smooth(aes(x = Distance, y = Median_R2))

Pm_decay |> ggplot() + geom_line(aes(x = Distance, y = Median_R2))

Pm_decay |> ggplot() + geom_col(aes(x = Distance, y = Median_R2))

Pm_decay |> ggplot() + geom_step(aes(x = Distance, y = Median_R2))

Pf_LD <- data.table::fread("Pf_LD.ld.gz")

Pf_LD <- Pf_LD |> dplyr::mutate(Distance = BP_B - BP_A)

Pf_decay <- Pf_LD |> dplyr::group_by(Distance) |> dplyr::summarise(Mean_R2 = mean(R2), Median_R2 = median(R2))

Pf_decay |> ggplot() + geom_smooth(aes(x = Distance, y = Mean_R2))

Pf_decay |> ggplot() + geom_line(aes(x = Distance, y = Mean_R2))

Pf_decay |> ggplot() + geom_col(aes(x = Distance, y = Mean_R2))

Pf_decay |> ggplot() + geom_step(aes(x = Distance, y = Mean_R2))

Pf_decay |> ggplot() + geom_smooth(aes(x = Distance, y = Median_R2))

Pf_decay |> ggplot() + geom_line(aes(x = Distance, y = Median_R2))

Pf_decay |> ggplot() + geom_col(aes(x = Distance, y = Median_R2))

Pf_decay |> ggplot() + geom_step(aes(x = Distance, y = Median_R2))

Pm_decay$Species <- "Pm"

Pf_decay$Species <- "Pf"

combined_decay <- rbind(Pm_decay, Pf_decay)

combined_decay_plot <- combined_decay |> ggplot() + geom_smooth(aes(x = Distance, y = Mean_R2, color = Species)) + scale_color_brewer(palette = "Dark2") + theme_classic() + ggtitle("LD Decay by Species") + ylab(Mean~R^2)

species_t <- t.test(combined_decay$Mean_R2 ~ combined_decay$Species)

inset_decay <- combined_decay_plot + theme(legend.position = "none") + scale_x_continuous(breaks = c(2, 4, 6, 8, 10), limits = c(0, 10))

library(patchwork)

LD_with_inset <- combined_decay_plot + inset_element(inset_decay, left = 0.5, bottom = 0.5, right = 1, top = 1) & theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20), plot.title = element_blank(), legend.title = element_text(size = 24), legend.text = element_text(size = 20))

ggsave("LD Decay by Species.png", LD_with_inset, dpi = 600, height = 10, width = 8, units = "in")

combined_decay_plot2 <- combined_decay |> ggplot() + geom_smooth(aes(x = Distance, y = Median_R2, color = Species)) + scale_color_brewer(palette = "Dark2") + theme_classic() + ggtitle("LD Decay by Species") + ylab(Median~R^2)

inset_decay2 <- combined_decay_plot2 + theme(legend.position = "none") + scale_x_continuous(breaks = c(2, 4, 6, 8, 10), limits = c(0, 10))

LD_with_inset2 <- combined_decay_plot2 + inset_element(inset_decay2, left = 0.5, bottom = 0.5, right = 1, top = 1) & theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20), plot.title = element_blank(), legend.title = element_text(size = 24), legend.text = element_text(size = 20))

ggsave("LD Decay by Species Median Smooth.png", LD_with_inset2, dpi = 600, height = 10, width = 8, units = "in")

combined_decay_plot3 <- combined_decay |> ggplot() + geom_line(aes(x = Distance, y = Mean_R2, color = Species)) + scale_color_brewer(palette = "Dark2") + theme_classic() + ggtitle("LD Decay by Species") + ylab(Mean~R^2)

inset_decay3 <- combined_decay_plot3 + theme(legend.position = "none") + scale_x_continuous(breaks = c(2, 4, 6, 8, 10), limits = c(0, 10))

LD_with_inset3 <- combined_decay_plot3 + inset_element(inset_decay3, left = 0.5, bottom = 0.5, right = 1, top = 1) & theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20), plot.title = element_blank(), legend.title = element_text(size = 24), legend.text = element_text(size = 20))

ggsave("LD Decay by Species Mean Line.png", LD_with_inset3, dpi = 600, height = 10, width = 8, units = "in")

combined_decay_plot4 <- combined_decay |> ggplot() + geom_line(aes(x = Distance, y = Median_R2, color = Species)) + scale_color_brewer(palette = "Dark2") + theme_classic() + ggtitle("LD Decay by Species") + ylab(Median~R^2)

inset_decay4 <- combined_decay_plot4 + theme(legend.position = "none") + scale_x_continuous(breaks = c(2, 4, 6, 8, 10), limits = c(0, 10))

LD_with_inset4 <- combined_decay_plot4 + inset_element(inset_decay4, left = 0.5, bottom = 0.5, right = 1, top = 1) & theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20), plot.title = element_blank(), legend.title = element_text(size = 24), legend.text = element_text(size = 20))

ggsave("LD Decay by Species Median Line.png", LD_with_inset4, dpi = 600, height = 10, width = 8, units = "in")

combined_decay_plot5 <- combined_decay |> ggplot() + geom_col(aes(x = Distance, y = Mean_R2, color = Species)) + scale_color_brewer(palette = "Dark2") + theme_classic() + ggtitle("LD Decay by Species") + ylab(Mean~R^2)

inset_decay5 <- combined_decay_plot5 + theme(legend.position = "none") + scale_x_continuous(breaks = c(2, 4, 6, 8, 10), limits = c(0, 10))

LD_with_inset5 <- combined_decay_plot5 + inset_element(inset_decay5, left = 0.5, bottom = 0.5, right = 1, top = 1) & theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20), plot.title = element_blank(), legend.title = element_text(size = 24), legend.text = element_text(size = 20))

ggsave("LD Decay by Species Mean Bar.png", LD_with_inset5, dpi = 600, height = 10, width = 8, units = "in")

combined_decay_plot6 <- combined_decay |> ggplot() + geom_col(aes(x = Distance, y = Median_R2, color = Species)) + scale_color_brewer(palette = "Dark2") + theme_classic() + ggtitle("LD Decay by Species") + ylab(Median~R^2)

inset_decay6 <- combined_decay_plot6 + theme(legend.position = "none") + scale_x_continuous(breaks = c(2, 4, 6, 8, 10), limits = c(0, 10))

LD_with_inset6 <- combined_decay_plot6 + inset_element(inset_decay6, left = 0.5, bottom = 0.5, right = 1, top = 1) & theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20), plot.title = element_blank(), legend.title = element_text(size = 24), legend.text = element_text(size = 20))

ggsave("LD Decay by Species Median Bar.png", LD_with_inset6, dpi = 600, height = 10, width = 8, units = "in")

combined_decay_plot7 <- combined_decay |> ggplot() + geom_step(aes(x = Distance, y = Mean_R2, color = Species)) + scale_color_brewer(palette = "Dark2") + theme_classic() + ggtitle("LD Decay by Species") + ylab(Mean~R^2)

inset_decay7 <- combined_decay_plot7 + theme(legend.position = "none") + scale_x_continuous(breaks = c(2, 4, 6, 8, 10), limits = c(0, 10))

LD_with_inset7 <- combined_decay_plot7 + inset_element(inset_decay7, left = 0.5, bottom = 0.5, right = 1, top = 1) & theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20), plot.title = element_blank(), legend.title = element_text(size = 24), legend.text = element_text(size = 20))

ggsave("LD Decay by Species Mean Step.png", LD_with_inset7, dpi = 600, height = 10, width = 8, units = "in")

combined_decay_plot8 <- combined_decay |> ggplot() + geom_step(aes(x = Distance, y = Median_R2, color = Species)) + scale_color_brewer(palette = "Dark2") + theme_classic() + ggtitle("LD Decay by Species") + ylab(Median~R^2)

inset_decay8 <- combined_decay_plot8 + theme(legend.position = "none") + scale_x_continuous(breaks = c(2, 4, 6, 8, 10), limits = c(0, 10))

LD_with_inset8 <- combined_decay_plot8 + inset_element(inset_decay8, left = 0.5, bottom = 0.5, right = 1, top = 1) & theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20), plot.title = element_blank(), legend.title = element_text(size = 24), legend.text = element_text(size = 20))

ggsave("LD Decay by Species Median Step.png", LD_with_inset8, dpi = 600, height = 10, width = 8, units = "in")

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
