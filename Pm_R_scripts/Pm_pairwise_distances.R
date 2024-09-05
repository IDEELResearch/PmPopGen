setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

inputVCF <- ("Pm_HC_missingness_filtered_first.vcf.gz")

Pm_distances <- fastreeR::vcf2dist(inputVCF)

Pm_distances_df <- Pm_distances |> as.matrix() |> reshape2::melt(varnames = c("Sample1", "Sample2")) |> subset(Sample1 != Sample2)

Pm_distances_df <- Pm_distances_df |> dplyr::mutate(prefix1 = dplyr::case_when(stringr::str_detect(Sample1, "^[:digit:]+\\_[:digit:]") ~ stringr::str_extract(Sample1, "^[:digit:]+"))) |> dplyr::mutate(prefix2 = dplyr::case_when(stringr::str_detect(Sample2, "^[:digit:]+\\_[:digit:]") ~ stringr::str_extract(Sample2, "^[:digit:]+"))) |> dplyr::mutate(paired = dplyr::case_when(prefix1 == prefix2 ~ TRUE,
                                                                                                                                                                                                                                                                                                                                                                                        .default = FALSE))

Pm_dist_density <- density(Pm_distances_df$value)

get_y <- approxfun(Pm_dist_density$x, Pm_dist_density$y)

Pm_distances_df <- Pm_distances_df |> dplyr::mutate(yval = get_y(value))

#Pm_distances_df <- Pm_distances_df |> dplyr::mutate(color = ifelse(paired == TRUE, "red", "black")) |> dplyr::mutate(alpha = ifelse(paired == TRUE, 1, 0.5))

cols <- c("TRUE" = "red", "FALSE" = "black")

alphas <- c("TRUE" = 1, "FALSE" = 0.8)

sizes <- c("TRUE" = 3, "FALSE" = 1)

Pm_dist <- Pm_distances_df |> ggplot(aes(x = value)) + theme_classic() + geom_histogram(aes(x = value, y = after_stat(density)), fill = "white", color = "black") + geom_density(color = "blue", linewidth = 1, linetype = 1) + geom_point(aes(x = value, y = yval, color = paired, alpha = paired, size = paired, fill = paired)) + scale_color_manual(values = cols) + scale_fill_manual(values = cols) + scale_alpha_manual(values = alphas) + scale_size_manual(values = sizes) + geom_text(stat = "bin", aes(y = after_stat(density), label = after_stat(count)), hjust = 2, angle = 90) + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.position = "none") + labs(x = "Pairwise Distance", y = "Density")

ggsave("Pm_pairwise_distances.svg", Pm_dist, dpi = 600, width = 10, height = 8, units = "in")

#manually moved labels and brought red dots to the front in Inkscape

Pm_distances_homo <- fastreeR::vcf2dist(inputVCF, ignoreHets = TRUE)

Pm_distances_homo_df <- Pm_distances_homo |> as.matrix() |> reshape2::melt(varnames = c("Sample1", "Sample2")) |> subset(Sample1 != Sample2)

Pm_distances_homo_df <- Pm_distances_homo_df |> dplyr::mutate(prefix1 = dplyr::case_when(stringr::str_detect(Sample1, "^[:digit:]+\\_[:digit:]") ~ stringr::str_extract(Sample1, "^[:digit:]+"))) |> dplyr::mutate(prefix2 = dplyr::case_when(stringr::str_detect(Sample2, "^[:digit:]+\\_[:digit:]") ~ stringr::str_extract(Sample2, "^[:digit:]+"))) |> dplyr::mutate(paired = dplyr::case_when(prefix1 == prefix2 ~ TRUE,
                                                                                                                                                                                                                                                                                                                                                                                        .default = FALSE))

Pm_dist_homo_density <- density(Pm_distances_homo_df$value)

get_y_homo <- approxfun(Pm_dist_homo_density$x, Pm_dist_homo_density$y)

Pm_distances_homo_df <- Pm_distances_homo_df |> dplyr::mutate(yval = get_y_homo(value))

#Pm_distances_df <- Pm_distances_df |> dplyr::mutate(color = ifelse(paired == TRUE, "red", "black")) |> dplyr::mutate(alpha = ifelse(paired == TRUE, 1, 0.5))

cols <- c("TRUE" = "red", "FALSE" = "black")

alphas <- c("TRUE" = 1, "FALSE" = 0.8)

sizes <- c("TRUE" = 3, "FALSE" = 1)

Pm_dist_homo <- Pm_distances_homo_df |> ggplot(aes(x = value)) + theme_classic() + geom_histogram(aes(x = value, y = after_stat(density)), fill = "white", color = "black") + geom_density(color = "blue", linewidth = 1, linetype = 1) + geom_point(aes(x = value, y = yval, color = paired, alpha = paired, size = paired, fill = paired)) + scale_color_manual(values = cols) + scale_fill_manual(values = cols) + scale_alpha_manual(values = alphas) + scale_size_manual(values = sizes) + geom_text(stat = "bin", aes(y = after_stat(density), label = after_stat(count)), hjust = 2, angle = 90) + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.position = "none") + labs(x = "Pairwise Distance", y = "Density")

ggsave("Pm_pairwise_distances_homo.svg", Pm_dist_homo, dpi = 600, width = 10, height = 8, units = "in")

#manually moved labels and brought red dots to the front in Inkscape
