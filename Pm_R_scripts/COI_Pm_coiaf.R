#devtools::install_github("bailey-lab/coiaf@v0.1.2")

library(vcfR)
library(coiaf)

setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun")

#bcftools view -e 'MAF[0]<0.05' -r PmUG01_01_v1,PmUG01_02_v1,PmUG01_03_v1,PmUG01_04_v1,PmUG01_05_v1,PmUG01_06_v1,PmUG01_07_v1,PmUG01_08_v1,PmUG01_09_v1,PmUG01_10_v1,PmUG01_11_v1,PmUG01_12_v1,PmUG01_13_v1,PmUG01_14_v1 -Oz -o Pm_all_exclude_MAF5.vcf.gz Pm_HC_missingness_filtered_first.vcf.gz


Pm_HC <- vcfR::read.vcfR("Pm_all_exclude_MAF5.vcf.gz")

# 1. create gt
gtmat <- vcfRmanip::gtmat012(Pm_HC) / 2
gtmat[is.na(gtmat)] <- -1
gtmat <- t(gtmat)
colnames(gtmat) <- paste0(Pm_HC@fix[, "CHROM"], "_", Pm_HC@fix[, "POS"])

# 2. create wsaf
# extract coverage and counts matrices
coverage <- t(vcfR::extract.gt(Pm_HC, element = "DP", as.numeric = T))
counts_raw <- t(vcfR::extract.gt(Pm_HC, element = "AD"))
counts <- vcfR::masplit(counts_raw, record = 1, sort = FALSE, decreasing = FALSE)
wsaf <- counts / coverage

wsaf_new <- wsaf
wsaf_new[(gtmat == 0)] <- 1
wsaf_new[(gtmat == 1)] <- 0

plaf <- colMeans(wsaf_new, na.rm = TRUE)

input_data <- purrr::map(seq_len(nrow(wsaf_new)), function(i) {
  tibble::tibble(wsmaf = wsaf_new[i, ], plmaf = plaf) |>
    tidyr::drop_na()
})

COI_coaif_discrete <- purrr::map_dbl(input_data, ~ optimize_coi(.x, data_type = "real"))

rounded_COI <- COI_coaif_discrete |> round()

rounded_labeled_COI <- rounded_COI |> cbind(wsaf_new[,0])

Pm_COI <- as.data.frame(rounded_labeled_COI)

#COI_coaif_continuous <- purrr::map(input_data, ~ compute_coi(.x, data_type = "real"))


Pm_COI <- Pm_COI |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(rownames(Pm_COI), "Gam_[:digit:]+") ~ "Nigeria",
                                                             stringr::str_detect(rownames(Pm_COI), "^[:digit:]+") ~ "DRC",
                                                             stringr::str_detect(rownames(Pm_COI), "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                             stringr::str_detect(rownames(Pm_COI), "Gam") ~ "Cameroon",
                                                             .default = "Tanzania"))

Pm_COI |> saveRDS("Pm_coiaf_COI.rds")

Pm_COI <- readRDS("Pm_coiaf_COI.rds")

library(viridis)
palette(turbo(4))
library(ggplot2)

# A function factory for getting integer y-axis values.
integer_breaks <- function(n = 10, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}


country_COI <- Pm_COI |>
  ggplot() +
  geom_violin(aes(x = Country, y = rounded_COI, color = Country)) +
  geom_jitter(aes(x = Country, y = rounded_COI, color = Country), 
              alpha = 0.3, size = 0.5, height=0.1) +
  theme_classic() +
  theme(legend.position = "none") + 
  ylim(0,10) +
  labs(y = "Complexity of Infection") +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 24)) +
  scale_y_continuous(breaks = integer_breaks())

ggsave("COI_results/coiaf_all_countries_violin.png", device = "png", dpi = 600, width = 6, height = 5, units = "in")

#monoclonals <- Pm_COI |> subset(rounded_COI == 1) |> tibble::rownames_to_column() |> dplyr::rename(Sample = rowname) |> dplyr::select(Sample) |> data.table::fwrite("monoclonals_redo.csv")

Pf_COI <- data.table::fread("COI_for_Kelly.csv") |> dplyr::rename(Sample = sample)

Pf_samples <- data.table::fread("Pf_COI_samples.txt", header = F) |> dplyr::rename(Sample = V1)

Pf_COI <- Pf_samples |> dplyr::left_join(Pf_COI) |> subset(is.na(mean_COI) == FALSE) |> dplyr::select(Sample, mean_COI, country)

#Pm_COI <- Pm_COI |> dplyr::select(Indiv, mean, Country)

Pf_COI$Species <- "Pf"

Pm_COI$Species <- "Pm"

Pf_COI <- Pf_COI |> dplyr::rename(Country = country)

Pm_COI <- Pm_COI |> tibble::rownames_to_column(var = "Sample") |> dplyr::rename(mean_COI = rounded_COI)

combined_COI <- rbind(Pf_COI, Pm_COI)

species_COI <- combined_COI |>
  ggplot() +
  geom_violin(aes(x = Species, y = mean_COI, color = Species)) +
  geom_jitter(aes(x = Species, y = mean_COI, color = Species), 
              alpha = 0.3, size = 0.5, height=0.1) +
  theme_classic() +
  theme(legend.position = "none") + 
  ylim(0,12) +
  labs(y = "Complexity of Infection") +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 24)) +
  scale_y_continuous(breaks = integer_breaks()) #+
#ggtitle(expression(paste(italic("P. malariae"), " Complexity of Infection")))

ggsave("COI_results/COI_species_violin_coiaf.png", device = "png", dpi = 600, width = 6, height = 5, units = "in")

species_anova <- aov(combined_COI$mean_COI ~ combined_COI$Species * combined_COI$Country)

summary(species_anova)

TukeyHSD(species_anova)

library(patchwork)

COI_figure <- country_COI + species_COI + plot_layout(axis_titles = "collect_y") + plot_annotation(tag_levels = "A")

ggsave("COI_figure.png", COI_figure, dpi = 600, width = 12, height = 5)

Pm_polyclonals <- Pm_COI |> subset(mean_COI > 1)

Pf_polyclonals <- Pf_COI |> subset(mean_COI > 1)

Pf_polyclonals |> dplyr::count(mean_COI)
