###############################################################
####################### ADMIXTURE #############################
###############################################################
#Description: Uses ADMIXTURE to calculate estimated number of 
#population clusters within sample pool

setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

#have to rename chromosomes first

system("bcftools annotate --rename-chrs Pm_chr_rename.txt -Oz -o ADMIXTURE_input.vcf.gz Pm_PCA.vcf.gz")

system("plink --vcf ADMIXTURE_input.vcf.gz --const-fid --allow-extra-chr --out ADMIXTURE")

#use admixture to evaluate cross-validation error among models using 1-10 population clusters
system("for K in 1 2 3 4 5 6 7 8 9 10; do admixture32 --cv --haploid="*" ADMIXTURE_input.bed $K | tee log${K}.out; done")

#extract cross-validation error values from each model
system("grep -h CV log*.out > admixture_CV.txt")

#extract and plot CV error to determine best model fit
CV_values <- data.table::fread("admixture_CV.txt")

CV_values <- CV_values |> dplyr::select(V3, V4) |> dplyr::rename(K = V3, CV = V4)

CV_values$K <- stringr::str_extract(CV_values$K, "[:digit:]+")

CV_values$K <- as.numeric(CV_values$K)

CV_plot <- CV_values |> ggplot(aes(x = K, y = CV)) + geom_line() + geom_point() +
  theme_classic() + scale_x_continuous(breaks = c(1:10))

ggsave("ADMIXTURE_CV.png", CV_plot, dpi = 600)

#for best model, plot sample assignment and cluster makeup
admixture_data <- data.table::fread("ADMIXTURE_input.2.Q")

admixture_samples <- data.table::fread("ADMIXTURE_input.fam")

admixture_data$Sample <- admixture_samples$V2

admixture_data <- admixture_data |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                   stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                   stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                   stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                   .default = "Tanzania"))

admixture_data <- admixture_data |> dplyr::arrange(V1) |> dplyr::rename(Pop1 = V1, Pop2 = V2)

admixture_data |> saveRDS("admixture_data.rds")

pure_samples <- admixture_data |> subset(Pop1 > 0.99 | Pop2 > 0.99) 

admixture_data <- admixture_data |> tidyr::pivot_longer(cols = c(Pop1, Pop2), names_to = "Population", values_to = "Percentage")

admixture_data <- admixture_data |> dplyr::mutate(Pop1_percentage = dplyr::case_when(Population == "Pop1" ~ Percentage)) |> dplyr::group_by(Sample) |> dplyr::arrange(Pop1_percentage)

admixture_data$Sample <- factor(admixture_data$Sample, levels = unique(admixture_data$Sample))

admixture_plot <- admixture_data |> ggplot(aes(fill = Population, y = Percentage, x = Sample, pattern = Country, pattern_key_scale_factor=0.5, pattern_density = 0.25)) +
  #geom_bar(position = "fill", stat = "identity") +
  #geom_text(aes(x = Sample, y = 0.95, label = scales::comma(All_Reads))) +
  #geom_text(aes(x = Sample, y = -0.05, label = scales::comma(Pm_Reads))) +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.text = element_text(size = 30)) +
  ggpattern::geom_bar_pattern(stat = "identity", color = "black", pattern_fill = "black") +
  ggpattern::scale_pattern_manual(values = c("stripe", "crosshatch", "circle", "none"), guide = guide_legend(override.aes = list(fill = "white"))) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2", direction = 1, guide = guide_legend(override.aes = list(pattern = "none"))) +
  #ylab("Percent of Reads") +
  #xlab("Sample (Sorted by Number of Pm Reads)") +
  ggtitle(expression(paste("ADMIXTURE Population Estimates")))

ggsave("ADMIXTURE_plot_wsaf_filtered.png", admixture_plot, dpi = 600, width = 12, height = 10, units = "in")

Cameroon_plot <- admixture_data |> 
  dplyr::filter(Country == "Cameroon") |>
  ggplot(aes(fill = Population, y = Percentage, x = Sample)) +
  geom_bar(position = "fill", stat = "identity") +
  #geom_text(aes(x = Sample, y = 0.95, label = scales::comma(All_Reads))) +
  #geom_text(aes(x = Sample, y = -0.05, label = scales::comma(Pm_Reads))) +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.text = element_text(size = 30)) +
  #ggpattern::geom_bar_pattern(stat = "identity", color = "black", pattern_fill = "black") +
  #ggpattern::scale_pattern_manual(values = c("stripe", "crosshatch", "circle", "none"), guide = guide_legend(override.aes = list(fill = "white"))) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2", direction = 1, guide = guide_legend(override.aes = list(pattern = "none"))) +
  #ylab("Percent of Reads") +
  #xlab("Sample (Sorted by Number of Pm Reads)") +
  ggtitle("Cameroon")


DRC_plot <- admixture_data |> 
  dplyr::filter(Country == "DRC") |>
  ggplot(aes(fill = Population, y = Percentage, x = Sample)) +
  geom_bar(position = "fill", stat = "identity") +
  #geom_text(aes(x = Sample, y = 0.95, label = scales::comma(All_Reads))) +
  #geom_text(aes(x = Sample, y = -0.05, label = scales::comma(Pm_Reads))) +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.text = element_text(size = 30)) +
  #ggpattern::geom_bar_pattern(stat = "identity", color = "black", pattern_fill = "black") +
  #ggpattern::scale_pattern_manual(values = c("stripe", "crosshatch", "circle", "none"), guide = guide_legend(override.aes = list(fill = "white"))) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2", direction = 1, guide = guide_legend(override.aes = list(pattern = "none"))) +
  #ylab("Percent of Reads") +
  #xlab("Sample (Sorted by Number of Pm Reads)") +
  ggtitle("DRC")


Nigeria_plot <- admixture_data |> 
  dplyr::filter(Country == "Nigeria") |>
  ggplot(aes(fill = Population, y = Percentage, x = Sample)) +
  geom_bar(position = "fill", stat = "identity") +
  #geom_text(aes(x = Sample, y = 0.95, label = scales::comma(All_Reads))) +
  #geom_text(aes(x = Sample, y = -0.05, label = scales::comma(Pm_Reads))) +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.text = element_text(size = 30)) +
  #ggpattern::geom_bar_pattern(stat = "identity", color = "black", pattern_fill = "black") +
  #ggpattern::scale_pattern_manual(values = c("stripe", "crosshatch", "circle", "none"), guide = guide_legend(override.aes = list(fill = "white"))) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2", direction = 1, guide = guide_legend(override.aes = list(pattern = "none"))) +
  #ylab("Percent of Reads") +
  #xlab("Sample (Sorted by Number of Pm Reads)") +
  ggtitle("Nigeria")


Tanzania_plot <- admixture_data |> 
  dplyr::filter(Country == "Tanzania") |>
  ggplot(aes(fill = Population, y = Percentage, x = Sample)) +
  geom_bar(position = "fill", stat = "identity") +
  #geom_text(aes(x = Sample, y = 0.95, label = scales::comma(All_Reads))) +
  #geom_text(aes(x = Sample, y = -0.05, label = scales::comma(Pm_Reads))) +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.text = element_text(size = 30)) +
  #ggpattern::geom_bar_pattern(stat = "identity", color = "black", pattern_fill = "black") +
  #ggpattern::scale_pattern_manual(values = c("stripe", "crosshatch", "circle", "none"), guide = guide_legend(override.aes = list(fill = "white"))) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2", direction = 1, guide = guide_legend(override.aes = list(pattern = "none"))) +
  #ylab("Percent of Reads") +
  #xlab("Sample (Sorted by Number of Pm Reads)") +
  ggtitle("Tanzania")

library(patchwork)

design <- "
ABBB
CDDD
"

combined_admixture <- Cameroon_plot + DRC_plot + Nigeria_plot + Tanzania_plot + plot_layout(design = design, guides = "collect", axes = "collect") 

ggsave("admixture_plot_revised.png", combined_admixture, dpi = 600, width = 7200, height = 6000, units = "px")

admixture_pop1 <- admixture_data |> subset(Pop1_percentage > 0.5)

admixture_pop2 <- admixture_data |> subset(Pop1_percentage < 0.5) |> dplyr::mutate(Population = "Pop2")

admixture_pops <- rbind(admixture_pop1, admixture_pop2)

dadi_admixture <- admixture_pops |> dplyr::select(Sample, Population)

dadi_admixture |> data.table::fwrite("dadi_admixture.txt", sep = "\t")

