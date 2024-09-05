setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

coverage_files <- list.files(path = "all_coverage/", pattern = "\\.txt")

setwd("all_coverage/")

coverage_tables <- lapply(coverage_files, data.table::fread)

names(coverage_tables) <- coverage_files

#coverage_tables <- lapply(coverage_tables, function(df) mutate(df, Sample = names(coverage_tables)))

coverage_tables <- purrr::imap(coverage_tables, ~dplyr::mutate(.x, Sample = .y))

all_coverage <- do.call(rbind, coverage_tables)

all_coverage <- all_coverage |> dplyr::rename(Contig = `#rname`)

coverage <- all_coverage |> dplyr::select(Contig, coverage, Sample)

coverage <- coverage |> subset(!stringr::str_detect(Contig, "archived|API|MIT"))

mean_coverage_by_contig <- coverage |> dplyr::group_by(Contig) |> dplyr::summarise(Mean = mean(coverage))


coverage_by_chromosome <- coverage |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                           stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                           stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                           stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                           .default = "Tanzania"))

#coverage_by_chromosome$Chromosome <- factor(coverage_by_chromosome$Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "API"))

#Tz_coverage_by_chromosome <- coverage_by_chromosome |> subset(Country == "Tanzania")

#determine what percent of genome all samples achieve 10X coverage at

coverage_thresholds <- coverage |> dplyr::mutate(`10X` = dplyr::case_when(coverage >= 10 ~ 1,
                                                                          .default = 0),
                                                 `20X` = dplyr::case_when(coverage >= 20 ~ 1,
                                                                          .default = 0),
                                                 `50X` = dplyr::case_when(coverage >= 50 ~ 1,
                                                                          .default = 0)) |>
  dplyr::group_by(Sample) |> dplyr::mutate(Prop_10X = sum(`10X`)/15, Prop_20X = sum(`20X`)/15, Prop_50X = sum(`50X`)/15)

sample_coverage_summary <- coverage_thresholds |> dplyr::group_by(Sample) |> dplyr::summarize(Chroms_10X = sum(`10X`), Chroms_20X = sum(`20X`), Chroms_50X = sum(`50X`))

Pm_wsaf_eigenvectors <- data.table::fread("../plink_WSAF_filtered.eigenvec")

final_samples <- vcfR::read.vcfR("../Pm_HC_missingness_filtered_first.vcf.gz") |> vcfR::vcfR2tidy()

variant_count <- final_samples$fix |> subset(!stringr::str_detect(CHROM, "archived")) |> subset(!stringr::str_detect(CHROM, "MIT")) |> subset(!stringr::str_detect(CHROM, "API")) |> subset(FILTER == "PASS")

final_samples <- final_samples$gt$Indiv |> unique() |> as.data.frame()

names(final_samples) <- "Sample"

sample_coverage_summary <- sample_coverage_summary |> dplyr::mutate(Sample = dplyr::case_when(stringr::str_detect(Sample, "Test") ~ stringr::str_extract(Sample, "^.*-Test_Pm_HC"),
                                                                                              .default = stringr::str_extract(Sample, "^.*-Pm_HC[:digit:]+")))

final_sample_coverage_summary <- dplyr::semi_join(sample_coverage_summary, final_samples) |> dplyr::mutate(Threshold = dplyr::case_when(Chroms_50X == 14 ~ "50",
                                                                                                                                                                      Chroms_20X == 14 ~ "20",
                                                                                                                                                                      Chroms_10X == 14 ~ "10",
                                                                                                                                                                      .default = "<10"))

below_10X_coverage <- final_sample_coverage_summary |> subset(Chroms_10X < 14)

unique(below_10X_coverage$Sample)

below_20X_coverage <- final_sample_coverage_summary |> subset(Chroms_20X < 14)

unique(below_20X_coverage$Sample)

below_50X_coverage <- final_sample_coverage_summary |> subset(Chroms_50X < 14)

unique(below_50X_coverage$Sample)

final_sample_coverage_summary |> writexl::write_xlsx("../final_sample_coverage_summary.xlsx")

sample_chrom_coverage_table <- final_sample_coverage_summary |> dplyr::count(Threshold)

library(viridis)
palette(turbo(15))

coverage_by_chromosome |> 
  ggplot() +
  geom_violin(aes(x = Contig, y = coverage, color = Contig, fill = Contig)) +
  #geom_jitter(aes(x = Chromosome, y = coverage, color = Chromosome), alpha = 0.3, size = 0.5, height=0.1) +
  theme_linedraw() +
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        plot.title = element_text(size = 28),
        axis.text.x = element_text(size = 20)) +
  #scale_y_continuous(labels = scales::percent) +
  #ylim(0,10) +
  labs(y = "Coverage", x = "Chromosome") +
  scale_x_discrete(labels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14))

ggsave("coverage_plot_violin.png", dpi = 600, width = 12, height = 10, units = "in")


coverage_by_chromosome |> 
  ggplot() +
  geom_boxplot(aes(x = Contig, y = coverage, color = Contig)) +
  #geom_jitter(aes(x = Chromosome, y = coverage, color = Chromosome), alpha = 0.3, size = 0.5, height=0.1) +
  theme_linedraw() +
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        plot.title = element_text(size = 30),
        axis.text.x = element_text(size = 30)) +
  #ylim(0,10) +
  labs(y = "Coverage", x = "Chromosome") +
  scale_x_discrete(labels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14))
  ggtitle(expression(paste("Sequencing Coverage across ", italic("P. malariae"), " Genome")))

ggsave("coverage_plot_boxplot.png", device = "png", dpi = 600, width = 12, height = 10, units = "in")

mean_coverage_by_chromosome <- coverage_by_chromosome |> dplyr::group_by(Contig) |> dplyr::summarise(Mean = mean(coverage))

all_chromosomes <- data.frame(Contig =c("All_Chromosomes"), Mean = c(mean(mean_coverage_by_chromosome$Mean)))

mean_coverage_by_chromosome <- rbind(mean_coverage_by_chromosome, all_chromosomes) |> dplyr::rename(Chromosome = Contig)

data.table::fwrite(mean_coverage_by_chromosome, "Mean Coverage.csv")

setwd("../all_refstats/")

refstats_files <- list.files(pattern = "\\.txt")

refstats_tables <- lapply(refstats_files, data.table::fread)

names(refstats_tables) <- refstats_files

refstats_tables <- purrr::imap(refstats_tables, ~dplyr::mutate(.x, Sample = .y))

#refstats_tables <- lapply(refstats_tables, function(df) mutate(df, percent_reads_assigned = (assignedReads/(sum(assignedReads)))*100))

all_refstats <- do.call(rbind, refstats_tables)

all_refstats <- all_refstats |> dplyr::rename(Species = `#name`)

percentages <- all_refstats |> tidyr::pivot_wider(names_from = Species, values_from = assignedReads, id_cols = c("Sample")) |> as.data.frame() |>
  mutate(All_Reads = Hs + Pf + Pm) |> mutate(Percent_Pm = (Pm / All_Reads) * 100)

percentages <- percentages |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                       stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                       stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                       stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                       .default = "Tanzania"))
#Tz_percentages <- percentages |> subset(Country == "Tanzania")

summary(percentages$Percent_Pm)

enrichment_plot <- percentages |> tidyr::pivot_longer(cols = c("Pm", "Hs", "Pf"), names_to = "read_type", values_to = "Count")

enrichment_plot <- enrichment_plot |> dplyr::group_by(Sample) |> dplyr::mutate(Pm_Reads = dplyr::case_when(read_type == "Pm" ~ Count))

Pm_reads_only <- enrichment_plot |> subset(read_type == "Pm")

summary(Pm_reads_only$Count)

summary(Pm_reads_only$Percent_Pm)

enrichment_plot |> ggplot(aes(fill = read_type, y = Count, x = factor(Sample, levels = unique(Sample[order(Pm_Reads, decreasing = TRUE)]), ordered = TRUE), angle = 90)) +
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
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "YlOrBr", direction = 1) +
  ylab("Percent of Reads") +
  xlab("Sample (Sorted by Number of Pm Reads)") +
  ggtitle(expression(paste("Proportion of ", italic("P. malariae "), "Reads"))) #+
  #annotate("text", x = 50, y = -0.05, label = "Pm Reads", angle = 90) +
  #annotate("text", x = 50, y = 0.95, label = "All Reads", angle = 90) +
  #coord_cartesian(xlim = c(0,49), clip = "off")

ggsave("enrichment_plot.png", dpi = 600, width = 12, height = 10, units = "in")

enrichment_plot <- enrichment_plot |> dplyr::mutate(Capture = dplyr::case_when(stringr::str_detect(Sample, "HC[:digit:]+") ~ stringr::str_extract(Sample, "HC[:digit:]+"),
                                                                               stringr::str_detect(Sample, "Test_Pm_HC") ~ "Test_HC"))
capture_efficiency_aov <- aov(Percent_Pm ~ Capture, data = enrichment_plot)

summary(capture_efficiency_aov)

TukeyHSD(capture_efficiency_aov)

enrichment_plot |> ggplot() +
  geom_boxplot(aes(x = Capture, y = Percent_Pm, color = Capture)) +
  theme_linedraw()

all_HC_metadata <- readxl::read_xlsx("../../Pm Samples for Hybrid Capture.xlsx") |> dplyr::select(Sample, `Pm Ct`, Source, `Sequence Prep`) 

#Tz_HC_metadata <- all_HC_metadata |> subset(Source != "DRC") |> subset(Source != "Gambia")

enrichment_by_parasitemia <- enrichment_plot |> dplyr::mutate(Sample = dplyr::case_when(stringr::str_detect(Sample,"-Pm.+") ~ stringr::str_remove(Sample,"-Pm.+"),
                                                                                        stringr::str_detect(Sample,"-Test_Pm.+") ~ stringr::str_remove(Sample,"-Test_Pm.+"))) |>
  subset(read_type == "Pm") |> dplyr::left_join(all_HC_metadata)

capture_efficiency_parasitemia <- glm(Percent_Pm ~ `Pm Ct`, data = enrichment_by_parasitemia)

capture_efficiency_parasitemia_LM <- lm(Percent_Pm ~ `Pm Ct`, data = enrichment_by_parasitemia)

summary(capture_efficiency_parasitemia)

summary(capture_efficiency_parasitemia_LM)

enrichment_by_parasitemia |> ggplot(aes(`Pm Ct`, Percent_Pm)) +
  geom_point() +
  theme_classic() +
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth")

