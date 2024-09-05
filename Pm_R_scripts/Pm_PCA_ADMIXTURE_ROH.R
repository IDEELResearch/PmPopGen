#bcftools annotate --rename-chrs Pm_chr_rename.txt -O z -o Pm_PCA.vcf.gz Pm_HC_missingness_filtered.vcf.gz
#bcftools index Pm_PCA.vcf.gz
#bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14 -O z -o Pm_PCA_chrs_only.vcf.gz Pm_PCA.vcf.gz
#bcftools index Pm_PCA_chrs_only.vcf.gz

setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")
#library("devtools")
library(gdsfmt)
library(SNPRelate)
library(SeqArray) 
library(VariantAnnotation)
library(vcfR)
library(tibble)
library(ggplot2)
library(viridis)

vcf.fn <- c("Pm_PCA_chrs_only.vcf.gz")

snpgdsVCF2GDS(vcf.fn, "PCA.gds")
genofile <- snpgdsOpen("PCA.gds")
snpgdsSummary("PCA.gds")

chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))

snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
snpset.id <- unlist(unname(snpset))

Pm_pca <-snpgdsPCA(genofile, snp.id=snpset.id)

study_table <- data.frame(sample.id = Pm_pca$sample.id)

study_table <- study_table |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                       stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                       stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                       .default = "Tanzania"))
#Tz_samples <- study_table |> subset(Country == "Tanzania")

#as.list(Tz_samples$sample.id) |> data.table::fwrite(file = "Tz_samples.txt", col.names = FALSE, sep = "\n")

#bcftools view -S Tz_samples.txt -O z -o Pm_PCA_Tz.vcf.gz Pm_PCA_chrs_only.vcf.gz

#bcftools index Pm_PCA_Tz.vcf.gz

#vcf.fn <- c("Pm_PCA_Tz.vcf.gz")

#snpgdsVCF2GDS(vcf.fn, "Tz_PCA.gds")
#genofile <- snpgdsOpen("Tz_PCA.gds")
#snpgdsSummary("Tz_PCA.gds")

#chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))

#snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
#snpset.id <- unlist(unname(snpset))

#Pm_pca <-snpgdsPCA(genofile, snp.id=snpset.id)

#study_table <- data.frame(sample.id = Pm_pca$sample.id)

#study_table <- study_table |> dplyr::mutate(Region = dplyr::case_when(stringr::str_detect(sample.id, "TransMIT") ~ "Pwani",
#                                                                      stringr::str_detect(sample.id, "DOCW") | stringr::str_detect(sample.id, "OMP") ~ "Dodoma",
#                                                                      stringr::str_detect(sample.id, "KG[A-Z]{2}") ~ "Kagera",
#                                                                      stringr::str_detect(sample.id, "KL[A-Z]{2}") ~ "Kilimanjaro",
#                                                                      stringr::str_detect(sample.id, "LUN-") ~ "Ruvuma",
#                                                                      stringr::str_detect(sample.id, "MAB-") | stringr::str_detect(sample.id, "MPY-") ~ "Tanga",
#                                                                      stringr::str_detect(sample.id, "MA[A-Z]{2}") ~ "Mara",
#                                                                      stringr::str_detect(sample.id, "MT[A-Z]{2}") ~ "Mtwara",
#                                                                      stringr::str_detect(sample.id, "MY[A-Z]{2}") ~ "Manyara",
#                                                                      stringr::str_detect(sample.id, "NJ[A-Z]{2}") ~ "Njombe",
#                                                                      stringr::str_detect(sample.id, "NYA-") ~ "Kigoma",
#                                                                      stringr::str_detect(sample.id, "SO[A-Z]{2}") ~ "Songwe",
#                                                                      stringr::str_detect(sample.id, "TB[A-Z]{2}") ~ "Tabora",
#                                                                      .default = "Geita"))

pc.percent <- Pm_pca$varprop*100

filtered_tab <- data.frame(sample.id = Pm_pca$sample.id,
                           EV1 = Pm_pca$eigenvect[,1],    # the first eigenvector
                           EV2 = Pm_pca$eigenvect[,2],    # the second eigenvector
                           stringsAsFactors = FALSE)

plot(filtered_tab$EV2, filtered_tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1",)

filtered_tab_site_info <- data.frame(sample.id = Pm_pca$sample.id,
                                     country = factor(study_table$Country),
                                     EV1 = Pm_pca$eigenvect[,1],    # the first eigenvector
                                     EV2 = Pm_pca$eigenvect[,2],    # the second eigenvector
                                     stringsAsFactors = FALSE)

filtered_tab_site_info |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab("PC1 (5.35%)") +
  ylab("PC2 (2.64%)") +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA.png", width = 15, height = 12, units = "in", dpi = 600)

Tz_Pm_summary <- study_table |> dplyr::summarise(Count = n(), .by = Region) |> data.table::fwrite("region_counts.csv")

Pm_summary <- study_table |> dplyr::summarise(Count = n(), .by = Country) |> data.table::fwrite("country_counts.csv")

#PLINK PCA

plink_eigenvalues <- data.table::fread("plink.eigenval")
plink_eigenvectors <- data.table::fread("plink.eigenvec")

plink_study_table <- data.frame(sample.id = plink_eigenvectors$FID)

plink_study_table <- plink_study_table |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(sample.id, "Gam") ~ "Nigeria",
                                                                       stringr::str_detect(sample.id, "^[:digit:]+") ~ "DRC",
                                                                       stringr::str_detect(sample.id, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                       .default = "Tanzania"))

plink_filtered_tab <- data.frame(sample.id = plink_eigenvectors$FID,
                           EV1 = plink_eigenvectors$PC1,    # the first eigenvector
                           EV2 = plink_eigenvectors$PC2,    # the second eigenvector
                           stringsAsFactors = FALSE)

plot(plink_filtered_tab$EV2, plink_filtered_tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1",)

plink_filtered_tab_site_info <- data.frame(sample.id = plink_eigenvectors$FID,
                                     country = factor(plink_study_table$Country),
                                     EV1 = plink_eigenvectors$PC1,    # the first eigenvector
                                     EV2 = plink_eigenvectors$PC2,    # the second eigenvector
                                     stringsAsFactors = FALSE)

plink_filtered_tab_site_info |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, color = factor(country)), size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab("PC1 (20.53%)") +
  ylab("PC2 (4.61%)") +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK.png", width = 15, height = 12, units = "in", dpi = 600)

Pm_summary_plink <- plink_study_table |> dplyr::summarise(Count = n(), .by = Country) |> data.table::fwrite("plink_country_counts.csv")

#redoing ADMIXTURE with MAF10 and WSAF10 filters

#have to rename chromosomes first

#bcftools annotate --rename-chrs Pm_chr_rename.txt -Oz -o ADMIXTURE_input.vcf.gz Pm_monoclonal_exclude_MAF10_WSAF10.vcf.gz

#replacing hets with major allele

#bcftools +setGT -Oz -o ADMIXTURE_input_nohets.vcf.gz ADMIXTURE_input.vcf.gz -- -t "b:AD>0" -n M

#$plink --vcf ADMIXTURE_input_nohets.vcf.gz --const-fid --allow-extra-chr --out ADMIXTURE_input_nohets

#for K in 1 2 3 4 5 6 7 8 9 10; do admixture32 --cv --haploid="*" ADMIXTURE_input_nohets.bed $K | tee haploid_nohets_log${K}.out; done

admixture_data <- data.table::fread("admixture_wsaf.2.Q")

admixture_samples <- data.table::fread("admixture_wsaf.fam")

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

#grep -h CV haploid_wsaf_filtered_log*.out > haploid_wsaf_CV.txt

CV_values <- data.table::fread("haploid_wsaf_CV.txt")

CV_values <- CV_values |> dplyr::select(V3, V4) |> dplyr::rename(K = V3, CV = V4)

CV_values$K <- stringr::str_extract(CV_values$K, "[:digit:]+")

CV_values$K <- as.numeric(CV_values$K)

CV_plot <- CV_values |> ggplot(aes(x = K, y = CV)) + geom_line() + geom_point() +
  theme_classic() + scale_x_continuous(breaks = c(1:10))

ggsave("ADMIXTURE_CV_wsaf_filtered.png", CV_plot, dpi = 600)

admixture_pop1 <- admixture_data |> subset(Pop1_percentage > 0.5)

admixture_pop2 <- admixture_data |> subset(Pop1_percentage < 0.5) |> dplyr::mutate(Population = "Pop2")

admixture_pops <- rbind(admixture_pop1, admixture_pop2)

dadi_admixture <- admixture_pops |> dplyr::select(Sample, Population)

dadi_admixture |> data.table::fwrite("dadi_admixture.txt", sep = "\t")

plinkadmix <- plink_filtered_tab_site_info |> dplyr::rename(Sample = sample.id) |> dplyr::left_join(admixture_pops)

plinkadmix |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, group = interaction(Country, Population), color = Country, shape = Population), size = 4) +
  #geom_point(aes(group = interaction(Country, Population))) +
  scale_shape_manual(values = c(17, 19), labels = c(1, 2)) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab("PC1 (20.53%)") +
  ylab("PC2 (4.61%)") +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK ADMIXTURE Pops.png", width = 15, height = 12, units = "in", dpi = 600)

######

setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

structure_data_2 <- data.table::fread("faststructure_2pops.2.meanQ")

admixture_samples <- data.table::fread("plink.fam")

structure_data_2$Sample <- admixture_samples$V1

structure_data_2 <- structure_data_2 |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                             stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                             stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                             stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                             .default = "Tanzania"))

structure_data_2 <- structure_data_2 |> dplyr::arrange(V1) |> dplyr::rename(Pop1 = V1, Pop2 = V2)

structure_data_2 <- structure_data_2 |> tidyr::pivot_longer(cols = c(Pop1, Pop2), names_to = "Population", values_to = "Percentage")

structure_data_2 <- structure_data_2 |> dplyr::mutate(Pop1_percentage = dplyr::case_when(Population == "Pop1" ~ Percentage)) |> dplyr::group_by(Sample) |> dplyr::arrange(Pop1_percentage)

structure_data_2$Sample <- factor(structure_data_2$Sample, levels = unique(structure_data_2$Sample))

structure_plot_2 <- structure_data_2 |> ggplot(aes(fill = Population, y = Percentage, x = Sample, pattern = Country, pattern_key_scale_factor=0.5, pattern_density = 0.25)) +
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
  ggtitle(expression(paste("fastStructure Population Estimates (Two Populations)")))

ggsave("structure_plot_2.png", structure_plot_2, dpi = 600, width = 12, height = 10, units = "in")

structure_2_pop1 <- structure_data_2 |> subset(Pop1_percentage > 0.5)

structure_2_pop2 <- structure_data_2 |> subset(Pop1_percentage < 0.5) |> dplyr::mutate(Population = "Pop2")

structure_2_pops <- rbind(structure_2_pop1, structure_2_pop2)

plinkstructure2 <- plink_filtered_tab_site_info |> dplyr::rename(Sample = sample.id) |> dplyr::left_join(structure_2_pops)

plinkstructure2 |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, group = interaction(Country, Population), color = Country, shape = Population), size = 4) +
  #geom_point(aes(group = interaction(Country, Population))) +
  scale_shape_manual(values = c(17, 19), labels = c(1, 2)) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab("PC1 (20.53%)") +
  ylab("PC2 (4.61%)") +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK fastStructure 2 Pops.png", width = 15, height = 12, units = "in", dpi = 600)

structure_data_3 <- data.table::fread("faststructure_3pops.3.meanQ")

admixture_samples <- data.table::fread("plink.fam")

structure_data_3$Sample <- admixture_samples$V1

structure_data_3 <- structure_data_3 |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                 stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                 stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                 stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                 .default = "Tanzania"))

structure_data_3 <- structure_data_3 |> dplyr::arrange(V1) |> dplyr::rename(Pop1 = V1, Pop2 = V2, Pop3 = V3)

structure_data_3 <- structure_data_3 |> tidyr::pivot_longer(cols = c(Pop1, Pop2, Pop3), names_to = "Population", values_to = "Percentage")

structure_data_3 <- structure_data_3 |> dplyr::mutate(Pop1_percentage = dplyr::case_when(Population == "Pop1" ~ Percentage)) |> dplyr::mutate(Pop2_percentage = dplyr::case_when(Population == "Pop2" ~ Percentage)) |> dplyr::mutate(Pop3_percentage = dplyr::case_when(Population == "Pop3" ~ Percentage)) |> dplyr::group_by(Sample) |> dplyr::mutate(major_pop = dplyr::case_when(Pop1_percentage > 0.5 ~ "Pop1",
                                                                                                                                                                                                                                                                                                                                                                                      Pop2_percentage > 0.5 ~ "Pop2",
                                                                                                                                                                                                                                                                                                                                                                                      Pop3_percentage > 0.5 ~ "Pop3")) |> dplyr::arrange(major_pop)

structure_data_3$Sample <- factor(structure_data_3$Sample, levels = unique(structure_data_3$Sample))

structure_plot_3 <- structure_data_3 |> ggplot(aes(fill = Population, y = Percentage, x = Sample, pattern = Country, pattern_key_scale_factor=0.5, pattern_density = 0.25)) +
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
  ggtitle(expression(paste("fastStructure Population Estimates (Three Populations)")))

ggsave("structure_plot_3.png", structure_plot_3, dpi = 600, width = 12, height = 10, units = "in")

structure_3_pop1 <- structure_data_3 |> subset(Pop1_percentage > 0.5)

structure_3_pop2 <- structure_data_3 |> subset(Pop2_percentage > 0.5) |> dplyr::mutate(Population = "Pop2")

structure_3_pop3 <- structure_data_3 |> subset(Pop3_percentage > 0.5) |> dplyr::mutate(Population = "Pop3")

structure_3_pops <- rbind(structure_3_pop1, structure_3_pop2, structure_3_pop3)

plinkstructure3 <- plink_filtered_tab_site_info |> dplyr::rename(Sample = sample.id) |> dplyr::left_join(structure_3_pops)

plinkstructure3 |> ggplot() + theme_bw() +
  geom_point(aes(x = EV1, y = EV2, group = interaction(Country, Population), color = Country, shape = Population), size = 4) +
  #geom_point(aes(group = interaction(Country, Population))) +
  scale_shape_manual(values = c(17, 19, 15), labels = c(1, 2, 3)) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  labs(color = "Country") +
  xlab("PC1 (20.53%)") +
  ylab("PC2 (4.61%)") +
  theme(axis.text.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size = 20)) +
  ggtitle(expression(paste("Principal Component Analysis of ", italic("P. malariae"), " Isolates")))

ggsave("Pm PCA PLINK fastStructure 3 Pops.png", width = 15, height = 12, units = "in", dpi = 600)

Pm_ROH <- data.table::fread("Pm_roh_region.txt") |> dplyr::rename(Sample = "[2]Sample", CHROM = "[3]Chromosome", START = "[4]Start", END = "[5]End", LENGTH = "[6]Length (bp)", NMARKERS = "[7]Number of markers", QUAL = "[8]Quality (average fwd-bwd phred score)")

Pm_ROH <- Pm_ROH |> subset(!stringr::str_detect(CHROM, "archived|MIT|API"))

ROH_n <- Pm_ROH |> dplyr::count(Sample)

Pm_ROH <- Pm_ROH |> dplyr::group_by(Sample) |> dplyr::mutate(TOTAL_LENGTH = sum(LENGTH))

ROH_length <- Pm_ROH |> dplyr::select(Sample, TOTAL_LENGTH) |> dplyr::distinct()

ROH_plot <- dplyr::left_join(ROH_n, ROH_length)

ROH_plot <- ROH_plot |> dplyr::mutate(length_mb = TOTAL_LENGTH/1000000)

ROH <- ROH_plot |> ggplot() + geom_point(aes(x = length_mb, y = n), size = 3) + theme_classic() + labs(x = "Total Length of ROH (Mb)", y = "Total n ROH per Sample") + scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25)) + scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 150)) + theme(axis.title = element_text(size = 24), axis.text = element_text(size = 20))

ggsave("Pm_ROH.png", ROH, dpi = 600, height = 8, width = 8, units = "in")

Pm_ROH <- dplyr::left_join(Pm_ROH, plinkadmix) |> dplyr::select(Sample, CHROM, POS, STATE, QUAL, Country, Population, Percentage)

#Pm_ROH <- Pm_ROH |> dplyr::mutate(STATE = dplyr::case_when(STATE == 0 ~ "HW",
#                                                           STATE == 1 ~ "AZ"))



Pm_ROH_by_sample <- Pm_ROH |> dplyr::group_by(Sample) |> dplyr::count(STATE) |> dplyr::mutate(Prop_AZ = n/sum(n)) |> subset(STATE == 1) |> dplyr::left_join(plinkadmix) |> dplyr::select(Sample, STATE, n, Prop_AZ, Country, Population) |> subset(is.na(Population) == FALSE)

Pm_ROH_aov <- aov(Pm_ROH_by_sample$Prop_AZ ~ Pm_ROH_by_sample$Country * Pm_ROH_by_sample$Population)

summary(Pm_ROH_aov)

TukeyHSD(Pm_ROH_aov)

Pm_ROH_by_sample |> ggplot() + geom_boxplot(aes(x = Population, y = Prop_AZ))

Pm_chrlen <- data.table::fread("../Pm_chrlen.txt") |> dplyr::select(V1, V3) |> dplyr::mutate(V1 = stringr::str_remove(V1, "PmUG01_")) |> dplyr::mutate(V1 = stringr::str_remove(V1, "_v1")) |> dplyr::mutate(V1 = paste0("chr", V1)) |> dplyr::rename(CHROM = V1, LENGTH = V3)

Pm_chrlen <- Pm_chrlen |> dplyr::mutate(cumsum = ave(LENGTH, FUN=cumsum), toadd = head(c(0, cumsum), -1)) 

Pm_chrlen <- Pm_chrlen |> dplyr::mutate(axis_breaks = (LENGTH/2) + toadd)

Pm_ROH <- Pm_ROH |> dplyr::mutate(CHROM = paste0("chr", stringr::str_replace(CHROM, "PmUG01_(\\d+)_v1", "\\1"))) |> dplyr::left_join(Pm_chrlen) |> dplyr::mutate(genome_pos = POS + toadd)

Pm_ROH_by_site <- Pm_ROH |> dplyr::group_by(CHROM, POS, Population) |> dplyr::count(STATE) |> dplyr::mutate(Prop_AZ = n/sum(n)) |> subset(STATE ==1) |> dplyr::left_join(Pm_chrlen) |> subset(is.na(Population) == FALSE)|> dplyr::mutate(genome_pos = POS + toadd)


Pm_ROH_by_site |> ggplot() + geom_boxplot(aes(x = Population, y = Prop_AZ))

Pm_ROH_plot <- Pm_ROH_by_site |> ggplot() + geom_line(aes(x = genome_pos, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + scale_x_continuous(name = "Chromosome", breaks = Pm_chrlen$axis_breaks, labels = c(1:14), minor_breaks = NULL) + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous")

ggsave("Pm_ROH_plot.png", Pm_ROH_plot, dpi = 600)

library(scales)

chr1plot <- Pm_ROH_by_site |> subset(CHROM == "chr01") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 1 Position") + scale_x_continuous(labels = label_comma())

chr2plot <- Pm_ROH_by_site |> subset(CHROM == "chr02") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 2 Position") + scale_x_continuous(labels = label_comma())

chr3plot <- Pm_ROH_by_site |> subset(CHROM == "chr03") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 3 Position") + scale_x_continuous(labels = label_comma())

chr4plot <- Pm_ROH_by_site |> subset(CHROM == "chr04") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 4 Position") + scale_x_continuous(labels = label_comma())

chr5plot <- Pm_ROH_by_site |> subset(CHROM == "chr05") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 5 Position") + scale_x_continuous(labels = label_comma())

chr6plot <- Pm_ROH_by_site |> subset(CHROM == "chr06") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 6 Position") + scale_x_continuous(labels = label_comma())

chr7plot <- Pm_ROH_by_site |> subset(CHROM == "chr07") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 7 Position") + scale_x_continuous(labels = label_comma())

chr8plot <- Pm_ROH_by_site |> subset(CHROM == "chr08") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 8 Position") + scale_x_continuous(labels = label_comma())

chr9plot <- Pm_ROH_by_site |> subset(CHROM == "chr09") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 9 Position") + scale_x_continuous(labels = label_comma())

chr10plot <- Pm_ROH_by_site |> subset(CHROM == "chr10") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 10 Position") + scale_x_continuous(labels = label_comma())

chr11plot <- Pm_ROH_by_site |> subset(CHROM == "chr11") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 11 Position") + scale_x_continuous(labels = label_comma())

chr12plot <- Pm_ROH_by_site |> subset(CHROM == "chr12") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 12 Position") + scale_x_continuous(labels = label_comma())

chr13plot <- Pm_ROH_by_site |> subset(CHROM == "chr13") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 13 Position") + scale_x_continuous(labels = label_comma())

chr14plot <- Pm_ROH_by_site |> subset(CHROM == "chr14") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Population)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 14 Position") + scale_x_continuous(labels = label_comma())

library(patchwork)

chr1plot + chr2plot + chr3plot + chr4plot + chr5plot + chr6plot + chr7plot + chr8plot + chr9plot + chr10plot + chr11plot + chr12plot + chr13plot + chr14plot + guide_area() + plot_layout(nrow = 4, guides = "collect") & theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20), legend.key.size = unit(10, "line"), legend.spacing.y = unit(1, "cm"))

ggsave("chromosome_ROH_plots.png", dpi = 600, width = 30, height = 30, units = "in")

Pm_ROH_by_site_country <- Pm_ROH |> dplyr::group_by(CHROM, POS, Country) |> dplyr::count(STATE) |> dplyr::mutate(Prop_AZ = n/sum(n)) |> subset(STATE ==1) |> dplyr::left_join(Pm_chrlen) |> subset(is.na(Country) == FALSE)|> dplyr::mutate(genome_pos = POS + toadd)

Pm_ROH_by_site_country |> ggplot() + geom_boxplot(aes(x = Country, y = Prop_AZ))

chr1plot_country<- Pm_ROH_by_site_country |> subset(CHROM == "chr01") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Country)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 1 Position") + scale_x_continuous(labels = label_comma())

chr2plot_country<- Pm_ROH_by_site_country |> subset(CHROM == "chr02") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Country)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 2 Position") + scale_x_continuous(labels = label_comma())

chr3plot_country<- Pm_ROH_by_site_country |> subset(CHROM == "chr03") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Country)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 3 Position") + scale_x_continuous(labels = label_comma())

chr4plot_country<- Pm_ROH_by_site_country |> subset(CHROM == "chr04") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Country)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 4 Position") + scale_x_continuous(labels = label_comma())

chr5plot_country<- Pm_ROH_by_site_country |> subset(CHROM == "chr05") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Country)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 5 Position") + scale_x_continuous(labels = label_comma())

chr6plot_country<- Pm_ROH_by_site_country |> subset(CHROM == "chr06") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Country)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 6 Position") + scale_x_continuous(labels = label_comma())

chr7plot_country<- Pm_ROH_by_site_country |> subset(CHROM == "chr07") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Country)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 7 Position") + scale_x_continuous(labels = label_comma())

chr8plot_country<- Pm_ROH_by_site_country |> subset(CHROM == "chr08") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Country)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 8 Position") + scale_x_continuous(labels = label_comma())

chr9plot_country<- Pm_ROH_by_site_country |> subset(CHROM == "chr09") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Country)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 9 Position") + scale_x_continuous(labels = label_comma())

chr10plot_country<- Pm_ROH_by_site_country |> subset(CHROM == "chr10") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Country)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 10 Position") + scale_x_continuous(labels = label_comma())

chr11plot_country<- Pm_ROH_by_site_country |> subset(CHROM == "chr11") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Country)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 11 Position") + scale_x_continuous(labels = label_comma())

chr12plot_country<- Pm_ROH_by_site_country |> subset(CHROM == "chr12") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Country)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 12 Position") + scale_x_continuous(labels = label_comma())

chr13plot_country<- Pm_ROH_by_site_country |> subset(CHROM == "chr13") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Country)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 13 Position") + scale_x_continuous(labels = label_comma())

chr14plot_country<- Pm_ROH_by_site_country |> subset(CHROM == "chr14") |> ggplot() + geom_line(aes(x = POS, y = Prop_AZ, color = Country)) + scale_color_brewer(palette = "Dark2") + theme_linedraw() + theme(panel.grid = element_blank(), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Proportion Autozygous") + xlab("Chromosome 14 Position") + scale_x_continuous(labels = label_comma())

chr1plot_country+ chr2plot_country+ chr3plot_country+ chr4plot_country+ chr5plot_country+ chr6plot_country+ chr7plot_country+ chr8plot_country+ chr9plot_country+ chr10plot_country+ chr11plot_country+ chr12plot_country+ chr13plot_country+ chr14plot_country+ guide_area() + plot_layout(nrow = 4, guides = "collect") & theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20), legend.key.size = unit(5, "line"))

ggsave("chromosome_ROH_plots_countries.png", dpi = 600, width = 32, height = 30, units = "in")
