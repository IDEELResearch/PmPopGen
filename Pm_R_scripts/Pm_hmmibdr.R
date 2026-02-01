###############################################################
###################### Pm_hmmibdr #############################
###############################################################
#Description: uses hidden Markov model to identify genomic segments 
#that are identical by descent

setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

#system("bcftools view -q 0.25 -Oz -o Pm_majors.vcf.gz Pm_monoclonals_wsaf_filtered.vcf.gz")
#system("bcftools view -e 'PLMAF<0.05' -Oz -o Pm_majors_PLMAF_filtered.vcf.gz Pm_majors.vcf.gz")

#subset to monoclonal samples
system("bcftools view -Oz -o Pm_monoclonals_wsaf_filtered.vcf.gz | vcfdo wsaf | -i 'FORMAT/WSAF =0 | FORMAT/WSAF =1' Pm_monoclonals_wsaf.vcf.gz")
system("bcftools view -q 0.25 -e 'PLMAF<0.05' -Oz -o Pm_IBD.vcf.gz Pm_monoclonals_wsaf_filtered.vcf.gz")

#prune intervals
system("bcftools +prune -m 0.25 -w 1000 -Oz -o Pm_IBD_pruned.vcf.gz Pm_IBD.vcf.gz")

#calculate IBD
system("python vcf2hmm.py Pm_IBD_pruned.vcf.gz Pm_IBD")

hmmibdr::hmm_ibd(input_file = "Pm_IBD_seq.txt", allele_freqs = "Pm_IBD_freq.txt", output_file = "Pm_hmmIBD")

#Plot IBD histogram across all samples

Pm_IBD <- data.table::fread("Pm_hmmIBD.hmm_fract.txt")

#plot IBD
library(ggplot2)
library(tidyverse)

main_hist <- Pm_IBD |>
  ggplot() +
  geom_histogram(aes(x=fract_sites_IBD, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "#d9d9d9") +
  xlab("IBD") + ylab("Frequency (%)") +
  theme_classic()

inset_hist <- Pm_IBD |>
  ggplot() +
  geom_histogram(aes(x=fract_sites_IBD, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "#d9d9d9") +
  xlab("IBD") + ylab("Frequency (%)") +
  theme_classic() +
  coord_cartesian(xlim = c(0.5,1), ylim = c(0,3)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

library(patchwork)

IBD_histogram <- main_hist + inset_element(inset_hist, left = 0.5, bottom = 0.4, right = 1, top = 1)

ggsave("IBD_histogram.png", dpi = 600, height = 6, width = 10, units = "in")

Pm_IBD <- Pm_IBD |> dplyr::mutate(country1 = dplyr::case_when(stringr::str_detect(sample1, "Gam_[:digit:]+") ~ "Nigeria",
                                                  stringr::str_detect(sample1, "^[:digit:]+") ~ "DRC",
                                                  stringr::str_detect(sample1, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                  stringr::str_detect(sample1, "Gam") ~ "Cameroon",
                                                  .default = "Tanzania"))

Pm_IBD <- Pm_IBD |> dplyr::mutate(country2 = dplyr::case_when(stringr::str_detect(sample2, "Gam_[:digit:]+") ~ "Nigeria",
                                                              stringr::str_detect(sample2, "^[:digit:]+") ~ "DRC",
                                                              stringr::str_detect(sample2, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                              stringr::str_detect(sample2, "Gam") ~ "Cameroon",
                                                              .default = "Tanzania"))

Pm_IBD_graph <- Pm_IBD |>
  tidygraph::as_tbl_graph(., directed = F) |>
  tidygraph::activate("nodes") |>
  dplyr::mutate(IBD_cluster = as.factor(tidygraph::group_louvain(weights = fract_sites_IBD)),
                Country = dplyr::case_when(stringr::str_detect(name, "Gam_[:digit:]+") ~ "Nigeria",
                                           stringr::str_detect(name, "^[:digit:]+") ~ "DRC",
                                           stringr::str_detect(name, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                           stringr::str_detect(name, "Gam") ~ "Cameroon",
                                           .default = "Tanzania"))

full_graph <- Pm_IBD_graph |>
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_link(aes(width = fract_sites_IBD,
                             color = fract_sites_IBD)) +
  ggraph::geom_node_point(aes(color = Country, shape = IBD_cluster),
                          size = 3) +
  ggraph::scale_edge_width_continuous(range = c(0, 1), guide = "none") +
  #ggraph::geom_node_text(aes(label = name), repel = T) +
  ggraph::scale_edge_color_viridis("IBD", values = c(0,1), option = "plasma", guide = guide_colorbar(available_aes =  "edge_colour")) +
  scale_color_brewer(palette = "Dark2") +
  ggraph::theme_graph() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), legend.position = "bottom") +
  ggtitle("All IBD Pairs")

ggsave("Full_IBD_network.png", full_graph, dpi = 600, height = 10, width = 10, units = "in")

IBD10_graph <- Pm_IBD_graph |>
  tidygraph::activate("edges") |>
  dplyr::filter(fract_sites_IBD >= 0.1) |>
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_link(aes(width = fract_sites_IBD,
                             color = fract_sites_IBD)) +
  ggraph::geom_node_point(aes(color = Country, shape = IBD_cluster),
                          size = 3) +
  ggraph::scale_edge_width_continuous(range = c(0, 1), guide = "none") +
  #ggraph::geom_node_text(aes(label = name), repel = T) +
  ggraph::scale_edge_color_viridis("IBD", values = c(0,1), option = "plasma", guide = "none") +
  scale_color_brewer(palette = "Dark2") +
  ggraph::theme_graph() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), legend.position = "bottom") +
    ggtitle("IBD \u2265 0.1")

ggsave("IBD10_network.png", IBD10_graph, dpi = 600, width = 10, height = 10, units = "in")

IBD25_graph <- Pm_IBD_graph |>
  tidygraph::activate("edges") |>
  dplyr::filter(fract_sites_IBD >= 0.25) |>
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_link(aes(width = fract_sites_IBD,
                             color = fract_sites_IBD)) +
  ggraph::geom_node_point(aes(color = Country, shape = IBD_cluster),
                          size = 3) +
  ggraph::scale_edge_width_continuous(range = c(0, 1), guide = "none") +
  #ggraph::geom_node_text(aes(label = name), repel = T) +
  ggraph::scale_edge_color_viridis("IBD", values = c(0,1), option = "plasma", guide = "none") +
  scale_color_brewer(palette = "Dark2") +
  ggraph::theme_graph() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), legend.position = "bottom") +
  ggtitle("IBD \u2265 0.25")

ggsave("IBD25_network.png", IBD25_graph, dpi = 600, width = 10, height = 10, units = "in")

IBD50_graph <- Pm_IBD_graph |>
  tidygraph::activate("edges") |>
  dplyr::filter(fract_sites_IBD >= 0.5) |>
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_link(aes(width = fract_sites_IBD,
                             color = fract_sites_IBD)) +
  ggraph::geom_node_point(aes(color = Country, shape = IBD_cluster),
                          size = 3) +
  ggraph::scale_edge_width_continuous(range = c(0, 1), guide = "none") +
  #ggraph::geom_node_text(aes(label = name), repel = T) +
  ggraph::scale_edge_color_viridis("IBD", values = c(0,1), option = "plasma", guide = "none") +
  scale_color_brewer(palette = "Dark2") +
  ggraph::theme_graph() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), legend.position = "bottom") +
  ggtitle("IBD \u2265 0.5")

ggsave("IBD50_network.png", IBD50_graph, dpi = 600, width = 10, height = 10, units = "in")

system("bcftools view -q 0.25 -e 'PLMAF<0.05' -Oz -o Pf_IBD.vcf.gz Pf_ortholog_samples_wsaf_filtered.vcf.gz")
system("bcftools +prune -m 0.25 -w 1000 -Oz -o Pf_IBD_pruned.vcf.gz Pf_IBD.vcf.gz")
system("python vcf2hmm.py Pf_IBD_pruned.vcf.gz Pf_IBD")

hmmibdr::hmm_ibd(input_file = "Pf_IBD_seq.txt", allele_freqs = "Pf_IBD_freq.txt", output_file = "Pf_hmmIBD")

Pf_IBD <- data.table::fread("Pf_hmmIBD.hmm_fract.txt")

library(ggplot2)
library(tidyverse)

main_hist_Pf <- Pf_IBD |>
  ggplot() +
  geom_histogram(aes(x=fract_sites_IBD, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "#d9d9d9") +
  xlab("IBD") + ylab("Frequency (%)") +
  theme_classic()

inset_hist_Pf <- Pf_IBD |>
  ggplot() +
  geom_histogram(aes(x=fract_sites_IBD, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "#d9d9d9") +
  xlab("IBD") + ylab("Frequency (%)") +
  theme_classic() +
  coord_cartesian(xlim = c(0.5,1), ylim = c(0,3)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))

library(patchwork)

IBD_histogram_Pf <- main_hist_Pf + inset_element(inset_hist_Pf, left = 0.5, bottom = 0.4, right = 1, top = 1)

ggsave("IBD_histogram_Pf.png", dpi = 600, height = 6, width = 10, units = "in")

Pf_metadata <- readxl::read_xlsx("Pf_ortholog_samples.xlsx")

#Pf_IBD <- Pf_IBD |> dplyr::mutate(sample1 = stringr::str_extract(sample1, "^[:print:]*_"), sample2 = stringr::str_extract(sample2, "^[:print:]*_")) |> dplyr::mutate(sample1 = stringr::str_remove(sample1, "_"), sample2 = stringr::str_remove(sample2, "_"))

Pf_IBD_metadata <- Pf_metadata |> dplyr::mutate(sample1 = Sample, sample2 = Sample)

Pf_IBD <- dplyr::left_join(Pf_IBD, Pf_IBD_metadata, by = "sample1") |> dplyr::mutate(country1 = Country, sample2 = sample2.x) |> dplyr::select(sample1, country1, sample2, N_informative_sites, discordance, log_p, N_fit_iteration, N_generation, N_state_transition, seq_shared_best_traj, fract_sites_IBD, fract_vit_sites_IBD) |> dplyr::left_join(Pf_IBD_metadata, by = "sample2") |> dplyr::mutate(country2 = Country, sample1 = sample1.x) |>  dplyr::select(sample1, sample2, country1, country2, N_informative_sites, discordance, log_p, N_fit_iteration, N_generation, N_state_transition, seq_shared_best_traj, fract_sites_IBD, fract_vit_sites_IBD)

Pf_IBD_graph <- Pf_IBD |>
  tidygraph::as_tbl_graph(., directed = F) |>
  tidygraph::activate("nodes") |>
  dplyr::mutate(IBD_cluster = as.factor(tidygraph::group_louvain(weights = fract_sites_IBD)), 
                Country = Pf_metadata$Country)
  

full_graph_Pf <- Pf_IBD_graph  |>
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_link(aes(width = fract_sites_IBD,
                             color = fract_sites_IBD)) +
  ggraph::geom_node_point(aes(color = Country, shape = IBD_cluster),
                          size = 3) +
  scale_shape_manual(values = c(19, 17, 15, 18, 25, 7, 9), guide = "none") +
  ggraph::scale_edge_width_continuous(range = c(0, 1), guide = "none") +
  #ggraph::geom_node_text(aes(label = name), repel = T) +
  ggraph::scale_edge_color_viridis("IBD", values = c(0,1), option = "plasma", guide = guide_colorbar(available_aes =  "edge_colour")) +
  scale_color_brewer(palette = "Dark2", labels = c("Cameroon", "DRC", "Tanzania")) +
  ggraph::theme_graph() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), legend.position = "bottom") +
  ggtitle("All IBD Pairs")


ggsave("Full_IBD_network_Pf.png", full_graph_Pf, dpi = 600, height = 10, width = 10, units = "in")

IBD10_graph_Pf <- Pf_IBD_graph |>
  tidygraph::activate("edges") |>
  dplyr::filter(fract_sites_IBD >= 0.1) |> 
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_link(aes(width = fract_sites_IBD,
                             color = fract_sites_IBD)) +
  ggraph::geom_node_point(aes(color = Country, shape = IBD_cluster),
                          size = 3) +
  scale_shape_manual(values = c(19, 17, 15, 18, 25, 7, 9), guide = "none") +
  ggraph::scale_edge_width_continuous(range = c(0, 1), guide = "none") +
  #ggraph::geom_node_text(aes(label = name), repel = T) +
  ggraph::scale_edge_color_viridis("IBD", values = c(0,1), option = "plasma", guide = "none") +
  scale_color_brewer(palette = "Dark2", labels = c("Cameroon", "DRC", "Tanzania")) +
  ggraph::theme_graph() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), legend.position = "bottom") +
  ggtitle("IBD \u2265 0.1")

ggsave("IBD10_network_Pf.png", IBD10_graph_Pf, dpi = 600, width = 10, height = 10, units = "in")

IBD25_graph_Pf <- Pf_IBD_graph |>
  tidygraph::activate("edges") |>
  dplyr::filter(fract_sites_IBD >= 0.25) |>
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_link(aes(width = fract_sites_IBD,
                             color = fract_sites_IBD)) +
  ggraph::geom_node_point(aes(color = Country, shape = IBD_cluster),
                          size = 3) +
  scale_shape_manual(values = c(19, 17, 15, 18, 25, 7, 9), guide = "none") +
  ggraph::scale_edge_width_continuous(range = c(0, 1), guide = "none") +
  #ggraph::geom_node_text(aes(label = name), repel = T) +
  ggraph::scale_edge_color_viridis("IBD", values = c(0,1), option = "plasma", guide = "none") +
  scale_color_brewer(palette = "Dark2", labels = c("Cameroon", "DRC", "Tanzania")) +
  ggraph::theme_graph() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), legend.position = "bottom") +
  ggtitle("IBD \u2265 0.25")

ggsave("IBD25_network_Pf.png", IBD25_graph_Pf, dpi = 600, width = 10, height = 10, units = "in")

IBD50_graph_Pf <- Pf_IBD_graph |>
  tidygraph::activate("edges") |>
  dplyr::filter(fract_sites_IBD >= 0.5) |>
  ggraph::ggraph(layout = 'kk') +
  ggraph::geom_edge_link(aes(width = fract_sites_IBD,
                             color = fract_sites_IBD)) +
  ggraph::geom_node_point(aes(color = Country, shape = IBD_cluster),
                          size = 3) +
  scale_shape_manual(values = c(19, 17, 15, 18, 25, 7, 9)) +
  ggraph::scale_edge_width_continuous(range = c(0, 1), guide = "none") +
  #ggraph::geom_node_text(aes(label = name), repel = T) +
  ggraph::scale_edge_color_viridis("IBD", values = c(0,1), option = "plasma", guide = "none") +
  scale_color_brewer(palette = "Dark2", labels = c("Cameroon", "DRC", "Tanzania")) +
  ggraph::theme_graph() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), legend.position = "bottom") +
  ggtitle("IBD \u2265 0.5")

ggsave("IBD50_network_Pf.png", IBD50_graph_Pf, dpi = 600, width = 10, height = 10, units = "in")

design <- "
AAAABBCC
AAAABBCC
AAAAEEDD
AAAAEEDD
"
Pm_IBD_plot <- full_graph + IBD10_graph + IBD25_graph + IBD50_graph + guide_area() + plot_layout(design = design, guides = "collect") + plot_annotation(title = expression(italic("P. malariae")), theme = theme(plot.title = element_text(face = "bold", size = 36)), tag_levels = "A") & theme(legend.text = element_text(size = 18), legend.key.size = unit(1, "cm"), legend.direction = "horizontal", legend.title = element_text(size = 20), plot.tag = element_text(face = "bold", size = 36))

ggsave("Pm_IBD_plot.png", Pm_IBD_plot, dpi = 600, width = 24, height = 12, units = "in")

Pf_IBD_plot <- full_graph_Pf + IBD10_graph_Pf + IBD25_graph_Pf + IBD50_graph_Pf + guide_area() + plot_layout(design = design, guides = "collect") + plot_annotation(title = expression(italic("P. falciparum")), theme = theme(plot.title = element_text(face = "bold", size = 36)), tag_levels = "A") & theme(legend.text = element_text(size = 18), legend.key.size = unit(1.5, "cm"), legend.direction = "horizontal", legend.title = element_text(size = 20), plot.tag = element_text(face = "bold", size = 36))

ggsave("Pf_IBD_plot.png", Pf_IBD_plot, dpi = 600, width = 24, height = 12, units = "in")

combined_IBD_plot <- Pm_IBD_plot / Pf_IBD_plot
 
ggsave("combined_IBD_plot.png", combined_IBD_plot, dpi = 600, width = 24, height = 24, units = "in")

