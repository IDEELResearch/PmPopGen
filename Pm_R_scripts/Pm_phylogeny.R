###############################################################
##################### Pm_phylogeny ############################
###############################################################
#Description: generates maximum likelihood phylogenetic trees of
#Pm samples using RaxML

setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

Pm_clusters <- fastreeR::vcf2clusters("Pm_IBD_pruned.vcf.gz")

Pm_newick <- Pm_clusters[2]

Pm_newick |> data.table::fwrite("Pm_newick.txt")

Pm_tree <- ape::read.tree("Pm_newick.txt")

plot(Pm_tree, direction = "down", cex = 0.5)
ape::add.scale.bar()
ape::axisPhylo(side = 2)

Pm_dist <- fastreeR::vcf2dist("Pm_IBD_pruned.vcf.gz")

stats_tree <- stats::hclust(Pm_dist)
plot(stats_tree, ann = FALSE, cex = 0.7)

system("python vcf2phylip.py -i Pm_IBD_pruned.vcf.gz")

system("raxml-ng/bin/raxml-ng --all --msa Pm_IBD_pruned.min4.phy --model LG+G8+F --tree pars{10} --bs-t
rees 200") #This will take a while - better to run it outside of R 

ml_tree <- ape::read.tree("Pm_IBD_pruned.min4.phy.raxml.bestTree")

ml_tree_collapsed <- ape::read.tree("Pm_IBD_pruned.min4.phy.raxml.bestTreeCollapsed")

plot(ml_tree, direction = "down", cex = 0.7)
ape::add.scale.bar()
ape::axisPhylo(side = 2)

plot(ml_tree_collapsed, cex = 0.7)

library(ggtree)
library(treeio)

ml_tibble <- as_tibble(ml_tree_collapsed)

ml_tibble <- ml_tibble |> dplyr::mutate(country = dplyr::case_when(stringr::str_detect(label, "Gam_[:digit:]+") ~ "Nigeria",
                                                                   stringr::str_detect(label, "^[:digit:]+") ~ "DRC",
                                                                   stringr::str_detect(label, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                   stringr::str_detect(label, "Gam") ~ "Cameroon",
                                                                   .default = "Tanzania"))

ml_tibble$sample <- ml_tibble$label

#ml_tibble$label <- ml_tibble$country

DAPC_clusters <- readRDS("DAPC_clusters.rds") |> dplyr::rename(cluster = "Pm_DAPC$grp", sample = Sample)

ml_tibble <- ml_tibble |> dplyr::left_join(DAPC_clusters)

admixture_data <- readRDS("admixture_data.rds") |> dplyr::mutate(Population = dplyr::case_when(Pop1 > 0.6 ~ "Population 1",
                                                                                               Pop2 > 0.6 ~ "Population 2",
                                                                                               .default = "Admixed"))

ml_tibble <- ml_tibble |> dplyr::left_join(admixture_data, by = dplyr::join_by(sample == Sample))

ml_treedata <- ml_tibble |> as.treedata()

library(ggnewscale)

color_palette <- viridis::viridis_pal(option = "turbo")(5)

ml_tree <- ml_treedata |> ggtree(ladderize = FALSE, aes(color = factor(interaction(isTip, country))), size = 2) + geom_treescale(fontsize = 8, linesize = 2) + 
  scale_color_manual(breaks = c("TRUE.Cameroon", "TRUE.DRC", "TRUE.Nigeria", "TRUE.Tanzania"), values = color_palette, name = "Country", labels = c("Cameroon", "DRC", "Nigeria", "Tanzania")) + geom_tiplab(aes(label = Population), show.legend = FALSE, size = 8, offset = 0.01) + new_scale_color() + 
  geom_tippoint(aes(color = cluster), size = 5) + scale_color_brewer(palette = "Set2", name = "DAPC Cluster") + theme(text = element_text(size = 24)) + expand_limits(x = 4)

ggsave("ml_tree.png", ml_tree, dpi = 600, width = 9600, height = 12000, units = "px")

API_tree_collapsed <- ape::read.tree("Pm_API.min4.phy.raxml.bestTreeCollapsed")

API_tibble <- as_tibble(API_tree_collapsed)

API_tibble <- API_tibble |> dplyr::mutate(country = dplyr::case_when(stringr::str_detect(label, "Gam_[:digit:]+") ~ "Nigeria",
                                                                     stringr::str_detect(label, "^[:digit:]+") ~ "DRC",
                                                                     stringr::str_detect(label, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                     stringr::str_detect(label, "Gam") ~ "Cameroon",
                                                                     .default = "Tanzania"))
API_treedata <- API_tibble |> as.treedata()

API_tree <- API_treedata |> ggtree(ladderize = FALSE, aes(color = factor(interaction(isTip, country))), size = 2) + geom_treescale(fontsize = 8, linesize = 2) + 
  scale_color_manual(breaks = c("TRUE.Cameroon", "TRUE.DRC", "TRUE.Nigeria", "TRUE.Tanzania"), values = color_palette, name = "Country", labels = c("Cameroon", "DRC", "Nigeria", "Tanzania")) + geom_tiplab(aes(label = country), show.legend = FALSE, size = 8, offset = 0.01) + new_scale_color() +
  geom_tippoint(aes(color = country), size = 5, show.legend = FALSE) + scale_color_manual(values = color_palette[c(1:3, 5)]) + theme(text = element_text(size = 24)) + expand_limits(x = 0.4)

ggsave("API_tree.png", API_tree, dpi = 600, width = 9600, height = 12000, units = "px")

MIT_tree_collapsed <- ape::read.tree("Pm_MIT.min4.phy.raxml.bestTreeCollapsed")

MIT_tibble <- as_tibble(MIT_tree_collapsed)

MIT_tibble <- MIT_tibble |> dplyr::mutate(country = dplyr::case_when(stringr::str_detect(label, "Gam_[:digit:]+") ~ "Nigeria",
                                                                     stringr::str_detect(label, "^[:digit:]+") ~ "DRC",
                                                                     stringr::str_detect(label, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                     stringr::str_detect(label, "Gam") ~ "Cameroon",
                                                                     .default = "Tanzania"))
MIT_tibble <- MIT_tibble |> dplyr::mutate(branch.length = dplyr::case_when(branch.length > 1 ~ 1,
                                                                           .default = branch.length))

MIT_treedata <- MIT_tibble |> as.treedata()

#library(ggbreak)

MIT_tree <- MIT_treedata |> ggtree(ladderize = FALSE, aes(color = factor(interaction(isTip, country))), size = 2) + geom_treescale(fontsize = 8, linesize = 2) + 
  scale_color_manual(breaks = c("TRUE.Cameroon", "TRUE.DRC", "TRUE.Nigeria", "TRUE.Tanzania"), values = color_palette, name = "Country", labels = c("Cameroon", "DRC", "Nigeria", "Tanzania")) + geom_tiplab(aes(label = country), show.legend = FALSE, size = 8, offset = 0.01) + new_scale_color() +
  geom_tippoint(aes(color = country), size = 5, show.legend = FALSE) + scale_color_manual(values = color_palette[c(1:3, 5)]) + theme(text = element_text(size = 24)) + expand_limits(x = 1.5)

ggsave("MIT_tree.png", MIT_tree, dpi = 600, width = 9600, height = 12000, units = "px")

library(patchwork)

organelle_trees <- API_tree + MIT_tree + plot_annotation(tag_levels = "A")

ggsave("organelle_trees.png", organelle_trees, dpi = 600, width = 19200, height = 12000, units = "px")

