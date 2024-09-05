setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

Pm_clusters <- fastreeR::vcf2clusters("Pm_monoclonals_wsaf.vcf.gz")

Pm_newick <- Pm_clusters[2]

Pm_newick |> data.table::fwrite("Pm_newick.txt")

Pm_tree <- ape::read.tree("Pm_newick.txt")

plot(Pm_tree, direction = "down", cex = 0.5)
ape::add.scale.bar()
ape::axisPhylo(side = 2)

Pm_dist <- fastreeR::vcf2dist("Pm_monoclonals_wsaf.vcf.gz")

stats_tree <- stats::hclust(Pm_dist)
plot(stats_tree, ann = FALSE, cex = 0.7)

Pm_clusters2 <- fastreeR::vcf2clusters("Pm_monoclonals_wsaf_filtered.vcf.gz")

Pm_newick2 <- Pm_clusters2[2]

Pm_newick2 |> data.table::fwrite("Pm_newick2.txt")

Pm_tree2 <- ape::read.tree("Pm_newick2.txt")

plot(Pm_tree2, direction = "down", cex = 0.7)
ape::add.scale.bar()
ape::axisPhylo(side = 2)

Pm_dist2 <- fastreeR::vcf2dist("Pm_monoclonals_wsaf_filtered.vcf.gz")

stats_tree2 <- stats::hclust(Pm_dist2)
plot(stats_tree2, ann = FALSE, cex = 0.7)

ml_tree <- ape::read.tree("Pm_monoclonals_wsaf_filtered.min4.phy.raxml.bestTree")

ml_tree_collapsed <- ape::read.tree("Pm_monoclonals_wsaf_filtered.min4.phy.raxml.bestTreeCollapsed")

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

ml_tibble$label <- ml_tibble$country

DAPC_clusters <- readRDS("DAPC_clusters.rds") |> dplyr::rename(cluster = "Pm_DAPC$grp", sample = Sample)

ml_tibble <- ml_tibble |> dplyr::left_join(DAPC_clusters)

ml_treedata <- ml_tibble |> as.treedata()

ml_treedata |> ggtree(ladderize = FALSE, aes(color = country)) + geom_treescale() + viridis::scale_color_viridis(discrete = TRUE, option = "turbo", name = "Country")
