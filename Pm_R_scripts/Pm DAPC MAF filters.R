setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

#MOST RECENT RUN WITH INTERMEDIATE WSAF REMOVED

Pm_wsaf_filtered <- vcfR::read.vcfR("Pm_monoclonals_wsaf_filtered.vcf.gz") |> vcfR::vcfR2genlight()

Pm_clust <- adegenet::find.clusters(Pm_wsaf_filtered, n.pca = 50, n.clust = 7)
#50 PCs, 7 clusters based on BIC

Pm_DAPC <- adegenet::dapc(Pm_wsaf_filtered, Pm_clust$grp, n.pca = 50, n.da = 6)

Pm_scatter <- ggplot(Pm_DAPC$ind.coord, aes(x = LD1, y = LD2, color = Pm_clust$grp)) +
  geom_vline(xintercept = 0, linewidth = 2) + geom_hline(yintercept = 0, linewidth = 2) +
  geom_point(size = 10, shape = 20) +
  scale_color_brewer(palette = "Set2", labels = c("Cluster 1 ", "Cluster 2 ", "Cluster 3 ", "Cluster 4 ", "Cluster 5 ", "Cluster 6 ", "Cluster 7 ")) +
  theme_void() + 
  theme(panel.background = element_rect(color = "black", linewidth = 2), legend.position = "inside", legend.position.inside = c(0.89, 0.84), 
        legend.background = element_rect(color = "black", linewidth = 1), legend.text = element_text(size = 20), legend.title = element_blank(), legend.key = element_blank())
  
#myCol <- RColorBrewer::brewer.pal(7, "Set2")

#svg(filename = "wsaf_filtered_DAPC_clusters.svg", width = 10, height = 10)

#Pm_scatter <- adegenet::scatter.dapc(Pm_DAPC, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=1,
#                                        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:7))

#dev.off()

#svg(filename = "wsaf_filtered_DAPC_discriminant.svg", width = 10, height = 10)

#Pm_pops <- adegenet::scatter.dapc(Pm_DAPC,1,2, col=myCol, bg="white",
#                                     scree.da=FALSE, legend=TRUE, solid=1)

#dev.off()

Pm_group_assignments <- Pm_DAPC$grp |> as.data.frame() |> tibble::rownames_to_column(var = "Sample")

Pm_group_assignments |> saveRDS("DAPC_clusters.rds")

Pm_group_assignments |> dplyr::count(Pm_DAPC$grp)

Pm_group_assignments <- Pm_group_assignments |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                               stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                               stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                               stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                               .default = "Tanzania"))

Pm_clusters_geog <- Pm_group_assignments |> ggplot() + geom_bar(aes(x = Pm_DAPC$grp, fill = Country)) + scale_fill_brewer(palette = "Dark2", name = "Country") + theme_classic() + xlab("Cluster") + ylab("Count") + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("wsaf_filtered_clusters_geography.png", Pm_clusters_geog, dpi = 600)

#Pm_scatter <- magick::image_read_svg("wsaf_filtered_DAPC_clusters.svg", width = 720) |> magick::image_ggplot(interpolate = TRUE) + theme(text = element_text(size = 16), plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"))

#Pm_pops <- magick::image_read_svg("wsaf_filtered_DAPC_discriminant.svg") |> magick::image_trim() |> magick::image_ggplot() + theme(text = element_text(size = 16), plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"))

Pm_wsaf_filtered <- Pm_scatter + plot_spacer() + Pm_clusters_geog + plot_annotation(tag_levels = "A") + plot_layout(widths = c(4, 0.5, 4))

ggsave("wsaf_filtered_DAPC.png", Pm_wsaf_filtered, dpi = 600, width = 20, height = 10, units = "in")


MAF5 <- vcfR::read.vcfR("Pm_monoclonal_exclude_MAF5.vcf.gz") |> vcfR::vcfR2genlight()

MAF5_clust <- adegenet::find.clusters(MAF5)
#50 PCs, 3 clusters based on BIC

MAF5_DAPC <- adegenet::dapc(MAF5, MAF5_clust$grp)

MAF5_scatter <- adegenet::scatter.dapc(MAF5_DAPC)

myCol <- RColorBrewer::brewer.pal(3, "Set2")

MAF5_scatter <- adegenet::scatter.dapc(MAF5_DAPC, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
                                     cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3))

#had to manually save as svg - can't figure out how to get into ggplot

MAF5_pops <- adegenet::scatter.dapc(MAF5_DAPC,1,1, col=myCol, bg="white",
                                  scree.da=FALSE, legend=TRUE, solid=.4)

#had to manually save as svg - can't figure out how to get into ggplot

MAF5_group_assignments <- MAF5_DAPC$grp |> as.data.frame() |> tibble::rownames_to_column(var = "Sample")

MAF5_group_assignments <- MAF5_group_assignments |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                   stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                   stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                   stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                   .default = "Tanzania"))

MAF5_clusters_geog <- MAF5_group_assignments |> ggplot() + geom_bar(aes(x = MAF5_DAPC$grp, fill = Country)) + scale_fill_brewer(palette = "Dark2", name = "Country") + theme_classic() + xlab("Cluster") + ylab("Count") + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("MAF5_clusters_geography.png", MAF5_clusters_geog, dpi = 600)

library(patchwork)
library(magick)

MAF5_scatter <- magick::image_read_svg("MAF5_DAPC_clusters.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()

MAF5_pops <- magick::image_read_svg("MAF5_DAPC_discriminant.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()

design <- "
AAAA##
AAAACC
BBBBCC
BBBB##
"

MAF5 <- MAF5_scatter + MAF5_pops + MAF5_clusters_geog + plot_layout(design = design, width = 1)

ggsave("MAF5_DAPC.png", MAF5, dpi = 600)

#redid this one incorporating a WSAF filter as well

MAF10 <- vcfR::read.vcfR("Pm_monoclonal_exclude_MAF10_WSAF10.vcf.gz") |> vcfR::vcfR2genlight()

MAF10_clust <- adegenet::find.clusters(MAF10, n.pca = 50)
#50 PCs, 3 clusters based on BIC

MAF10_DAPC <- adegenet::dapc(MAF10, MAF10_clust$grp, n.pca = 50, n.da = 2)

myCol <- RColorBrewer::brewer.pal(3, "Set2")

MAF10_scatter <- adegenet::scatter.dapc(MAF10_DAPC, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
                                       cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3))

MAF10_pops <- adegenet::scatter.dapc(MAF10_DAPC,1,1, col=myCol, bg="white",
                                    scree.da=FALSE, legend=TRUE, solid=.4)

MAF10_group_assignments <- MAF10_DAPC$grp |> as.data.frame() |> tibble::rownames_to_column(var = "Sample")

MAF10_group_assignments <- MAF10_group_assignments |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                             stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                             stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                             stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                             .default = "Tanzania"))

MAF10_clusters_geog <- MAF10_group_assignments |> ggplot() + geom_bar(aes(x = MAF10_DAPC$grp, fill = Country)) + scale_fill_brewer(palette = "Dark2", name = "Country") + theme_classic() + xlab("Cluster") + ylab("Count") + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("MAF10_clusters_geography.png", MAF10_clusters_geog, dpi = 600)

MAF10_scatter <- magick::image_read_svg("MAF10_DAPC_clusters.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()

MAF10_pops <- magick::image_read_svg("MAF10_DAPC_discriminant.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()


MAF10 <- MAF10_scatter + MAF10_pops + MAF10_clusters_geog + plot_layout(design = design, width = 1)

ggsave("MAF10_DAPC.png", MAF10, dpi = 600)

dadi_DAPC <- MAF10_group_assignments |> dplyr::rename(Group = "MAF10_DAPC$grp" ) |> dplyr::select(Sample, Group)

dadi_DAPC |> data.table::fwrite("dadi_DAPC.txt", sep = "\t")

MAF15 <- vcfR::read.vcfR("Pm_monoclonal_exclude_MAF15.vcf.gz") |> vcfR::vcfR2genlight()

MAF15_clust <- adegenet::find.clusters(MAF15)
#50 PCs, 3 clusters based on BIC

MAF15_DAPC <- adegenet::dapc(MAF15, MAF15_clust$grp)

myCol <- RColorBrewer::brewer.pal(3, "Set2")

svg(filename = "MAF15_DAPC_clusters.svg")

MAF15_scatter <- adegenet::scatter.dapc(MAF15_DAPC, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
                                        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3))

dev.off()

svg(filename = "MAF15_DAPC_discriminant.svg")

MAF15_pops <- adegenet::scatter.dapc(MAF15_DAPC,1,1, col=myCol, bg="white",
                                     scree.da=FALSE, legend=TRUE, solid=.4)

dev.off()

MAF15_group_assignments <- MAF15_DAPC$grp |> as.data.frame() |> tibble::rownames_to_column(var = "Sample")

MAF15_group_assignments <- MAF15_group_assignments |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                               stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                               stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                               stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                               .default = "Tanzania"))

MAF15_clusters_geog <- MAF15_group_assignments |> ggplot() + geom_bar(aes(x = MAF15_DAPC$grp, fill = Country)) + scale_fill_brewer(palette = "Dark2", name = "Country") + theme_classic() + xlab("Cluster") + ylab("Count") + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("MAF15_clusters_geography.png", MAF15_clusters_geog, dpi = 600)

MAF15_scatter <- magick::image_read_svg("MAF15_DAPC_clusters.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()

MAF15_pops <- magick::image_read_svg("MAF15_DAPC_discriminant.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()


MAF15 <- MAF15_scatter + MAF15_pops + MAF15_clusters_geog + plot_layout(design = design, width = 1)

ggsave("MAF15_DAPC.png", MAF15, dpi = 600)

MAF20 <- vcfR::read.vcfR("Pm_monoclonal_exclude_MAF20.vcf.gz") |> vcfR::vcfR2genlight()

MAF20_clust <- adegenet::find.clusters(MAF20)
#50 PCs, 3 clusters based on BIC

MAF20_DAPC <- adegenet::dapc(MAF20, MAF20_clust$grp)

myCol <- RColorBrewer::brewer.pal(3, "Set2")

svg(filename = "MAF20_DAPC_clusters.svg")

MAF20_scatter <- adegenet::scatter.dapc(MAF20_DAPC, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
                                        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3))

dev.off()

svg(filename = "MAF20_DAPC_discriminant.svg")

MAF20_pops <- adegenet::scatter.dapc(MAF20_DAPC,1,1, col=myCol, bg="white",
                                     scree.da=FALSE, legend=TRUE, solid=.4)

dev.off()

MAF20_group_assignments <- MAF20_DAPC$grp |> as.data.frame() |> tibble::rownames_to_column(var = "Sample")

MAF20_group_assignments <- MAF20_group_assignments |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                               stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                               stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                               stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                               .default = "Tanzania"))

MAF20_clusters_geog <- MAF20_group_assignments |> ggplot() + geom_bar(aes(x = MAF20_DAPC$grp, fill = Country)) + scale_fill_brewer(palette = "Dark2", name = "Country") + theme_classic() + xlab("Cluster") + ylab("Count") + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("MAF20_clusters_geography.png", MAF20_clusters_geog, dpi = 600)

MAF20_scatter <- magick::image_read_svg("MAF20_DAPC_clusters.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()

MAF20_pops <- magick::image_read_svg("MAF20_DAPC_discriminant.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()


MAF20 <- MAF20_scatter + MAF20_pops + MAF20_clusters_geog + plot_layout(design = design, width = 1)

ggsave("MAF20_DAPC.png", MAF20, dpi = 600)

MAF25 <- vcfR::read.vcfR("Pm_monoclonal_exclude_MAF25.vcf.gz") |> vcfR::vcfR2genlight()

MAF25_clust <- adegenet::find.clusters(MAF25)
#50 PCs, 3 clusters based on BIC

MAF25_DAPC <- adegenet::dapc(MAF25, MAF25_clust$grp, n.pca = 50, n.da = 2)

myCol <- RColorBrewer::brewer.pal(3, "Set2")

svg(filename = "MAF25_DAPC_clusters.svg")

MAF25_scatter <- adegenet::scatter.dapc(MAF25_DAPC, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
                                        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3))

dev.off()

svg(filename = "MAF25_DAPC_discriminant.svg")

MAF25_pops <- adegenet::scatter.dapc(MAF25_DAPC,1,1, col=myCol, bg="white",
                                     scree.da=FALSE, legend=TRUE, solid=.4)

dev.off()

MAF25_group_assignments <- MAF25_DAPC$grp |> as.data.frame() |> tibble::rownames_to_column(var = "Sample")

MAF25_group_assignments <- MAF25_group_assignments |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                               stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                               stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                               stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                               .default = "Tanzania"))

MAF25_clusters_geog <- MAF25_group_assignments |> ggplot() + geom_bar(aes(x = MAF25_DAPC$grp, fill = Country)) + scale_fill_brewer(palette = "Dark2", name = "Country") + theme_classic() + xlab("Cluster") + ylab("Count") + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("MAF25_clusters_geography.png", MAF25_clusters_geog, dpi = 600)

MAF25_scatter <- magick::image_read_svg("MAF25_DAPC_clusters.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()

MAF25_pops <- magick::image_read_svg("MAF25_DAPC_discriminant.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()


MAF25 <- MAF25_scatter + MAF25_pops + MAF25_clusters_geog + plot_layout(design = design, width = 1)

ggsave("MAF25_DAPC.png", MAF25, dpi = 600)

MAF30 <- vcfR::read.vcfR("Pm_monoclonal_exclude_MAF30.vcf.gz") |> vcfR::vcfR2genlight()

MAF30_clust <- adegenet::find.clusters(MAF30, n.pca = 50)
#50 PCs, 3 clusters based on BIC

MAF30_DAPC <- adegenet::dapc(MAF30, MAF30_clust$grp, n.pca = 50, n.da = 2)

myCol <- RColorBrewer::brewer.pal(3, "Set2")

svg(filename = "MAF30_DAPC_clusters.svg")

MAF30_scatter <- adegenet::scatter.dapc(MAF30_DAPC, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
                                        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3))

dev.off()

svg(filename = "MAF30_DAPC_discriminant.svg")

MAF30_pops <- adegenet::scatter.dapc(MAF30_DAPC,1,1, col=myCol, bg="white",
                                     scree.da=FALSE, legend=TRUE, solid=.4)

dev.off()

MAF30_group_assignments <- MAF30_DAPC$grp |> as.data.frame() |> tibble::rownames_to_column(var = "Sample")

MAF30_group_assignments <- MAF30_group_assignments |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                               stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                               stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                               stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                               .default = "Tanzania"))

MAF30_clusters_geog <- MAF30_group_assignments |> ggplot() + geom_bar(aes(x = MAF30_DAPC$grp, fill = Country)) + scale_fill_brewer(palette = "Dark2", name = "Country") + theme_classic() + xlab("Cluster") + ylab("Count") + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("MAF30_clusters_geography.png", MAF30_clusters_geog, dpi = 600)

MAF30_scatter <- magick::image_read_svg("MAF30_DAPC_clusters.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()

MAF30_pops <- magick::image_read_svg("MAF30_DAPC_discriminant.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()


MAF30 <- MAF30_scatter + MAF30_pops + MAF30_clusters_geog + plot_layout(design = design, width = 1)

ggsave("MAF30_DAPC.png", MAF30, dpi = 600)

MAF35 <- vcfR::read.vcfR("Pm_monoclonal_exclude_MAF35.vcf.gz") |> vcfR::vcfR2genlight()

MAF35_clust <- adegenet::find.clusters(MAF35, n.pca = 50)
#50 PCs, 3 clusters based on BIC

MAF35_DAPC <- adegenet::dapc(MAF35, MAF35_clust$grp, n.pca = 50, n.da = 2)

myCol <- RColorBrewer::brewer.pal(3, "Set2")

svg(filename = "MAF35_DAPC_clusters.svg")

MAF35_scatter <- adegenet::scatter.dapc(MAF35_DAPC, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
                                        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3))

dev.off()

svg(filename = "MAF35_DAPC_discriminant.svg")

MAF35_pops <- adegenet::scatter.dapc(MAF35_DAPC,1,1, col=myCol, bg="white",
                                     scree.da=FALSE, legend=TRUE, solid=.4)

dev.off()

MAF35_group_assignments <- MAF35_DAPC$grp |> as.data.frame() |> tibble::rownames_to_column(var = "Sample")

MAF35_group_assignments <- MAF35_group_assignments |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                               stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                               stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                               stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                               .default = "Tanzania"))

MAF35_clusters_geog <- MAF35_group_assignments |> ggplot() + geom_bar(aes(x = MAF35_DAPC$grp, fill = Country)) + scale_fill_brewer(palette = "Dark2", name = "Country") + theme_classic() + xlab("Cluster") + ylab("Count") + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("MAF35_clusters_geography.png", MAF35_clusters_geog, dpi = 600)

MAF35_scatter <- magick::image_read_svg("MAF35_DAPC_clusters.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()

MAF35_pops <- magick::image_read_svg("MAF35_DAPC_discriminant.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()


MAF35 <- MAF35_scatter + MAF35_pops + MAF35_clusters_geog + plot_layout(design = design, width = 1)

ggsave("MAF35_DAPC.png", MAF35, dpi = 600)

MAF40 <- vcfR::read.vcfR("Pm_monoclonal_exclude_MAF40.vcf.gz") |> vcfR::vcfR2genlight()

MAF40_clust <- adegenet::find.clusters(MAF40, n.pca = 50)
#50 PCs, 3 clusters based on BIC

MAF40_DAPC <- adegenet::dapc(MAF40, MAF40_clust$grp, n.pca = 50, n.da = 2)

myCol <- RColorBrewer::brewer.pal(3, "Set2")

svg(filename = "MAF40_DAPC_clusters.svg")

MAF40_scatter <- adegenet::scatter.dapc(MAF40_DAPC, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
                                        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3))

dev.off()

svg(filename = "MAF40_DAPC_discriminant.svg")

MAF40_pops <- adegenet::scatter.dapc(MAF40_DAPC,1,1, col=myCol, bg="white",
                                     scree.da=FALSE, legend=TRUE, solid=.4)

dev.off()

MAF40_group_assignments <- MAF40_DAPC$grp |> as.data.frame() |> tibble::rownames_to_column(var = "Sample")

MAF40_group_assignments <- MAF40_group_assignments |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                               stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                               stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                               stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                               .default = "Tanzania"))

MAF40_clusters_geog <- MAF40_group_assignments |> ggplot() + geom_bar(aes(x = MAF40_DAPC$grp, fill = Country)) + scale_fill_brewer(palette = "Dark2", name = "Country") + theme_classic() + xlab("Cluster") + ylab("Count") + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("MAF40_clusters_geography.png", MAF40_clusters_geog, dpi = 600)

MAF40_scatter <- magick::image_read_svg("MAF40_DAPC_clusters.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()

MAF40_pops <- magick::image_read_svg("MAF40_DAPC_discriminant.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()


MAF40 <- MAF40_scatter + MAF40_pops + MAF40_clusters_geog + plot_layout(design = design, width = 1)

ggsave("MAF40_DAPC.png", MAF40, dpi = 600)

MAF45 <- vcfR::read.vcfR("Pm_monoclonal_exclude_MAF45.vcf.gz") |> vcfR::vcfR2genlight()

MAF45_clust <- adegenet::find.clusters(MAF45, n.pca = 50)
#50 PCs, 3 clusters based on BIC

MAF45_DAPC <- adegenet::dapc(MAF45, MAF45_clust$grp, n.pca = 50, n.da = 2)

myCol <- RColorBrewer::brewer.pal(3, "Set2")

svg(filename = "MAF45_DAPC_clusters.svg")

MAF45_scatter <- adegenet::scatter.dapc(MAF45_DAPC, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
                                        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3))

dev.off()

svg(filename = "MAF45_DAPC_discriminant.svg")

MAF45_pops <- adegenet::scatter.dapc(MAF45_DAPC,1,1, col=myCol, bg="white",
                                     scree.da=FALSE, legend=TRUE, solid=.4)

dev.off()

MAF45_group_assignments <- MAF45_DAPC$grp |> as.data.frame() |> tibble::rownames_to_column(var = "Sample")

MAF45_group_assignments <- MAF45_group_assignments |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                               stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                               stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                               stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                               .default = "Tanzania"))

MAF45_clusters_geog <- MAF45_group_assignments |> ggplot() + geom_bar(aes(x = MAF45_DAPC$grp, fill = Country)) + scale_fill_brewer(palette = "Dark2", name = "Country") + theme_classic() + xlab("Cluster") + ylab("Count") + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("MAF45_clusters_geography.png", MAF45_clusters_geog, dpi = 600)

MAF45_scatter <- magick::image_read_svg("MAF45_DAPC_clusters.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()

MAF45_pops <- magick::image_read_svg("MAF45_DAPC_discriminant.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()


MAF45 <- MAF45_scatter + MAF45_pops + MAF45_clusters_geog + plot_layout(design = design, width = 1)

ggsave("MAF45_DAPC.png", MAF45, dpi = 600)

MAF49 <- vcfR::read.vcfR("Pm_monoclonal_exclude_MAF50.vcf.gz") |> vcfR::vcfR2genlight()

MAF49_clust <- adegenet::find.clusters(MAF49, n.pca = 50)
#50 PCs, 3 clusters based on BIC

MAF49_DAPC <- adegenet::dapc(MAF49, MAF49_clust$grp, n.pca = 50, n.da = 2)

myCol <- RColorBrewer::brewer.pal(3, "Set2")

svg(filename = "MAF49_DAPC_clusters.svg")

MAF49_scatter <- adegenet::scatter.dapc(MAF49_DAPC, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
                                        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3))

dev.off()

svg(filename = "MAF49_DAPC_discriminant.svg")

MAF49_pops <- adegenet::scatter.dapc(MAF49_DAPC,1,1, col=myCol, bg="white",
                                     scree.da=FALSE, legend=TRUE, solid=.4)

dev.off()

MAF49_group_assignments <- MAF49_DAPC$grp |> as.data.frame() |> tibble::rownames_to_column(var = "Sample")

MAF49_group_assignments <- MAF49_group_assignments |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                               stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                               stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                               stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                               .default = "Tanzania"))

MAF49_clusters_geog <- MAF49_group_assignments |> ggplot() + geom_bar(aes(x = MAF49_DAPC$grp, fill = Country)) + scale_fill_brewer(palette = "Dark2", name = "Country") + theme_classic() + xlab("Cluster") + ylab("Count") + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("MAF49_clusters_geography.png", MAF49_clusters_geog, dpi = 600)

MAF49_scatter <- magick::image_read_svg("MAF49_DAPC_clusters.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()

MAF49_pops <- magick::image_read_svg("MAF49_DAPC_discriminant.svg", width = 4767) |> magick::image_trim() |> magick::image_ggplot()


MAF49 <- MAF49_scatter + MAF49_pops + MAF49_clusters_geog + plot_layout(design = design, width = 1)

ggsave("MAF49_DAPC.png", MAF49, dpi = 600)