setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

Pm_VCF <- vcfR::read.vcfR("Pm_PCA.vcf.gz") |> vcfR::vcfR2genlight()

Pm_clust <- adegenet::find.clusters(Pm_VCF, n.pca = 50, n.clust = 3)
#50 PCs, 3 clusters based on BIC

Pm_DAPC <- adegenet::dapc(Pm_VCF, Pm_clust$grp, n.pca = 50, n.da = 2)

Pm_scatter <- ggplot(Pm_DAPC$ind.coord, aes(x = LD1, y = LD2, color = Pm_clust$grp)) +
  geom_vline(xintercept = 0, linewidth = 2) + geom_hline(yintercept = 0, linewidth = 2) +
  geom_point(size = 10, shape = 20) +
  scale_color_brewer(palette = "Set2", labels = c("Cluster 1 ", "Cluster 2 ", "Cluster 3 ", "Cluster 4 ", "Cluster 5 ", "Cluster 6 ", "Cluster 7 ")) +
  theme_void() + 
  theme(panel.background = element_rect(color = "black", linewidth = 2), legend.position = "inside", legend.position.inside = c(0.89, 0.84), 
        legend.background = element_rect(color = "black", linewidth = 1), legend.text = element_text(size = 20), legend.title = element_blank(), legend.key = element_blank())
  
Pm_group_assignments <- Pm_DAPC$grp |> as.data.frame() |> tibble::rownames_to_column(var = "Sample")

Pm_group_assignments |> saveRDS("DAPC_clusters.rds")

Pm_group_assignments |> dplyr::count(Pm_DAPC$grp)

Pm_group_assignments <- Pm_group_assignments |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                               stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                               stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                               stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                               .default = "Tanzania"))

Pm_clusters_geog <- Pm_group_assignments |> ggplot() + geom_bar(aes(x = Pm_DAPC$grp, fill = Country)) + scale_fill_brewer(palette = "Dark2", name = "Country") + theme_classic() + xlab("Cluster") + ylab("Count") + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("Pm_clusters_geography.png", Pm_clusters_geog, dpi = 600)

#Pm_scatter <- magick::image_read_svg("wsaf_filtered_DAPC_clusters.svg", width = 720) |> magick::image_ggplot(interpolate = TRUE) + theme(text = element_text(size = 16), plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"))

#Pm_pops <- magick::image_read_svg("wsaf_filtered_DAPC_discriminant.svg") |> magick::image_trim() |> magick::image_ggplot() + theme(text = element_text(size = 16), plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"))

Pm_DAPC_plot <- Pm_scatter + plot_spacer() + Pm_clusters_geog + plot_annotation(tag_levels = "A") + plot_layout(widths = c(4, 0.5, 4))

ggsave("Pm_DAPC.png", Pm_DAPC_plot, dpi = 600, width = 20, height = 10, units = "in")
