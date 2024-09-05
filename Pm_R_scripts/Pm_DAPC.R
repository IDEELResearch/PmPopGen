setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun")

Pm_SNPs <- vcfR::read.vcfR("Pm_PCA_chrs_only.vcf.gz") |> vcfR::vcfR2genlight()

Pm_clust <- adegenet::find.clusters(Pm_SNPs)
#I selected 50 PCs to retain and 3 clusters - the BIC value indicates that 3 clusters is the best choice

#Pm_clust2 <- adegenet::find.clusters(Pm_SNPs, n.pca = 40, n.clust = 3)

Pm_DAPC <- adegenet::dapc(Pm_SNPs, Pm_clust$grp)
#this uses the groups detected in the previous step - I retained 40 PCs and all discriminant functions (there are only 2)

Pm_scatter <- adegenet::scatter.dapc(Pm_DAPC)

myCol <- RColorBrewer::brewer.pal(3, "Set2")

Pm_scatter <- adegenet::scatter.dapc(Pm_DAPC, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3))

Pm_pops <- adegenet::scatter.dapc(Pm_DAPC,1,1, col=myCol, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)

group_assignments <- Pm_DAPC$grp |> as.data.frame() |> tibble::rownames_to_column(var = "Sample")

group_assignments <- group_assignments |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                                      stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                                      stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                                      stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                                      .default = "Tanzania"))

Pm_clusters_geog <- group_assignments |> ggplot() + geom_bar(aes(x = Pm_DAPC$grp, fill = Country)) + scale_fill_brewer(palette = "Dark2", name = "Country") + theme_classic() + xlab("Cluster") + ylab("Count") + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("Pm_clusters_geography.png", Pm_clusters_geog, dpi = 600)
