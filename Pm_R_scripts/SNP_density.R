setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/")

#bears::gff2bed(gffpath = "PlasmoDB-67_PmalariaeUG01.gff")

#bears::gff2bed(gffpath = "PlasmoDB-67_Pfalciparum3D7.gff")

Pm_mask <- rtracklayer::import.bed("rerun/Pm_mask_merged.bed")

Pf_core <- rtracklayer::import.bed("rerun/pf3d7_core.bed")

PmGFF <- ape::read.gff("PlasmoDB-67_PmalariaeUG01.gff") |> subset(!stringr::str_detect(seqid, "_archived_|MIT|API"))

PfGFF <- ape::read.gff("PlasmoDB-67_Pfalciparum3D7.gff") |> subset(!stringr::str_detect(seqid, "_archived_|MIT|API"))

PmGFF_genes <- PmGFF |> subset(stringr::str_detect(type, "_gene")) |> GenomicRanges::makeGRangesFromDataFrame() |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> rtracklayer::export.bed("Pm_genes.bed", ignore.strand = TRUE)

#PmGFF_genes <- data.table::fread("Pm_genes.bed")

PfGFF_genes <- PfGFF |> subset(stringr::str_detect(type, "_gene")) |> GenomicRanges::makeGRangesFromDataFrame() |> IRanges::subsetByOverlaps(Pf_core) |> rtracklayer::export.bed("Pf_genes.bed", ignore.strand = TRUE)

PmGFF_exons <- PmGFF |> subset(type == "exon") |> GenomicRanges::makeGRangesFromDataFrame() |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> rtracklayer::export.bed("Pm_exons.bed", ignore.strand = TRUE)

PfGFF_exons <- PfGFF |> subset(type == "exon") |> GenomicRanges::makeGRangesFromDataFrame() |> IRanges::subsetByOverlaps(Pf_core) |> rtracklayer::export.bed("Pf_exons.bed", ignore.strand = TRUE)

PmGFF_CDS <- PmGFF |> subset(type == "CDS") |> GenomicRanges::makeGRangesFromDataFrame() |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> rtracklayer::export.bed("Pm_CDS.bed", ignore.strand = TRUE)

PfGFF_CDS <- PfGFF |> subset(type == "CDS") |> GenomicRanges::makeGRangesFromDataFrame() |> IRanges::subsetByOverlaps(Pf_core) |> rtracklayer::export.bed("Pf_CDS.bed", ignore.strand = TRUE)

PmBED <- rtracklayer::import.bed("PlasmoDB-67_PmalariaeUG01.bed")

PfBED <- rtracklayer::import.bed("PlasmoDB-67_Pfalciparum3D7.bed")

Pm_chroms <- unique(PmGFF$seqid) |> as.vector()

Pm_introns <- GenomicFeatures::makeTxDbFromGFF("PlasmoDB-67_PmalariaeUG01.gff") |> GenomicFeatures::intronsByTranscript() |> unlist() |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> plyranges::filter(seqnames %in% Pm_chroms) |> rtracklayer::export.bed("Pm_introns.bed", ignore.strand = TRUE)

Pf_introns <- GenomicFeatures::makeTxDbFromGFF("PlasmoDB-67_Pfalciparum3D7.gff") |> GenomicFeatures::intronsByTranscript() |> unlist() |> IRanges::subsetByOverlaps(Pf_core) |> rtracklayer::export.bed("Pf_introns.bed")

#make intergenic bed file in bash with bedtools complement for the complete bed file

Pm_intergenic <- rtracklayer::import.bed("rerun/Pm_intergenic.bed")

Pf_intergenic <- rtracklayer::import.bed("rerun/Pf_intergenic.bed")

#vcftools --SNPdensity 1000 --gzvcf "/work/users/z/p/zpopkinh/Pm_rerun/Pf_VCFs/ortholog_samples/Pf_ortholog_samples.vcf.gz" --out Pf_SNPdensity_wholegenome

setwd("rerun/SNPdensity/")

####redoing this by subsetting the whole genome files because the math is not mathing

Pm_whole <- data.table::fread("Pm_SNPdensity_wholegenome.snpden") |> subset(!stringr::str_detect(CHROM, "_archived_|MIT|API")) |> dplyr::rename(START = BIN_START) |> dplyr::mutate(END = START + 999) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

Pf_whole <- data.table::fread("Pf_SNPdensity_wholegenome.snpden") |> subset(!stringr::str_detect(CHROM, "_archived_|MIT|API")) |> dplyr::rename(START = BIN_START) |> dplyr::mutate(END = START + 999) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pf_core) |> GenomicRanges::as.data.frame()

Pm_CDS <- Pm_whole |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(PmGFF_CDS) |> GenomicRanges::as.data.frame()

Pf_CDS <- Pf_whole |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(PfGFF_CDS) |> GenomicRanges::as.data.frame()

Pm_genes <- Pm_whole |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(PmGFF_genes) |> GenomicRanges::as.data.frame()

Pf_genes <- Pf_whole |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(PfGFF_genes) |> GenomicRanges::as.data.frame()

Pm_exons <- Pm_whole |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(PmGFF_exons)|> GenomicRanges::as.data.frame()
 
Pf_exons <- Pf_whole |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(PfGFF_exons) |> GenomicRanges::as.data.frame()

Pm_introns <- Pm_whole |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_introns) |> GenomicRanges::as.data.frame()

Pf_introns <- Pf_whole |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pf_introns) |> GenomicRanges::as.data.frame()

Pm_intergenic <- Pm_whole |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_intergenic) |> GenomicRanges::as.data.frame()

Pf_intergenic <- Pf_whole |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pf_intergenic) |> GenomicRanges::as.data.frame()

SNP_densities <- list(Pm_whole, Pf_whole, Pm_CDS, Pf_CDS, Pm_genes, Pf_genes, Pm_exons, Pf_exons, Pm_introns, Pf_introns, Pm_intergenic, Pf_intergenic)

library(foreach)

#SNP_densities <- foreach(i = SNP_densities) %do% {GenomicRanges::as.data.frame(i)}
  
names(SNP_densities) <- c("Pm_whole", "Pf_whole", "Pm_CDS", "Pf_CDS", "Pm_genes", "Pf_genes", "Pm_exons", "Pf_exons", "Pm_introns", "Pf_introns", "Pm_intergenic", "Pf_intergenic")
  
SNP_densities <- purrr::imap(SNP_densities, ~dplyr::mutate(.x, Species = dplyr::case_when(stringr::str_detect(.y, "Pf") ~ "Pf",
                                                                                     stringr::str_detect(.y, "Pm") ~ "Pm")))

SNP_densities <- purrr::imap(SNP_densities, ~dplyr::mutate(.x, Type = dplyr::case_when(stringr::str_detect(.y, "CDS") ~ "CDS",
                                                                                       stringr::str_detect(.y, "gene") ~ "Gene",
                                                                                       stringr::str_detect(.y, "exon") ~ "Exon",
                                                                                       stringr::str_detect(.y, "intron") ~ "Intron",
                                                                                       stringr::str_detect(.y, "intergenic") ~ "Intergenic",
                                                                                       stringr::str_detect(.y, "whole") ~ "Whole")))

SNP_densities <- purrr::imap(SNP_densities, ~dplyr::select(.x, seqnames, start, end, width, SNP_COUNT, VARIANTS.KB, Species, Type))

list2env(SNP_densities,envir=.GlobalEnv)

all_SNP_densities <- rbind(Pf_whole, Pm_whole, Pf_intergenic, Pm_intergenic, Pm_CDS, Pm_exons, Pm_genes, Pm_introns, Pf_CDS, Pf_exons, Pf_genes, Pf_intergenic, Pf_introns)

library(ggplot2)

library(RColorBrewer)
getPalette <- brewer.pal(6, "Set2")

SNP_density_boxplot <- all_SNP_densities |> ggplot() +
  geom_boxplot(aes(x = Species, y = SNP_COUNT, fill = Type))+ theme_linedraw() +
  scale_fill_manual(values = getPalette, name = "Genomic Region", labels = c ("CDS", "Exons", "Genes", "Intergenic Regions", "Introns", "Whole Genome")) +
  scale_x_discrete(labels = c(expression(italic("P. falciparum"), italic("P. malariae")))) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 20), legend.title = element_text(size = 24), legend.text = element_text(size = 20), panel.grid.major.x = element_blank()) + labs(y = "SNP Density (SNPs/Kb)")

Pm_zoom_boxplot <- all_SNP_densities |> subset(Species == "Pm") |> ggplot() +
  geom_boxplot(aes(x = Species, y = SNP_COUNT, fill = Type))+ theme_linedraw() +
  scale_fill_manual(values = getPalette, name = "Genomic Region", labels = c ("CDS", "Exons", "Genes", "Intergenic Regions", "Introns", "Whole Genome")) + ylim(0, 25) +
  scale_x_discrete(labels = c(expression(italic("P. falciparum"), italic("P. malariae")))) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.x = element_blank(), axis.text.y = element_text(size = 16), legend.position = "none", panel.grid.major.x = element_blank(), axis.ticks.x = element_blank(), plot.background = element_rect(color = " black", size = 1)) + labs(y = "SNP Density (SNPs/Kb)")

boxplot_with_zoom <- cowplot::ggdraw() +
  cowplot::draw_plot(SNP_density_boxplot, x = 0, y = 0, width = 1, height = 1, scale = 1) +
  cowplot::draw_plot(Pm_zoom_boxplot, x = 0.415, y= 0.4, width = 0.35, height = 0.4)

ggsave("../SNP_density_boxplot.png", boxplot_with_zoom, width = 16, height = 8, units = "in", dpi = 600)

density_summary <- all_SNP_densities |> dplyr::group_by(Species, Type) |> dplyr::summarize(median_density = median(SNP_COUNT), average_density = median(SNP_COUNT), max_density = max(SNP_COUNT), min_density = min(SNP_COUNT))

density_summary |> writexl::write_xlsx("../SNP_density_summary.xlsx")

aggregated_densities <- all_SNP_densities |> dplyr::group_by(Species, Type) |> dplyr::summarize(aggregated_density = sum(SNP_COUNT)/sum(width) * 1000)

SNP_density_barplot <- aggregated_densities |> ggplot(aes(x = Species, y = aggregated_density, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1)) + theme_linedraw() +
  scale_fill_manual(values = getPalette, name = "Genomic Region", labels = c ("CDS", "Exons", "Genes", "Intergenic Regions", "Introns", "Whole Genome")) +
  scale_x_discrete(labels = c(expression(italic("P. falciparum"), italic("P. malariae")))) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 24), axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 20), legend.title = element_text(size = 24), legend.text = element_text(size = 20), panel.grid.major.x = element_blank()) + labs(y = "SNP Density (SNPs/Kb)")

Pm_zoom <- aggregated_densities |> subset(Species == "Pm") |> ggplot(aes(x = Type, y = aggregated_density, fill = Type)) +
  geom_bar(stat = "identity") + theme_linedraw() +
  scale_fill_manual(values = getPalette, name = "Genomic Region", labels = c ("CDS", "Exons", "Genes", "Intergenic Regions", "Introns", "Whole Genome")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 18), axis.text.x = element_blank(), axis.text.y = element_text(size = 16), legend.position = "none", panel.grid.major.x = element_blank(), axis.ticks.x = element_blank(), plot.background = element_rect(color = "black", size = 1)) + labs(y = "SNP Density (SNPs/Kb)")
  
barplot_with_zoom <- cowplot::ggdraw() +
  cowplot::draw_plot(SNP_density_barplot, x = 0, y = 0, width = 1, height = 1, scale = 1) +
  cowplot::draw_plot(Pm_zoom, x = 0.44, y= 0.2, width = 0.35, height = 0.4)

ggsave("../SNP_density_barplot.png", barplot_with_zoom, width = 16, height = 8, units = "in", dpi = 600)
