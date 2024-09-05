setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

Pm_genes_Tajima <- data.table::fread("Pm_genes_Tajima.txt") 

#Pm_genes_Tajima|> dplyr::arrange(desc(abs(TajimaD))) |> writexl::write_xlsx("Pm_genes_Tajima.xlsx")

nonneutral_Pm_genes <- Pm_genes_Tajima |> subset(abs(TajimaD) > 2)

top_Pm_genes <- Pm_genes_Tajima |> dplyr::arrange(desc(abs(TajimaD))) |> dplyr::slice_head(prop = 0.005)
  
Pm_exons_Tajima <- data.table::fread("Pm_exons_Tajima.txt")

Pm_exons_Tajima|> dplyr::arrange(desc(abs(TajimaD))) |> writexl::write_xlsx("Pm_exons_Tajima.xlsx")

nonneutral_Pm_exons <- Pm_exons_Tajima |> subset(abs(TajimaD) > 2)

top_Pm_exons <- Pm_exons_Tajima |> dplyr::arrange(desc(abs(TajimaD))) |> dplyr::slice_head(prop = 0.005)

masked_orthologs <- readxl::read_xlsx("Pf-Pm_masked_orthologs.xlsx")

#Pm_masked_orthos <- masked_orthologs |> dplyr::select(Group_ID, Pm_CHROM, Pm_START, Pm_END, Pm_LENGTH, Pm_ortho) |> dplyr::rename(CHROM = Pm_CHROM, START = Pm_START, END = Pm_END) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE)

PmGFF_genes <- ape::read.gff("../PlasmoDB-67_PmalariaeUG01.gff") |> subset(!stringr::str_detect(seqid, "_archived_|MIT|API")) |> subset(stringr::str_detect(type, "_gene")) |> tidyr::separate_wider_delim(attributes, delim = ";", names = c("ID", "Name", "description", "biotype"), too_few = "align_start", too_many = "debug") |> 
  dplyr::mutate(ID = stringr::str_remove(ID, "ID="), biotype = dplyr::case_when(stringr::str_detect(biotype, "ebi_biotype=") ~ stringr::str_remove(biotype, "ebi_biotype="),
                                                                                stringr::str_detect(description, "ebi_biotype=") ~ stringr::str_remove(description, "ebi_biotype=")), description = dplyr::case_when(stringr::str_detect(description, "description=") ~ stringr::str_remove(description, "description="),
                                                                                    stringr::str_detect(Name, "description=") ~ stringr::str_remove(Name, "description=")), Name = dplyr::case_when(stringr::str_detect(Name, "Name=") ~ stringr::str_remove(Name, "Name="),
                                                                             .default = "")) |> dplyr::select(seqid, source, type, start, end, strand, ID, Name, description, biotype) |> GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

top_genes <- top_Pm_genes |> dplyr::rename(START = BIN_START, END = BIN_END) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE)

top_exons <- top_Pm_exons |> dplyr::rename(START = BIN_START, END = BIN_END) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE)

gene_intervals <- plyranges::join_overlap_intersect(top_genes, PmGFF_genes) |> GenomicRanges::as.data.frame()

exon_intervals <- plyranges::join_overlap_intersect(top_exons, PmGFF_genes) |> GenomicRanges::as.data.frame()

#unique_genes <- unique(gene_intervals$ID)

unique_genes_metadata <- gene_intervals |> dplyr::group_by(ID) |> dplyr::slice_max(abs(TajimaD)) |> dplyr::distinct(ID, .keep_all = TRUE)

unique_exons_metadata <- exon_intervals |> dplyr::group_by(ID) |> dplyr::slice_max(abs(TajimaD)) |> dplyr::distinct(ID, .keep_all = TRUE)

unique_genes_metadata <- unique_genes_metadata |> dplyr::mutate(Overlap = dplyr::case_when(ID %in% unique_exons_metadata$ID ~ "Yes",
                                                                                           .default = "No"))
unique_exons_metadata <- unique_exons_metadata |> dplyr::mutate(Overlap = dplyr::case_when(ID %in% unique_genes_metadata$ID ~ "Yes",
                                                                                           .default = "No"))

overlapping_genes <- unique_genes_metadata |> subset(Overlap == "Yes") |> dplyr::select(seqnames, TajimaD, ID, Name, description, biotype) |> dplyr::rename(CHROM = seqnames, gene_Tajima = TajimaD)

overlapping_exons <- unique_exons_metadata |> subset(Overlap == "Yes") |> dplyr::select(seqnames, TajimaD, ID, Name, description, biotype) |> dplyr::rename(CHROM = seqnames, exon_Tajima = TajimaD)

overlaps <- dplyr::left_join(overlapping_genes, overlapping_exons) |> 
  dplyr::mutate(diff_Tajima = dplyr::case_when(gene_Tajima != exon_Tajima ~ "Yes",
                                               .default = "No"))

overlaps |> data.table::fwrite("gene_exon_overlapping_sig_Tajima.txt", sep = "\t")

nonoverlapping_genes <- unique_genes_metadata |> subset(Overlap == "No") |> dplyr::select(seqnames, TajimaD, ID, Name, description, biotype) |> dplyr::rename(CHROM = seqnames, gene_Tajima = TajimaD)

nonoverlapping_genes |> data.table::fwrite("gene_overlapping_sig_Tajima.txt", sep = "\t")

nonoverlapping_exons <- unique_exons_metadata |> subset(Overlap == "No") |> dplyr::select(seqnames, TajimaD, ID, Name, description, biotype) |> dplyr::rename(CHROM = seqnames, exon_Tajima = TajimaD)

nonoverlapping_exons |> data.table::fwrite("exons_overlapping_sig_Tajima.txt", sep = "\t")

nonneutral_genes <- nonneutral_Pm_genes |> dplyr::rename(START = BIN_START, END = BIN_END) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE)

nonneutral_gene_intervals <- plyranges::join_overlap_intersect(nonneutral_genes, PmGFF_genes) |> GenomicRanges::as.data.frame()

nonneutral_genes_metadata <- nonneutral_gene_intervals |> dplyr::group_by(ID) |> dplyr::slice_max(abs(TajimaD)) |> dplyr::distinct(ID, .keep_all = TRUE)

Tajima_genes <- Pm_genes_Tajima |> dplyr::rename(START = BIN_START, END = BIN_END) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE)

Tajima_gene_intervals<- plyranges::join_overlap_intersect(Tajima_genes, PmGFF_genes) |> GenomicRanges::as.data.frame()

Tajima_genes_metadata <- Tajima_gene_intervals |> dplyr::group_by(ID) |> dplyr::slice_max(abs(TajimaD)) |> dplyr::distinct(ID, .keep_all = TRUE)

Tajima_genes_metadata |> writexl::write_xlsx("Pm_genes_Tajima.xlsx")

#unique_exons <- unique(exon_intervals$ID)

#gene_exon_overlap <- intersect(unique_genes, unique_exons) |> as.data.frame()

#gene_exon_overlap |> data.table::fwrite("gene_exon_overlapping_sig_Tajima.txt", col.names =F, sep = "\n")

#genes_only <- setdiff(unique_genes, unique_exons) |> as.data.frame()

#genes_only |> data.table::fwrite("gene_overlapping_sig_Tajima.txt", col.names =F, sep = "\n")

#exons_only <- setdiff(unique_exons, unique_genes) |> as.data.frame()

#exons_only |> data.table::fwrite("exons_overlapping_sig_Tajima.txt", col.names =F, sep = "\n")

Pm_chrlen <- data.table::fread("../Pm_chrlen.txt") |> dplyr::select(V1, V3) |> dplyr::mutate(V1 = stringr::str_remove(V1, "PmUG01_")) |> dplyr::mutate(V1 = stringr::str_remove(V1, "_v1")) |> dplyr::mutate(V1 = paste0("chr", V1)) |> dplyr::rename(CHROM = V1, LENGTH = V3)

Pm_chrlen <- Pm_chrlen |> dplyr::mutate(cumsum = ave(LENGTH, FUN=cumsum), toadd = head(c(0, cumsum), -1)) 

setwd("nSL/")

nsl_filenames <- list.files(pattern = "\\d+.res.nsl.out.100bins")

nsl <- lapply(nsl_filenames, read.table)

names(nsl) <- nsl_filenames |> stringr::str_remove(".res.nsl.out.100bins.norm")

nsl <- purrr::imap(nsl, ~dplyr::mutate(.x, V9 = paste0("chr", .y)))

nsl_big <- dplyr::bind_rows(nsl) |> dplyr::rename(locusID = V1, POS = V2, Freq1 = V3, sl1 = V4, sl0 = V5, raw_nSL = V6, normalized_nSL = V7, CHROM = V9) |> dplyr::select(locusID, CHROM, POS, normalized_nSL)

nsl_big <- nsl_big |> dplyr::left_join(Pm_chrlen) |> dplyr::mutate(genome_pos = POS + toadd)
  
top_nsl <- nsl_big |> dplyr::arrange(desc(abs(normalized_nSL))) |> dplyr::slice_head(prop = 0.005)

top_nsl_grange <- top_nsl |> dplyr::mutate(seqnames = paste0("PmUG01_", stringr::str_remove(CHROM, "chr"), "_v1")) |> dplyr::rename(start = POS) |> dplyr::mutate(end = start + 1) |> dplyr::select(seqnames, start, end, normalized_nSL, genome_pos) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE)

nsl_intervals <- plyranges::join_overlap_intersect(top_nsl_grange, PmGFF_genes) |> GenomicRanges::as.data.frame()

unique_nsl_metadata <- nsl_intervals |> dplyr::group_by(ID) |> dplyr::slice_max(abs(normalized_nSL)) |> dplyr::distinct(ID, .keep_all = TRUE)

unique_nsl_metadata |> data.table::fwrite("top_nsl.txt", sep = "\t")

Pm_chrlen <- Pm_chrlen |> dplyr::mutate(axis_breaks = (LENGTH/2) + toadd)

getPalette <- brewer.pal(6, "Paired")

nsl_MAF5_plot <- nsl_big |> ggplot() + geom_point(aes(x = genome_pos, y = abs(normalized_nSL), color = CHROM)) + discrete_scale("color", "custom", function(n){getPalette[1:2][1:n%%2+1]}) + geom_point(data = top_nsl, aes(x = genome_pos, y = abs(normalized_nSL)), color = getPalette[6]) + theme_linedraw() + scale_x_continuous(name = "Chromosome", breaks = Pm_chrlen$axis_breaks, labels = c(1:14), minor_breaks = NULL) + ggrepel::geom_label_repel(data = unique_nsl_metadata, min.segment.length = 0, aes(x = genome_pos, y = abs(normalized_nSL), label = Name, segment.color = getPalette[6], size = 16), nudge_x = 1000, nudge_y = 0.05, color = getPalette[6]) + theme(panel.grid = element_blank(), legend.position = "none", axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab(expression(italic("nS")[L]))

ggsave("../nSL_MAF5.png", nsl_MAF5_plot, width = 15, height = 5, units = "in", dpi = 600)

nsl_filenames_nomaf <- list.files(pattern = "*no_maf.res.nsl.out.100bins*")

nsl_nomaf <- lapply(nsl_filenames_nomaf, read.table)

names(nsl_nomaf) <- nsl_filenames |> stringr::str_remove("_no_maf.res.nsl.out.100bins.norm")

nsl_nomaf <- purrr::imap(nsl_nomaf, ~dplyr::mutate(.x, V9 = paste0("chr", .y)))

nsl_big_nomaf <- dplyr::bind_rows(nsl_nomaf) |> dplyr::rename(locusID = V1, POS = V2, Freq1 = V3, sl1 = V4, sl0 = V5, raw_nSL = V6, normalized_nSL = V7, CHROM = V9) |> dplyr::select(locusID, CHROM, POS, normalized_nSL) |> dplyr::mutate(CHROM = stringr::str_remove(CHROM, ".res.nsl.out.100bins.norm"))

nsl_big_nomaf <- nsl_big_nomaf |> dplyr::left_join(Pm_chrlen) |> dplyr::mutate(genome_pos = POS + toadd)

top_nsl_nomaf <- nsl_big_nomaf |> dplyr::arrange(desc(abs(normalized_nSL))) |> dplyr::slice_head(prop = 0.005)

top_nsl_nomaf_grange <- top_nsl_nomaf |> dplyr::mutate(seqnames = paste0("PmUG01_", stringr::str_remove(CHROM, "chr"), "_v1")) |> dplyr::rename(start = POS) |> dplyr::mutate(end = start + 1) |> dplyr::select(seqnames, start, end, normalized_nSL, genome_pos) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE)

nsl_intervals_nomaf <- plyranges::join_overlap_intersect(top_nsl_nomaf_grange, PmGFF_genes) |> GenomicRanges::as.data.frame()

unique_nsl_metadata_nomaf <- nsl_intervals_nomaf |> dplyr::group_by(ID) |> dplyr::slice_max(abs(normalized_nSL)) |> dplyr::distinct(ID, .keep_all = TRUE)

unique_nsl_metadata_nomaf |> data.table::fwrite("top_nsl_nomaf.txt", sep = "\t")

nsl_noMAF_plot <- nsl_big_nomaf |> ggplot() + geom_point(aes(x = genome_pos, y = abs(normalized_nSL), color = CHROM)) + discrete_scale("color", "custom", function(n){getPalette[1:2][1:n%%2+1]}) + geom_point(data = top_nsl_nomaf, aes(x = genome_pos, y = abs(normalized_nSL)), color = getPalette[6]) + theme_linedraw() + scale_x_continuous(name = "Chromosome", breaks = Pm_chrlen$axis_breaks, labels = c(1:14), minor_breaks = NULL) + ggrepel::geom_label_repel(data = unique_nsl_metadata_nomaf, min.segment.length = 0, aes(x = genome_pos, y = abs(normalized_nSL), label = Name, segment.color = getPalette[6], size = 16), nudge_x = 1000, nudge_y = 0.05, color = getPalette[6], max.overlaps = 20) + theme(panel.grid = element_blank(), legend.position = "none", axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab(expression(italic("nS")[L]))

ggsave("../nSL_no_MAF.png", nsl_noMAF_plot, width = 15, height = 5, units = "in", dpi = 600)

library(patchwork)

nSL_plots <- nsl_MAF5_plot + nsl_noMAF_plot + plot_annotation(tag_levels = "A") + plot_layout(ncol = 1, nrow = 2)

ggsave("../nSL_plots.png", nSL_plots, width = 15, height = 10, units = "in", dpi = 600)

big_Tajima_gene <- Pm_genes_Tajima |> dplyr::mutate(CHROM = paste0("chr", stringr::str_replace(CHROM, "PmUG01_(\\d+)_v1", "\\1"))) |> 
  dplyr::mutate(POS = BIN_END - 150) |> dplyr::left_join(Pm_chrlen) |> dplyr::mutate(genome_pos = POS + toadd)

top_gene_highlight <- gene_intervals |> dplyr::mutate(CHROM = paste0("chr", stringr::str_replace(seqnames, "PmUG01_(\\d+)_v1", "\\1"))) |> 
  dplyr::mutate(POS = end - 150) |> dplyr::left_join(Pm_chrlen) |> dplyr::mutate(genome_pos = POS + toadd) |> dplyr::distinct(ID, .keep_all = TRUE)

Tajima_gene_plot <- big_Tajima_gene |> ggplot() + geom_point(aes(x = genome_pos, y = TajimaD, color = CHROM)) + discrete_scale("color", "custom", function(n){getPalette[1:2][1:n%%2+1]}) + geom_point(data = top_gene_highlight, aes(x = genome_pos, y = TajimaD), color = getPalette[6]) + theme_linedraw() + scale_x_continuous(name = "Chromosome", breaks = Pm_chrlen$axis_breaks, labels = c(1:14), minor_breaks = NULL) + ggrepel::geom_label_repel(data = top_gene_highlight, min.segment.length = 0, aes(x = genome_pos, y = TajimaD, label = Name, segment.color = getPalette[6], size = 16), nudge_x = 1000, nudge_y = 0.05, color = getPalette[6], max.overlaps = 20) + theme(panel.grid = element_blank(), legend.position = "none", axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Tajima's D")

ggsave("../Tajima_genes.png", Tajima_gene_plot, width = 15, height = 10, units = "in", dpi = 600)

big_Tajima_exon <- Pm_exons_Tajima |> dplyr::mutate(CHROM = paste0("chr", stringr::str_replace(CHROM, "PmUG01_(\\d+)_v1", "\\1"))) |> 
  dplyr::mutate(POS = BIN_END - 150) |> dplyr::left_join(Pm_chrlen) |> dplyr::mutate(genome_pos = POS + toadd)

top_exon_highlight <- exon_intervals |> dplyr::mutate(CHROM = paste0("chr", stringr::str_replace(seqnames, "PmUG01_(\\d+)_v1", "\\1"))) |> 
  dplyr::mutate(POS = end - 150) |> dplyr::left_join(Pm_chrlen) |> dplyr::mutate(genome_pos = POS + toadd) |> dplyr::distinct(ID, .keep_all = TRUE)

Tajima_exon_plot <- big_Tajima_exon |> ggplot() + geom_point(aes(x = genome_pos, y = TajimaD, color = CHROM)) + discrete_scale("color", "custom", function(n){getPalette[1:2][1:n%%2+1]}) + geom_point(data = top_exon_highlight, aes(x = genome_pos, y = TajimaD), color = getPalette[6]) + theme_linedraw() + scale_x_continuous(name = "Chromosome", breaks = Pm_chrlen$axis_breaks, labels = c(1:14), minor_breaks = NULL) + ggrepel::geom_label_repel(data = top_exon_highlight, min.segment.length = 0, aes(x = genome_pos, y = TajimaD, label = Name, segment.color = getPalette[6], size = 16), nudge_x = 1000, nudge_y = 0.05, color = getPalette[6], max.overlaps = 20) + theme(panel.grid = element_blank(), legend.position = "none", axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Tajima's D")

ggsave("../Tajima_exons.png", Tajima_exon_plot, width = 15, height = 10, units = "in", dpi = 600)

Tajima_plots <- Tajima_gene_plot + Tajima_exon_plot + plot_annotation(tag_levels = "A") + plot_layout(ncol = 1, nrow = 2)

ggsave("../Tajima_plot.png", Tajima_plots, width = 15, height = 20, units = "in", dpi = 600)

gene_list <- c("MSP1", "AMA1", "CSP", "TRAP", "P25", "P48/45", "MDR1", "MDR2", "CRT", "Kelch13", "PPPK-DHPS", "DHFR-TS", "MRP1", "MRP2", "LSA1")

ID_list <- c("PmUG01_10046900")

targeted_genes <- PmGFF_genes |> subset(Name %in% gene_list | ID %in% ID_list)|> GenomicRanges::as.data.frame() |> dplyr::rename(name = Name) |> GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

targeted_genes |> rtracklayer::export.bed("target_genes.bed")

target_Tajima <- data.table::fread("target_Tajima.txt") 


PfGFF_genes <- ape::read.gff("../PlasmoDB-67_Pfalciparum3D7.gff") |> subset(!stringr::str_detect(seqid, "_archived_|MIT|API")) |> subset(stringr::str_detect(type, "_gene")) |> tidyr::separate_wider_delim(attributes, delim = ";", names = c("ID", "Name", "description", "biotype"), too_few = "align_start", too_many = "debug") |> 
  dplyr::mutate(ID = stringr::str_remove(ID, "ID="), biotype = dplyr::case_when(stringr::str_detect(biotype, "ebi_biotype=") ~ stringr::str_remove(biotype, "ebi_biotype="),
                                                                                stringr::str_detect(description, "ebi_biotype=") ~ stringr::str_remove(description, "ebi_biotype=")), description = dplyr::case_when(stringr::str_detect(description, "description=") ~ stringr::str_remove(description, "description="),
                                                                                                                                                                                                                     stringr::str_detect(Name, "description=") ~ stringr::str_remove(Name, "description=")), Name = dplyr::case_when(stringr::str_detect(Name, "Name=") ~ stringr::str_remove(Name, "Name="),
                                                                                                                                                                                                                                                                                                                                     .default = "")) |> dplyr::select(seqid, source, type, start, end, strand, ID, Name, description, biotype) |> GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

Pf_targets <- PfGFF_genes |> subset(Name %in% gene_list | ID %in% ID_list)|> GenomicRanges::as.data.frame() |> dplyr::rename(name = Name) |> GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

Pf_targets |> rtracklayer::export.bed("Pf_target_genes.bed")

targets <- target_Tajima |> dplyr::rename(START = BIN_START, END = BIN_END) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE)

target_intervals <- plyranges::join_overlap_intersect(targets, PmGFF_genes) |> GenomicRanges::as.data.frame()

target_metadata <- target_intervals |> dplyr::group_by(ID) |> dplyr::slice_max(abs(TajimaD)) |> dplyr::distinct(ID, .keep_all = TRUE)

big_target <- target_Tajima |> dplyr::mutate(CHROM = paste0("chr", stringr::str_replace(CHROM, "PmUG01_(\\d+)_v1", "\\1"))) |> 
  dplyr::mutate(POS = BIN_END - 150) |> dplyr::left_join(Pm_chrlen) |> dplyr::mutate(genome_pos = POS + toadd)

target_highlight <- target_intervals |> dplyr::mutate(CHROM = paste0("chr", stringr::str_replace(seqnames, "PmUG01_(\\d+)_v1", "\\1"))) |> 
  dplyr::mutate(POS = end - 150) |> dplyr::left_join(Pm_chrlen) |> dplyr::mutate(genome_pos = POS + toadd) |> dplyr::group_by(ID) |> dplyr::slice_max(abs(TajimaD)) |> dplyr::distinct(ID, .keep_all = TRUE) |> dplyr::mutate(Name = dplyr::case_when(ID == "PmUG01_10046900" ~ "LSA1",
                                                                                                                                                                                                                                                      .default = Name))

target_Tajima_plot <- big_target |> ggplot() + geom_point(aes(x = genome_pos, y = TajimaD, color = CHROM)) + discrete_scale("color", "custom", function(n){getPalette[1:2][1:n%%2+1]}) + geom_point(data = target_highlight, aes(x = genome_pos, y = TajimaD), color = getPalette[6]) + theme_linedraw() + scale_x_continuous(name = "Chromosome", breaks = Pm_chrlen$axis_breaks, labels = c(1:14), minor_breaks = NULL) + ggrepel::geom_label_repel(data = target_highlight, min.segment.length = 0, aes(x = genome_pos, y = TajimaD, label = Name, segment.color = getPalette[6], size = 16), nudge_x = 1000, nudge_y = 0.05, color = getPalette[6], max.overlaps = 20) + theme(panel.grid = element_blank(), legend.position = "none", axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab("Tajima's D")

ggsave("Tajima_targets.png", target_Tajima_plot, width = 15, height = 10, units = "in", dpi = 600)

nsl_MAF5_grange <- nsl_big |> dplyr::mutate(seqnames = paste0("PmUG01_", stringr::str_remove(CHROM, "chr"), "_v1")) |> dplyr::rename(start = POS) |> dplyr::mutate(end = start + 1) |> dplyr::select(seqnames, start, end, normalized_nSL, genome_pos) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE)

nsl_MAF5_intervals <- plyranges::join_overlap_intersect(nsl_MAF5_grange, targeted_genes) |> GenomicRanges::as.data.frame()

nsl_MAF5_target_highlight <- nsl_MAF5_intervals |> dplyr::group_by(ID) |> dplyr::slice_max(abs(normalized_nSL)) |> dplyr::distinct(ID, .keep_all = TRUE)

nsl_MAF5_target_plot <- nsl_big |> ggplot() + geom_point(aes(x = genome_pos, y = abs(normalized_nSL), color = CHROM)) + discrete_scale("color", "custom", function(n){getPalette[1:2][1:n%%2+1]}) + geom_point(data =nsl_MAF5_intervals, aes(x = genome_pos, y = abs(normalized_nSL)), color = getPalette[6]) + theme_linedraw() + scale_x_continuous(name = "Chromosome", breaks = Pm_chrlen$axis_breaks, labels = c(1:14), minor_breaks = NULL) + ggrepel::geom_label_repel(data = nsl_MAF5_target_highlight, min.segment.length = 0, aes(x = genome_pos, y = abs(normalized_nSL), label = Name, segment.color = getPalette[6], size = 16), nudge_x = 1000, nudge_y = 0.05, color = getPalette[6], max.overlaps = 20) + theme(panel.grid = element_blank(), legend.position = "none", axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab(expression(italic("nS")[L]))

ggsave("nSL_targets_MAF5.png", nsl_MAF5_target_plot, width = 15, height = 5, units = "in", dpi = 600)

nsl_noMAF_grange <- nsl_big_nomaf |> dplyr::mutate(seqnames = paste0("PmUG01_", stringr::str_remove(CHROM, "chr"), "_v1")) |> dplyr::rename(start = POS) |> dplyr::mutate(end = start + 1) |> dplyr::select(seqnames, start, end, normalized_nSL, genome_pos) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE)

nsl_noMAF_intervals <- plyranges::join_overlap_intersect(nsl_noMAF_grange, targeted_genes) |> GenomicRanges::as.data.frame()

nsl_noMAF_target_highlight <- nsl_noMAF_intervals |> dplyr::group_by(ID) |> dplyr::slice_max(abs(normalized_nSL)) |> dplyr::distinct(ID, .keep_all = TRUE)

nsl_noMAF_target_plot <- nsl_big |> ggplot() + geom_point(aes(x = genome_pos, y = abs(normalized_nSL), color = CHROM)) + discrete_scale("color", "custom", function(n){getPalette[1:2][1:n%%2+1]}) + geom_point(data =nsl_noMAF_intervals, aes(x = genome_pos, y = abs(normalized_nSL)), color = getPalette[6]) + theme_linedraw() + scale_x_continuous(name = "Chromosome", breaks = Pm_chrlen$axis_breaks, labels = c(1:14), minor_breaks = NULL) + ggrepel::geom_label_repel(data = nsl_noMAF_target_highlight, min.segment.length = 0, aes(x = genome_pos, y = abs(normalized_nSL), label = Name, segment.color = getPalette[6], size = 16), nudge_x = 1000, nudge_y = 0.05, color = getPalette[6], max.overlaps = 20) + theme(panel.grid = element_blank(), legend.position = "none", axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.ticks.length = unit(0.5,"cm")) + ylab(expression(italic("nS")[L]))

ggsave("nSL_targets_noMAF.png", nsl_noMAF_target_plot, width = 15, height = 5, units = "in", dpi = 600)

PmGFF_CDS <- ape::read.gff("../PlasmoDB-67_PmalariaeUG01.gff") |> subset(!stringr::str_detect(seqid, "_archived_|MIT|API")) |> subset(stringr::str_detect(type, "CDS")) |> tidyr::separate_wider_delim(attributes, delim = ";", names = c("CDS_ID", "Parent", "gene_id", "protein_source"), too_few = "align_start", too_many = "debug") |> 
  dplyr::mutate(CDS_ID = stringr::str_remove(CDS_ID, "ID="), gene_id = stringr::str_remove(gene_id, "gene_id=")) |> dplyr::select(seqid, source, type, start, end, strand, CDS_ID, gene_id) |> GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

target_IDs <- unique(targeted_genes$ID)

target_names <- unique(targeted_genes$name)

target_info <- cbind(target_IDs, target_names) |> as.data.frame()

by <- dplyr::join_by(gene_id == target_IDs)

targeted_CDS <- PmGFF_CDS |> subset(gene_id %in% target_IDs)|> GenomicRanges::as.data.frame() |> dplyr::left_join(target_info, by) |> dplyr::rename(name = target_names) |> GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

targeted_CDS |> rtracklayer::export.bed("target_CDS.bed")

PfGFF_CDS <- ape::read.gff("../PlasmoDB-67_Pfalciparum3D7.gff") |> subset(!stringr::str_detect(seqid, "_archived_|MIT|API")) |> subset(stringr::str_detect(type, "CDS")) |> tidyr::separate_wider_delim(attributes, delim = ";", names = c("CDS_ID", "Parent", "gene_id", "protein_source"), too_few = "align_start", too_many = "debug") |> 
  dplyr::mutate(CDS_ID = stringr::str_remove(CDS_ID, "ID="), gene_id = stringr::str_remove(gene_id, "gene_id=")) |> dplyr::select(seqid, source, type, start, end, strand, CDS_ID, gene_id) |> GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

Pf_target_IDs <- unique(Pf_targets$ID)

Pf_target_names <- unique(Pf_targets$name)

Pf_target_info <- cbind(Pf_target_IDs, Pf_target_names) |> as.data.frame() |> dplyr::rename(gene_id = Pf_target_IDs, name = Pf_target_names)

Pf_targeted_CDS <- PfGFF_CDS |> subset(gene_id %in% Pf_target_IDs)|> GenomicRanges::as.data.frame() |> dplyr::left_join(Pf_target_info) |> GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

Pf_targeted_CDS |> rtracklayer::export.bed("Pf_target_CDS.bed")
