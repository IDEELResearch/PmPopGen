setwd("/work/users/z/p/zpopkinh/Pm_rerun/Picard_output/")

library("ggcoverage")

#system("conda activate dadi_r")

#system("for i in *dedupped.bam; do bamCoverage -b ${i} -o ${i%.bam}.bamCoverage; done")

bams <- list.files(pattern = "*dedupped.bam$")

bamCoverage <- list.files(pattern = "*.bamCoverage")

library(foreach)

library(iterators)

chrs <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "API", "MIT") |> as.character()

#k <- foreach(i = 1:length(chrs)) %do% print(paste0("PmUG01_", i, "_v1"))

chroms <- paste0("PmUG01_", formatC(1:length(chrs), width = 2, format = "d", flag = "0"), "_v1") #has to be encoded like this to keep the leading zeroes

#system("module load samtools")

#system("for i in *dedupped.bam; do samtools index $i; done")

#trial_bams <- bams[1:5]

#trial_coverage <- bamCoverage[1:5]

#trial_run <- foreach(i = trial_bams, j = trial_coverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_01_v1:1-10,000,000"))

chr1_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_01_v1:1-10,000,000"))

chr1_tracks |> saveRDS("chr1_coverage.rds")
                                                                       
chr2_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_02_v1:1-10,000,000"))

chr2_tracks |> saveRDS("chr2_coverage.rds")

chr3_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_03_v1:1-10,000,000"))

chr3_tracks |> saveRDS("chr3_coverage.rds")

chr4_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_04_v1:1-10,000,000"))

chr4_tracks |> saveRDS("chr4_coverage.rds")

chr5_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_05_v1:1-10,000,000"))

chr5_tracks |> saveRDS("chr5_coverage.rds")

chr6_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_06_v1:1-10,000,000"))

chr6_tracks |> saveRDS("chr6_coverage.rds")

chr7_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_07_v1:1-10,000,000"))

chr7_tracks |> saveRDS("chr7_coverage.rds")

chr8_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_08_v1:1-10,000,000"))

chr8_tracks |> saveRDS("chr8_coverage.rds")

chr9_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_09_v1:1-10,000,000"))

chr9_tracks |> saveRDS("chr9_coverage.rds")

chr10_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_10_v1:1-10,000,000"))

chr10_tracks |> saveRDS("chr10_coverage.rds")

chr11_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_11_v1:1-10,000,000"))

chr11_tracks |> saveRDS("chr11_coverage.rds")

chr12_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_12_v1:1-10,000,000"))

chr12_tracks |> saveRDS("chr12_coverage.rds")

chr13_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_13_v1:1-10,000,000"))

chr13_tracks |> saveRDS("chr13_coverage.rds")

chr14_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_14_v1:1-10,000,000"))

chr14_tracks |> saveRDS("chr14_coverage.rds")

API_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_API_v1:1-10,000,000"))

API_tracks |> saveRDS("API_coverage.rds")

MIT_tracks <- foreach(i = bams, j = bamCoverage) %do% LoadTrackFile(track.file = i, bamcoverage.path = j, norm.method = "None", region = c("PmUG01_MIT_v1:1-10,000,000"))

MIT_tracks |> saveRDS("MIT_coverage.rds")

all_tracks <- foreach(i = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "API", "MIT")) %do% readRDS(paste0(i,"_coverage.rds"))

Pm_mask <- rtracklayer::import.bed("Pm_mask_merged.bed")

chr1_tracks <- do.call(rbind, all_tracks[[1]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

chr2_tracks <- do.call(rbind, all_tracks[[2]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame() 

chr3_tracks <- do.call(rbind, all_tracks[[3]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

chr4_tracks <- do.call(rbind, all_tracks[[4]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

chr5_tracks <- do.call(rbind, all_tracks[[5]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

chr6_tracks <- do.call(rbind, all_tracks[[6]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

chr7_tracks <- do.call(rbind, all_tracks[[7]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

chr8_tracks <- do.call(rbind, all_tracks[[8]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

chr9_tracks <- do.call(rbind, all_tracks[[9]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

chr10_tracks <- do.call(rbind, all_tracks[[10]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

chr11_tracks <- do.call(rbind, all_tracks[[11]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

chr12_tracks <- do.call(rbind, all_tracks[[12]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

chr13_tracks <- do.call(rbind, all_tracks[[13]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

chr14_tracks <- do.call(rbind, all_tracks[[14]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

API_tracks <- do.call(rbind, all_tracks[[15]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

MIT_tracks <- do.call(rbind, all_tracks[[16]]) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

#foreach(i = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "API", "MIT")) %do% paste0(i, "_tracks") saveRDS(paste0(i,"_tracks.rds"))

snp_density <- data.table::fread("SNPdensity/Pm_SNPdensity_wholegenome.snpden")

snp_density <- snp_density |> dplyr::rename(START = BIN_START) |> dplyr::mutate(END = START + 999) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE) |> IRanges::subsetByOverlaps(Pm_mask, invert = TRUE) |> GenomicRanges::as.data.frame()

chr_palette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))

chr_palette2 <- chr_palette(16)

chr1_snps <- snp_density |> subset(seqnames == "PmUG01_01_v1")

chr1_median <- chr1_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

library(ggplot2)

chr1_coverage_plot <- chr1_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[1]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

chr1_snp_plot <- chr1_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[1]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Chromosome 1")

library(patchwork)

chr1_plot <- chr1_snp_plot / chr1_coverage_plot

chr2_snps <- snp_density |> subset(seqnames == "PmUG01_02_v1")

chr2_median <- chr2_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

chr2_coverage_plot <- chr2_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[2]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

chr2_snp_plot <- chr2_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[2]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Chromosome 2")

chr2_plot <- chr2_snp_plot / chr2_coverage_plot

chr3_snps <- snp_density |> subset(seqnames == "PmUG01_03_v1")

chr3_median <- chr3_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

chr3_coverage_plot <- chr3_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[3]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

chr3_snp_plot <- chr3_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[3]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Chromosome 3")

chr3_plot <- chr3_snp_plot / chr3_coverage_plot

chr4_snps <- snp_density |> subset(seqnames == "PmUG01_04_v1")

chr4_median <- chr4_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

chr4_coverage_plot <- chr4_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[4]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

chr4_snp_plot <- chr4_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[4]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Chromosome 4")

chr4_plot <- chr4_snp_plot / chr4_coverage_plot

chr5_snps <- snp_density |> subset(seqnames == "PmUG01_05_v1")

chr5_median <- chr5_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

chr5_coverage_plot <- chr5_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[5]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

chr5_snp_plot <- chr5_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[5]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Chromosome 5")

chr5_plot <- chr5_snp_plot / chr5_coverage_plot

chr6_snps <- snp_density |> subset(seqnames == "PmUG01_06_v1")

chr6_median <- chr6_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

chr6_coverage_plot <- chr6_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[6]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

chr6_snp_plot <- chr6_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[6]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Chromosome 6")

chr6_plot <- chr6_snp_plot / chr6_coverage_plot

chr7_snps <- snp_density |> subset(seqnames == "PmUG01_07_v1")

chr7_median <- chr7_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

chr7_coverage_plot <- chr7_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[7]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

chr7_snp_plot <- chr7_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[7]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Chromosome 7")

chr7_plot <- chr7_snp_plot / chr7_coverage_plot

chr8_snps <- snp_density |> subset(seqnames == "PmUG01_08_v1")

chr8_median <- chr8_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

chr8_coverage_plot <- chr8_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[8]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

chr8_snp_plot <- chr8_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[8]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Chromosome 8")

chr8_plot <- chr8_snp_plot / chr8_coverage_plot

chr9_snps <- snp_density |> subset(seqnames == "PmUG01_09_v1")

chr9_median <- chr9_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

chr9_coverage_plot <- chr9_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[9]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

chr9_snp_plot <- chr9_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[9]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Chromosome 9")

chr9_plot <- chr9_snp_plot / chr9_coverage_plot

chr10_snps <- snp_density |> subset(seqnames == "PmUG01_10_v1")

chr10_median <- chr10_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

chr10_coverage_plot <- chr10_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[10]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

chr10_snp_plot <- chr10_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[10]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Chromosome 10")

chr10_plot <- chr10_snp_plot / chr10_coverage_plot

chr11_snps <- snp_density |> subset(seqnames == "PmUG01_11_v1")

chr11_median <- chr11_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

chr11_coverage_plot <- chr11_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[11]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

chr11_snp_plot <- chr11_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[11]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Chromosome 11")

chr11_plot <- chr11_snp_plot / chr11_coverage_plot

chr12_snps <- snp_density |> subset(seqnames == "PmUG01_12_v1")

chr12_median <- chr12_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

chr12_coverage_plot <- chr12_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[12]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

chr12_snp_plot <- chr12_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[12]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Chromosome 12")

chr12_plot <- chr12_snp_plot / chr12_coverage_plot

chr13_snps <- snp_density |> subset(seqnames == "PmUG01_13_v1")

chr13_median <- chr13_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

chr13_coverage_plot <- chr13_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[13]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

chr13_snp_plot <- chr13_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[13]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Chromosome 13")

chr13_plot <- chr13_snp_plot / chr13_coverage_plot

chr14_snps <- snp_density |> subset(seqnames == "PmUG01_14_v1")

chr14_median <- chr14_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

chr14_coverage_plot <- chr14_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[14]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

chr14_snp_plot <- chr14_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[14]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Chromosome 14")

chr14_plot <- chr14_snp_plot / chr14_coverage_plot

API_snps <- snp_density |> subset(seqnames == "PmUG01_API_v1")

API_median <- API_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

API_coverage_plot <- API_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[15]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

API_snp_plot <- API_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[15]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Apicoplast")

API_plot <- API_snp_plot / API_coverage_plot

MIT_snps <- snp_density |> subset(seqnames == "PmUG01_MIT_v1")

MIT_median <- MIT_tracks |> dplyr::group_by(seqnames, start, end) |> dplyr::summarize(median_coverage = median(score))

MIT_coverage_plot <- MIT_median |> ggplot() + geom_col(aes(x = start, y = median_coverage), color = chr_palette2[16]) + theme_classic() + labs(x = "Position (bp)", y = "Median Coverage") + scale_x_continuous(labels = scales::comma)

MIT_snp_plot <- MIT_snps |> ggplot() + geom_col(aes(x = start, y = SNP_COUNT), color = chr_palette2[16]) + theme_classic() + labs(x = "Position (bp)", y = "SNP Count") + scale_x_continuous(labels = scales::comma) + ggtitle("Mitochondrion")

MIT_plot <- MIT_snp_plot / MIT_coverage_plot

all_plots <- (chr1_plot | chr2_plot | chr3_plot | chr4_plot | chr5_plot | chr6_plot | chr7_plot | chr8_plot) / (chr9_plot | chr10_plot | chr11_plot | chr12_plot | chr13_plot | chr14_plot | API_plot | MIT_plot) & theme()# + plot_layout(ncol = 8, nrow = 4) & theme(text = element_text(size = 30))

ggsave("all_plots.png", dpi = 600, height = 15, width = 30, units = "in")
