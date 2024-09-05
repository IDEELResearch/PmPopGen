setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun/")

all_orthologs <- readxl::read_xlsx("../Pfal-Pmal_1-1_orthologs.xlsx", sheet = "Final set")

filtered_orthologs <- all_orthologs |> subset(!stringr::str_detect(`Pmal chrom`, "archived")) |> subset(!stringr::str_detect(`Pmal chrom`, "API")) |> subset(!stringr::str_detect(`Pmal chrom`, "MIT"))

Pm_mask <- rtracklayer::import.bed("Pm_mask_merged.bed")

Pm_coverage <- rtracklayer::import.bed("Pm_60percentcov5_merged.bed")

Pf_coverage <- rtracklayer::import.bed("Pf_60percentcov5_merged.bed")

Pf_core <- rtracklayer::import.bed("pf3d7_core.bed")

#filtered_orthologs <- filtered_orthologs |> dplyr::mutate(Pm_mask_overlap = dplyr::case_when(`Pmal chrom` == Pm_mask$CHROM && `Pmal start` <= Pm_mask$STOP && Pm_mask$START <= `Pmal stop` ~ "YES"))

#rtracklayer::export.bed(all_orthologs, "orthologs.bed")

Pm_orthologs <- filtered_orthologs |> dplyr::select(`Pmal chrom`, `Pmal start`, `Pmal stop`, "Group ID", `Pmal ortholog`) |> dplyr::rename(CHROM = `Pmal chrom`, START = `Pmal start`, STOP = `Pmal stop`) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE)

Pf_orthologs <- filtered_orthologs |> dplyr::select(`Pfal chrom`, `Pfal start`, `Pfal stop`, "Group ID", `Pfal ortholog`) |> dplyr::rename(CHROM = `Pfal chrom`, START = `Pfal start`, STOP = `Pfal stop`) |> GenomicRanges::makeGRangesFromDataFrame(na.rm = TRUE, keep.extra.columns = TRUE)

Pm_masked <- IRanges::subsetByOverlaps(Pm_orthologs, Pm_mask, invert = TRUE) |> IRanges::subsetByOverlaps(Pm_coverage)

Pf_masked <- IRanges::subsetByOverlaps(Pf_orthologs, Pf_core) |> IRanges::subsetByOverlaps(Pf_coverage)

Pm_masked <- Pm_masked |> GenomicRanges::as.data.frame()

Pf_masked <- Pf_masked |> GenomicRanges::as.data.frame()

combined_masked <- dplyr::left_join(Pm_masked, Pf_masked, by = "Group.ID") #|> writexl::write_xlsx("Pf-Pm_masked_orthologs.xlsx")

Pm_masked_bed <- combined_masked |> dplyr::select(seqnames.x, start.x, end.x) |> dplyr::rename(CHROM = seqnames.x, START = start.x, STOP = end.x) |> GenomicRanges::makeGRangesFromDataFrame() |> rtracklayer::export.bed("Pm_masked_orthologs.bed")

Pf_masked_bed <- combined_masked |> dplyr::select(seqnames.y, start.y, end.y) |> dplyr::rename(CHROM = seqnames.y, START = start.y, STOP = end.y) |> GenomicRanges::makeGRangesFromDataFrame() |> rtracklayer::export.bed("Pf_masked_orthologs.bed")
