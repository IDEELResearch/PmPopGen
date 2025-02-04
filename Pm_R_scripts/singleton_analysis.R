setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun")

system("vk calc sample_hom_gt Pm_PCA.vcf.gz > Pm_PCA_sample_hom_gt.txt")

system("vk calc sample_hom_gt Pm_IBD_pruned.vcf.gz > Pm_IBD_pruned_sample_hom_gt.txt")

system("vk calc sample_hom_gt Pm_HC_missingness_filtered_first.vcf.gz > Pm_HC_missingness_filtered_first_pruned_sample_hom_gt.txt")

initial_genotypes <- data.table::fread("Pm_HC_missingness_filtered_first_pruned_sample_hom_gt.txt")

initial_singletons <- initial_genotypes |> dplyr::filter(freq_of_gt == 1) |> dplyr::select(sample, n_gt_at_freq) |> dplyr::rename(all_singletons = n_gt_at_freq)

initial_non_singletons <- initial_genotypes |> dplyr::filter(freq_of_gt > 1) |> dplyr::group_by(sample) |> dplyr::summarize(all_non_singletons = sum(n_gt_at_freq))

initial_df <- dplyr::left_join(initial_singletons, initial_non_singletons) |> dplyr::group_by(sample) |>
  dplyr::mutate(prop_singletons = (all_singletons/sum(all_singletons + all_non_singletons)), perc_singletons = prop_singletons * 100)

initial_summary <- summary(initial_df$perc_singletons)

pruned_IBD_genotypes <- data.table::fread("Pm_IBD_pruned_sample_hom_gt.txt")

pruned_IBD_singletons <- pruned_IBD_genotypes |> dplyr::filter(freq_of_gt == 1) |> dplyr::select(sample, n_gt_at_freq) |> dplyr::rename(all_singletons = n_gt_at_freq)

pruned_IBD_non_singletons <- pruned_IBD_genotypes |> dplyr::filter(freq_of_gt > 1) |> dplyr::group_by(sample) |> dplyr::summarize(all_non_singletons = sum(n_gt_at_freq))

pruned_IBD_df <- dplyr::left_join(pruned_IBD_singletons, pruned_IBD_non_singletons) |> dplyr::group_by(sample) |>
  dplyr::mutate(prop_singletons = (all_singletons/sum(all_singletons + all_non_singletons)), perc_singletons = prop_singletons * 100)

pruned_IBD_summary <- summary(pruned_IBD_df$perc_singletons)

PCA_genotypes <- data.table::fread("Pm_PCA_sample_hom_gt.txt")

PCA_singletons <- PCA_genotypes |> dplyr::filter(freq_of_gt == 1) |> dplyr::select(sample, n_gt_at_freq) |> dplyr::rename(all_singletons = n_gt_at_freq)

PCA_non_singletons <- PCA_genotypes |> dplyr::filter(freq_of_gt > 1) |> dplyr::group_by(sample) |> dplyr::summarize(all_non_singletons = sum(n_gt_at_freq))

PCA_df <- dplyr::left_join(PCA_singletons, PCA_non_singletons) |> dplyr::group_by(sample) |>
  dplyr::mutate(prop_singletons = (all_singletons/sum(all_singletons + all_non_singletons)), perc_singletons = prop_singletons * 100)

PCA_summary <- summary(PCA_df$perc_singletons)
