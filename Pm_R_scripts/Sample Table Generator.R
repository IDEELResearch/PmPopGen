setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun")

Pm_HC <- vcfR::read.vcfR("Pm_HC_missingness_filtered_first.vcf.gz") |> vcfR::vcfR2tidy()

sample_list <- Pm_HC$gt$Indiv |> unique() |> as.data.frame()

names(sample_list) <- "Sample"

sample_list <- sample_list |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(Sample, "Gam_[:digit:]+") ~ "Nigeria",
                                                                       stringr::str_detect(Sample, "^[:digit:]+") ~ "DRC",
                                                                       stringr::str_detect(Sample, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                                       stringr::str_detect(Sample, "Gam") ~ "Cameroon",
                                                                       .default = "Tanzania"))

country_totals <- sample_list |> dplyr::count(Country)

Tz_samples <- sample_list |> subset(Country == "Tanzania")

Tz_samples <- Tz_samples |> dplyr::mutate(Region = dplyr::case_when(stringr::str_detect(Sample, "Tran") ~ "Pwani",
                                                                    stringr::str_detect(Sample, "DOCW") | stringr::str_detect(Sample, "OMP") ~ "Dodoma",
                                                                    stringr::str_detect(Sample, "KG[A-Z]{2}") ~ "Kagera",
                                                                    stringr::str_detect(Sample, "KL[A-Z]{2}") ~ "Kilimanjaro",
                                                                    stringr::str_detect(Sample, "LUN-") ~ "Ruvuma",
                                                                    stringr::str_detect(Sample, "MAB-") | stringr::str_detect(Sample, "MPY-") ~ "Tanga",
                                                                    stringr::str_detect(Sample, "MA[A-Z]{2}") ~ "Mara",
                                                                    stringr::str_detect(Sample, "MT[A-Z]{2}") ~ "Mtwara",
                                                                    stringr::str_detect(Sample, "MY[A-Z]{2}") ~ "Manyara",
                                                                    stringr::str_detect(Sample, "NJ[A-Z]{2}") ~ "Njombe",
                                                                    stringr::str_detect(Sample, "NYA-") ~ "Kigoma",
                                                                    stringr::str_detect(Sample, "SO[A-Z]{2}") ~ "Songwe",
                                                                    stringr::str_detect(Sample, "TB[A-Z]{2}") ~ "Tabora",
                                                                    .default = "Geita"))

region_totals <- Tz_samples |> dplyr::count(Region)

region_totals <- region_totals |> dplyr::mutate(Type = dplyr::case_when(Region %in% c("Dodoma", "Kagera", "Kilimanjaro", "Manyara", "Mara", "Mtwara", "Njombe", "Songwe", "Tabora") ~ "HF",
                                                                        Region %in% c("Kigoma", "Ruvuma", "Tanga") ~ "CS",
                                                                        .default = "non-MSMT"))

type_totals <- region_totals |> dplyr::group_by(Type) |> dplyr::summarize(Total = sum(n))
