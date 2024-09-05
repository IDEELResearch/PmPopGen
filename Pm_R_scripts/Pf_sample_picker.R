setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun")

Pf_samples <- readxl::read_xlsx("Pf7_samples_for_zach.xlsx")

Oscar_samples <- data.table::fread("samples.txt", header = F) |> dplyr::rename(Sample = V1)

Oscar_samples$Available <- "Yes"

Pf_samples <- Pf_samples |> subset(`Exclusion reason` == "Analysis_set")

callable_Pf <- Pf_samples |> subset(`% callable` >= 85)

#unique(callable_Pf$Country)

relevant_Pf <- callable_Pf |> subset(Country == "Democratic Republic of the Congo" | Country == "Tanzania" | Country == "Cameroon")

relevant_Pf <- relevant_Pf |> dplyr::left_join(Oscar_samples)

relevant_Pf <- relevant_Pf |> subset(Available == "Yes")

Pf_COI <- data.table::fread("COI_for_Kelly.csv")

Pf_monoclonals <- Pf_COI |> subset(mean_COI == 1) |> dplyr::rename(Sample = sample)

relevant_Pf <- relevant_Pf |> dplyr::semi_join(Pf_monoclonals)

#unique(relevant_Pf$`Admin level 1`)

Pf_regions <- relevant_Pf |> dplyr::group_by(Country) |> dplyr::count(`Admin level 1`)

monoclonals <- data.table::fread("monoclonals.csv", header = F)

monoclonals <- monoclonals |> dplyr::mutate(Country = dplyr::case_when(stringr::str_detect(V1, "Gam_[:digit:]+") ~ "Nigeria",
                                                             stringr::str_detect(V1, "^[:digit:]+") ~ "DRC",
                                                             stringr::str_detect(V1, "P[A-Z][:digit:]{3}") ~ "Cameroon",
                                                             stringr::str_detect(V1, "Gam") ~ "Cameroon",
                                                             .default = "Tanzania"))

#country_counts <- data.table::fread("../country_counts.csv")

country_counts <- monoclonals |> dplyr::count(Country)

Tz_counts <- monoclonals |> subset(Country == "Tanzania") |> dplyr::mutate(Region = dplyr::case_when(stringr::str_detect(V1, "DOCW") | stringr::str_detect(V1, "DOMP") ~ "Dodoma",
                                                                                                     stringr::str_detect(V1, "^KG") ~ "Kagera",
                                                                                                     stringr::str_detect(V1, "^KL") ~ "Kilimanjaro",
                                                                                                     stringr::str_detect(V1, "MABU") | stringr::str_detect(V1, "MAMS") | stringr::str_detect(V1, "MARO") | stringr::str_detect(V1, "MATA") ~ "Mara",
                                                                                                     stringr::str_detect(V1, "^MT") ~ "Mtwara",
                                                                                                     stringr::str_detect(V1, "^MY") ~ "Manyara",
                                                                                                     stringr::str_detect(V1, "^NJ") ~ "Njombe",
                                                                                                     stringr::str_detect(V1, "^SO") ~ "Songwe",
                                                                                                     stringr::str_detect(V1, "^TB") ~ "Tabora",
                                                                                                     stringr::str_detect(V1, "^Mq") | stringr::str_detect(V1, "^Trans") ~ "Pwani",
                                                                                                     .default = "Geita"))

#Tz_region_counts <- data.table::fread("../../Pm_Po_Sample_Locations_ALL.csv") |> subset(Species == "Pm" & Country == "Tanzania")

#So, essentially I'm just going to try and get as close as possible. Songwe and Mara don't border any regions with Pf samples, so for Songwe I'll sub Kigoma and for Mara I'll sub Kagera

Tz_counts <- Tz_counts |> dplyr::mutate(Closest_Region = dplyr::case_when(Region == "Pwani" ~ "Tanga",
                                                                       Region == "Dodoma" ~ "Morogoro",
                                                                       Region == "Geita" ~ "Kagera",
                                                                       Region == "Kagera" ~ "Kagera",
                                                                       Region == "Kigoma" ~ "Kigoma",
                                                                       Region == "Kilimanjaro" ~ "Tanga",
                                                                       Region == "Manyara" ~ "Tanga",
                                                                       Region == "Mara" ~ "Kagera",
                                                                       Region == "Mtwara" ~ "Lindi",
                                                                       Region == "Njombe" ~ "Morogoro",
                                                                       Region == "Ruvuma" ~ "Morogoro",
                                                                       Region == "Songwe" ~ "Kigoma",
                                                                       Region == "Tabora" ~ "Kigoma",
                                                                       Region == "Tanga" ~ "Tanga"))

Tz_Pf_counts <- Tz_counts |> dplyr::group_by(Closest_Region) |> dplyr::select(V1, Closest_Region) |> dplyr::count(Closest_Region) 

relevant_Pf_no_Tz <- relevant_Pf |> subset(Country == "Cameroon" | Country == "Democratic Republic of the Congo")

Pf_ortholog_samples <- relevant_Pf_no_Tz |> splitstackshape::stratified("Country", c(`Democratic Republic of the Congo` = 16, Cameroon = 10))

relevant_Pf_Tz_only <- relevant_Pf |> subset(Country == "Tanzania")

Tz_Pf_samples <- relevant_Pf_Tz_only |> splitstackshape::stratified("Admin level 1", c(Tanga = 11, Morogoro = 7, Kagera = 26, Lindi = 1))

Pf_ortholog_samples <- rbind(Pf_ortholog_samples, Tz_Pf_samples)

relevant_Pf$Sample |> as.list() |> data.table::fwrite("Pf_COI_samples.txt")

Pf_ortholog_samples$Sample |> as.list() |> data.table::fwrite("Pf_ortholog_samples.txt", sep = "\n")

relevant_Pf |> writexl::write_xlsx("Pf_COI_samples.xlsx")

Pf_ortholog_samples |> writexl::write_xlsx("Pf_ortholog_samples.xlsx")

#redoing this to extract all monoclonal samples from relevant countries

Pf_COI <- data.table::fread("COI_for_Kelly.csv") |> dplyr::rename(Sample = sample)

Pf_samples <- relevant_Pf |> dplyr::left_join(Pf_COI) |> subset(mean_COI == 1) |> dplyr::select(Sample)

Pf_samples |> data.table::fwrite("Pf_ortholog_samples.txt", sep = "\n", col.names = F)

country_count <- relevant_Pf |> dplyr::count(Country)
