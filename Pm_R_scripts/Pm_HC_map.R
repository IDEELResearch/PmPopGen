#load required libraries for base maps
library(tidyverse)
library(rnaturalearth) 
library(rnaturalearthdata)
library(sf)
library(ggspatial)
library(colorspace)
library(dplyr)
library(ggpattern)
library(cowplot)

setwd("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/Pm and Po Sequencing/Twist Pm/rerun")

sample_metadata <- readxl::read_xlsx("Attributes_Upload_Microbe.xlsx", skip = 12) |> dplyr::select(`*sample_name`, `*geo_loc_name`)

colnames(sample_metadata) <- c("Sample", "Location")

map_data <- sample_metadata |> dplyr::group_by(Location) |> dplyr::count()

map_data <- map_data |> tidyr::separate_wider_delim(Location, ": ", names = c("Country", "Region", "Subregion", "Locality"), too_few = "align_start")

map_data |> writexl::write_xlsx("Pm_sample_locations.xlsx")

#Added coordinates then re-imported

map_data <- readxl::read_xlsx("Pm_sample_locations.xlsx")

map_data <- map_data |> tidyr::separate_wider_delim(Coords, ",", names = c("Lat", "Long"))

map_data$Lat <- as.numeric(map_data$Lat)

map_data$Long <- as.numeric(map_data$Long)

map_data$Region <- map_data$Region |> stringr::str_remove(" State") #fixes an issue specific to Cross River State in Nigeria

map_data$Region <- map_data$Region |> stringr::str_replace("Adamawa", "Adamaoua") #fixes an issue where code matches Adamawa in Nigeria instead of Adamawa in Cameroon

africa <- ne_countries(scale="medium", type = "sovereignty", continent = "Africa", returnclass = "sf")

study_countries <- geodata::gadm(country = c("Cameroon", "Democratic Republic of the Congo", "Nigeria", "United Republic of Tanzania"), level = 0, path = "study_countries") |> tidyterra::as_sf()

study_regions <- geodata::gadm(country = c("CMR", "COD", "NGA", "TZA"), level = 1, path = "study_regions") |> tidyterra::as_sf()

study_regions2 <- study_regions |> dplyr::semi_join(map_data, by = dplyr::join_by(NAME_1 == Region))

study_regions3 <- study_regions |> dplyr::semi_join(map_data, by = dplyr::join_by(VARNAME_1 == Region))

study_regions <- rbind(study_regions2, study_regions3)

water<-st_read("C:/Users/zpopkinh/OneDrive - University of North Carolina at Chapel Hill/MSMT non-falciparum Projects/Africa_waterbody.shp")

Pm_map <- ggplot() + 
  geom_sf(data=africa, fill="gray90")+
  geom_sf(data = study_countries, fill="gray96", lwd=1) +
  geom_sf(data = study_regions, fill="gray96") +
  geom_sf(data = water, fill="lightblue") +
  annotate("label", x = 12.2, y = 5, label = "Cameroon", color="grey30", size=2.3, fontface="bold", angle = 30, fill = "gray96", label.size = NA) +
  annotate("label", x = 23.8, y = -3, label = "DRC", color="grey30", size=5 , fontface="bold", fill = "gray96", label.size = NA) +
  annotate("label", x = 8, y = 9.6, label = "Nigeria", color="grey30", size=5 , fontface="bold", fill = "gray96", label.size = NA) +
  annotate("label", x = 34.5, y = -6, label = "Tanzania", color="grey30", size=3 , fontface="bold", angle = 30, fill = "gray96", label.size = NA) +
 geom_point(data=map_data, aes(x=Long, y=Lat, size=n), shape=21, alpha=1, fill = "red", position = "jitter")+
  scale_size_continuous(breaks = c(1, 2, 5, 9)) +
  labs(size = "# of Samples") +
  #ggtitle(expression(paste(italic("P. malariae "), "Sequencing Sample Locations"))) +
  ylim(-33, 37) +
  theme_minimal() + theme(panel.background = element_rect(fill = 'lightblue', color = "black"), plot.background = element_rect(fill = "white", color = "white"),
                          axis.title = element_blank())

ggsave("Pm_sample_collections.png", dpi=600, width=10, height=10, units = "in")

knitr::plot_crop("Pm_sample_collections.png")
