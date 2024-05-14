#script for dynamic rasterization of LegacyVegetation data
rm(list = ls())

#necessary packages
packages <- c("raster","tidyverse","rnaturalearth","sf","pangaear","data.table")

lapply(packages,
       require,
       character.only =TRUE,
       quietly = TRUE)

#settings
x_resolution <- 5 #in degrees
y_resolution <- 5 #in degrees
temporal_resolution <- 500 #in yrs 14000 should be divisible by temp resolution without remainder!
taxon <- "forest" #taxon/PFT coverage to grid

#load reconstructions for all three NH continents
asia <- "963610/files/composition_forest_REVEALS_Asia.csv"
north_america <- "963612/files/composition_forest_REVEALS_North_America.csv"
europe <- "963611/files/composition_forest_REVEALS_Europe.csv"

links <- c(asia, north_america, europe)

get_data <- function(link){
  link <- paste0("https://download.pangaea.de/dataset/",link)
  necessary_columns <- c("Dataset_ID", "Age_mean [yrs BP]", "Longitude", "Latitude", taxon)
  df <- data.table::fread(link,
                          select = necessary_columns)
  return(df)
}

NH_timeseries <- lapply(links,
                        get_data) %>% 
  bind_rows()

NH_timeseries %>% 
  as.data.frame() %>% 
  mutate(Age_slice = cut(`Age_mean [yrs BP]`,
                         seq(0,14000,temporal_resolution),
                         include.lowest = TRUE,
                         labels = paste(seq(temporal_resolution,14000,temporal_resolution),"to",
                                        seq(0,14000 - temporal_resolution,temporal_resolution))))%>% 
  arrange(desc(Age_slice)) %>% 
  group_by(Dataset_ID,Age_slice,Longitude,Latitude) %>% 
  summarize(forest = mean(get(taxon)))
