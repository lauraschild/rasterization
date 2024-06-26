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

NH_sliced <- NH_timeseries %>% 
  as.data.frame() %>% 
  mutate(Age_slice = cut(`Age_mean [yrs BP]`,
                         seq(0,14000,temporal_resolution),
                         include.lowest = TRUE,
                         labels = paste(seq(temporal_resolution,14000,temporal_resolution),"to",
                                        seq(0,14000 - temporal_resolution,temporal_resolution),"BP")))%>% 
  arrange(desc(Age_slice)) %>% 
  group_by(Dataset_ID,Age_slice) %>% 
  summarize(forest = mean(get(taxon))) %>% 
  merge(read.csv("meta_data.csv"),
        by = "Dataset_ID")

#grid the data
#create grid
grid <- raster(nrow = 90/y_resolution,
               ncol = 180/x_resolution,
               ymn = 0,
               ymx = 90)
values(grid)  <- 1:(90/y_resolution*180/x_resolution)

#get record coordinates
coords <- NH_sliced %>% 
  distinct(Dataset_ID,Longitude, Latitude) %>% 
  st_as_sf(coords = c("Longitude","Latitude"),
           crs = 4326)

#extract grid cell assignments
cells <- raster::extract(grid,
                         coords) %>% 
  cbind(coords) %>% 
  rename(cell = ".") %>% 
  st_drop_geometry()

NH_grid <-NH_sliced %>% 
  merge(cells,
        by= "Dataset_ID") %>% 
  mutate(validity = ifelse(record_type == "large_lake",
                           2,
                           1)) %>% 
  group_by(Age_slice, cell) %>% 
  summarize(forest = mean(forest),
            validity = sum(validity),
            validity = (validity >=2)) %>% 
  filter(!(is.na(Age_slice)))

#plot gridded data
world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")

grid %>% 
  rasterToPolygons() %>% 
  st_as_sf() %>% 
  rename(cell = layer) %>% 
  merge(NH_grid,
        by = "cell",
        all.y = TRUE) %>% 
  ggplot()+
  geom_sf(data = world,
          col = "darkgrey",
          fill = NA)+
  geom_sf(aes(fill = forest,
              col = validity))+
  coord_sf(ylim = c(0,90))+
  scale_color_manual(values = c("red",NA),
                     limits = c(FALSE),
                     labels = c("not enough records"))+
  facet_wrap(.~Age_slice)+
  labs(col = "",
       fill = "forest cover")+
  theme_void()+
  scale_fill_viridis_c(direction = -1)+
  theme(legend.position = "bottom")

#save plot
ggsave(paste0("gridded_overview_",x_resolution,".png"),
       width = 10,
       height = 5)

#save data as shapefile
grid %>% 
  rasterToPolygons() %>% 
  st_as_sf() %>% 
  rename(cell = layer) %>% 
  merge(NH_grid,
        by = "cell",
        all.y = TRUE) %>% 
  write_sf(paste0("gridded_forest_",x_resolution,".shp"))
