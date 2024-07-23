#script for dynamic rasterization of LegacyVegetation data
rm(list = ls())

#necessary packages
packages <- c("raster","tidyverse","rnaturalearth","sf","pangaear","data.table","LegacyData")

lapply(packages,
       require,
       character.only =TRUE,
       quietly = TRUE)

#settings
x_resolution <- 5 #in degrees
y_resolution <- 5 #in degrees
temporal_resolution <- 500 #in yrs (14000 should be divisible by temp resolution without remainder!)

#load reconstructions for all three NH continents
# asia <- "12793806/files/Asia_complete_revised.csv?download=1"
# north_america <- "12793806/files/North_America_complete_revised.csv?download=1"
# europe <- "12793806/files/Europe_complete_revised.csv?download=1"
asia <- "Asia_complete_revised_SD.csv"
europe <- "Europe_complete_revised_SD.csv"
north_america <- "North_America_complete_revised_SD.csv"

links <- c(asia, north_america, europe)

# get_data <- function(link){
#   path <- "https://zenodo.org/records/"
#   data.table::fread(paste0(path,link)) %>% 
#     as.data.frame() %>% 
#     return()
# }

get_data <- function(link){
  path <- "//dmawi.de/potsdam/data/bioing/user/lschild/optimize_REVEALS/output/reconstructions/for_upload/"
  data.table::fread(paste0(path,link)) %>% 
    as.data.frame() %>% 
    return()
}

NH_timeseries <- lapply(links,
                        get_data) %>% 
  bind_rows()

NH_timeseries[is.na(NH_timeseries)] <- 0

delta_method <- function(v){
  v <- v[!(is.na(v))]
  n <- length(v)
  v^2 %>% 
    sum() %>% 
    sqrt() %>% 
    '/'(n) %>% 
    return()
}

  NH_sliced <- NH_timeseries %>% 
    as.data.frame() %>% 
    mutate(Age_slice = cut(`Age_mean [yrs BP]`,
                           seq(-500,14000,temporal_resolution),
                           include.lowest = TRUE,
                           labels = paste(seq(0,14000,temporal_resolution),"to",
                                          seq(-500,14000 - temporal_resolution,temporal_resolution))))%>% 
    arrange(desc(Age_slice)) %>% 
    dplyr::select(-contains("yrs BP")) %>% 
    group_by(Dataset_ID,Age_slice,Longitude,Latitude,Continent,valid_as_site) %>% 
    summarise(across(contains("mean of cover"), mean),
              across(contains("sd of cover"), delta_method)) %>% 
    ungroup()


#grid the data
#create grid

grid_sizes <- rev(c(1,2,5,10))

for(size in grid_sizes){
  y_resolution <- size
  x_resolution <- size
  
  
  grid <- raster(nrow = 90/y_resolution,
                 ncol = 180/x_resolution,
                 ymn = 0,
                 ymx = 90)
  values(grid)  <- 1:(90/y_resolution*180/x_resolution)
  
  #get record coordinates
  coords <- NH_sliced %>% 
    ungroup() %>% 
    distinct(Dataset_ID,Longitude, Latitude) %>% 
    st_as_sf(coords = c("Longitude","Latitude"),
             crs = 4326)
  
  #extract grid cell assignments
  cells <- raster::extract(grid,
                           coords) %>% 
    cbind(coords$Dataset_ID) %>% 
    as.data.frame() 
  
  names(cells) <- c("cell","Dataset_ID")
  
  NH_grid <- NH_sliced %>% 
    ungroup() %>% 
    merge(cells,
          by= "Dataset_ID") %>% 
    mutate(validity = ifelse(valid_as_site == TRUE,
                             2,
                             1)) %>% 
    dplyr::select(-valid_as_site, -Longitude, -Latitude, -Dataset_ID,-Continent) %>% 
    filter(!(is.na(Age_slice))) %>% 
    group_by(Age_slice, cell) %>% 
    mutate(validity = sum(validity)) %>% 
    summarise(validity = mean(validity),
              validity = validity > 1,
              across(contains("mean of cover"), mean),
              across(contains("sd of cover"), delta_method)) 
  
  #plot gridded data
  world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")
  
  # grid %>% 
  #   rasterToPolygons() %>% 
  #   st_as_sf() %>% 
  #   rename(cell = layer) %>% 
  #   merge(NH_grid,
  #         by = "cell",
  #         all.y = TRUE) %>% 
  #   mutate(validity = validity >1) %>% 
  #   filter(validity) %>% 
  #   ggplot()+
  #   geom_sf(data = world,
  #           col = "darkgrey",
  #           fill = NA)+
  #   geom_sf(aes(fill = forest_mean,
  #               col = validity))+
  #   coord_sf(ylim = c(0,90))+
  #   scale_color_manual(values = c("red",NA),
  #                      limits = c(FALSE,TRUE),
  #                      labels = c("not enough records"),
  #                      na.value = NA)+
  #   facet_wrap(.~Age_slice)+
  #   labs(col = "",
  #        fill = "forest cover")+
  #   theme_void()+
  #   scale_fill_viridis_c(direction = -1)+
  #   theme(legend.position = "bottom")
  # 
  # #save plot
  # ggsave(paste0("gridded_overview_",x_resolution,".png"),
  #        width = 10,
  #        height = 5)
  

  #save data
  merged_grid <-grid %>% 
    rasterToPolygons() %>% 
    st_as_sf() %>% 
    rename(cell = layer) %>% 
    merge(NH_grid,
          by = "cell",
          all.y = TRUE)
  
  coordinates <-merged_grid %>% 
    st_coordinates() %>% 
    as.data.frame() %>% 
    group_by(L2) %>% 
    summarize(xmin = min(X),
              xmax = max(X),
              ymin = min(Y),
              ymax = max(Y)) %>% 
    dplyr::select(-L2)
  
  #save empty grid to pair with later
  grid %>% 
    rasterToPolygons() %>% 
    st_as_sf() %>% 
    rename(cell = layer) %>% 
    write_sf(paste0("empty_grid_",
                    size,
                    ".shp"))
    #save as csv
    merged_grid %>% 
      st_drop_geometry() %>% 
      cbind(coordinates) %>% 
      data.table::fwrite(paste0("SD_gridded_forest_",
                                x_resolution,
                                ".csv"))
  
 
}

