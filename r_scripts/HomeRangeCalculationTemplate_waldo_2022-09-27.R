#   Script Details                                                          ####

# Author: Waldemar Ortiz-Calo

# Date:2022-09-06 

# Purpose: This code is designed to model animal home ranges from the 3 study 
# populations. 

###############################################################################
#   Library / Functions / Data                                              ####

#      Library                                                              ####
library(ctmm)
library(adehabitatHR)
library(tidyverse)
library(move)
library(mapview)
library(sf)
library(foreach)
library(doParallel)
library(lubridate)

#      Functions                                                            ####

#      Data                                                                 ####

#        [Importing Data]                                                   ####

# Deer Data
deer_df <- read_csv("1.DataManagement/CleanData/deer_all_clean.csv") %>% 
  rename(individual.local.identifier = id,
         timestamp = t,
         location.long = x,
         location.lat = y)

#        [Adding Season classification]                                     ####

# Spring:  1 March – 31 May (3 mo.)
# Summer:  1 June – 30 August (3 mo.)
# Fall:  1 September – 30 November (3 mo.)
# Winter:  1 December – 28 February (3 mo.)

# Start
time_partition_start <- expand_grid(year = (deer_df$timestamp %>% year() %>% unique() %>% .[1]-1) : 
                                      (deer_df$timestamp %>% year() %>% unique() %>% .[length(.)]+1),
                                    month = c("03","06","09","12"),
                                    day = "01") %>% 
  mutate(season = ifelse(month == "03", "spring",
                         ifelse(month == "06", "summer",
                                ifelse(month == "09", "fall",
                                       ifelse(month == "12", "winter",NA))))) %>% 
  mutate(date = paste(year,month,day,sep = "-"),.keep = "unused",.before = 1) %>% 
  mutate(date = paste0(date, " 00:00:00"))

# End
time_partition_end <- expand_grid(year = (deer_df$timestamp %>% year() %>% unique() %>% .[1]-1) : 
                                    (deer_df$timestamp %>% year() %>% unique() %>% .[length(.)]+1),
                                  month = c("05","08","11","02"))

time_partition_end <- time_partition_end %>% 
  mutate(day = ifelse(time_partition_end$month == "05","31",
                      ifelse(time_partition_end$month == "08","31",
                             ifelse(time_partition_end$month == "11","30",
                                    ifelse(time_partition_end$month == "02","28",NA))))) %>% 
  mutate(year = ifelse(.$month == "02",year + 1, year)) %>% 
  mutate(date = paste(year,month,day,sep = "-"),.keep = "unused") %>% 
  mutate(date = paste0(date, " 23:59:59"))

# Creating Time Partition dataframe
time_partitions <- data.frame(start = time_partition_start$date,
                              end = time_partition_end$date,
                              season = time_partition_start$season)

ints <- list()

for (i in 1:nrow(time_partitions)) {
  ints[[i]] <- interval(time_partitions[i,1], time_partitions[i,2])
}

# Mutating dataframe and adding season
deer_df <- deer_df %>% 
  mutate(season = case_when(
    timestamp %within% ints[[1]] ~ time_partitions$season %>% unique() %>% .[1],
    timestamp %within% ints[[2]] ~ time_partitions$season %>% unique() %>% .[2],
    timestamp %within% ints[[3]] ~ time_partitions$season %>% unique() %>% .[3],
    timestamp %within% ints[[4]] ~ time_partitions$season %>% unique() %>% .[4],
    timestamp %within% ints[[5]] ~ time_partitions$season %>% unique() %>% .[1],
    timestamp %within% ints[[6]] ~ time_partitions$season %>% unique() %>% .[2],
    timestamp %within% ints[[7]] ~ time_partitions$season %>% unique() %>% .[3],
    timestamp %within% ints[[8]] ~ time_partitions$season %>% unique() %>% .[4],
    timestamp %within% ints[[9]] ~ time_partitions$season %>% unique() %>% .[1],
    timestamp %within% ints[[10]] ~ time_partitions$season %>% unique() %>% .[2],
    timestamp %within% ints[[11]] ~ time_partitions$season %>% unique() %>% .[3],
    timestamp %within% ints[[12]] ~ time_partitions$season %>% unique() %>% .[4],
    timestamp %within% ints[[13]] ~ time_partitions$season %>% unique() %>% .[1],
    timestamp %within% ints[[14]] ~ time_partitions$season %>% unique() %>% .[2],
    timestamp %within% ints[[15]] ~ time_partitions$season %>% unique() %>% .[3],
    timestamp %within% ints[[16]] ~ time_partitions$season %>% unique() %>% .[4],
    timestamp %within% ints[[17]] ~ time_partitions$season %>% unique() %>% .[1],
    timestamp %within% ints[[18]] ~ time_partitions$season %>% unique() %>% .[2],
    timestamp %within% ints[[19]] ~ time_partitions$season %>% unique() %>% .[3],
    timestamp %within% ints[[20]] ~ time_partitions$season %>% unique() %>% .[4],
    timestamp %within% ints[[21]] ~ time_partitions$season %>% unique() %>% .[1],
    timestamp %within% ints[[22]] ~ time_partitions$season %>% unique() %>% .[2],
    timestamp %within% ints[[23]] ~ time_partitions$season %>% unique() %>% .[3],
    timestamp %within% ints[[24]] ~ time_partitions$season %>% unique() %>% .[4],
    timestamp %within% ints[[25]] ~ time_partitions$season %>% unique() %>% .[1],
    timestamp %within% ints[[26]] ~ time_partitions$season %>% unique() %>% .[2],
    timestamp %within% ints[[27]] ~ time_partitions$season %>% unique() %>% .[3],
    timestamp %within% ints[[28]] ~ time_partitions$season %>% unique() %>% .[4],
    timestamp %within% ints[[29]] ~ time_partitions$season %>% unique() %>% .[1],
    timestamp %within% ints[[30]] ~ time_partitions$season %>% unique() %>% .[2],
    timestamp %within% ints[[31]] ~ time_partitions$season %>% unique() %>% .[3],
    timestamp %within% ints[[32]] ~ time_partitions$season %>% unique() %>% .[4],
    TRUE ~ "NA"),.after = timestamp)

# Fixing leap years
deer_df$season <- ifelse(month(deer_df$timestamp) == 2 & day(deer_df$timestamp) == 29, "winter",deer_df$season)

# Exporting csv file
write_csv(deer_df,
          file = "1.DataManagement/RawData/csv_files/deer_all_revised.csv",
          append = FALSE)

# Removing unnecessary objects and clearing RAM
rm(list = c("time_partition_start","time_partition_end","time_partitions","ints","i"))
gc()

###############################################################################
#   [North]                                                                 ####
#      [Data Prep]                                                          ####

# Creating Telemetry Object for the Data 

data <- deer_df %>% 
  filter(site == "North")

data_aggregated <- data

telemetry_object <- move(x=data$location.long, 
                         y=data$location.lat, time=as.POSIXct(data$timestamp, format="%Y-%m-%d %H:%M:%OS", tz="UTC"), 
                         proj=CRS("+init=epsg:5070"),
                         data=data, 
                         animal=data$individual.local.identifier) %>% 
  as.telemetry()


#      [Aggregated Home Ranges]                                             ####
#        [All Locations]                                                    ####
#           [Creating Dataframe]                                            ####

# Creating SpatialPointsDataFrame
spdf <- data_aggregated
coordinates(spdf) <- ~location.long + location.lat
proj4string(spdf) <- CRS("+init=epsg:5070")

#           [MCP]                                                           ####

# Creating MCP 
mcp <- mcp(spdf, percent=95, 
           unin = c("m"),
           unout = c("km2")) %>% 
  st_as_sf() %>% 
  st_write("1.DataManagement/HomeRangePolygons/north/AggregatedPolygons/mcp_north.shp",
           append=FALSE)

# Verifying Polygon

# mapview(mcp)

#           [KDE]                                                           ####

# Creating KDE
kde <- kernelUD(spdf, h = "href")
kde_UD <- getverticeshr(kde, 95)%>% 
  st_as_sf() %>% 
  st_write("1.DataManagement/HomeRangePolygons/north/AggregatedPolygons/kde_north.shp",
           append=FALSE)

# Verifying Polygon

# mapview(kde_UD)

#        [Aggregated - Seasonal Home Ranges]                                ####
#           [HR Calculation]                                                ####

# Creating SpatialPointsDataFrame
spdf <- data_aggregated
coordinates(spdf) <- ~location.long + location.lat
proj4string(spdf) <- CRS("+init=epsg:5070")

# Subsetting by season
seasons <- unique(spdf$season)

# Creating home ranges and exporting them
for (i in 1:length(seasons)) {
  
  print(i)
  
  # Subsetting by season
  spdf_season <- subset(spdf, spdf$season == seasons[i])
  
  # Creating MCP 
  mcp <- mcp(spdf_season, percent= 95, 
             unin = c("m"),
             unout = c("km2")) %>% 
    st_as_sf() %>% 
    st_write(paste0("1.DataManagement/HomeRangePolygons/north/AggregatedPolygons/mcp_north_",seasons[i],".shp"),
             append=FALSE)
  
  # Creating KDE
  kde <- kernelUD(spdf_season, h = "href",grid = 500, extent = 5)
  kde_UD <- getverticeshr(kde, 95)%>% 
    st_as_sf() %>% 
    st_write(paste0("1.DataManagement/HomeRangePolygons/north/AggregatedPolygons/kde_north_",seasons[i],".shp"),
             append=FALSE)
}

#        [Seasonal Home Ranges - By Sex]                                    ####

# Creating SpatialPointsDataFrame
spdf <- data_aggregated
coordinates(spdf) <- ~location.long + location.lat
proj4string(spdf) <- CRS("+init=epsg:5070")

# Subsetting parameters
seasons <- unique(spdf$season)
sex <- unique(spdf$sex)

for (i in 1:length(seasons)) {
  
  # Progress
  print(seasons[i])
  
  for(j in 1:length(sex)){
    
    # Progress
    print(sex[j])
    
    # Subsetting by season
    spdf_season <- spdf %>% subset(season == seasons[i]) %>% 
      subset(sex == sex[j])
    
    # Creating MCP 
    mcp <- mcp(spdf_season, percent= 95, 
               unin = c("m"),
               unout = c("km2")) %>% 
      st_as_sf() %>% 
      st_write(paste0("1.DataManagement/HomeRangePolygons/north/AggregatedPolygons/mcp_north_",sex[j],"_",seasons[i],".shp"),
               append=FALSE)
    
    # Creating KDE
    kde <- kernelUD(spdf_season, h = "href",grid = 500, extent = 5)
    kde_UD <- getverticeshr(kde, 95)%>% 
      st_as_sf() %>% 
      st_write(paste0("1.DataManagement/HomeRangePolygons/north/AggregatedPolygons/kde_north_",sex[j],"_",seasons[i],".shp"),
               append=FALSE)
  }
}

#      [Individual Home Ranges - All Locations]                             ####
#        [Setting Up Cluster for Parallel Computing]                        ####

# Setting the Settings for the Cluster
myCluster <- makeCluster(4,
                         type = "PSOCK") 

# Registering 
registerDoParallel(myCluster)

# Exporting Packages
clusterEvalQ(myCluster,
             {
               library(ctmm)
               library(tidyverse)
             })

# Exporting data to clusters
clusterExport(cl = myCluster, varlist = c("telemetry_object"), envir = environment())

#        [Home Range Calculation and Export]                                ####

ctmm_modeltypes <- foreach(i = 1:length(telemetry_object), 
                           .errorhandling="pass",
                           .combine = bind_rows) %dopar% {
                             
                             # Subsetting and Individual 
                             indiv <- telemetry_object[[i]]
                             
                             # Creating CTMM fit objects for the smoothing
                             M.IID <- ctmm.fit(indiv) 
                             GUESS <- ctmm.guess(indiv,interactive=FALSE) 
                             M.OUF <- ctmm.fit(indiv,GUESS) 
                             
                             # Home Range Polygons
                             KDE <- akde(indiv,M.IID) 
                             AKDE <- akde(indiv,M.OUF) 
                             
                             # Exporting
                             writeShapefile(object = KDE,
                                            folder = "1.DataManagement/HomeRangePolygons/north/IndividualPolygons",
                                            file=paste0(AKDE@info$identity,"_kde_all"),
                                            overwrite = T)
                             
                             writeShapefile(object = AKDE,
                                            folder = "1.DataManagement/HomeRangePolygons/north/IndividualPolygons",
                                            file=paste0(AKDE@info$identity,"_akde_all"),
                                            overwrite = T)
                             
                             # Model Details
                             model_types <- ctmm.select(indiv, CTMM=GUESS, verbose=TRUE) %>% 
                               summary(fitted.mods) %>% 
                               as.data.frame() %>% 
                               add_column(model = row.names(.), .before = 1) %>% 
                               add_column(ID = indiv@info[[1]], .before = 1)
                             
                             rownames(model_types)<- c()
                             
                             return(model_types)
                           }

# Exporting dataframe with Model info.
ctmm_modeltypes %>% rename(delta_AIC = ΔAICc,
                           delta_RMSPE = "ΔRMSPE (m)") %>% 
  write_csv(file = "1.DataManagement/HomeRangePolygons/north/IndividualPolygons/north_ctmm_models.csv",
            append = F)

#        [Closing Cluster]                                                  ####

# Stopping Cluster
stopCluster(myCluster)
###############################################################################
#   [South]                                                                 ####
#      [Data Prep]                                                          ####

# Creating Telemetry Object for the Data 

data <- deer_df %>% 
  filter(site == "South")

data_aggregated <- data

telemetry_object <- move(x=data$location.long, 
                         y=data$location.lat, time=as.POSIXct(data$timestamp, format="%Y-%m-%d %H:%M:%OS", tz="UTC"), 
                         proj=CRS("+init=epsg:5070"),
                         data=data, 
                         animal=data$individual.local.identifier) %>% 
  as.telemetry()


#      [Aggregated Home Ranges]                                             ####
#        [All Locations]                                                    ####
#           [Creating Dataframe]                                            ####

# Creating SpatialPointsDataFrame
spdf <- data_aggregated
coordinates(spdf) <- ~location.long + location.lat
proj4string(spdf) <- CRS("+init=epsg:5070")

#           [MCP]                                                           ####

# Creating MCP 
mcp <- mcp(spdf, percent=95, 
           unin = c("m"),
           unout = c("km2")) %>% 
  st_as_sf() %>% 
  st_write("1.DataManagement/HomeRangePolygons/south/AggregatedPolygons/mcp_south.shp",
           append=FALSE)

# Verifying Polygon

# mapview(mcp)

#           [KDE]                                                           ####

# Creating KDE
kde <- kernelUD(spdf, h = "href")
kde_UD <- getverticeshr(kde, 95)%>% 
  st_as_sf() %>% 
  st_write("1.DataManagement/HomeRangePolygons/south/AggregatedPolygons/kde_south.shp",
           append=FALSE)

# Verifying Polygon

# mapview(kde_UD)

#        [Aggregated - Seasonal Home Ranges]                                ####
#           [HR Calculation]                                                ####

# Creating SpatialPointsDataFrame
spdf <- data_aggregated
coordinates(spdf) <- ~location.long + location.lat
proj4string(spdf) <- CRS("+init=epsg:5070")

# Subsetting by season
seasons <- unique(spdf$season)

# Creating home ranges and exporting them
for (i in 1:length(seasons)) {
  
  print(i)
  
  # Subsetting by season
  spdf_season <- subset(spdf, spdf$season == seasons[i])
  
  # Creating MCP 
  mcp <- mcp(spdf_season, percent= 95, 
             unin = c("m"),
             unout = c("km2")) %>% 
    st_as_sf() %>% 
    st_write(paste0("1.DataManagement/HomeRangePolygons/south/AggregatedPolygons/mcp_south_",seasons[i],".shp"),
             append=FALSE)
  
  # Creating KDE
  kde <- kernelUD(spdf_season, h = "href",grid = 500, extent = 5)
  kde_UD <- getverticeshr(kde, 95)%>% 
    st_as_sf() %>% 
    st_write(paste0("1.DataManagement/HomeRangePolygons/south/AggregatedPolygons/kde_south_",seasons[i],".shp"),
             append=FALSE)
}


#        [Seasonal Home Ranges - By Sex]                                    ####

# Creating SpatialPointsDataFrame
spdf <- data_aggregated
coordinates(spdf) <- ~location.long + location.lat
proj4string(spdf) <- CRS("+init=epsg:5070")

# Subsetting parameters
seasons <- unique(spdf$season)
sex <- unique(spdf$sex)

for (i in 1:length(seasons)) {
  
  # Progress
  print(seasons[i])
  
  for(j in 1:length(sex)){
    
    # Progress
    print(sex[j])
    
    # Subsetting by season
    spdf_season <- spdf %>% subset(season == seasons[i]) %>% 
      subset(sex == sex[j])
    
    # Creating MCP 
    mcp <- mcp(spdf_season, percent= 95, 
               unin = c("m"),
               unout = c("km2")) %>% 
      st_as_sf() %>% 
      st_write(paste0("1.DataManagement/HomeRangePolygons/south/AggregatedPolygons/mcp_south_",sex[j],"_",seasons[i],".shp"),
               append=FALSE)
    
    # Creating KDE
    kde <- kernelUD(spdf_season, h = "href",grid = 500, extent = 5)
    kde_UD <- getverticeshr(kde, 95)%>% 
      st_as_sf() %>% 
      st_write(paste0("1.DataManagement/HomeRangePolygons/south/AggregatedPolygons/kde_south_",sex[j],"_",seasons[i],".shp"),
               append=FALSE)
  }
}

#      [Individual Home Ranges - All Locations]                             ####
#        [Setting Up Cluster for Parallel Computing]                        ####

# Setting the Settings for the Cluster
myCluster <- makeCluster(5,
                         type = "PSOCK") 

# Registering 
registerDoParallel(myCluster)

# Exporting Packages
clusterEvalQ(myCluster,
             {
               library(ctmm)
               library(tidyverse)
             })

# Exporting data to clusters
clusterExport(cl = myCluster, varlist = c("telemetry_object"), envir = environment())

#        [Home Range Calculation and Export]                                ####

ctmm_modeltypes <- foreach(i = 1:length(telemetry_object), 
                           .errorhandling="pass",
                           .combine = bind_rows) %dopar% {
                             
                             # Subsetting and Individual 
                             indiv <- telemetry_object[[i]]
                             
                             # Creating CTMM fit objects for the smoothing
                             M.IID <- ctmm.fit(indiv) 
                             GUESS <- ctmm.guess(indiv,interactive=FALSE) 
                             M.OUF <- ctmm.fit(indiv,GUESS) 
                             
                             # Home Range Polygons
                             KDE <- akde(indiv,M.IID) 
                             AKDE <- akde(indiv,M.OUF) 
                             
                             # Exporting
                             writeShapefile(object = KDE,
                                            folder = "1.DataManagement/HomeRangePolygons/south/IndividualPolygons",
                                            file=paste0(AKDE@info$identity,"_kde_all"),
                                            overwrite = T)
                             
                             writeShapefile(object = AKDE,
                                            folder = "1.DataManagement/HomeRangePolygons/south/IndividualPolygons",
                                            file=paste0(AKDE@info$identity,"_akde_all"),
                                            overwrite = T)
                             
                             # Model Details
                             model_types <- ctmm.select(indiv, CTMM=GUESS, verbose=TRUE) %>% 
                               summary(fitted.mods) %>% 
                               as.data.frame() %>% 
                               add_column(model = row.names(.), .before = 1) %>% 
                               add_column(ID = indiv@info[[1]], .before = 1)
                             
                             rownames(model_types)<- c()
                             
                             return(model_types)
                           }

# Exporting dataframe with Model info.
ctmm_modeltypes %>% rename(delta_AIC = ΔAICc,
                           delta_RMSPE = "ΔRMSPE (m)") %>% 
  write_csv(file = "1.DataManagement/HomeRangePolygons/south/IndividualPolygons/south_ctmm_models.csv",
            append = F)

#        [Closing Cluster]                                                  ####

# Stopping Cluster
stopCluster(myCluster)
###############################################################################
#   [Southeast]                                                             ####
#      [Data Prep]                                                          ####

# Creating Telemetry Object for the Data 

data <- deer_df %>% 
  filter(site == "CroplandStudy")

data_aggregated <- data

telemetry_object <- move(x=data$location.long, 
                         y=data$location.lat, time=as.POSIXct(data$timestamp, format="%Y-%m-%d %H:%M:%OS", tz="UTC"), 
                         proj=CRS("+init=epsg:5070"),
                         data=data, 
                         animal=data$individual.local.identifier) %>% 
  as.telemetry()


#      [Aggregated Home Ranges]                                             ####
#        [All Locations]                                                    ####

#           [Creating Dataframe]                                            ####

# Creating SpatialPointsDataFrame
spdf <- data_aggregated
coordinates(spdf) <- ~location.long + location.lat
proj4string(spdf) <- CRS("+init=epsg:5070")

#           [MCP]                                                           ####

# Creating MCP 
mcp <- mcp(spdf, percent=95, 
           unin = c("m"),
           unout = c("km2")) %>% 
  st_as_sf() %>% 
  st_write("1.DataManagement/HomeRangePolygons/southeast/AggregatedPolygons/mcp_southeast.shp",
           append=FALSE)

# Verifying Polygon

# mapview(mcp)

#           [KDE]                                                           ####

# Creating KDE
kde <- kernelUD(spdf, h = "href")
kde_UD <- getverticeshr(kde, 95)%>% 
  st_as_sf() %>% 
  st_write("1.DataManagement/HomeRangePolygons/southeast/AggregatedPolygons/kde_southeast.shp",
           append=FALSE)

# Verifying Polygon

# mapview(kde_UD)

#        [Aggregated - Seasonal Home Ranges]                                ####
#           [HR Calculation]                                                ####

# Creating SpatialPointsDataFrame
spdf <- data_aggregated
coordinates(spdf) <- ~location.long + location.lat
proj4string(spdf) <- CRS("+init=epsg:5070")

# Subsetting by season
seasons <- unique(spdf$season)

# Creating home ranges and exporting them
for (i in 1:length(seasons)) {
  
  print(i)
  
  # Subsetting by season
  spdf_season <- subset(spdf, spdf$season == seasons[i])
  
  # Creating MCP 
  mcp <- mcp(spdf_season, percent= 95, 
             unin = c("m"),
             unout = c("km2")) %>% 
    st_as_sf() %>% 
    st_write(paste0("1.DataManagement/HomeRangePolygons/southeast/AggregatedPolygons/mcp_southeast_",seasons[i],".shp"),
             append=FALSE)
  
  # Creating KDE
  kde <- kernelUD(spdf_season, h = "href",grid = 500, extent = 5)
  kde_UD <- getverticeshr(kde, 95)%>% 
    st_as_sf() %>% 
    st_write(paste0("1.DataManagement/HomeRangePolygons/southeast/AggregatedPolygons/kde_southeast_",seasons[i],".shp"),
             append=FALSE)
}

#        [Seasonal Home Ranges - By Sex]                                    ####

# Creating SpatialPointsDataFrame
spdf <- data_aggregated
coordinates(spdf) <- ~location.long + location.lat
proj4string(spdf) <- CRS("+init=epsg:5070")

# Subsetting parameters
seasons <- unique(spdf$season)
sex <- unique(spdf$sex)

for (i in 1:length(seasons)) {
  
  # Progress
  print(seasons[i])
  
  for(j in 1:length(sex)){
    
    # Progress
    print(sex[j])
    
    # Subsetting by season
    spdf_season <- spdf %>% subset(season == seasons[i]) %>% 
      subset(sex == sex[j])
    
    # Creating MCP 
    mcp <- mcp(spdf_season, percent= 95, 
               unin = c("m"),
               unout = c("km2")) %>% 
      st_as_sf() %>% 
      st_write(paste0("1.DataManagement/HomeRangePolygons/southeast/AggregatedPolygons/mcp_southeast_",sex[j],"_",seasons[i],".shp"),
               append=FALSE)
    
    # Creating KDE
    kde <- kernelUD(spdf_season, h = "href",grid = 500, extent = 5)
    kde_UD <- getverticeshr(kde, 95)%>% 
      st_as_sf() %>% 
      st_write(paste0("1.DataManagement/HomeRangePolygons/southeast/AggregatedPolygons/kde_southeast_",sex[j],"_",seasons[i],".shp"),
               append=FALSE)
  }
}

#      [Individual Home Ranges - All Locations]                             ####
#        [Setting Up Cluster for Parallel Computing]                        ####

# Setting the Settings for the Cluster
myCluster <- makeCluster(5,
                         type = "PSOCK") 

# Registering 
registerDoParallel(myCluster)

# Exporting Packages
clusterEvalQ(myCluster,
             {
               library(ctmm)
               library(tidyverse)
             })

# Exporting data to clusters
clusterExport(cl = myCluster, varlist = c("telemetry_object"), envir = environment())

#        [Home Range Calculation and Export]                                ####

ctmm_modeltypes <- foreach(i = 1:length(telemetry_object), 
                           .errorhandling="pass",
                           .combine = bind_rows) %dopar% {
                             
                             # Subsetting and Individual 
                             indiv <- telemetry_object[[i]]
                             
                             # Creating CTMM fit objects for the smoothing
                             M.IID <- ctmm.fit(indiv) 
                             GUESS <- ctmm.guess(indiv,interactive=FALSE) 
                             M.OUF <- ctmm.fit(indiv,GUESS) 
                             
                             # Home Range Polygons
                             KDE <- akde(indiv,M.IID) 
                             AKDE <- akde(indiv,M.OUF) 
                             
                             # Exporting
                             writeShapefile(object = KDE,
                                            folder = "1.DataManagement/HomeRangePolygons/southeast/IndividualPolygons",
                                            file=paste0(AKDE@info$identity,"_kde_all"),
                                            overwrite = T)
                             
                             writeShapefile(object = AKDE,
                                            folder = "1.DataManagement/HomeRangePolygons/southeast/IndividualPolygons",
                                            file=paste0(AKDE@info$identity,"_akde_all"),
                                            overwrite = T)
                             
                             # Model Details
                             model_types <- ctmm.select(indiv, CTMM=GUESS, verbose=TRUE) %>% 
                               summary(fitted.mods) %>% 
                               as.data.frame() %>% 
                               add_column(model = row.names(.), .before = 1) %>% 
                               add_column(ID = indiv@info[[1]], .before = 1)
                             
                             rownames(model_types)<- c()
                             
                             return(model_types)
                           }

# Exporting dataframe with Model info.
ctmm_modeltypes %>% rename(delta_AIC = ΔAICc,
                           delta_RMSPE = "ΔRMSPE (m)") %>% 
  write_csv(file = "1.DataManagement/HomeRangePolygons/southeast/IndividualPolygons/southeast_ctmm_models.csv",
            append = F)

#        [Closing Cluster]                                                  ####

# Stopping Cluster
stopCluster(myCluster)
###############################################################################