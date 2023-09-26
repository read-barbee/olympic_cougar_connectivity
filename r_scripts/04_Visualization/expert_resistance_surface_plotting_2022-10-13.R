library(sf)
library(raster)
library(terra)
library(mapview)
library(tidyverse)
library(ggmap)

exp <-rast("/Users/tb201494/Desktop/Cougar Cost-Weighted Distance1.tif", na)
  

exp_reproj <- project(x = exp, 
                      y = "epsg:4326")
plot(exp_reproj)

exp_breaks <- c(0, 1, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 774936.5)

labels = c(0, 1, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 774936.5)
  
mapview(raster(exp_reproj),
        at= exp_breaks,
        trim = TRUE,
        map.types= "Esri.WorldImagery")
