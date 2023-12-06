library(sf)
library(raster)
library(terra)
library(mapview)
library(tidyverse)
library(ggmap)

exp <-rast("fitted_models/Model_Outputs/expert_opinion/Cougar Cost-Weighted Distance1.tif")
  
exp_ag <- terra::aggregate(exp, fact = 10)

exp_reproj <- project(x = exp_ag, 
                      y = "epsg:5070")
plot(exp_reproj)

exp_breaks <- c(0, 1, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 774936.5)

labels = c(0, 1, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 774936.5)
  
mapview(raster(exp_reproj),
        at= exp_breaks,
        trim = TRUE,
        map.types= "Esri.WorldImagery")

pred_vals <- terra::values(exp_reproj)

breaks <- quantile(pred_vals, probs = 0:10/10, na.rm = T)

bins <- cut(pred_vals, breaks, include.lowest=TRUE)

binned <- terra::classify(exp_reproj, rcl=as.vector(breaks))

terra::plot(binned)

levels(binned) <- 1:10


#viridis magma
ggplot()+
  tidyterra::geom_spatraster(data=binned, mapping=aes()) +
  #geom_sf(data=red_deer_used, aes(geometry=geometry))+
  #coord_sf(datum = st_crs(5070))+
  scale_fill_manual(values = viridis::magma(10), na.value = NA, guide = guide_legend(reverse = TRUE), na.translate=FALSE)  #rev() to reverse pallete

mapview::mapview(raster::raster(binned))