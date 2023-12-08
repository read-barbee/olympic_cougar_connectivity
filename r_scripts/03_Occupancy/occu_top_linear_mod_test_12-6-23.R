#### Global model all subsets selection occupancy annual covariates####

# Author: Read Barbee

# Date:2023-11-26

# Purpose:


################################ Libraries #################################
library(tidyverse)
library(sf)
library(terra)
library(ubms)
library(beepr)
library(DataExplorer)


#########################################################################
##
##  1. Import stacked occupancy data
##
######################################################################

#activity/detection data
occ_dat <- read_csv("data/Camera_Data/master/ocp_onp_occ_dat_annual_covs_14_day_period_11-21-23.csv")

#########################################################################
##
##  2. Make objects for ubms
##
##########################################################################
occ_dat_scaled <- occ_dat %>% 
  select(!contains("usfs")) %>% 
  mutate(across(c(elevation:dist_all_roads_annual), scale)) %>% 
  mutate(across(c(elevation:dist_all_roads_annual), as.numeric))


#remove station rows with missing covariate values
complete_cases <- occ_dat_scaled %>% select(elevation:dist_all_roads_annual) %>% 
  complete.cases()

occ_dat_complete <- occ_dat_scaled %>% 
  mutate(comp = complete_cases) %>% 
  #unite("cell_id_year", cell_id, year, remove = FALSE) %>% 
  #mutate(cell_id_year = as.factor(cell_id_year)) %>% 
  mutate(cell_id = as.factor(cell_id)) %>% 
  filter(comp==TRUE)


#construct umf stack

nsite <- nrow(occ_dat_complete)
y <- occ_dat_complete %>% 
  dplyr::select(contains("detections")) %>% 
  as.matrix()

# Number of surveys detected per site
summary(rowSums(y, na.rm = TRUE))	

eff <- occ_dat_complete %>%
  dplyr::select(contains("cam")) %>%
  as.matrix()

bait <- occ_dat_complete %>%
  dplyr::select(contains("bait")) %>%
  as.matrix()

snare <- occ_dat_complete %>%
  dplyr::select(contains("snare")) %>%
  as.matrix()

covs_scaled <- occ_dat_complete %>% 
  select(cell_id, year, elevation:dist_all_roads_annual)

umf_stack <- unmarkedFrameOccu(y = y, 
                               siteCovs = covs_scaled,
                               obsCovs = list(eff = eff,
                                              bait = bait,
                                              snare = snare)) #, survey = surveyID))
head(umf_stack)

#########################################################################
##
##  2. Fit top model
##
##########################################################################

mod <- stan_occu(~ scale(eff) + scale(bait) + scale(snare) ~ aspect_northness + evi_annual + dist_minor_roads_annual + perc_nonveg_annual + tri + hii_annual + mtpi + distance_water + (1|cell_id), data = umf_stack, chains = 3, iter = 10000, cores = 8)

traceplot(mod, pars = c("sigma_state"), inc_warmup = TRUE)



#########################################################################
##
##  2. Project probability surface from top model
##
##########################################################################

#import_cov_stack
cov_stack_pred <- rast("/Users/tb201494/Desktop/1km_buffer/cov_stack_pred_30m_11-30-2023.tif")

#extract fixed effect beta values from fitted_model
betas <- summary(mod, submodel = "state")$mean

#ext <- ext(-2153256.05769231, -1938727.21153846, 2837586.92307692, 3152175.86538462)

#cov_stack_pred2 <- terra::crop(cov_stack_pred, ext)


#manual predictions
preds <- exp( betas[1] +
              (betas[2] * cov_stack_pred$aspect_northness) + 
              (betas[3] * cov_stack_pred$evi_annual) +
              (betas[4] * cov_stack_pred$dist_minor_roads_annual) +
              (betas[5] * cov_stack_pred$perc_nonveg_annual) +
              (betas[6] * cov_stack_pred$tri) +
              (betas[7] * cov_stack_pred$hii_annual) +
              (betas[8] * cov_stack_pred$mtpi) +
              (betas[9] * cov_stack_pred$distance_water)) /
      (1 + exp(betas[1] +
              (betas[2] * cov_stack_pred$aspect_northness) + 
              (betas[3] * cov_stack_pred$evi_annual) +
              (betas[4] * cov_stack_pred$dist_minor_roads_annual) +
              (betas[5] * cov_stack_pred$perc_nonveg_annual) +
              (betas[6] * cov_stack_pred$tri) +
              (betas[7] * cov_stack_pred$hii_annual) +
              (betas[8] * cov_stack_pred$mtpi) +
              (betas[9] * cov_stack_pred$distance_water)))

terra::plot(preds)

#~ scale(eff) + scale(bait) + scale(snare) ~ aspect_northness + evi_annual + dist_minor_roads_annual + perc_nonveg_annual + tri + hii_annual + mtpi + distance_water + (1|cell_id)
#writeRaster(preds, "occu_prob_use_surface_12-6-23.tif")


#########################################################################
##
##  2. Make Circuitscape resistance file
##
##########################################################################
# From Koen et al. 2014:
#  We found that using nodes that were randomly located around the perimeter of the buffered 
#  study area was less biased by node placement than randomly selecting nodes within the study 
#  area. We also found that a buffer of ≥ 20% of the study area width was sufficient to remove the 
#  effects of node placement on current density.

# Circuitscape wants values ranging from 1 to 100. Use stretch to increase range but maintain 
#  relative values. Plot should look exactly the same as original.

#load study area polygon
op_poly <- st_read("data/Habitat_Covariates/study_area_polys/ocp_study_area_poly_wa_only_10-30-23.shp")

tmp_rast <- rast(vals = 0, extent = ext(resistance), crs = "EPSG:5070", resolution = 30)

op_ocean_mask <- mask(tmp_rast, op_poly, updatevalue = 1)


# #split the multipolygon into disjoing polygons
# op_poly_split <- st_cast(op_poly, "POLYGON")
# 
# #select only the main polygon to get a continuous boundary
# op_poly_main <-op_poly_split[1,]

#calculate resistance from predicted probability of use surface from model
resistance <- 1 - preds
terra::plot(resistance)

# resistance <- rast("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/fitted_models/Model_Outputs/occu_resistance_surface_12-6-23.tif")

#stretch from 0 to 1 scale to 1 to 100 scale
resist_stretch <- stretch(resistance, minv = 1, maxv = 100)
plot(resist_stretch)

#set oceans as low resistance values
resist_ocean_masked <- mask(resist_stretch, op_ocean_mask, maskvalue = 1, updatevalue = 1)

#set inland water bodies as high resistance but not as complete barriers
resist_ocean_masked2 <- resist_ocean_masked
resist_ocean_masked2[is.na(resist_ocean_masked2)] <- 99

#set no data value
NAflag(resist_ocean_masked2)<- -9999


#########################################################################
##
##  2. Make Circuitscape Node file
##
##########################################################################

#around raster border

bbox <- st_bbox(resist_ocean_masked2)

bbox_coords <- matrix(c(bbox$xmin, bbox$ymin,
                        bbox$xmax, bbox$ymin,
                        bbox$xmax, bbox$ymax,
                        bbox$xmin, bbox$ymax,
                        bbox$xmin, bbox$ymin), 
                      ncol = 2, byrow = TRUE)


ext_poly <- st_polygon(list(bbox_coords))

#convert the polygon to a linestring
boundary <- st_boundary(ext_poly)
#sample points along the linestring
points <-  st_line_sample(boundary, n = 50, type = "regular") %>% st_as_sf(crs = 5070)

mapview::mapview(raster::raster(resist_ocean_masked2)) + mapview::mapview(points)


nodes <- st_coordinates(points) %>% as_tibble() %>% select(-L1) %>% mutate(ID = 1:nrow(.), .before=X)


#write_csv(nodes, "circuitscape_nodes_50_12-07-23.csv")




#following land shape
#convert the polygon to a linestring
# boundary <- st_boundary(op_poly_main)
# #sample points along the linestring
# points <-  st_line_sample(boundary, n = 100, type = "regular")
# 
# nodes <- st_coordinates(points) %>% as_tibble() %>% select(-L1) %>% mutate(ID = 1:nrow(.), .before=X)
# 

#write_csv(nodes, "circuitscape_nodes_100_12-07-23.csv")

#########################################################################
##
##  2. Extend raster bounds to make sure they cover nodes
##
##########################################################################


r_extend <- st_buffer(ext_poly, 5000) %>% st_bbox() %>% ext()

r_final <- terra::extend(resist_ocean_masked2, r_extend, fill = 1)

mapview::mapview(raster::raster(r_final)) + mapview::mapview(points)

writeRaster(r_final, "occu_resistance_surface_12-7-23.tif")
