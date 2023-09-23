
library(tidyverse)
library(terra)
library(sf)
library(glmmTMB)
library(GGally)

#########################################################################
##
## 1. Import and format step data
##
##########################################################################
#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv")


#set all negative elevations to 0 and filter dispersal tracks for post dispersal event
steps <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad)) %>% 
  filter(is.na(dispersing) | dispersing==TRUE)

#collapse landcover categories
steps <- steps %>% mutate(land_cover_usfs_lumped = fct_collapse(land_cover_usfs, trees = c("trees", "tall_trees_shrubs", "gfh_tree_mix", "tree_shrub_mix",  "barren_tree_mix"),
                                                                shrubs = c("tall_shrubs", "shrubs", "gfh_shrub_mix", "barren_shrub_mix"),
                                                                gfh = c("gfh", "barren_gfh_mix"),
                                                                barren = c("barren_impervious", "snow_ice")),
                          land_use_usfs_lumped = fct_collapse(land_use_usfs, agriculture = c("agriculture", "rangeland_pasture")), .after = land_cover_usfs)

#scale and dummify covariates
steps_scaled <- steps %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric)) %>% 
  select(-c(land_cover_usfs, land_use_usfs)) %>% 
  DataExplorer::dummify(select=c("land_use_usfs_lumped", "land_cover_usfs_lumped", "season", "hunting_season", "calving_season"))

cov_stack <- terra::rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2.tif")

names(cov_stack) <- c("tree_cover_hansen",
                      "gpp",
                      "infra_hii",
                      "landuse_hii",
                      "land_cover_usfs",
                      "land_use_usfs",
                      "npp",
                      "popdens_hii",
                      "power_hii",
                      "precip",
                      "rails_hii",
                      "roads_hii",
                      "elevation",
                      "slope",
                      "aspect",
                      "tri",
                      "tpi",
                      "perc_tree_cover",
                      "perc_nontree_veg",
                      "perc_nonveg",
                      "ndvi",
                      "evi",
                      "dist_water")

cov_stack<- terra::aggregate(cov_stack, fact=10, fun="mean", cores=5)

cov_stack$aspect_rad <- (pi*cov_stack$aspect)/180
cov_stack$northing <- cos(cov_stack$aspect_rad)
cov_stack$easting <- sin(cov_stack$aspect_rad)

cov_stack2 <- terra::scale(cov_stack)


global3 <- glmmTMB(case_ ~ -1 +
                     #fixed effects
                     tree_cover_hansen +
                     #gpp +
                     northing +
                     easting + 
                     perc_tree_cover +
                     npp +
                     landuse_hii + 
                     ndvi +  I(ndvi^2) +
                     popdens_hii + I(popdens_hii^2) +
                     rails_hii +
                     infra_hii +
                     #stratum-based intercept
                     (1|step_id_) +
                     #random slopes
                     (0 +  tree_cover_hansen | animal_id) +
                     (0 +  northing | animal_id) +
                     (0 +  easting | animal_id) +
                     (0 +  perc_tree_cover | animal_id) +
                     (0 +  npp | animal_id) +
                     (0 +  landuse_hii | animal_id) +
                     (0 +  ndvi | animal_id) +
                     (0 +  popdens_hii | animal_id) +
                     (0 +  rails_hii | animal_id) +
                     (0 +  infra_hii | animal_id),
                   family=poisson,
                   data = steps_scaled,
                   doFit=FALSE); 
global3$parameters$theta[1] <- log(1e3)
global3$mapArg <- list(theta=factor(c(NA, 1:10)))

system.time(fit <- fitTMB(global3))#; beep("fanfare")




#viridis magma
ggplot()+
  tidyterra::geom_spatraster(data=binned, mapping=aes()) +
  #geom_sf(data=red_deer_used, aes(geometry=geometry))+
  #coord_sf(datum = st_crs(5070))+
  scale_fill_manual(values = rev(viridis::magma(10)), na.value = NA, guide = guide_legend(reverse = TRUE), na.translate=FALSE)  #rev() to reverse pallete

mapview::mapview(raster::raster(binned))

sjPlot::plot_model(fit)

sjPlot::plot_model(fit, type="re", terms= "ndvi")







