#### iSSF Simulation ####

# Author: Read Barbee

# Date:2023-07-26 

# Purpose:


 #install.packages("https://cran.r-project.org/src/contrib/Archive/amt/amt_0.1.7.tar.gz", repos=NULL, type="source")


library(tidyverse)
library(amt)

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
  DataExplorer::dummify(select=c("land_use_usfs_lumped", "land_cover_usfs_lumped", "season", "hunting_season", "calving_season")) #%>% 
#select(case_, animal_id, step_id_, gpp:calving_season_yes) %>% 

al_dat <- steps %>% filter(animal_id=="Al")
al_dat_scaled <- steps_scaled %>% filter(animal_id=="Al")

# g_fit <- fit_issf(case_ ~
#                     ndvi + I(ndvi^2) +
#                     northing +
#                     easting +
#                     land_cover_usfs_lumped_shrubs +
#                     land_cover_usfs_lumped_trees +
#                     land_cover_usfs_lumped_water +
#                     tree_cover_hansen +
#                     popdens_hii + I(popdens_hii^2) +
#                     slope +
#                     sl_ +
#                     log_sl_ +
#                     cos_ta_ +
#                     strata(step_id_), data=al_dat_scaled, model=TRUE
#                     )

g_fit2 <- fit_issf(case_ ~
                    tree_cover_hansen +
                    gpp +
                    northing +
                    easting +
                    perc_nontree_veg +
                    tpi +
                    popdens_hii + #I(popdens_hii^2) +
                    landuse_hii +
                    ndvi + #I(ndvi^2) +
                    infra_hii +
                    sl_ +
                    log_sl_ +
                    cos_ta_ +
                    strata(step_id_), data=al_dat_scaled, model=TRUE
)


library(survival)
g_fit3 <- clogit(case_ ~
                     tree_cover_hansen +
                     gpp +
                     northing +
                     easting +
                     perc_nontree_veg +
                     tpi +
                     popdens_hii + #I(popdens_hii^2) +
                     landuse_hii +
                     ndvi + #I(ndvi^2) +
                     infra_hii +
                     sl_ +
                     log_sl_ +
                     cos_ta_ +
                     strata(step_id_), data=al_dat_scaled, model=TRUE
)


cov_stack <- raster::brick("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2.tif", crs="+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs")

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

cov_stack$aspect_rad <- (pi*cov_stack$aspect)/180
cov_stack$northing <- cos(cov_stack$aspect_rad)
cov_stack$easting <- sin(cov_stack$aspect_rad)

#old version of amt (working)

al_dat <- na.omit(al_dat)

dist <- fit_distr(al_dat$sl_, dist_name = "gamma")

shape <- dist$params$shape
scale <- dist$params$scale

#cov_stack$ndvi2 = raster::scale(cov_stack[["ndvi"]]^2)


#cov_stack2 <-raster::as.raster(cov_stack)

movement_kernel <- amt::movement_kernel(shape=shape, scale=scale, template = cov_stack) 

raster::crs(movement_kernel) <- 5070

#resources layer has to be a raster not a spatraster
habitat_kernel <- amt::habitat_kernel(coef=list(tree_cover_hansen = coef(g_fit2)["tree_cover_hansen"],
                                                gpp = coef(g_fit2)["gpp"],
                                                northing = coef(g_fit2)["northing"],
                                                easting = coef(g_fit2)["easting"],
                                                perc_nontree_veg = coef(g_fit2)["perc_nontree_veg"],
                                                tpi = coef(g_fit2)["tpi"],
                                                popdens_hii = coef(g_fit2)["popdens_hii"],
                                                landuse_hii = coef(g_fit2)["landuse_hii"],
                                                ndvi = coef(g_fit2)["ndvi"],
                                                infra_hii = coef(g_fit2)["infra_hii"]),
                                      resources=cov_stack, exp=FALSE)

coefs <-  coef(g_fit2)[1:10]

hk_layers <- list()
for(i in 1:3){ #length(coefs)
  name <- names(coefs)[i]
  hk_layers[[i]] <- cov_stack[[name]]
}



test <- raster::predict(cov_stack, g_fit3)

raster::crs(habitat_kernel) <- 5070

#takes about 11 minutes for one covariate. Not working. just looks blank
system.time(sim_ud <- amt::simulate_ud(movement_kernel = movement_kernel,
                           habitat_kernel = habitat_kernel,
                           start = as.numeric(steps[1,c("x1_", "y1_")]),
                           n=1e7))

raster::plot(sim_ud)

#new version of amt (not working)
#amt::get_max_dist(g_fit2, p=0.99)

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

start <- make_start(x=c(-2046788, 2987782), ta_ = -0.04204607, dt = hours(2), time = ymd_hms("2021-03-18 14:00:00"), crs=5070)

rd <- amt::redistribution_kernel(x=g_fit2, start=start, map = cov_stack, max.dist = 2000)


