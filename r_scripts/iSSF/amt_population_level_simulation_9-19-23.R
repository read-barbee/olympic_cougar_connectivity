#### amt Population Level Simulation ####

# Author: Read Barbee

# Date:2023-09-19

# Purpose:


################################ libraries #################################
library(tidyverse)
library(glmmTMB)
library(amt)
#library(MuMIn)
#library(INLA)
library(beepr)
library(sf)
library(terra)


#import top model fit
top_mod <- readRDS("muff_top_global_fit_res_9-19-23.rds")

########################################################################
##
## 1. Import and format step data
##
##########################################################################

#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv") 

#locs <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_7-11-2023.csv")

#%>% mutate(ndvi = ndvi*0.0001)

#set all negative elevations to 0 and filter dispersal tracks for post dispersal event
steps_scaled <- steps %>% 
  filter(dispersal_status=="resident") %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad, dispersing:disp_qual)) %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric))


########################################################################
##
## 2. Import covariate data
##
##########################################################################

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

cov_stack$northing <- cos((pi*cov_stack$aspect)/180)
cov_stack$easting <- sin((pi*cov_stack$aspect)/180)
cov_stack$station_id = as.factor("LEKT_Station1")


########################################################################
##
## 2. Set starting points and make redistribution kernel
##
##########################################################################

top_estimates <- summary(top_mod)$coefficients$cond %>% as.data.frame() %>% rownames_to_column("term") %>% select(term, Estimate)

# coefs <- vector()
# for(i in 1:nrow(top_estimates)){
#   tm <- top_estimates$term[i] %>% 
#     str_replace(coll("^2)"), "2") %>% 
#     str_remove(coll("I("))
#   
#   coefs[i] <- paste0(tm, " = ", top_estimates$Estimate[i])
# }

coefs <- vector()
names <- vector()
for(i in 1:nrow(top_estimates)){
  coefs[i] <- top_estimates$Estimate[i]
  names[i] <- top_estimates$term[i] %>% 
    str_replace(coll("^2)"), "2") %>% 
    str_remove(coll("I(")) %>% 
    str_c("_end")
}



names(coefs) <- names

#linear terms only
coefs_linear <- coefs[str_detect(names, "2", negate = TRUE)] 

# cov_stack2 <- cov_stack[[names(cov_stack) %in% names(coefs_linear)]]
# coefs_linear <- coefs_linear[names(cov_stack2)]

cov_stack2 <- terra::aggregate(cov_stack, fact=10, cores=5)



#plot constructred gamma distribution
x <- seq(0.1, 10, 0.1)
plot(x, dgamma(x, shape = 2, scale = 2), type = "l")


mod <- make_issf_model(coefs = c(coefs_linear, sl_ = median(steps$sl_, na.rm=T), ta_ = median(steps$ta_, na.rm=T)),
                       sl = fit_distr(steps$sl_, "gamma", na.rm = TRUE),
                       ta =fit_distr(steps$ta_, "vonmises", na.rm = TRUE))


start <- make_start(x= c(-2114054, 3007263),
                    #time = ymd_hms("2022-04-05 05:00:35"),
                    ta = 0,
                    dt = hours(2),
                    crs = 5070)


k1 <- amt::redistribution_kernel(mod, map = cov_stack2, start = start)

s1 <- simulate_path(k1, n.steps = 1000)

terra::plot(cov_stack2$tree_cover_hansen)
lines(s1$x_, s1$y_, col = "red")

# 
# k1 <- redistribution_kernel(mod, map = cov_stack, start = start,
#                             landscape = "continuous", tolerance.outside = 0.2, 
#                             n.control = 1e4)

bbox <- terra::ext(cov_stack2)

study_area <- terra::vect(bbox)

# Ensure that start nodes are drawn from within start_zone and save matrix of the start nodes to use for each iteration:

#Import water body polygons
water_polys <- st_read("data/Habitat_Covariates/washington_water_polygons/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation.shp") %>% st_transform(crs = 5070)

#filter water body polygons to only include permanent water bodies
water_polys_filtered <- water_polys %>% filter(WB_PERIOD_ =="PER" ) #|WB_PERIOD_ =="INT"

water_polys_cropped <- water_polys_filtered %>%
  sf::st_union() %>%
  sf::st_convex_hull() %>%
  sf::st_intersection(water_polys_filtered) %>%
  st_collection_extract('POLYGON') %>%
  st_union() %>% 
  st_sf()


generate_start_nodes <- function (start_zone, barriers)
{
  x1 <- start_zone$xmin
  x2 <- start_zone$xmax
  y1 <- start_zone$ymin
  y2 <- start_zone$ymax
  
  x_rand <- runif(1, x1, x2)
  y_rand <- runif(1, y1, y2)
    
  return(cbind(x_rand, y_rand))
}

################################ Graveyard #################################
#cov_stack$roads_hii2 <- (cov_stack$roads_hii^2)


#reordering doesn't seem to make a difference
#na.omit doesn't help
#converting to raster stack doesn't help

#cov_stack3 <- na.omit(cov_stack2)

#check to make sure point is within raster area--it is
#terra::plot(cov_stack2[[1]])
#terra::points(terra::vect("POINT (-2114054 3007263)", crs="+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"))
#reducing cov_stack_resolution doesn't work either

# mod <- make_issf_model(coefs = c(coefs_linear, sl_ = median(steps$sl_), ta_ = median(steps$ta_), strata = steps$step_id_),
#                        sl = make_gamma_distr(shape = 2, scale = 2, vcov = NULL),
#                        ta =make_vonmises_distr(kappa = 4, vcov = NULL))

#cov_stack2_scaled <- terra::scale(cov_stack2)

# cov_stack2$sl <- 1
# cov_stack2$ta <- 0

#cov_stack2_rast <- raster::stack(cov_stack2)