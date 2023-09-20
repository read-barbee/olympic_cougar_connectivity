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


#import top model fit
top_mod <- readRDS("muff_top_global_fit_res_9-19-23.rds")

########################################################################
##
## 1. Import and format step data
##
##########################################################################

#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv") 

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
    str_remove(coll("I("))
}

names(coefs) <- names


mod <- make_issf_model(coefs = coefs,
                       sl = make_gamma_distr(shape = 1, scale = 1, vcov = NULL),
                       ta =make_vonmises_distr(kappa = 1, vcov = NULL))



start <- make_start(x= c(-2114054, 3007263),
                    time = ymd_hms("2022-04-05 05:00:35"),
                    ta = 0,
                    dt = hours(2),
                    crs = 5070)


#stops working here. RESUME HERE ########

k1 <- redistribution_kernel(mod, map = cov_stack, start = start)

# 
# k1 <- redistribution_kernel(mod, map = cov_stack, start = start,
#                             landscape = "continuous", tolerance.outside = 0.2, 
#                             n.control = 1e4)







