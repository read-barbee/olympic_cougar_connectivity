#### Resident - Disperser Global Model Validation ####

# Author: Read Barbee

# Date:2023-09-26 

# Purpose:


################################ libraries #################################
library(tidyverse)
library(glmmTMB)
library(beepr)
library(terra)

#load CV function
source("r_scripts/02_iSSF/08_cross_validation/cv_function_glmm_tmb_8-2-23.R")


########################################################################
##
## 1. Import location data and top model fit
##
##########################################################################

#import top model fit
top_mod <- readRDS("fitted_models/muff_top_global_fit_res_9-19-23.rds")


#import screened locations and split into residents and dispersers
locs <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_9-26-2023.csv")

#import steps for model predictions
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv") 

residents <- locs %>% 
  filter(dispersal_status =="resident") %>% 
  select(-c(disp_date_nsd:dispersing))


dispersers <- locs %>% 
  filter(dispersal_status =="disperser") %>% 
  filter(dispersing == TRUE)


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

#add aspect transformations
cov_stack$northing <- cos((pi*cov_stack$aspect)/180)
cov_stack$easting <- sin((pi*cov_stack$aspect)/180)
cov_stack$step_id_ = sample(steps$step_id_, ncell(cov_stack), replace = TRUE)
cov_stack$animal_id = sample(steps$animal_id, ncell(cov_stack), replace = TRUE)



#resample cov_stack to lower resolution for mapping
cov_stack2 <- terra::aggregate(cov_stack, fact=10, cores=5)

########################################################################
##
## 3. Make predictions from top model
##
##########################################################################

betas <- fixef(top_mod)$cond

reff <- ranef(top_mod)

#clean up names of quadratic terms
new_names <- vector()
for(i in 1:length(betas)){
  new_names[i] <- names(betas[i]) %>% 
    str_replace(coll("^2)"), "2") %>% 
    str_remove(coll("I("))
}

names(betas) <- new_names

#add covariate layers for quadratic terms
cov_stack2$roads_hii2 <- cov_stack2$roads_hii^2
cov_stack2$perc_tree_cover2 <- cov_stack2$perc_tree_cover^2
cov_stack2$ndvi2 <- cov_stack2$ndvi^2

#reorder raster layers to match order of terms in model
cov_stack_pred <- cov_stack2[[names(betas)]]
cov_stack_pred$step_id_ = sample(steps$step_id_, ncell(cov_stack_pred), replace = TRUE)
cov_stack_pred$animal_id = sample(steps$animal_id, ncell(cov_stack_pred), replace = TRUE)



#about 1.5 min with resampled raster. Not working though. predicts all values as 1
system.time(test <- predict(top_mod, 
                 newdata = cov_stack_pred,
                 re.form = ~0,
                 allow.new.levels = TRUE,
                 type = "link"))

proj <- cov_stack2

values(proj) <- test

cov_stack2 <- scale(cov_stack2)

#manual predictions
preds <- (betas[1] * cov_stack2$roads_hii) + (betas[2] * cov_stack2$roads_hii^2) + 
  (betas[3] * cov_stack2$popdens_hii) +
  (betas[4] * cov_stack2$infra_hii) +
  (betas[5] * cov_stack2$tpi) +
  (betas[6] * cov_stack2$npp) +
  (betas[7] * cov_stack2$perc_tree_cover) + (betas[8] * cov_stack2$perc_tree_cover^2) +
  (betas[9] * cov_stack2$ndvi) + (betas[10] * cov_stack2$ndvi^2) +
  (betas[11] * cov_stack2$tree_cover_hansen) +
  (betas[12] * cov_stack2$northing) +
  (betas[13] * cov_stack2$easting) +
  (betas[14] * cov_stack2$tri) +
  (betas[15] * cov_stack2$precip)


plot(preds)






