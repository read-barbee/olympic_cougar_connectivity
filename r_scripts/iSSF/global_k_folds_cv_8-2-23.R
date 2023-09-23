#### Title ####

# Author: Read Barbee

# Date:2023-08-02 

# Purpose:


################################ libraries #################################
library(tidyverse)
library(glmmTMB)
library(beepr)

#load CV function
source("r_scripts/cv_function_glmm_tmb_8-2-23.R")

#########################################################################
##
## 1. Import and format step data
##
##########################################################################
#no imputation
#steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv")

#imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_imputed_7-12-2023.csv")

#set all negative elevations to 0 and filter dispersal tracks for post dispersal event
steps <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad)) %>% 
  filter(is.na(dispersing) | dispersing==TRUE) %>% 
  select(-c(disp_qual:dispersing))

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


#########################################################################
##
## 2. Import covariate data
##
##########################################################################

cov_stack <- terra::rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2.tif")

cov_stack <- terra::aggregate(cov_stack, fact=10, cores=5)

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

cov_stack_scaled <- terra::scale(cov_stack)

#########################################################################
##
## 3. Fit global model
##
##########################################################################

#steps_scaled_no_na <- steps_scaled %>% na.omit()

form <- as.formula(case_ ~ -1 +
                     #fixed effects
                     tree_cover_hansen +
                     ndvi + #I(ndvi^2) +
                     elevation +
                     dist_water +
                     popdens_hii +
                     roads_hii +
                     #stratum-based intercept
                     (1|step_id_) +
                     #random slopes
                     (0 + tree_cover_hansen | animal_id) +
                     (0 + ndvi | animal_id) +
                     (0 + elevation| animal_id) +
                     (0 + dist_water | animal_id) +
                     (0 + popdens_hii | animal_id) +
                     (0 + roads_hii | animal_id))


#Fit global glmmTMB model: takes about 1.7 hours
global <- glmmTMB(form,
                  family=poisson,
                  data = steps_scaled,
                  doFit=FALSE); 
global$parameters$theta[1] <- log(1e3)
global$mapArg <- list(theta=factor(c(NA, 1:6)))

system.time(fit <- fitTMB(global)) #; beep("fanfare")

plot <- sjPlot::plot_model(fit, show.values = TRUE)

plotly::ggplotly(plot)

#########################################################################
##
## 3. K-folds cross_validation
##
##########################################################################

#takes about 30 minutes for 4 folds and 10 random slopes
system.time(k_folds <- k_fold_tmb(dat = steps_scaled,
                      form = form,
                      cov_stack = cov_stack_scaled,
                      n_folds = 4))

cv_scores <- k_folds[1:4] %>% unlist()

mean_score <- mean(cv_scores)


#form1: -0.7797619

#form2: 0.5333333

#form3: 0.3454545

#form4: 0.6030303





# form1 <- as.formula(case_ ~ -1 +
#                      #fixed effects
#                      tree_cover_hansen +
#                      gpp +
#                      northing +
#                      easting +
#                      perc_tree_cover +
#                      landuse_hii +
#                      ndvi +
#                      popdens_hii +
#                      rails_hii +
#                      infra_hii +
#                      #stratum-based intercept
#                      (1|step_id_) +
#                      #random slopes
#                      (0 + tree_cover_hansen | animal_id) +
#                      (0 + gpp | animal_id) +
#                      (0 + northing| animal_id) +
#                      (0 + easting | animal_id) +
#                      (0 + perc_tree_cover | animal_id) +
#                      (0 + landuse_hii | animal_id) +
#                      (0 + ndvi| animal_id) +
#                      (0 + popdens_hii | animal_id) +
#                      (0 + rails_hii | animal_id) +
#                      (0 + infra_hii | animal_id))



# form2 <- as.formula(case_ ~ -1 +
#                      #fixed effects
#                      tree_cover_hansen +
#                      ndvi +
#                      northing +
#                      easting +
#                      popdens_hii +
#                      slope +
#                      #stratum-based intercept
#                      (1|step_id_) +
#                      #random slopes
#                      (0 + tree_cover_hansen | animal_id) +
#                      (0 + ndvi | animal_id) +
#                      (0 + northing| animal_id) +
#                      (0 + easting | animal_id) +
#                      (0 + popdens_hii | animal_id) +
#                      (0 + slope | animal_id))


# form3 <- as.formula(case_ ~ -1 +
#                      #fixed effects
#                      tree_cover_hansen +
#                      ndvi + I(ndvi^2) +
#                      popdens_hii + 
#                      slope +
#                      tpi +
#                      #stratum-based intercept
#                      (1|step_id_) +
#                      #random slopes
#                      (0 + tree_cover_hansen | animal_id) +
#                      (0 + ndvi | animal_id) +
#                      (0 + popdens_hii | animal_id) +
#                      (0 + slope | animal_id) +
#                      (0 + tpi | animal_id))

#form4 <- as.formula(case_ ~ -1 +
#fixed effects
# tree_cover_hansen +
#   ndvi + #I(ndvi^2) +
#   elevation +
#   dist_water +
#   popdens_hii +
#   roads_hii +
#   #stratum-based intercept
#   (1|step_id_) +
#   #random slopes
#   (0 + tree_cover_hansen | animal_id) +
#   (0 + ndvi | animal_id) +
#   (0 + elevation| animal_id) +
#   (0 + dist_water | animal_id) +
#   (0 + popdens_hii | animal_id) +
#   (0 + roads_hii | animal_id))
