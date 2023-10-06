#### Resident - Disperser Global Model Validation ####

# Author: Read Barbee

# Date:2023-09-26 

# Purpose:


################################ libraries #################################
library(tidyverse)
library(glmmTMB)
library(beepr)
library(terra)
library(sf)

#load CV function
source("r_scripts/02_iSSF/08_cross_validation/cv_function_glmm_tmb_8-2-23.R")


########################################################################
##
## 1. Import location data and top model fit
##
##########################################################################

# #import top model fit
# top_mod <- readRDS("fitted_models/muff_top_global_fit_res_9-19-23.rds")


#import screened locations and split into residents and dispersers
locs <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_10-02-2023.csv")

#import steps for model predictions
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv") 

residents <- locs %>% 
  filter(dispersal_status =="resident") %>% 
  select(-c(disp_date_nsd:dispersing))


dispersers <- locs %>% 
  filter(dispersal_status =="disperser") %>% 
  filter(dispersing == TRUE)


#########################################################################
##
## 2. Import covariate data
##
##########################################################################

cov_stack <- terra::rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2_asp.tif")

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
                      "dist_water",
                      "northing",
                      "easting")


cov_stack_scaled <- terra::scale(cov_stack)

#########################################################################
##
## 3. Fit global model
##
##########################################################################

#steps_scaled_no_na <- steps_scaled %>% na.omit()

# form <- as.formula(case_ ~ -1 +
#                      #fixed effects
#                      tree_cover_hansen +
#                      ndvi + #I(ndvi^2) +
#                      tpi +
#                      northing +
#                      tri +
#                      perc_tree_cover +
#                      roads_hii +
#                      popdens_hii +
#                      landuse_hii +
#                      infra_hii +
#                      #stratum-based intercept
#                      (1|step_id_) +
#                      #random slopes
#                      (0 + tree_cover_hansen | animal_id) +
#                      (0 + ndvi | animal_id) +
#                      (0 + tpi| animal_id) +
#                      (0 + northing | animal_id) +
#                      (0 + tri | animal_id) +
#                      (0 + perc_tree_cover | animal_id) +
#                      (0 + roads_hii| animal_id) +
#                      (0 + popdens_hii | animal_id) +
#                      (0 + landuse_hii | animal_id) +
#                      (0 + infra_hii | animal_id))


#Fit global glmmTMB model: takes about 10 min for 10 linear terms
# global <- glmmTMB(form,
#                   family=poisson,
#                   data = steps_scaled,
#                   doFit=FALSE); 
# global$parameters$theta[1] <- log(1e3)
# global$mapArg <- list(theta=factor(c(NA, 1:10)))
# 
# system.time(global_fit <- fitTMB(global)) ; beep("fanfare")

#saveRDS(global_fit, "fitted_models/global_fit_residents_linear_tmb_10-02-2023.rds")

#global_fit <- readRDS("fitted_models/global_fit_residents_linear_tmb_10-02-2023.rds")

#plot <- sjPlot::plot_model(fit, show.values = TRUE, transform = NULL)

#plotly::ggplotly(plot)


########################################################################
##
## 3. Make predictions from top model
##
##########################################################################

betas <- fixef(global_fit)$cond

reff <- ranef(global_fit)

# #clean up names of quadratic terms
# new_names <- vector()
# for(i in 1:length(betas)){
#   new_names[i] <- names(betas[i]) %>% 
#     str_replace(coll("^2)"), "2") %>% 
#     str_remove(coll("I("))
# }

#names(betas) <- new_names

#add covariate layers for quadratic terms
# cov_stack2$roads_hii2 <- cov_stack2$roads_hii^2
# cov_stack2$perc_tree_cover2 <- cov_stack2$perc_tree_cover^2
# cov_stack2$ndvi2 <- cov_stack2$ndvi^2

#reorder raster layers to match order of terms in model
cov_stack_pred <- cov_stack_scaled[[names(betas)]]
cov_stack_pred$step_id_ = sample(steps$step_id_, ncell(cov_stack_pred), replace = TRUE)
cov_stack_pred$animal_id = sample(steps$animal_id, ncell(cov_stack_pred), replace = TRUE)




#manual predictions
preds <- exp((betas[1] * cov_stack_pred$tree_cover_hansen) + 
               (betas[2] * cov_stack_pred$ndvi) + 
               (betas[3] * cov_stack_pred$tpi) + 
               (betas[4] * cov_stack_pred$northing) + 
               (betas[5] * cov_stack_pred$tri) + 
               (betas[6] * cov_stack_pred$perc_tree_cover) +
               (betas[7] * cov_stack_pred$roads_hii) + 
               (betas[8] * cov_stack_pred$popdens_hii) +
               (betas[9] * cov_stack_pred$landuse_hii) + 
               (betas[10] * cov_stack_pred$infra_hii))/
  (1 + exp((betas[1] * cov_stack_pred$tree_cover_hansen) + 
             (betas[2] * cov_stack_pred$ndvi) + 
             (betas[3] * cov_stack_pred$tpi) + 
             (betas[4] * cov_stack_pred$northing) + 
             (betas[5] * cov_stack_pred$tri) + 
             (betas[6] * cov_stack_pred$perc_tree_cover) +
             (betas[7] * cov_stack_pred$roads_hii) + 
             (betas[8] * cov_stack_pred$popdens_hii) +
             (betas[9] * cov_stack_pred$landuse_hii) + 
             (betas[10] * cov_stack_pred$infra_hii)))


plot(preds)



########################################################################
##
## 3. Bin predictions
##
##########################################################################

pred_vals <- terra::values(preds)

breaks <- quantile(pred_vals, probs = 0:10/10, na.rm = T)

bins <- cut(pred_vals, breaks, include.lowest=TRUE)

binned <- terra::classify(preds, rcl=as.vector(breaks))

plot(binned)

levels(binned) <- 1:10


########################################################################
##
## 3. Internal validation (residents)
##
##########################################################################
eval_res <- residents %>% st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070)

loc_preds<- terra::extract(binned, eval_res)[,2] %>% as.numeric() %>% na.omit() 


#plot proportion of test points in each bin. Currently storing all of the same plot for some reason.
ggplot() +
  geom_bar(aes(x=loc_preds, y = after_stat(prop))) 

n_bins <- length(unique(loc_preds))

#calculate pearson correlation between the bin number and the proportion of used locations in each bin
cor_cv_res <- cor.test(x = 1:n_bins, y = as.numeric(table(loc_preds))/(length(loc_preds)/n_bins), method = "spearman", exact = FALSE)$estimate


########################################################################
##
## 3. External validation (dispersers)
##
##########################################################################
eval_disp <- dispersers %>% st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070)

loc_preds<- terra::extract(binned, eval_disp)[,2] %>% as.numeric() %>% na.omit() 


#plot proportion of test points in each bin. Currently storing all of the same plot for some reason.
ggplot() +
  geom_bar(aes(x=loc_preds, y = after_stat(prop))) 

n_bins <- length(unique(loc_preds))

#calculate pearson correlation between the bin number and the proportion of used locations in each bin
cor_cv_disp <- cor.test(x = 1:n_bins, y = as.numeric(table(loc_preds))/(length(loc_preds)/n_bins), method = "spearman", exact = FALSE)$estimate


