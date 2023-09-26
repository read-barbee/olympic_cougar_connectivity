#### Individual Predictive iSSF Models ####

# Author: Read Barbee

# Date:2023-07-05 

# Purpose: Develop predictive isSF models for each individual mountain lion

library(tidyverse)
library(amt)
library(MuMIn)
library(glmmTMB)

#########################################################################
##
## 1. Import and format step data
##
##########################################################################

steps <- read_csv("data/Location_Data/Steps/6h_steps_unscaled_7-03-2023.csv")

#set all negative elevations to 0
steps <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad))

#########################################################################
##
## 2. Reduce land cover and land use categories
##
##########################################################################
steps <- steps %>% mutate(land_cover_usfs_lumped = fct_collapse(land_cover_usfs, trees = c("trees", "tall_trees_shrubs", "gfh_tree_mix", "tree_shrub_mix",  "barren_tree_mix"),
                                                                shrubs = c("tall_shrubs", "shrubs", "gfh_shrub_mix", "barren_shrub_mix"),
                                                                gfh = c("gfh", "barren_gfh_mix"),
                                                                barren = c("barren_impervious", "snow_ice")),
                          land_use_usfs_lumped = fct_collapse(land_use_usfs, agriculture = c("agriculture", "rangeland_pasture")), .after = land_cover_usfs)

#plot_bar(steps %>% select(gpp:calving_season))

#########################################################################
##
## 3. Scale continuous covariates for model comparison and correlation analysis
##
##########################################################################

steps_scaled <- steps %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric))

#plot_histogram(steps_scaled %>% select(sl_, ta_, gpp:power_hii))

#plot_boxplot(steps_scaled %>% select(case_, sl_, ta_, gpp:power_hii), by = "case_")

#########################################################################
##
## 4. Dredging
##
##########################################################################

steps_scaled2 <- steps_scaled %>% select(-c(land_cover_usfs, land_use_usfs)) %>% na.omit()

#%>% dummify()

al_dat <- steps_scaled2 %>% filter(animal_id == "Al") %>% na.omit()

library(survival)

#global model--has convergence issues
full_mod_issf <-clogit(formula = case_ ~ 
                       #prey 
                       npp + dist_water +
                       #cover
                       perc_tree_cover + #I(land_cover_usfs_lumped) +
                       #terrain
                       tpi +
                       #humans
                       popdens_hii + infra_hii + roads_hii + landuse_hii +
                       #seasons
                       season +
                       #steps
                       sl_ + log(sl_) + cos(ta_) +
                       #interactions
                       sl_:tpi +
                       sl_:perc_tree_cover +
                       sl_:infra_hii +
                       sl_:season +
                       #strata
                       strata(step_id_),
                       data = steps_scaled2,
                       na.action = "na.fail")

#works for clogit but not for fit_issf
mod_sel_test <- MASS::stepAIC(full_mod_issf, direction = "backward")



full_mod_muff <-glmmTMB(case_ ~ 
                          #prey 
                          npp + dist_water +
                          #cover
                          perc_tree_cover + 
                          #I(land_cover_usfs_lumped) +
                          #terrain
                          tpi +
                          #humans
                          popdens_hii + infra_hii + roads_hii + landuse_hii +
                          #seasons
                          #season +
                          #quadratics
                          I(npp^2) + I(dist_water^2) + I(perc_tree_cover^2) + I(tpi^2) + I(roads_hii^2) +
                          #random intercept for individual
                          (1|animal_id),
                          family =poisson, 
                          data = steps_scaled2, 
                          doFit = FALSE,
                          na.action = "na.fail")

full_mod_muff$parameters$theta[1] <-log(1e3)
full_mod_muff$mapArg <-list(theta=factor(1))
fit <- glmmTMB::fitTMB(full_mod_muff)

# library(future)
# plan(multisession)
# options(future.workers = 8)
sjPlot::plot_model(fit)

sjPlot::plot_model(fit, type="re",
                   transform=NULL)

#not working
dredge_full <- MuMIn::dredge(fit, rank=BIC, fixed = "cond((Int))")




###### Muff round 2

#global model--has convergence issues
full_mod <- glmmTMB(case_ ~ gpp + npp + ndvi + evi + tree_cover_hansen + perc_tree_cover + perc_nontree_veg + perc_nonveg + land_cover_usfs_lumped + land_use_usfs_lumped + precip + dist_water + elevation + slope + northing + easting + tri + tpi + roads_hii + popdens_hii + landuse_hii + infra_hii + rails_hii + power_hii + season + hunting_season + calving_season + (1|animal_id), family = poisson, data = steps_scaled2, doFit = FALSE, na.action = "na.fail")


full_mod$parameters$theta[1] <-log(1e3)
full_mod$mapArg <-list(theta=factor(1))
full_mod_fit <- glmmTMB::fitTMB(full_mod)


dredge_full <- dredge(full_mod_fit, rank=BIC) #not currently working




