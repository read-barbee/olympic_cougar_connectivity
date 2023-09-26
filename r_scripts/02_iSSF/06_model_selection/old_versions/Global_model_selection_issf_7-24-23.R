#### iSSF global stepwise model seleciton ####

# Author: Read Barbee

# Date:2023-07-24 

# Purpose:


################################ libraries #################################
library(tidyverse)
library(survival)
library(beepr)

#########################################################################
##
## 1. Import and format step data
##
##########################################################################
#no imputation
#steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv")

#with imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_imputed_7-12-2023.csv")

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

DataExplorer::plot_missing(steps_scaled_no_na)

steps_scaled_no_na <- steps_scaled %>% 
  select(-c(disp_qual, disp_date_nsd, dispersing)) %>% 
  na.omit()

#########################################################################
##
## 2. Fit global model and auto-dredge (not working)
##
##########################################################################

#fits in 9 seconds with the approximate method
system.time(global <- survival::clogit(case_ ~
                             ndvi + I(ndvi^2) +
                             tree_cover_hansen +
                             gpp +
                             northing +
                             land_cover_usfs_lumped_shrubs +
                             perc_tree_cover + I(perc_tree_cover^2) +
                             popdens_hii + I(popdens_hii^2) +
                             land_cover_usfs_lumped_trees +
                             perc_nontree_veg +
                             land_cover_usfs_lumped_water +
                             land_use_usfs_lumped_developed +
                             land_use_usfs_lumped_non_forest_wetland +
                             land_cover_usfs_lumped_barren +
                             landuse_hii + I(landuse_hii^2) +
                             perc_nonveg +
                             tpi +
                             easting + I(easting^2) +
                             dist_water + I(dist_water^2) +
                             slope + I(slope^2) +
                             infra_hii +
                             rails_hii +
                             roads_hii + I(roads_hii^2) +
                             sl_ +
                             log_sl_ +
                             cos_ta_ +
                             strata(step_id_), cluster = animal_id, method = "approximate", data = steps_scaled_no_na))


system.time(global2 <- survival::clogit(case_ ~
                                         tree_cover_hansen +
                                         gpp +
                                         northing +
                                         easting + 
                                         perc_nontree_veg +
                                         tpi +
                                         popdens_hii + 
                                         landuse_hii + 
                                         ndvi +
                                         infra_hii +
                                         sl_ +
                                         log_sl_ +
                                         cos_ta_ +
                                         strata(step_id_), cluster = animal_id, method = "approximate", data = steps_scaled_no_na))


# global_linear <- survival::clogit(case_ ~
#                              ndvi + 
#                              tree_cover_hansen +
#                              gpp +
#                              northing +
#                              land_cover_usfs_lumped_shrubs +
#                              perc_tree_cover + 
#                              popdens_hii + 
#                              land_cover_usfs_lumped_trees +
#                              perc_nontree_veg +
#                              land_cover_usfs_lumped_water +
#                              land_use_usfs_lumped_developed +
#                              land_use_usfs_lumped_non_forest_wetland +
#                              land_cover_usfs_lumped_barren +
#                              landuse_hii + 
#                              perc_nonveg +
#                              tpi +
#                              easting + 
#                              dist_water + 
#                              slope + 
#                              infra_hii +
#                              rails_hii +
#                              roads_hii +
#                              sl_ +
#                              log_sl_ +
#                              cos_ta_ +
#                              strata(step_id_), cluster = animal_id, method = "approximate", data = steps_scaled_no_na)

#backward stepwise selection (takes about 36 minutes)
#system.time(stepwise <- MASS::stepAIC(global2, direction= "backward"))

summary <- summary(stepwise)

step_df <- summary$coefficients %>% as.data.frame()

#write_csv(step_df, "issf_backward_stepwise_7-24-23.csv")

#backward stepwise with linear terms only---doesn't seem to eliminate anything
system.time(stepwise2 <- MASS::stepAIC(global2, direction = "backward", trace=0, scope = list(lower = ~ sl_ + log_sl_ + cos_ta_ + strata(step_id_))))



#########################################################################
##
## 2. Manual stepwise (no imp)
##
##########################################################################

#removing terms
global2 <- survival::clogit(case_ ~
                              ndvi + I(ndvi^2) +
                              tree_cover_hansen +
                              gpp +
                              northing +
                              land_cover_usfs_lumped_shrubs +
                              perc_tree_cover + I(perc_tree_cover^2) +
                              popdens_hii + I(popdens_hii^2) +
                              land_cover_usfs_lumped_trees +
                              perc_nontree_veg +
                              land_cover_usfs_lumped_water +
                              #land_use_usfs_lumped_developed +
                              #land_use_usfs_lumped_non_forest_wetland +
                              land_cover_usfs_lumped_barren +
                              #landuse_hii + I(landuse_hii^2) +
                              #perc_nonveg +
                              #tpi +
                              easting + #I(easting^2) +
                              #dist_water + I(dist_water^2) +
                              slope + #I(slope^2) +
                              #infra_hii +
                              rails_hii +
                              #roads_hii + I(roads_hii^2) +
                              sl_ +
                              log_sl_ +
                              cos_ta_ +
                              strata(step_id_), cluster = animal_id, method = "approximate", data = steps_scaled_no_na)



#########################################################################
##
## 2. Manual stepwise (imp)
##
##########################################################################

global_imp <- survival::clogit(case_ ~
                                 gpp +
                                 tree_cover_hansen +
                                 ndvi + I(ndvi^2) +
                                 land_cover_usfs_lumped_water +
                                 northing +
                                 perc_tree_cover + I(perc_tree_cover^2) +
                                 land_cover_usfs_lumped_trees +
                                 land_cover_usfs_lumped_shrubs +
                                 perc_nontree_veg + I(perc_nontree_veg^2) +
                                 popdens_hii + I(popdens_hii^2) +
                                 landuse_hii +
                                 easting + I(easting^2) +
                                 tri +
                                 land_use_usfs_lumped_developed +
                                 infra_hii +
                                 land_use_usfs_lumped_non_forest_wetland +
                                 tpi +
                                 rails_hii +
                                 roads_hii +
                                 dist_water + I(dist_water^2) +
                                 perc_nonveg + I(perc_nonveg^2) +
                                 land_cover_usfs_lumped_barren +
                                 sl_ +
                                 log_sl_ +
                                 cos_ta_ +
                                 strata(step_id_), cluster = animal_id, method = "approximate", data = steps_scaled_no_na)


#remove terms
global_imp2 <- survival::clogit(case_ ~
                                 #gpp +
                                 tree_cover_hansen +
                                 ndvi + I(ndvi^2) +
                                 land_cover_usfs_lumped_water +
                                 northing +
                                 #perc_tree_cover + I(perc_tree_cover^2) +
                                 land_cover_usfs_lumped_trees +
                                 land_cover_usfs_lumped_shrubs +
                                 perc_nontree_veg + I(perc_nontree_veg^2) +
                                 popdens_hii + I(popdens_hii^2) +
                                 #landuse_hii +
                                 #easting + I(easting^2) +
                                 tri +
                                 land_use_usfs_lumped_developed +
                                 #infra_hii +
                                 land_use_usfs_lumped_non_forest_wetland +
                                 #tpi +
                                 #rails_hii +
                                 #roads_hii +
                                 dist_water + I(dist_water^2) +
                                 perc_nonveg + I(perc_nonveg^2) +
                                 #land_cover_usfs_lumped_barren +
                                 sl_ +
                                 log_sl_ +
                                 cos_ta_ +
                                 strata(step_id_), cluster = animal_id, method = "approximate", data = steps_scaled_no_na)
