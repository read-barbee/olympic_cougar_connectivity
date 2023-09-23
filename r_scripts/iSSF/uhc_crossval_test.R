#### UHC Plots ####

# Author: Read Barbee

# Date:2023-07-25 

# Purpose:


############################### libraries #################################
library(tidyverse)
library(amt)
library(beepr)

#########################################################################
##
## 1. Import and format step data
##
##########################################################################
#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv")

#with imputation
#steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_imputed_7-12-2023.csv")

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

#ataExplorer::plot_missing(steps_scaled_no_na)

steps_scaled_no_na <- steps_scaled %>% 
  select(-c(disp_qual, disp_date_nsd, dispersing)) %>% 
  na.omit()

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
#########################################################################
##
## 2. Split data into training and test datasets
##
##########################################################################

steps_nested <- steps_scaled_no_na %>% 
  nest_by(animal_id)

#Split into train (80%) and test (20%)
set.seed(1)
steps_nested$train <- rbinom(n = nrow(steps_nested),
                      size = 1, prob = 0.8)
train <- steps_nested %>% filter(train==1) %>% unnest(cols=c(data))  #[case_t$train == 1, ]
test <- steps_nested %>%  
  filter(train==0) %>% 
  unnest(cols=c(data)) %>% 
  dplyr::select(-unique_step) %>% 
  mutate(animal_id = as.factor(animal_id),
         sex= as.factor(sex),
         dispersal_status=as.factor(dispersal_status)) %>% 
  mutate_at(vars(contains("land_")), as.factor)
  #[case_t$train == 0, ]

fit <- fit_issf(case_ ~
                  slope +
                  northing +
                  easting +
                  ndvi + I(ndvi^2) +
                  gpp +
                  perc_nontree_veg +
                  tree_cover_hansen +
                  perc_tree_cover + I(perc_tree_cover^2) +
                  land_cover_usfs_lumped_shrubs +
                  land_cover_usfs_lumped_trees +
                  land_cover_usfs_lumped_barren +
                  land_cover_usfs_lumped_water +
                  popdens_hii + I(popdens_hii^2) +
                  rails_hii +
                  sl_ +
                  log_sl_ +
                  cos_ta_ +
                  strata(step_id_), data = train, model = TRUE)


uhc_prep <- prep_uhc(object = fit, test_dat = test,
                           n_samp = 20, verbose = TRUE)


dev.new()
par(mfrow = c(4, 5))

plot(uhc_prep)
