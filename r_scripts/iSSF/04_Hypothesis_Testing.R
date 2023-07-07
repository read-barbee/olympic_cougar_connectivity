#### iSSF Hypothesis Testing ####

# Author: Read Barbee

# Date:2023-07-07 

# Purpose: Examine log-RSS responses of each individual to each covariate and look for differences between sexes and dispersal groups


################################ Libraries #################################
library(tidyverse)
library(terra)
library(amt)
library(janitor)
library(GGally)

#########################################################################
##
## 0. Helper Functions
##
##########################################################################
run_global <- function (dat){
  library(survival)
  mod <- survival::clogit(formula = case_ ~ 
                                        #prey 
                                        evi + I(evi^2) +
                                        dist_water + I(dist_water^2) +
                                        #cover
                                        perc_tree_cover + I(perc_tree_cover^2) +
                                        #terrain
                                        tpi + I(tpi^2) +
                                        #humans
                                        popdens_hii + I(popdens_hii^2) +
                                        roads_hii + I(roads_hii^2) +
                                        infra_hii + I(infra_hii^2) +
                                        #steps
                                        sl_ + log(sl_) + cos(ta_) +
                                        #strata
                                        strata(step_id_),
                                      data = dat,
                                      na.action = "na.omit")
  return(mod)
}

set_df_s1 <- function ()
{
  s1 <- data.frame(sl_ = 100, log_sl_ = log(100), cos_ta_ = 1,
                   ruggedness_e = seq(from = -2, to = 2, length.out = 200),
                   d2forestedge_e = 0, densforestedge_e = 0, ndvi_e = 0,
                   densriparian_e = 0, densbuildings_e = 0, d2core_e = 0,
                   ruggedness_s = 0, d2forestedge_s = 0, densforestedge_s = 0,
                   ndvi_s = 0, densriparian_s = 0, densbuildings_s = 0,
                   d2core_s = 0, `(sl_ + log_sl_ + cos_ta_)` = 0)
  s1$ruggedness_e <- 0
  return(s1)
}

#function to calculate log_rss object for elevation for each individual incorporating quadratic terms, sl and ta
l_rss <- function(dat, indiv, curr_param){
  indiv_dat <- dat %>% 
    na.omit() %>% 
    filter(animal_id == indiv) %>% 
    unnest(cols=c(steps))
  
  #data frame varying elevation from min value to max value encountered by Al, holding all other covariates at the mean
  s1 <- data.frame(
    evi <-seq(from = -2, to =2, length.out = 200),
    
    evi2 <-mean(indiv_dat$evi^2),
    
    dist_water <- mean(indiv_dat$dist_water),
    
    dist_water2 <-mean(indiv_dat$dist_water^2),
    
    perc_tree_cover <- mean(indiv_dat$perc_tree_cover),
    
    perc_tree_cover2 <-mean(indiv_dat$perc_tree_cover^2),
    
    tpi <- mean(indiv_dat$tpi),
    
    tpi2 <-mean(indiv_dat$tpi^2),
    
    popdens_hii <- mean(indiv_dat$popdens_hii),
    
    popdens_hii2 <-mean(indiv_dat$popdens_hii^2),
    
    roads_hii <- mean(indiv_dat$roads_hii),
    
    roads_hii2 <-mean(indiv_dat$roads_hii^2),
    
    infra_hii <-mean(indiv_dat$infra_hii),
    
    infra_hii2 <-mean(indiv_dat$infra_hii^2),
    
    sl <- mean(indiv_dat$sl_),
    
    log_sl <- mean(indiv_dat$log_sl_),
    
    cos_ta <- mean(indiv_dat$cos_ta_)
  ) %>% 
    rename(evi = 1,
           evi2 = 2,
           dist_water = 3,
           dist_water2 = 4,
           perc_tree_cover = 5,
           perc_tree_cover2 = 6,
           tpi = 7,
           tpi2 = 8,
           popdens_hii = 9,
           popdens_hii2 = 10,
           roads_hii = 11,
           roads_hii2 = 12,
           infra_hii = 13,
           infra_hii2 = 14,
           sl= 15,
           log_sl = 16,
           cos_ta = 17
    )
  
  s1$evi <- mean(indiv_dat$evi)
  
  #### Continue Editing Function Here #######
  
  #data frame with means of all covariates encountered by Al
  s2 <- data.frame(
    elev_end <-mean(indiv_dat$elev_end),
    
    elev_end2 <-mean(indiv_dat$elev_end^2),
    
    ndvi_end <- mean(indiv_dat$ndvi_end),
    
    ndvi_end2 <-mean(indiv_dat$ndvi_end^2),
    
    dist_water_end <- mean(indiv_dat$dist_water_end),
    
    dist_water_end2 <-mean(indiv_dat$dist_water_end^2),
    
    roads_hii_end <- mean(indiv_dat$roads_hii_end),
    
    roads_hii_end2 <-mean(indiv_dat$roads_hii_end^2),
    
    forest_end <- mean(indiv_dat$forest_end),
    
    forest_end2 <-mean(indiv_dat$forest_end^2),
    
    landuse_hii_end <- mean(indiv_dat$landuse_hii_end),
    
    landuse_hii_end2 <-mean(indiv_dat$landuse_hii_end^2),
    
    sl <- mean(indiv_dat$sl_),
    
    log_sl <- mean(log(indiv_dat$sl_)),
    
    cos_ta <- mean(cos(indiv_dat$ta_))
  ) %>% 
    rename(elev_end = 1,
           elev_end2 = 2,
           ndvi_end = 3,
           ndvi_end2 = 4,
           dist_water_end = 5,
           dist_water_end2 = 6,
           roads_hii_end = 7,
           roads_hii_end2 = 8,
           forest_end = 9,
           forest_end2 = 10,
           landuse_hii_end = 11,
           landuse_hii_end2 = 12,
           sl_= 13,
           log_sl = 14,
           ta_ = 15
    )
  
  if (curr_param == "elev_end"){
    s1$elev_end <- seq(from = min(indiv_dat$elev_end, na.rm=T), to = max(indiv_dat$elev_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "elev_end2"){
    s1$elev_end2 <- seq(from = min((indiv_dat$elev_end)^2, na.rm=T), to = max((indiv_dat$elev_end)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "ndvi_end"){
    s1$ndvi_end <- seq(from = min(indiv_dat$ndvi_end, na.rm=T), to = max(indiv_dat$ndvi_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "ndvi_end2"){
    s1$ndvi_end2 <- seq(from = min((indiv_dat$ndvi_end)^2, na.rm=T), to = max((indiv_dat$ndvi_end)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "forest_end"){
    s1$forest_end <- seq(from = min(indiv_dat$forest_end, na.rm=T), to = max(indiv_dat$forest_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "forest_end2"){
    s1$forest_end2 <- seq(from = min((indiv_dat$forest_end)^2, na.rm=T), to = max((indiv_dat$forest_end)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "dist_water_end"){
    s1$dist_water_end <- seq(from = min(indiv_dat$dist_water_end, na.rm=T), to = max(indiv_dat$dist_water_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "dist_water_end2"){
    s1$dist_water_end2 <- seq(from = min((indiv_dat$dist_water_end)^2, na.rm=T), to = max((indiv_dat$dist_water_end)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "roads_hii_end"){
    s1$roads_hii_end <- seq(from = min(indiv_dat$roads_hii_end, na.rm=T), to = max(indiv_dat$roads_hii_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "roads_hii_end2"){
    s1$roads_hii_end2 <- seq(from = min((indiv_dat$roads_hii_end)^2, na.rm=T), to = max((indiv_dat$roads_hii_end)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "landuse_hii_end"){
    s1$landuse_hii_end <- seq(from = min(indiv_dat$landuse_hii_end, na.rm=T), to = max(indiv_dat$landuse_hii_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "landuse_hii_end2"){
    s1$landuse_hii_end2 <- seq(from = min((indiv_dat$landuse_hii_end)^2, na.rm=T), to = max((indiv_dat$landuse_hii_end)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "sl_"){
    s1$sl <- seq(from = min(indiv_dat$sl_, na.rm=T), to = max(indiv_dat$sl_, na.rm=T), length.out = 200)
  }
  if (curr_param == "log_sl"){
    s1$log_sl <- seq(from = min(indiv_dat$log_sl, na.rm=T), to = max(indiv_dat$log_sl, na.rm=T), length.out = 200)
  }
  if (curr_param == "ta_"){
    s1$ta_ <- seq(from = min(indiv_dat$ta_, na.rm=T), to = max(indiv_dat$ta_, na.rm=T), length.out = 200)
  }
  
  indiv_dat_nested <- dat %>% 
    filter(animal_id == indiv)
  
  mod <- indiv_dat_nested$fit2[[1]]
  
  ### Working. variable names have to be the same across all data frames and model
  l_rss_indiv <- amt::log_rss(mod, s1, s2, ci = "se", ci_level = 0.95)
  
  return(l_rss_indiv$df)
}

#########################################################################
##
## 1. Import and format step data
##
##########################################################################

steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_7-07-2023.csv")

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
## 3. Scale continuous covariates for model comparison and analysis
##
##########################################################################

steps_scaled <- steps %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric)) %>% 
  nest_by(animal_id)

#plot_histogram(steps_scaled %>% select(sl_, ta_, gpp:power_hii))

#plot_boxplot(steps_scaled %>% select(case_, sl_, ta_, gpp:power_hii), by = "case_")

#########################################################################
##
## 1. Fit global iSSF model to each individual
##
##########################################################################


global_fits <- steps_scaled %>%  
  pull(data) %>% 
  map(run_global)

steps_scaled$global_fit <- global_fits


#########################################################################
##
## 4. Calculate log-RSS and classify individual responses to each covariate
##
##########################################################################















