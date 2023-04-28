#### OCP iSSF Module_03: Explore Inidivual Variation in Selection ####

# Author: Read Barbee

# Date:2023-04-28 

# Inputs:
#   •	GPS2: Unscaled amt data frame for fitting iSSFs
# 
# Outputs:
#   •	Figures
# •	Fitted Models for each individual
# 
# Steps
# •	Import scaled and unscaled amt data frame for fitting iSSFs
# •	Fit global iSSF to each individual
# •	Explore summary stats and figures
# •	Calculate log-RSS and CI for range of each covariate for each individual
# •	Compare proportions of individuals overall and in each demographic category with positive, negative or neutral response to each covariate



################################ Libraries #################################
library(tidyverse)
library(terra)
library(lubridate)
library(amt)
library(janitor)

################################ User-Defined Parameters #################################


############################### Import amt-formatted Data #################################
# Import SSF data for 100 individual mountain lions (November 2022)
steps_raw <- read_csv("6h_steps_cov_4-28-2023.csv")



################################ Prepare data to fit models #################################

#nest and scale steps and remove locations with missing data
steps_scaled_nested <- steps_raw %>% 
  mutate(across(elev_start:landuse_hii_end, scale)) %>% 
  mutate(across(elev_start:landuse_hii_end, as.numeric))%>% 
  na.omit() %>% 
  nest_by(animal_id, sex, dispersal_status) %>% 
  rename(steps=data)

#unscaled model data
steps_unscaled_nested <- steps_raw %>% 
  na.omit() %>% 
  nest_by(animal_id, sex, dispersal_status) %>% 
  rename(steps=data)

################################ Fit iSSF to each individual #################################

#the piping workflow isn't working for some reason
# indiv_issfs_global <- steps_scaled_nested %>% dplyr::mutate(fit = map(steps, ~ amt::fit_issf(., case_ ~ elev_end + ndvi_end + dist_water_end + roads_hii_end + forest_end + landuse_hii_end +  strata(unique_step))))

steps_scaled_nested$fit <-  map(steps_scaled_nested$steps, ~ amt::fit_issf(., case_ ~ elev_end + ndvi_end + dist_water_end + roads_hii_end + forest_end + landuse_hii_end  + strata(step_id_), model = TRUE))

#+ sl_ + log(sl_) --not converging for some reason

#inspect fitted object
steps_scaled_nested$fit

#inspect model fit for first individual
steps_scaled_nested$fit[[1]]$model



################################ Calculate log_RSS for one individual #################################

# log_rss(amt_steps_scaled$fit[[1]], x1 = )


al <- steps_raw %>% 
  na.omit() %>% 
  filter(animal_id =="Al")

al2 <- steps_scaled %>% 
  na.omit() %>% 
  filter(animal_id =="Al")

#data frame varying elevation from min value to max value encountered by Al, holding all other covariates at the mean
s1 <- data.frame(
elev_end <- seq(from = min(al2$elev_end), to = max(al2$elev_end), length.out = 200),

ndvi_end <- mean(al2$ndvi_end),

dist_water_end <- mean(al2$dist_water_end),

roads_hii_end <- mean(al2$roads_hii_end),

forest_end <- mean(al2$forest_end),

landuse_hii_end <- mean(al2$landuse_hii_end)

#sl <- mean(al$sl_),
  
#log_sl <- log(sl)
)

#data frame with means of all covariates encountered by Al
s2 <- data.frame(
  elev_end <- mean(al2$elev_end),
  
  ndvi_end <- mean(al2$ndvi_end),
  
  dist_water_end <- mean(al2$dist_water_end),
  
  roads_hii_end <- mean(al2$roads_hii_end),
  
  forest_end <- mean(al2$forest_end),
  
  landuse_hii_end <- mean(al2$landuse_hii_end)
  
  #sl <- mean(al$sl_),
  
  #log_sl <- log(sl)
)


### Not working. RESUME HERE
l_rss_al <- amt::log_rss(steps_scaled_nested$fit[[1]], s1, s2)


