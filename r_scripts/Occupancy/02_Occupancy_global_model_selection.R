#### Global Occupancy Model Selection ####

# Author: Read Barbee

# Date:2023-09-19

# Purpose:


################################ libraries #################################
library(tidyverse)
library(unmarked)
library(ubms)
library(camtrapR)
library(beepr)
library(doParallel)
library(DataExplorer)

#########################################################################
##
##  Specify Model Parameters
##
##########################################################################

params <- c("popdens_hii",
            "easting",
            "npp",
            "roads_hii",
            "precip",
            "dist_water",
            "perc_tree_cover",
            "infra_hii")



quad_params <- c("dist_water",
                 "roads_hii")


quad_terms <- vector()

for(i in 1:length(quad_params)){
  quad_terms[i] <- paste0("I(",quad_params[i],"^2)")
}


#########################################################################
##
## 1. Import and format step data
##
##########################################################################

occ_dat <- read_csv("data/Camera_Data/master/ocp_occ_dat_9-18-23.csv") %>% 
  mutate(across(station_id_year:year, as.factor)) %>% 
  mutate(aspect_rad = (pi*aspect)/180, .after=aspect) %>%
  mutate(northing = cos(aspect_rad),
         easting = sin(aspect_rad), .after=aspect_rad) %>% 
  rename(aspect_deg = aspect) %>% 
  select(-c(aspect_deg, aspect_rad, land_cover_usfs, land_use_usfs))

#scale covariates
occ_dat_scaled <- occ_dat %>% 
  mutate(across(tree_cover_hansen:dist_water, scale)) %>% 
  mutate(across(tree_cover_hansen:dist_water, as.numeric))


#remove station rows with missing covariate values
complete_cases <- occ_dat_scaled %>% select(tree_cover_hansen:dist_water) %>% complete.cases()

occ_dat_complete <- occ_dat_scaled %>% 
  mutate(comp = complete_cases) %>% 
  filter(comp==TRUE)


#convert to unmarked dataframe
umf <- unmarkedFrameOccu(y = occ_dat_complete %>% select(d_1:d_365),
                         siteCovs = occ_dat_complete %>% select(station_id, tree_cover_hansen:dist_water))

#########################################################################
##
## 2. Fit global model without quadratics
##
##########################################################################

#Fit global glmmTMB model: takes about ~ 15 min for 15 covs
run_global <- function (dat){
  #library(survival)  survival::clogit
  form <- as.formula(paste0("~1 ~", #remove standard intercept to be replace with stratum-based intercept
                            #fixed effects
                            paste(c(params, quad_terms), collapse = " + "), "+",
                            #random intercept (strata)
                            "(1|station_id)"))
  
  
  fit <- stan_occu(form, data=dat, chains=3, iter=1000)
  return(fit)
}


system.time(global_fit <- run_global(umf))


summary(global_fit)

#saveRDS(global_fit, "occu_global_fit_9-19-23.rds")

#global_fit <- readRDS("occu_global_fit_9-19-23.rds")

traceplot(global_fit)
