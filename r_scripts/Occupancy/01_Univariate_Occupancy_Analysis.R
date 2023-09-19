#### Univriate Occupancy Analysis####

# Author: Read Barbee

# Date:2023-09-18

# Purpose:

#Note: need to run with more iterations


################################ libraries #################################
library(tidyverse)
library(unmarked)
library(ubms)
library(camtrapR)
library(beepr)
library(doParallel)

#########################################################################
##
##  Specify Model Parameters
##
##########################################################################

occ_covs <- c("roads_hii",
              "popdens_hii",
              "landuse_hii",
              "infra_hii",
              "rails_hii",
              "power_hii",
              "tpi",
              "npp",
              "perc_tree_cover",
              "perc_nontree_veg",
              "perc_nonveg",
              "ndvi",
              "evi",
              "tree_cover_hansen",
              "gpp",
              "northing",
              "easting",
              "slope",
              "tri",
              "elevation",
              "precip",
              "dist_water")


det_covs <- NULL


#########################################################################
##
##  1. Import stacked occupancy data
##
##########################################################################


occ_dat <- read_csv("data/Camera_Data/master/ocp_occ_dat_9-18-23.csv") %>% 
  mutate(across(station_id_year:year, as.factor)) %>% 
  mutate(aspect_rad = (pi*aspect)/180, .after=aspect) %>%
  mutate(northing = cos(aspect_rad),
         easting = sin(aspect_rad), .after=aspect_rad) %>% 
  rename(aspect_deg = aspect) %>% 
  select(-c(aspect_deg, aspect_rad))

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
##  2. Fit linear univariate stacked occupancy models
##
##########################################################################

#register doParallel backend. Note:the cores argument implies forking with doMC backend, but you can check with the getDoParName() function
doParallel::registerDoParallel(cores = 9) 


#### NULL MODEL ###

null_fit <- stan_occu(~1 ~ (1|station_id), data=umf, chains=3, iter=100)

traceplot(null_fit)

null_fit_list <- list(null_fit)

names(null_fit_list) <- "null"


#### UNIVARIATE MODELS ###

cov_names <- occ_covs

#parallelized for loop: fits in ~ 9 minutes with 3 chains and 100 iterations
uni_fits <- list()
system.time(uni_fits <- foreach (i = 1:length(cov_names), .packages=c("ubms"))%dopar%{
  
  cov <- cov_names[i]
  
  form <- as.formula(paste0("~1 ~", #remove standard intercept to be replace with stratum-based intercept
                            #fixed effects
                            cov, "+",
                            #random intercept (strata)
                            "(1|station_id)"))
  
  
  stan_occu(form, data=umf, chains=3, iter=100)
  
}); beep("fanfare")


names(uni_fits) <- cov_names


# fit_list <- fitList(c(uni_fits, null_fit))
# 
# mod_sel <- modSel(fit_list, nullmod="null")

#use loo function for leave one out cross-validation


#########################################################################
##
##  3. Fit quadratic univariate stacked occupancy models
##
##########################################################################

#parallelized for loop: fits in ~ 10 minutes with 3 chains and 100 iterations
uni_fits_quad <- list()
system.time(uni_fits_quad <- foreach (i = 1:length(cov_names), .packages=c("ubms"))%dopar%{
  
  cov <- cov_names[i]
  
  form <- as.formula(paste0("~1 ~", #remove standard intercept to be replace with stratum-based intercept
                            #fixed effects
                            cov, "+", "I(", cov,  "^2) + ",
                            #random intercept (strata)
                            "(1|station_id)"))
  
  
  stan_occu(form, data=umf, chains=3, iter=100)
  
}); beep("fanfare")

#rename qudaratic covariates with "2" at the end
quad_names <- vector()
for(i in 1:length(cov_names)){
  old_name <- cov_names[i]
  quad_names[i] <- paste0(old_name, "2")
}

names(uni_fits_quad) <- quad_names


#########################################################################
##
##  4. Create model summary table
##
##########################################################################

uni_fits_all <- c(uni_fits, uni_fits_quad, null_fit_list)

fit_list <- fitList(uni_fits_all)

mod_sel <- modSel(fit_list, nullmod="null") %>% 
  rownames_to_column(var = "model")



estimates <- list()
for(i in 1:length(uni_fits_all)){
  estimates[[i]] <- summary(uni_fits_all[[i]], "state") %>% 
    rownames_to_column(var = "term") %>% 
    filter(str_detect(term, "Intercept")== FALSE & str_detect(term, "station_id")== FALSE) %>% 
    mutate(model = names(uni_fits_all)[i], .before = term)
}

estimates <- bind_rows(estimates)


occu_uni_summ <- mod_sel %>% left_join(estimates, by=join_by(model)) %>% relocate(elpd:weight, .after=Rhat)



#remove redundant linear terms from quadratic models
occu_uni_summ_filt <- occu_uni_summ %>% filter(str_detect(term, "I") ==TRUE | str_detect(model, "2")==FALSE)

#write_csv(occu_uni_summ_filt, "feature_selection/Occupancy/occu_uni_summary_quad_residents_9-18-23.csv")







