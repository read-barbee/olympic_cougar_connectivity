#### Global model selection all subsets glmmTMB ####

# Author: Read Barbee

# Date:2023-09-28 

# Purpose: Feature selection for global ssf model using INLA


################################ Libraries #################################
library(tidyverse)
library(glmmTMB)
library(beepr)


#########################################################################
##
## 1. Import and format step data
##
##########################################################################

#choose whether to use steps without imputation, imputation but no rerouting, or iimputation and rerouting

#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv")

#imputed and rerouted
#steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_imputed_7-12-2023.csv")


#set all negative elevations to 0 and filter for residents only
steps <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad, disp_qual:dispersing)) %>% 
  filter(dispersal_status =="resident")

#########################################################################
##
## 4. Scale continuous covariates for model comparison and correlation analysis
##
##########################################################################

#DON'T Scale step length and turn angle

steps_scaled <- steps %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric))

#plot_histogram(steps_scaled %>% select(sl_, ta_, gpp:power_hii))
#plot_boxplot(steps_scaled %>% select(case_, sl_, ta_, gpp:power_hii), by = "case_")

#########################################################################
##
## 3. Manual dredge (INLA)
##
##########################################################################

#List of parameters to compute all subset models for
params <- c("tree_cover_hansen",
            "ndvi",
            "tpi",
            "northing", 
            "tri", 
            "perc_tree_cover",
            "roads_hii",
            "popdens_hii",
            "landuse_hii",
            #"rails_hii",
            "infra_hii")

#List of quadratic parameters if necessary 
# quad_params <- c("roads_hii",
#                  "tpi",
#                  "npp",
#                  "perc_tree_cover",
#                  "ndvi",
#                  "tree_cover_hansen",
#                  "northing",
#                  "easting",
#                  "tri",
#                  "precip")


#########################################################################
##
## 3. Manual dredge (glmmTMB)
##
##########################################################################


vars_tmb <- vector()
for(i in 1:length(params)){
  vars_tmb[i] <- paste0(params[i]," + (0 + ", params[i], " | animal_id)")
}


all_comb_tmb <- do.call("c", lapply(seq_along(vars_tmb), function(i) combn(vars_tmb, i, FUN = list)))


#~ 10 min for 2,097,151 combinations of 21 covariates
forms <- list()

#Takes about 30 hours for 1023 combinations of 10 covariates
for (i in 1:length(all_comb_tmb)){
  var_i <- all_comb_tmb[[i]]
  forms[[i]] <- as.formula(paste("case_",  paste("-1", "(1|step_id_)", paste(var_i, collapse="+"), sep="+"), sep="~"))
}

mods_tmb <- list()
aic_list <- list()
bic_list <- list()
system.time(for (i in 1:length(forms)){
  form <- forms[[i]]
  mod <- glmmTMB(form, family=poisson, data=steps_scaled, doFit=FALSE)
  mod$parameters$theta[1] <- log(1e3)
  map_length <- 1:(length(mod$parameters$theta)-1)
  mod$mapArg <- list(theta=factor(c(NA, map_length)))
  
  fit <- fitTMB(mod)
  
  mods_tmb[[i]] <- form
  aic_list[[i]] <- AIC(fit)
  bic_list[[i]] <- BIC(fit)
  
  print(paste0(i, "/", length(all_comb_tmb)))
})

#dredge selection table
dredge_table_tmb <- tibble(model = unlist(as.character(mods_tmb)),
                           aic = unlist(aic_list),
                           bic = unlist(bic_list)) %>% 
  mutate(model = str_remove(model, coll("case_ ~ -1 + (1 | step_id_) + "))) 
  

# write_csv(dredge_table_tmb, "feature_selection/all_subsets_model_selection_residents_tmb_9-30-23.csv")

