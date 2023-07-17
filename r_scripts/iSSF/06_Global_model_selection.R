#### All Subsets Model Selection ####

# Author: Read Barbee

# Date:2023-07-17 

# Purpose:

################################ Libraries #################################
library(tidyverse)
library(amt)

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


#set all negative elevations to 0 and filter dispersal tracks for post dispersal event
steps <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad)) %>% 
  filter(is.na(dispersing) | dispersing==TRUE)


#########################################################################
##
## 3. Reduce land cover and land use categories
##
##########################################################################
steps <- steps %>% mutate(land_cover_usfs_lumped = fct_collapse(land_cover_usfs, trees = c("trees", "tall_trees_shrubs", "gfh_tree_mix", "tree_shrub_mix",  "barren_tree_mix"),
                                                                shrubs = c("tall_shrubs", "shrubs", "gfh_shrub_mix", "barren_shrub_mix"),
                                                                gfh = c("gfh", "barren_gfh_mix"),
                                                                barren = c("barren_impervious", "snow_ice")),
                          land_use_usfs_lumped = fct_collapse(land_use_usfs, agriculture = c("agriculture", "rangeland_pasture")), .after = land_cover_usfs)

plot_bar(steps %>% select(gpp:calving_season))


#########################################################################
##
## 4. Scale continuous covariates for model comparison and correlation analysis
##
##########################################################################

steps_scaled <- steps %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric))

plot_histogram(steps_scaled %>% select(sl_, ta_, gpp:power_hii))


#########################################################################
##
## 4. Fit all model subsets (Dredge)
##
##########################################################################


vars <- c("evi",
          "dist_water",
          "perc_tree_cover",
          "tip",
          "popdens_hii",
          "roads_hii",
          "infra_hii")



all_comb <- do.call("c", lapply(seq_along(vars), function(i) combn(vars, i, FUN = list)))
response <- "case_"
mods <- list()
aic_list <- list()
## Took 32273 seconds (~9hrs) on Legion for 11 covariates (2047 combinations)
system.time(
  for(i in 1:length(all_comb)){
    var_i <- all_comb[[i]]
    form <- as.formula(paste(response, paste(paste(var_i, " strata(step_id_)", collapse="+"), sep="+"), sep="~"))
    
    mods[[i]] <- fit_issf(form, data=steps_scaled)
    
    aic_list[[i]] <- mods[[i]]$aic$aic
    
    print(i/length(all_comb))
  }
)








