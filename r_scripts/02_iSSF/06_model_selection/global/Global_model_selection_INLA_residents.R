#### Global model selection all subsets INLA ####

# Author: Read Barbee

# Date:2023-09-28 

# Purpose: Feature selection for global ssf model using INLA


################################ Libraries #################################
library(tidyverse)
#library(INLA)
library(glmmTMB)
library(beepr)


################################ User defined Parameters ###########################

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

#set means of prior distributions for fixed effects means and precisions. from https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40#inla-1

mean.beta = 0
prec.beta = 1e-4

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

#get number of individuals from dataframe
n_indiv = steps_scaled %>% distinct(animal_id) %>% count() %>% pull()

#create variable names for INLA
vars <- vector()
for(i in 1:length(params)){
  vars[i] <- paste0(params[i]," + f(", paste0("id", i), ", ", params[i], ", ", "values = ", paste0("1:",n_indiv), ", ", 'model = "iid", ', "hyper = list(theta = list(initial = log(1), fixed = FALSE, ", 'prior = "pc.prec", ', "param = c(1, .05))))")
}

#make case numeric instead of logical for INLA
dat <- steps_scaled %>% mutate(case_ = as.numeric(case_))

#create separate animal_id columns for each random effect
for(i in 1:length(vars)){
  name <- as.name(paste0("id", i))
  dat[[name]] <- as.numeric(factor(dat$animal_id))
}

#~7 min for model 18 to run without priors specified

#create list of all combinations of variables
all_comb <- do.call("c", lapply(seq_along(vars), function(i) combn(vars, i, FUN = list)))
response <- "case_"

#initialize lists to store model formulas and waic scores
mods <- list()
waic_list <- list()
## Took 32273 seconds (~9hrs) on Legion for 11 covariates (2047 combinations)
system.time(
  for(i in 1:length(all_comb)){
    var_i <- all_comb[[i]]
    form <- as.formula(paste(response, paste("-1", "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T)))", paste(var_i, collapse="+"), sep="+"), sep="~"))
    
    mods[[i]] <- inla(form, family ="Poisson", data=dat,
                      #control.fixed = list(
                        #mean = mean.beta,
                        #prec = list(default = prec.beta)),
                      control.compute = list(waic=TRUE)) #dic = TRUE,
    
    waic_list[[i]] <- mods[[i]]$waic$waic
    
    print(paste0(i, "/", length(all_comb)))
  }
)





# write_csv(dredge_table_tmb, "feature_selection/all_subsets_model_selection_residents_tmb_9-30-23.csv")

