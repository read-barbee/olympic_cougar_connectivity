#### OCP iSSF Module_04: Univariate Analysis ####

# Author: Read Barbee

# Date:2023-09-27 

# Purpose: Univariate feature selection for step selection models


################################ Libraries #################################
library(tidyverse)
# library(terra)
library(amt)
# library(janitor)
library(DataExplorer)
library(GGally)
library(sjPlot)
library(doParallel)
library(beepr)
library(INLA)

## imputation just seems to make confidence intervals wider and estimates worse

#register doParallel backend. Note:the cores argument implies forking with doMC backend, but you can check with the getDoParName() function
doParallel::registerDoParallel(cores = 9) 

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


#%>% update_columns(c("sex", "dispersal_status", "case_", "land_cover_usfs", "land_use_usfs", "season", "hunting_season", "calving_season"), as.factor) %>% 


#########################################################################
##
## 2. Check covariate distributions
##
##########################################################################

#data overview
introduce(steps)
plot_intro(steps)

#check for missing values
plot_missing(steps)


#continuous histograms
plot_histogram(steps %>% select(sl_, ta_, gpp:power_hii))

#categorical bar plots
plot_bar(steps %>% select(gpp:calving_season))
plot_bar(steps %>% select(sex, dispersal_status, gpp:calving_season), by="dispersal_status")

#full report with all DataExplorer metrics
#create_report(steps, y="case_")

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

#DON'T Scale step length and turn angle

steps_scaled <- steps %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric))

#plot_histogram(steps_scaled %>% select(sl_, ta_, gpp:power_hii))

#plot_boxplot(steps_scaled %>% select(case_, sl_, ta_, gpp:power_hii), by = "case_")

#########################################################################
##
## 5. Identify highly correlated/redundant covariates
##
##########################################################################

#Correlation plot
plot_correlation(steps_scaled %>% 
                   select(sl_, ta_, gpp:calving_season) %>% 
                   select_if(is.numeric) %>% 
                   na.omit())

#select covariates of interest and dummy code factors
dummy <- steps_scaled %>% select(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii) %>% 
 #select(-c(land_cover_usfs, land_use_usfs)) %>% 
  na.omit() #%>% 
  #dummify()

#create a correlation matrix
cor_matrix <- cor(dummy, use="pairwise.complete.obs")

# #identify correlation values above 0.3
# high_cor <- caret::findCorrelation(cor_matrix, cutoff = .3)
# 
# #convert identified indices to covariate names
# dummy %>% select(!!!high_cor) %>% names()

threshold <- 0.5

high_cor_pairs <- which(abs(cor_matrix) > threshold, arr.ind = TRUE)
high_cor_pairs <- high_cor_pairs[high_cor_pairs[, 1] != high_cor_pairs[, 2], ]

var1 <- vector()
var2 <- vector()
correlation <- vector()

for (i in 1:nrow(high_cor_pairs)) {
  var1[i] <- rownames(cor_matrix)[high_cor_pairs[i, 1]]
  var2[i] <- colnames(cor_matrix)[high_cor_pairs[i, 2]]
  correlation[i] <- cor_matrix[high_cor_pairs[i, 1], high_cor_pairs[i, 2]]
  #print(paste("Variables:", var1, "and", var2, "- Correlation:", correlation))
}

#ordered list of correlated covariates
cor_dat <- tibble(variables = paste0(var1, "_", var2), corr = correlation) %>% 
  #filter(corr < 0.977) %>% 
  arrange(desc(corr))

#write_csv(cor_dat, "feature_selection/pairwise_cov_corr_t5_no_imp_residents_9-18-23.csv")

#the GGAlly package corrplot is more informative for numeric variables
ggcorr(steps_scaled %>% 
         select(case_, sl_, ta_, gpp:calving_season), label = TRUE)

#VERY slow
#ggpairs(steps_scaled %>% 
#select(case_, gpp:elevation), label = TRUE)



#########################################################################
##
## 6. Try out PCA for feature reduction (optional)
##
##########################################################################

#plot_prcomp(steps_scaled %>% na.omit() %>% select(gpp:calving_season), nrow=1, ncol=2, parallel = TRUE)

#PC1: Slope, elevation, tri, and some veg
#PC2: Veg, landuse and pop density
#PC3: Seasonality

#environmental factors most important, followed by human factors


#########################################################################
##
## 7. Univariate Muff models (INLA)
##
##########################################################################
#subset dataframe to only include animal_id, case_, and all predictors
covs <- steps_scaled %>% 
  select(case_, animal_id, step_id_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii) %>% 
  mutate(case_ = as.numeric(case_),
         animal_id = as.numeric(as.factor(animal_id)))


#### NULL MODEL ###
run_null <- function(dat){
  #convert response from categorical to numeric
  dat <- dat %>% mutate(case_ = as.numeric(case_))
  
  form <- case_ ~  f(step_id_, 
                     model = "iid",
                     hyper = list(theta = list(initial = log(1e-6), fixed = TRUE)))
  
  fit <- inla(form, family = "poisson", data = dat, control.compute=list(dic=TRUE, 
                                                                         cpo=TRUE, 
                                                                         waic=TRUE, 
                                                                         return.marginals.predictor=TRUE))
  
  return(fit)
}


null_fit <- run_null(covs)

null_fit <- list(null_fit)

names(null_fit) <- "null"


#### LINEAR MODELS ###

cov_names <- covs %>% select(-c(case_:step_id_)) %>% names()

n_indiv = covs %>% distinct(animal_id) %>% count() %>% pull()

#inla doesn't work in parallel loops. Parallelizes automatically

uni_fits <- list()
for (i  in 1:length(cov_names)) {
  
  cov <- cov_names[i]
  
  #construct the model formula with linear terms only
  form <- as.formula(paste("case_ ~ ",  
                           #fixed effects
                           cov, "+",
                           #random intercept (strata)
                           "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6), fixed=T))) + ",
                           #random slopes
                           paste0("f(animal_id, ", cov, 
                                  ", values = ", 
                                  paste0("1:",n_indiv), 
                                  ", model = \"iid\", hyper = list(theta = list(initial = log(1),",
                                  "fixed = FALSE, prior = \"pc.prec\", param = c(1, .05))))")))
  
  
 uni_fits[[i]] <- inla(form, family = "poisson", data = covs, control.compute=list(dic=TRUE, 
                                                                                 cpo=TRUE, 
                                                                                 waic=TRUE, 
                                                                                 return.marginals.predictor=TRUE))
  
} ; beep("fanfare")


names(uni_fits) <- cov_names


#### QUADRATIC MODELS ###

#takes about 3 minutes to fit 22 models; ~26 minutes when calculating model selection statistics
uni_fits_quad <- list()
system.time(for (i  in 1:length(cov_names)){
  
  cov <- cov_names[i]
  
  form <- as.formula(paste("case_ ~ ",  
                           #fixed effects
                           cov, "+", "I(", cov,  "^2) + ",
                           #random intercept (strata)
                           "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6), fixed=T))) + ",
                           #random slopes
                           paste0("f(animal_id, ", cov, 
                                  ", values = ", 
                                  paste0("1:",n_indiv), 
                                  ", model = \"iid\", hyper = list(theta = list(initial = log(1),",
                                  "fixed = FALSE, prior = \"pc.prec\", param = c(1, .05))))")))
  
  uni_fits_quad[[i]] <- inla(form, family = "poisson", data = covs, control.compute=list(dic=TRUE, 
                                                                                         cpo=TRUE, 
                                                                                         waic=TRUE, 
                                                                                         return.marginals.predictor=TRUE))
  
}); beep("fanfare")


#rename qudaratic covariates with "2" at the end
quad_names <- vector()
for(i in 1:length(cov_names)){
  old_name <- cov_names[i]
  quad_names[i] <- paste0(old_name, "2")
}

names(uni_fits_quad) <- quad_names


#combine linear and quadratic fits
uni_fits_all <- c(uni_fits, uni_fits_quad, null_fit)


#########################################################################
##
## 8. Create model summary table
##
##########################################################################

mod_names_all <- names(uni_fits_all)



mod_table <- list()
for (i in 1:length(uni_fits_all)){
  name <- mod_names_all[i]
  
  if(str_detect(name, "2")){
    tmp1 <-  uni_fits_all[[i]]$summary.fixed %>% 
      rownames_to_column() %>%  
      as_tibble() %>% filter(str_detect(rowname, "2")==T)
  } else if (name=="null"){
    tmp1 <-  uni_fits_all[[i]]$summary.fixed %>% 
      rownames_to_column() %>%  
      as_tibble()
  } else{
    tmp1 <-  uni_fits_all[[i]]$summary.fixed %>% 
      rownames_to_column() %>%  
      as_tibble() %>% filter(str_detect(rowname, coll("("))==F)
  }
  
  tmp2 <-  tmp1 %>% 
    dplyr::mutate(term = name, .before = mean) %>% 
    dplyr::mutate(tau_mean = uni_fits_all[[i]]$summary.hyperpar$mean,
                  tau_sd = uni_fits_all[[i]]$summary.hyperpar$sd,
                  waic = uni_fits_all[[i]]$waic$waic,
                  p.eff_waic = uni_fits_all[[i]]$waic$p.eff,
                  dic = uni_fits_all[[i]]$dic$dic,
                  p.eff_dic = uni_fits_all[[i]]$dic$p.eff,
                  mean_cpo = mean(uni_fits_all[[i]]$cpo$cpo),
                  mean_pit = mean(uni_fits_all[[i]]$cpo$pit),
                  mlik = uni_fits_all[[i]]$mlik[1])
  
  if(name=="null"){
    mod_table[[i]] <- tmp2 %>% dplyr::mutate(tau_mean = NA, tau_sd = NA) %>% 
      dplyr::select(term, mean, everything(), -rowname)
  } else{
    mod_table[[i]] <-tmp2 %>% 
      dplyr::select(term, mean, everything(), -rowname)
    }
    
  
  print(i)
}

mod_table <- bind_rows(mod_table)


#write_csv(mod_table, "feature_selection/uni_muff_summary_no_imp_residents_inla_9-27-23.csv")
# write_csv(muff_mod_summ_filt, "feature_selection/uni_muff_summary_quad_mods_no_imp_residents_9-18-23.csv")



#remove redundant linear terms from quadratic models
#muff_mod_summ_filt <- muff_mod_summ %>% filter(str_detect(term, "I") ==TRUE | str_detect(name, "2")==FALSE)






