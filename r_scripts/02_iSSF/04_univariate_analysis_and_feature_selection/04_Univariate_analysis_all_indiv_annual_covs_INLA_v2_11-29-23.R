#### OCP iSSF Module_04: Univariate Analysis ####

# Author: Read Barbee

# Date:2023-09-27 

#updated: 2023-11-29

# Purpose: Univariate feature selection for step selection models

#added random effects for quadratic models 11/29/23

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
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_annual_cov_11-30-2023.csv")

#imputed and rerouted
#steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_imputed_7-12-2023.csv")


#round and define factors
steps <- steps %>% 
  mutate(across(c(perc_tree_cov_annual, perc_nonveg_annual, roads_hii_annual), round)) %>% 
  mutate(across(c(popdens_hii_annual, 
                  infra_hii_annual, 
                  #roads_hii_annual, 
                  hii_annual), function(x){return(round(x, digits=2))})) %>% 
  mutate(across(c(#animal_id, 
                  sex, 
                  dispersal_status, 
                  year, 
                  #case_, 
                  tod_end_, 
                  intersects_water, 
                  intersects_study_area, 
                  land_cover_usfs_annual, 
                  land_use_usfs_annual, 
                  land_use_change_usfs_annual, 
                  dispersing, 
                  disp_qual), as.factor)) 

#remove steps when dispersers were not dispersing and drop NA factor values
steps_filt <- steps %>% 
  filter(dispersal_status=="resident" | dispersing ==TRUE) %>% 
  mutate(land_cover_usfs_annual = case_when(land_cover_usfs_annual == 15 ~ NA,
                                            .default = land_cover_usfs_annual),
         land_use_usfs_annual = case_when(land_use_usfs_annual == 7 ~ NA,
                                          .default = land_use_usfs_annual),
         land_use_change_usfs_annual = case_when(land_use_change_usfs_annual == 5 ~ NA,
                                                 .default = land_use_change_usfs_annual)) %>% 
  select(-c(disp_qual:dispersing))


#inspect steps
# steps_filt %>% 
#   #filter(intersects_water==TRUE) %>% 
#   filter(intersects_study_area==FALSE) %>% 
#   #filter(case_==FALSE) %>%
#   sf::st_as_sf(coords=c("x2_", "y2_"), crs=5070) %>% 
#   mapview::mapview()

#########################################################################
##
## 2. Check covariate distributions
##
##########################################################################

#data overview
introduce(steps_filt)
plot_intro(steps_filt)

#check for missing values
plot_missing(steps_filt)


#continuous histograms all
plot_histogram(steps_filt %>% select(sl_, ta_, elevation:dens_all_roads_annual))

#continuous histograms used
plot_histogram(steps_filt %>% select(case_, sl_, ta_, elevation:dens_all_roads_annual) %>% filter(case_==TRUE))

#continuous histograms unused
plot_histogram(steps_filt %>% select(case_, sl_, ta_, elevation:dens_all_roads_annual) %>% filter(case_==FALSE))

#categorical bar plots
plot_bar(steps_filt %>% select(land_cover_usfs_annual, land_use_usfs_annual, land_use_change_usfs_annual))

plot_bar(steps_filt %>% select(land_cover_usfs_annual, land_use_usfs_annual, land_use_change_usfs_annual, dispersal_status), by="dispersal_status")

#full report with all DataExplorer metrics
#create_report(steps, y="case_")

#########################################################################
##
## 3. Reduce land cover and land use categories
##
##########################################################################
# steps <- steps %>% mutate(land_cover_usfs_lumped = fct_collapse(land_cover_usfs, trees = c("trees", "tall_trees_shrubs", "gfh_tree_mix", "tree_shrub_mix",  "barren_tree_mix"),
#                                                                 shrubs = c("tall_shrubs", "shrubs", "gfh_shrub_mix", "barren_shrub_mix"),
#                                                                 gfh = c("gfh", "barren_gfh_mix"),
#                                                                 barren = c("barren_impervious", "snow_ice")),
#                           land_use_usfs_lumped = fct_collapse(land_use_usfs, agriculture = c("agriculture", "rangeland_pasture")), .after = land_cover_usfs)
# 
# plot_bar(steps %>% select(gpp:calving_season))


#########################################################################
##
## 4. Scale continuous covariates for model comparison and correlation analysis
##
##########################################################################

#DON'T Scale step length and turn angle

steps_scaled <- steps_filt %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dens_all_roads_annual), scale)) %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dens_all_roads_annual), as.numeric))

#plot_histogram(steps_scaled %>% select(sl_, ta_, gpp:power_hii))

#plot_boxplot(steps_scaled %>% select(case_, sl_, ta_, gpp:power_hii), by = "case_")

#########################################################################
##
## 5. Identify highly correlated/redundant covariates
##
##########################################################################

#Correlation plot
plot_correlation(steps_scaled %>% 
                   select(sl_, ta_, elevation:dens_all_roads_annual) %>% 
                   select_if(is.numeric) %>% 
                   na.omit())

#select covariates of interest and dummy code factors
cont <- steps_scaled %>% 
  select(elevation:dens_all_roads_annual, -aspect) %>% 
  dummify(select = c("land_cover_usfs_annual", "land_use_usfs_annual", "land_use_change_usfs_annual")) %>%
  mutate(across(elevation:dens_all_roads_annual, as.numeric)) %>% 
  select(!contains("NA"))
  

#create a correlation matrix
cor_matrix <- cor(cont, use="pairwise.complete.obs", method = "pearson")

# Set correlation threshold
threshold <- 0.0

high_cor_pairs <- which(abs(cor_matrix) >= threshold, arr.ind = TRUE)
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

rows_to_remove <- seq(from = 2, to = nrow(cor_dat), by = 2)
cor_dat_no_dupes <- cor_dat[-rows_to_remove, ]

#write_csv(cor_dat_no_dupes, "feature_selection/pairwise_cov_corr_no_imp_all_indiv_annual_cov_11-30-23.csv")

#the GGAlly package corrplot is more informative for numeric variables
ggcorr(steps_scaled %>% 
         select(elevation:dist_all_roads_annual), label = TRUE)




#########################################################################
##
## 6. Try out PCA for feature reduction (optional)
##
##########################################################################

#plot_prcomp(steps_scaled %>% na.omit() %>% select(elevation:dist_all_roads_annual), nrow=1, ncol=2, parallel = TRUE)


#########################################################################
##
## 7. Univariate Muff models (INLA)
##
##########################################################################
#subset dataframe to only include animal_id, case_, and all predictors
covs <- steps_scaled %>% 
  select(case_, animal_id, step_id_, elevation:dens_all_roads_annual, -aspect) %>%
  mutate(across(c(case_, land_cover_usfs_annual:land_use_change_usfs_annual), as.numeric)) %>% 
  mutate(animal_id = as.numeric(as.factor(animal_id)))


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

covs <- covs %>% mutate(animal_id2 = animal_id, .after=animal_id)

# ~ 11 min
uni_fits <- list()
system.time(for (i  in 1:length(cov_names)) { 
  
  cov <- cov_names[i]
  
  #construct the model formula with linear terms only
  form <- as.formula(paste("case_ ~ -1 + ",  
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
  
  print(paste0(i, "/", length(cov_names)))
 
}) ; beep("fanfare")


names(uni_fits) <- cov_names


#### CATEGORICAL MODELS ###

covs_cat <- covs %>% 
  select(case_, animal_id, step_id_, land_cover_usfs_annual:land_use_change_usfs_annual) %>%
  mutate(across(land_cover_usfs_annual:land_use_change_usfs_annual, as.factor)) %>% 
  dummify(select = c("land_cover_usfs_annual", "land_use_usfs_annual", "land_use_change_usfs_annual")) %>% 
  select(!contains("NA")) %>% 
  rename(land_cover_usfs_trees = land_cover_usfs_annual_1,
         land_cover_usfs_tall_tree_shrub = land_cover_usfs_annual_2,
         land_cover_usfs_tree_shrub = land_cover_usfs_annual_3,
         land_cover_usfs_gfh_tree = land_cover_usfs_annual_4,
         land_cover_usfs_barren_tree = land_cover_usfs_annual_5,
         land_cover_usfs_tall_shrub = land_cover_usfs_annual_6,
         land_cover_usfs_shrubs = land_cover_usfs_annual_7,
         land_cover_usfs_gfh_shrub = land_cover_usfs_annual_8,
         land_cover_usfs_barren_shrub = land_cover_usfs_annual_9,
         land_cover_usfs_gfh = land_cover_usfs_annual_10,
         land_use_usfs_agr = land_use_usfs_annual_1,
         land_use_usfs_dev = land_use_usfs_annual_2,
         land_use_usfs_forest = land_use_usfs_annual_3,
         land_use_usfs_non_forest_wet = land_use_usfs_annual_4,
         land_use_usfs_other = land_use_usfs_annual_5,
         land_use_usfs_rangeland = land_use_usfs_annual_6,
         land_use_change_stable = land_use_change_usfs_annual_1,
         land_use_change_slow_loss = land_use_change_usfs_annual_2,
         land_use_change_fast_loss = land_use_change_usfs_annual_3,
         land_use_change_gain = land_use_change_usfs_annual_4)
  

cov_names_cat <- covs_cat %>% select(-c(case_:step_id_)) %>% names()

# ~12 minutes
uni_fits_cat <- list()
system.time(for (i  in 1:length(cov_names_cat)) {
  
  cov <- cov_names_cat[i]
  
  #construct the model formula with linear terms only
  form <- as.formula(paste("case_ ~ -1 + ",  
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
  
  
  uni_fits_cat[[i]] <- inla(form, family = "poisson", data = covs_cat, control.compute=list(dic=TRUE, 
                                                                                    cpo=TRUE, 
                                                                                    waic=TRUE, 
                                                                                    return.marginals.predictor=TRUE))
  
  print(paste0(i, "/", length(cov_names_cat)))
  
}) ; beep("fanfare")


names(uni_fits_cat) <- cov_names_cat


#### QUADRATIC MODELS ###

#takes about 13 minutes to fit 25 models;
cov_names_quad <- covs %>% select(-c(case_:step_id_, land_cover_usfs_annual:land_use_change_usfs_annual)) %>% names()

uni_fits_quad <- list()
system.time(for (i  in 1:length(cov_names_quad)){
  
  cov <- cov_names_quad[i]
  
  form <- as.formula(paste("case_ ~ -1 + ",  
                           #fixed effects
                           cov, "+", "I(",cov,"^2) + ",
                           #random intercept (strata)
                           "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6), fixed=T))) + ",
                           #random slopes
                           paste0("f(animal_id, ", cov, 
                                  ", values = ", 
                                  paste0("1:",n_indiv), 
                                  ", model = \"iid\", hyper = list(theta = list(initial = log(1),",
                                  "fixed = FALSE, prior = \"pc.prec\", param = c(1, .05)))) + "),
                           paste0("f(animal_id2, ", "I(", cov,  "^2)", 
                                  ", values = ", 
                                  paste0("1:",n_indiv), 
                                  ", model = \"iid\", hyper = list(theta = list(initial = log(1),",
                                  "fixed = FALSE, prior = \"pc.prec\", param = c(1, .05))))")))
  
  uni_fits_quad[[i]] <- inla(form, family = "poisson", data = covs, control.compute=list(dic=TRUE, 
                                                                                         cpo=TRUE, 
                                                                                         waic=TRUE, 
                                                                                         return.marginals.predictor=TRUE))
  
  print(paste0(i, "/", length(cov_names_quad)))
  
}); beep("fanfare")


#rename qudaratic covariates with "2" at the end
quad_names <- vector()
for(i in 1:length(cov_names_quad)){
  old_name <- cov_names_quad[i]
  quad_names[i] <- paste0(old_name, "2")
}

names(uni_fits_quad) <- quad_names



#combine linear and quadratic fits
uni_fits_all <- c(uni_fits, uni_fits_cat, uni_fits_quad, null_fit)


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
    tau_mean = uni_fits_all[[i]]$summary.hyperpar$mean[2]
    tau_sd = uni_fits_all[[i]]$summary.hyperpar$sd[2]
  } else if (name=="null"){
    tmp1 <-  uni_fits_all[[i]]$summary.fixed %>% 
      rownames_to_column() %>%  
      as_tibble()
    tau_mean = NA
    tau_sd = NA
  } else{
    tmp1 <-  uni_fits_all[[i]]$summary.fixed %>% 
      rownames_to_column() %>%  
      as_tibble() %>% filter(str_detect(rowname, coll("("))==F)
    
    tau_mean = uni_fits_all[[i]]$summary.hyperpar$mean
    tau_sd = uni_fits_all[[i]]$summary.hyperpar$sd
  }
  
  tmp2 <-  tmp1 %>% 
    dplyr::mutate(term = name, .before = mean) %>% 
    dplyr::mutate(tau_mean = tau_mean,
                  tau_sd = tau_sd,
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

mod_table_filt <- mod_table %>% filter(!(term %in% c("land_cover_usfs_annual", "land_use_usfs_annual", "land_use_change_usfs_annual", "land_cover_usfs_annual2", "land_use_usfs_annual2", "land_use_change_usfs_annual2")))


#write_csv(mod_table_filt "feature_selection/uni_muff_fits_no_imp_all_indiv_annual_covs_inla_11-30-23.csv")

#remove redundant linear terms from quadratic models
#muff_mod_summ_filt <- muff_mod_summ %>% filter(str_detect(term, "I") ==TRUE | str_detect(name, "2")==FALSE)

# write_csv(muff_mod_summ_filt, "feature_selection/uni_muff_summary_quad_mods_no_imp_residents_9-18-23.csv")
