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
library(DataExplorer)

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
##  2. Check correlation of covariates
##
##########################################################################

c_dat <- occ_dat_scaled %>% select(tree_cover_hansen:dist_water)

#check for missing values
plot_missing(c_dat)

#continuous histograms
plot_histogram(c_dat)

#create a correlation matrix
cor_matrix <- cor(c_dat, use="pairwise.complete.obs")

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

#write_csv(cor_dat, "feature_selection/Occupancy/pairwise_cov_corr_t5_9-19-23.csv")

GGally::ggcorr(c_dat, label = TRUE)

#########################################################################
##
##  3. Fit linear univariate stacked occupancy models
##
##########################################################################

#register doParallel backend. Note:the cores argument implies forking with doMC backend, but you can check with the getDoParName() function
doParallel::registerDoParallel(cores = 9) 


#### NULL MODEL ###

null_fit <- stan_occu(~1 ~ (1|station_id), data=umf, chains=3, iter=1000)

traceplot(null_fit)

traceplot(null_fit, pars = "sigma_state")

null_fit_list <- list(null_fit)

names(null_fit_list) <- "null"


#### UNIVARIATE MODELS ###

cov_names <- occ_covs

#parallelized for loop: fits in ~ 9 minutes with 3 chains and 100 iterations; ~27 min with 1000 iterations
uni_fits <- list()
system.time(uni_fits <- foreach (i = 1:length(cov_names), .packages=c("ubms"))%dopar%{
  
  cov <- cov_names[i]
  
  form <- as.formula(paste0("~1 ~", #remove standard intercept to be replace with stratum-based intercept
                            #fixed effects
                            cov, "+",
                            #random intercept (strata)
                            "(1|station_id)"))
  
  
  stan_occu(form, data=umf, chains=3, iter=1000)
  
}); beep("fanfare")


names(uni_fits) <- cov_names


# fit_list <- fitList(c(uni_fits, null_fit))
# 
# mod_sel <- modSel(fit_list, nullmod="null")

#use loo function for leave one out cross-validation


#########################################################################
##
##  4. Fit quadratic univariate stacked occupancy models
##
##########################################################################

#parallelized for loop: fits in ~ 10 minutes with 3 chains and 100 iterations; 35 min for 1000 iterations
uni_fits_quad <- list()
system.time(uni_fits_quad <- foreach (i = 1:length(cov_names), .packages=c("ubms"))%dopar%{
  
  cov <- cov_names[i]
  
  form <- as.formula(paste0("~1 ~", #remove standard intercept to be replace with stratum-based intercept
                            #fixed effects
                            cov, "+", "I(", cov,  "^2) + ",
                            #random intercept (strata)
                            "(1|station_id)"))
  
  
  stan_occu(form, data=umf, chains=3, iter=1000)
  
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
##  5. Create model summary table
##
##########################################################################

uni_fits_all <- c(uni_fits, uni_fits_quad, null_fit_list)

saveRDS(uni_fits_all, "occu_uni_fits_9-19-23.rds")

uni_fits_all <- readRDS("occu_uni_fits_9-19-23.rds")

fit_list <- fitList(uni_fits_all)

mod_sel <- modSel(fit_list, nullmod="null") %>% 
  rownames_to_column(var = "model")

test <- waic(uni_fits_all[[1]])

waic(uni_fits_all[[1]])$estimates %>% as.data.frame() %>% rownames_to_column(var = "metric") %>% pivot_wider(names_from = "metric", values_from = c("Estimate", "SE"))


estimates <- list()
waic <- list()
for(i in 1:length(uni_fits_all)){
  estimates[[i]] <- summary(uni_fits_all[[i]], "state") %>% 
    rownames_to_column(var = "term") %>% 
    filter(str_detect(term, "Intercept")== FALSE & str_detect(term, "station_id")== FALSE) %>% 
    mutate(model = names(uni_fits_all)[i], .before = term)
  
  waic[[i]] <- waic(uni_fits_all[[i]])$estimates %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "metric") %>% 
    pivot_wider(names_from = "metric", values_from = c("Estimate", "SE")) %>% 
    mutate(model = names(uni_fits_all)[i], .before = Estimate_elpd_waic)
}

estimates <- bind_rows(estimates)
waic <- bind_rows(waic)


occu_uni_summ <- mod_sel %>% left_join(estimates, by=join_by(model)) %>% relocate(elpd:weight, .after=Rhat)



#remove redundant linear terms from quadratic models
occu_uni_summ_filt <- occu_uni_summ %>% 
  filter(str_detect(term, "I") ==TRUE | str_detect(model, "2")==FALSE) %>% 
  left_join(waic, by=join_by(model))



write_csv(occu_uni_summ_filt, "feature_selection/Occupancy/occu_uni_summary_quad_residents_9-18-23.csv")







