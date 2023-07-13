#### OCP iSSF Module_04: Univariate Analysis ####

# Author: Read Barbee

# Date:2023-07-12 

# Purpose: Univariate feature selection for step selection models


################################ Libraries #################################
library(tidyverse)
# library(terra)
library(amt)
# library(janitor)
library(DataExplorer)
library(GGally)
library(glmmTMB)
library(sjPlot)
# #library(performance)
# library(caret)
# library(MuMIn)
library(randomForest)
library(doParallel)
library(beepr)


#########################################################################
##
## 1. Import and format step data
##
##########################################################################

#choose whether to use steps without imputation, imputation but no rerouting, or iimputation and rerouting

#no imputation
#steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv")

#imputed and rerouted
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_imputed_7-12-2023.csv")


#set all negative elevations to 0 and filter dispersal tracks for post dispersal event
steps <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad)) %>% 
  filter(is.na(dispersing) | dispersing==TRUE)


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

steps_scaled <- steps %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric))

plot_histogram(steps_scaled %>% select(sl_, ta_, gpp:power_hii))

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
dummy <- steps_scaled %>% select(sl_, ta_, gpp:calving_season) %>% 
  select(-c(land_cover_usfs, land_use_usfs)) %>% 
  na.omit() %>% 
  dummify()

#create a correlation matrix
cor_matrix <- cor(dummy)

# #identify correlation values above 0.3
# high_cor <- caret::findCorrelation(cor_matrix, cutoff = .3)
# 
# #convert identified indices to covariate names
# dummy %>% select(!!!high_cor) %>% names()

threshold <- 0.3

high_cor_pairs <- which(cor_matrix > threshold, arr.ind = TRUE)
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

#write_csv(cor_dat, "feature_selection/pairwise_cov_corr_imp_7-12-23.csv")

#the GGAlly package corrplot is more informative for numeric variables
ggcorr(steps_scaled %>% 
         select(case_, sl_, ta_, gpp:calving_season), label = TRUE)

#VERY slow
#ggpairs(steps_scaled %>% 
#select(case_, gpp:elevation), label = TRUE)



# #not sure if this is useful or not
# cor.prob <- function(X, dfr = nrow(X) - 2) {
#   R <- cor(X, use="complete.obs")
#   above <- row(R) < col(R)
#   r2 <- R[above]^2
#   Fstat <- r2 * dfr / (1 - r2)
#   R[above] <- 1 - pf(Fstat, 1, dfr)
#   Rstar = ifelse(R[above]<0.05, "***", "NS")
#   R[above]=paste(R[above],Rstar)
#   R
# }
# 
# cor_prob <- steps_scaled %>% 
#   select(sl_, ta_, gpp:calving_season) %>%
#   select(where(is.numeric)) %>% 
#   cor.prob()



#########################################################################
##
## 6. Try out PCA for feature reduction (optional)
##
##########################################################################

plot_prcomp(steps_scaled %>% na.omit() %>% select(gpp:calving_season), nrow=1, ncol=2, parallel = TRUE)

#PC1: Slope, elevation, tri, and some veg
#PC2: Veg, landuse and pop density
#PC3: Seasonality

#environmental factors most important, followed by human factors


#########################################################################
##
## 7. Univariate Muff models
##
##########################################################################
#subset dataframe to only include animal_id, case_, and all predictors
covs <- steps_scaled %>% 
  select(-c(land_cover_usfs, land_use_usfs)) %>% 
  dummify(select = c("land_cover_usfs_lumped", 
                     "land_use_usfs_lumped", 
                     "season",
                     "hunting_season",
                     "calving_season")) %>% 
  select(case_, animal_id, step_id_, gpp:calving_season_yes) %>% 
  select(-c(sex:disp_qual))

#function to fit a univariate Muff model with random effect for individual for each covariate
uni_fit <- function(cov){
  form <- as.formula(paste0("case_ ~ ", 
                            #fixed effects
                            cov, "+", 
                            #random intercept (strata)
                            "(1|step_id_) +",
                            #random slopes
                            "(0 +", cov, "| animal_id )"))
  
  uni_mod <- glmmTMB(form, family =poisson, data = covs, doFit = FALSE)
  uni_mod$parameters$theta[1] <-log(1e3)
  uni_mod$mapArg <-list(theta=factor(c(NA, 1)))
  fit <- glmmTMB::fitTMB(uni_mod)
  
  return(fit)
}

# library(furrr)
# plan(multisession, workers = 8)

cov_names <- covs %>% select(-c(case_:step_id_)) %>% names()

#map the function across all covariates. Too big for future_map(). ~ 47 min
system.time(uni_fits <- map(cov_names, uni_fit, .progress=TRUE)); beep("fanfare")

#save model fits--huge file and takes forever. Probably not worth it.
#save(uni_fits, file = "fitted_models/muff_uni_fits_imp_7-13-23.RData")

 # The same operation as above but in a for loop
#initialize clusters for parallel computing. parallel doesn't seem to work though
# n.clusters = 8
# my.cluster <- parallel::makeCluster(n.clusters, type = "PSOCK")
# doParallel::registerDoParallel(cl = my.cluster, cores = ncores)
# 
# cov_names <- covs %>% select(-c(case_:step_id_)) %>% names()
# 
# uni_fits <- list()
# 
# system.time(uni_fits <- foreach (i = 1:length(cov_names), .packages=c("glmmTMB"))%dopar%{
#   
#   cov <- cov_names[i]
#   
#   form <- as.formula(paste0("case_ ~ ", 
#                        #fixed effects
#                        cov, "+", 
#                        #random intercept (strata)
#                        "(1|step_id_) +",
#                        #random slopes
#                        "(0 +", cov, "| animal_id )"))
#   
#   
#   uni_mod <- glmmTMB(form, family =poisson, data = covs, doFit = FALSE)
#   
#   uni_mod$parameters$theta[1] <-log(1e3)
#   uni_mod$mapArg <-list(theta=factor(c(NA, 1)))
#   uni_fits[[i]] <- glmmTMB::fitTMB(uni_mod)
#   
# }); beep("fanfare")
# 
# parallel::stopCluster(cl = my.cluster)

#muff_uni_table <- sjPlot::tab_model(uni_fits)

#make table of AIC values for all univariate models
# muff_aic <- AICcmodavg::aictab(uni_fits, second.ord = FALSE)
# muff_bic <- AICcmodavg::bictab(uni_fits, second.ord = FALSE)


uni_fits2 <- uni_fits[3:length(uni_fits)]

### Create model summary table of coefficient estimates and AIC/BIC values
test_summ <- list()

for (i in 1:length(uni_fits2)){
  name <- names(uni_fits2)[[i]]
  aic<- broom.mixed::glance(uni_fits2[[i]])$AIC
  bic <- broom.mixed::glance(uni_fits2[[i]])$BIC
  ll <- broom.mixed::glance(uni_fits2[[i]])$logLik
  
  test_summ[[i]] <- broom.mixed::tidy(uni_fits2[[i]]) %>% 
    filter(!(term %in% c("(Intercept)", "sd__(Intercept)"))) %>% 
    #summarize(mean = mean(estimate), se = plotrix::std.error(estimate)) %>% 
    select(estimate:p.value) %>% 
    mutate(AIC = aic,
           BIC = bic,
           LogLik = ll,
           name = name) %>% 
    select(name, everything())
}

muff_mod_summ <- bind_rows(test_summ) %>% filter(!grepl("usfs", name))


#look at covariate estimates by individual
sjPlot::plot_model(uni_fits$elevation, type="re",
                   transform=NULL,
                   axis.title = "Coefficient Estimates (Untransformed)",
                   title = "Responses to Elevation by Individual")

#not implemented for glmmTMB yet
#performance::check_model(uni_fits$elevation)

#########################################################################
##
## 8. Univariate SSF models
##
##########################################################################

covs2 <- steps_scaled %>% 
  select(case_, step_id_, gpp:calving_season)
#function to fit a univariate issf model for each covariate
uni_fit_issf <- function(cov){
  uni_mod <- fit_issf(data=covs2, case_ ~ cov + strata(step_id_))
  return(uni_mod)
}

#map the function across all covariates
uni_fits_issf <- map(covs2, uni_fit_issf)


uni_fits_issf2 <- uni_fits_issf[3:length(uni_fits_issf)]

#create model summary table of coefficient estimates and AIC/BIC values
tidy_list <- list()
mod_summ_issf <- list()

for (i in 1:length(uni_fits_issf2)){
  Name <- names(uni_fits_issf2)[[i]]
  tidy_list[[i]] <- broom::tidy(uni_fits_issf2[[i]]$model) %>% 
    mutate(name = Name, .before=term)
  
  mod_summ_issf[[i]] <- broom::glance(uni_fits_issf2[[i]]$model) %>% 
    select(r.squared:BIC, statistic.wald, p.value.wald) %>% 
    mutate(name = Name, .before=r.squared)
}

tidy_list <- bind_rows(tidy_list) 
mod_summ_issf <- bind_rows(mod_summ_issf)

issf_summary_tab <- tidy_list %>% left_join(mod_summ_issf, by= join_by(name)) %>% 
  mutate(name= case_when(term!="cov" ~ paste0(name, "_", term), .default = name)) %>% select(-term)

#write_csv(issf_summary_tab, "feature_selection/uni_issf_summary_imp_7-12-23.csv")
