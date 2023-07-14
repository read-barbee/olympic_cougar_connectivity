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

## imputation just seems to make confidence intervals wider and estimates worse


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

#control = glmmTMBControl(parallel = 9), in glmmTMB doesn't appear to be supported in mac

#test <- uni_fit("gpp")

# library(furrr)
# plan(multisession, workers = 8)

cov_names <- covs %>% select(-c(case_:step_id_)) %>% names()

#map the function across all covariates. Too big for future_map(). ~ 47 min
system.time(uni_fits <- map(cov_names, uni_fit, .progress=TRUE)); beep("fanfare")

names(uni_fits) <- cov_names

#save model fits--huge file and takes forever. Probably not worth it.
#save(uni_fits, file = "fitted_models/muff_uni_fits_imp_7-13-23.RData")

#muff_uni_table <- sjPlot::tab_model(uni_fits)

#make table of AIC values for all univariate models--not working
#muff_aic <- AICcmodavg::aictab(uni_fits, second.ord = FALSE)
# muff_bic <- AICcmodavg::bictab(uni_fits, second.ord = FALSE)


#uni_fits2 <- uni_fits[3:length(uni_fits)]

### Create model summary table of coefficient estimates and AIC/BIC values
test_summ <- list()

for (i in 1:length(uni_fits)){
  name <- cov_names[[i]]
  aic<- broom.mixed::glance(uni_fits[[i]])$AIC
  bic <- broom.mixed::glance(uni_fits[[i]])$BIC
  ll <- broom.mixed::glance(uni_fits[[i]])$logLik
  
  test_summ[[i]] <- broom.mixed::tidy(uni_fits[[i]]) %>% 
    filter(effect=="fixed" & !(term %in% c("(Intercept)", "sd__(Intercept)"))) %>% 
    #summarize(mean = mean(estimate), se = plotrix::std.error(estimate)) %>% 
    select(estimate:p.value) %>% 
    mutate(AIC = aic,
           BIC = bic,
           LogLik = ll,
           name = name) %>% 
    select(name, everything())
}

muff_mod_summ <- bind_rows(test_summ)# %>% filter(!grepl("usfs", name))

#write_csv(muff_mod_summ, "feature_selection/uni_muff_summary_no_imp_7-13-23.csv")

# TPI, roads, popdens, rails, and power didn’t converge (no AIC or loglik vals)

#look at covariate estimates by individual
sjPlot::plot_model(uni_fits$evi, type="re",
                   transform=NULL,
                   axis.title = "Coefficient Estimates (Untransformed)",
                   title = "Responses to EVI by Individual")

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


#########################################################################
##
## 8. Random Forest
##
##########################################################################

# #Use unscaled steps to get actual values in partial dependence plots.
# covs3 <- steps %>% 
#   #select(animal_id, gpp:calving_season, -c(land_cover_usfs, land_use_usfs)) %>% 
#   dummify() %>% 
#   dplyr::mutate(case_=as.factor(covs$case_))

covs3 <- steps %>% 
  select(animal_id, gpp:calving_season, -c(land_cover_usfs, land_use_usfs)) %>% 
  dummify() %>% 
  dplyr::mutate(case_=as.factor(steps$case_))

covs3 <- covs3 %>% na.omit()

rf_fits <- list()

n_trials <- 5

#fit multiple rf models sampling the data differently
for (i in 1:n_trials){
  data_set_size <- floor(nrow(covs3)*0.80)
  index <- sample(1:nrow(covs3), size = 5000)
  training <- covs3[index,]
  testing <- covs3[-index,]
  
  rf_fits[[i]] <- randomForest(case_ ~ ., data = training, mtry = 4, ntree = 501, importance = TRUE, type="classification")
}

#combine all fitted rf models
rf_combined <- randomForest::combine(rf_fits[[1]], rf_fits[[2]] , rf_fits[[3]], rf_fits[[4]], rf_fits[[5]])

#create clean dataframe
rf_combined_imp <- as_tibble(rf_combined$importance) %>% 
  janitor::clean_names() %>% 
  mutate(cov = rownames(rf_combined$importance), .before = false)

#write_csv(rf_combined_imp, "rf_combined_imp_7-6-23.csv")

#basic plot
plot(rf_fits[[1]])

varImpPlot(rf_combined)

#partial dependence plots n=40 variables
dev.new(width = 8, height = 8)
#par(mar = c(2, 2, 2, 2))
par(mfrow = c(4, 5))
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = evi, main= "EVI")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = elevation, main= "Elevation")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = perc_nontree_veg, main = "Perc_nontree_veg")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = ndvi, main= "NDVI")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = npp, main="NPP")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = gpp, main = "GPP")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = perc_nonveg, main = "Perc_nonveg")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = precip, main = "Precip")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = popdens_hii, main = "Popdens_hii")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = perc_tree_cover, main = "Perc_tree_cover")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = tri, main = "TRI")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = tpi, main = "TPI")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = northing, main = "Northing")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = slope, main = "Slope")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = roads_hii, main = "Roads_hii")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = tree_cover_hansen, main = "Tree_cover_hansen")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = easting, main = "Easting")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = dist_water, main = "Dist_water")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = power_hii, main = "Power_hii")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = landuse_hii, main = "Landuse_hii")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = rails_hii, main = "Rails_hii")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = infra_hii, main = "Infra_hii")

partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = land_cover_usfs_lumped_trees, main = "Landcover_trees")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = land_cover_usfs_lumped_shrubs, main = "Landcover_shrubs")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = land_cover_usfs_lumped_gfh, main = "Landcover_gfh")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = land_cover_usfs_lumped_barren, main = "Landcover_barren")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = land_cover_usfs_lumped_water, main = "Landcover_water")

partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = land_use_usfs_lumped_forest, main = "Landuse_forest")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = land_use_usfs_lumped_developed, main = "Landuse_developed")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = land_use_usfs_lumped_agriculture, main = "Landuse_agriculture")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = land_use_usfs_lumped_non_forest_wetland, main = "Landuse_wetland")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = land_use_usfs_lumped_other, main = "Landuse_other")

partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = season_dry, main = "Season_dry")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = season_wet, main = "Season_wet")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = hunting_season_yes, main = "Hunting_season_yes")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = hunting_season_no, main = "Hunting_season_no")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = calving_season_no, main = "Calving_season_no")
partialPlot(rf_combined, pred.data = as.data.frame(training), x.var = calving_season_yes, main = "Calving_season_yes")


################################ GRAVEYARD #################################

# Muff uni fits in a for loop
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