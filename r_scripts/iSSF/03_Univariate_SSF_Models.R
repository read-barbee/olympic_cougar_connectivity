#### OCP iSSF Module_03: Univariate SSF Models ####

# Author: Read Barbee

# Date:2023-06-02 

# Purpose: Fit univariate SSF models to full dataset with random effect for individual to explore individual covariate effects

# Inputs:
#   Unscaled amt step data frame for fitting iSSFs
#
# Outputs:
#   •	univariate regression plots
#   •	univariate regression tables

#Steps
# •	Import amt step data




################################ Libraries #################################
library(tidyverse)
library(terra)
library(lubridate)
library(amt)
library(janitor)
library(DataExplorer)
library(GGally)
library(glmmTMB)
library(sjPlot)
library(performance)
library(caret)
library(MuMIn)
library(randomForest)


############################### Import Steps #################################

steps <- read_csv("data/Location_Data/Steps/6h_steps_unscaled_6-02-2023.csv")

steps <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad))

# update_columns(sex, dispersal_Status, case_, land_cover_usfs, land_use_usfs, season, hunting_season, calving_season, as.factor)


################################ Check covariate distributions #################################

#data overview
introduce(steps)
plot_intro(steps)

#check for missing values
plot_missing(steps)


#continuous histograms
plot_histogram(steps %>% select(sl_, ta_, gpp:power_hii))

#categorical bar plots
plot_bar(steps %>% select(gpp:calving_season))
plot_bar(steps %>% select(sex, gpp:calving_season), by="sex")

#full report with all DataExplorer metrics
#create_report(steps, y="case_")


################################ Scale continuous covariates #################################

steps_scaled <- steps %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric))

plot_histogram(steps_scaled %>% select(sl_, ta_, gpp:power_hii))

#plot_boxplot(steps_scaled %>% select(case_, sl_, ta_, gpp:power_hii), by = "case_")


################################ Identify highly correlated/redundant covariates #################################

#Correlation plot
plot_correlation(steps_scaled %>% 
                   select(sl_, ta_, gpp:calving_season))

#select covariates of interest and dummy code factors
dummy <- steps %>% select(sl_, ta_, gpp:calving_season) %>% 
  dummify()

#create a correlation matrix
cor_matrix <- cor(dummy)

#identify correlation values above 0.3
high_cor <- caret::findCorrelation(cor_matrix, cutoff = .3)

#convert identified indices to covariate names
dummy %>% select(!!!high_cor) %>% names()

#the GGAlly package corrplot is more informative for numeric variables
ggcorr(steps_scaled %>% 
         select(case_, sl_, ta_, gpp:calving_season), label = TRUE)

#VERY slow
ggpairs(steps_scaled %>% 
         select(case_, gpp:elevation), label = TRUE)



#not sure if this is useful or not
cor.prob <- function(X, dfr = nrow(X) - 2) {
  R <- cor(X, use="complete.obs")
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr / (1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  Rstar = ifelse(R[above]<0.05, "***", "NS")
  R[above]=paste(R[above],Rstar)
  R
}

cor_prob <- steps_scaled %>% 
  select(sl_, ta_, gpp:calving_season) %>%
  select(where(is.numeric)) %>% 
  cor.prob()



################################ Try out PCA for feature reduction #################################

plot_prcomp(steps_scaled %>% select(gpp:calving_season), nrow=1, ncol=2, parallel = TRUE)

#PC1: Slope, elevation, tri, and some veg
#PC2: Veg, landuse and pop density
#PC3: Seasonality


################################ Univariate Muff Models #################################
#subset dataframe to only include animal_id, case_, and all predictors
covs <- steps_scaled %>% 
  select(case_, animal_id, gpp:calving_season)

#function to fit a univariate Muff model with random effect for individual for each covariate
uni_fit <- function(cov){
  uni_mod <- glmmTMB(case_ ~ cov + (1|animal_id), family =
                       poisson, data = covs, doFit = FALSE)
  uni_mod$parameters$theta[1] <-log(1e3)
  uni_mod$mapArg <-list(theta=factor(1))
  fit <- glmmTMB::fitTMB(uni_mod)
  
  return(fit)
}

#map the function across all covariates
uni_fits <- map(covs, uni_fit)

#make table of AIC values for all univariate models
AICcmodavg::aictab(uni_fits, second.ord = FALSE)
AICcmodavg::bictab(uni_fits, second.ord = FALSE)

#top performing covariates for each element
#ndvi
#popdens_hii
#tree_cover_hansen
#northing

sjPlot::plot_model(elev_uni_fit, type="re",
                   transform=NULL,
                   axis.title = "Coefficient Estimates (Untransformed)",
                   title = "Responses to Elevation by Individual")

check_model(elev_uni_fit)



################################ Univariate SSF Models #################################

covs2 <- steps %>% 
  select(case_, step_id_, gpp:calving_season)
#function to fit a univariate issf model for each covariate
uni_fit_issf <- function(cov){
  uni_mod <- fit_issf(data=covs2, case_ ~ cov + strata(step_id_))
  return(uni_mod)
}

#map the function across all covariates
uni_fits_issf <- map(covs2, uni_fit_issf)

#make table of AIC values for all univariate models
AICcmodavg::aictab(uni_fits_issf, second.ord = FALSE)
AICcmodavg::bictab(uni_fits_issf, second.ord = FALSE)

map(uni_fits_issf, AIC) %>% 
  unlist() %>% 
  enframe(name = "model", value = "AIC")


################################ Dredging #################################

#Variables to keep based on univariate analysis:
#land_cover_usfs
#pop_dens_hii


#global model--has convergence issues
full_mod <- glmmTMB(case_ ~ gpp + npp + ndvi + evi + tree_cover_hansen + perc_tree_cover + perc_nontree_veg + perc_nonveg + land_cover_usfs + precip + dist_water + elevation + slope + northing + easting + tri +tpi + land_use_usfs + roads_hii + popdens_hii + landuse_hii + infra_hii + rails_hii + power_hii + season + hunting_season + calving_season + (1|animal_id), family = poisson, data = steps_scaled, doFit = FALSE)


full_mod$parameters$theta[1] <-log(1e3)
full_mod$mapArg <-list(theta=factor(1))
full_mod_fit <- glmmTMB::fitTMB(full_mod)


dredge_full <- dredge(full_mod, rank=BIC)


################################ Random Forest #################################

#not currently working. Dataset too large. works for around 10,000 points
covs2 <- steps %>% 
  select(animal_id, gpp:calving_season) %>% 
  dummify() %>% 
  mutate(case_=as.factor(covs$case_))

data_set_size <- floor(nrow(covs2)*0.80)
index <- sample(1:nrow(covs2), size = 5000)
training <- covs2[index,]
testing <- covs2[-index,]

rf <- randomForest(case_ ~ ., data = training, mtry = 4, ntree = 501, importance = TRUE)
rf
#basic plot
plot(rf)

#variable importance plot
importance(rf)
varImpPlot(rf)

#partial dependence plots
partialPlot(rf, pred.data = as.data.frame(training), x.var = ndvi)
partialPlot(rf, pred.data = as.data.frame(training), x.var = elevation)
partialPlot(rf, pred.data = as.data.frame(training), x.var = tree_cover_hansen)
partialPlot(rf, pred.data = as.data.frame(training), x.var = perc_tree_cover)
partialPlot(rf, pred.data = as.data.frame(training), x.var = roads_hii)


################################ XGBoost #################################
library(xgboost)
library(cvms)
library(caTools)

covs3 <- covs2 %>% select(-animal_id)

set.seed(42)

sample_split <- sample.split(Y = covs3$case_, SplitRatio = 0.7)
train_set <- subset(x = covs3, sample_split == TRUE)
test_set <- subset(x = covs3, sample_split == FALSE)

y_train <- as.integer(train_set$case_) - 1
y_test <- as.integer(test_set$case_) - 1
X_train <- train_set %>% select(-case_)
X_test <- test_set %>% select(-case_)


xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
xgb_params <- list(
  booster = "gbtree",
  eta = 0.01,
  max_depth = 8,
  gamma = 4,
  subsample = 0.75,
  colsample_bytree = 1,
  objective = "multi:softprob",
  eval_metric = "mlogloss",
  num_class = length(levels(covs2$case_))
)

xgb_model <- xgb.train(
  params = xgb_params,
  data = xgb_train,
  nrounds = 10, #5000
  verbose = 1
)
xgb_model

#variable importance scores
xgb_importance <- xgb.importance(model = xgb_model)

xgb.ggplot.importance(xgb_importance)




