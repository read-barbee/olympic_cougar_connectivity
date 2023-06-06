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


############################### Import Steps #################################

steps <- read_csv("data/Location_Data/Steps/6h_steps_unscaled_6-02-2023.csv")


steps <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation))

# update_columns(sex, dispersal_Status, case_, land_cover_usfs, land_use_usfs, season, hunting_season, calving_season, as.factor)


################################ Check covariate distributions #################################

#data overview
introduce(steps)
plot_intro(steps)

#check for missing values
plot_missing(steps)


#continuous histograms
plot_histogram(steps %>% select(sl_, ta_, gpp:power_hii, -c(aspect_deg, aspect_rad)))

#categorical bar plots
plot_bar(steps %>% select(gpp:calving_season))
plot_bar(steps %>% select(sex, gpp:calving_season), by="sex")

plot_correlation(steps %>% select(sl_, ta_, gpp:power_hii, -c(aspect_deg, aspect_rad)))


plot_prcomp(steps %>% select(gpp:calving_season))

plot_prcomp(steps %>% select(gpp:evi))

#create_report(steps, y="case_")

#caret::findCorrelation(steps)


################################ Muff Method #################################
#subset dataframe to only include animal_id, case_, and all predictors
covs <- steps %>% select(case_, animal_id, gpp:calving_season)

#function to fit a univariate Muff model with random effect for individual for each covariate
uni_fit <- function(cov){
  #cov_name <- names(steps %>% select(cov))
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
AICcmodavg::aictab(uni_fits)



sjPlot::plot_model(elev_uni_fit, type="re",
                   transform=NULL,
                   axis.title = "Coefficient Estimates (Untransformed)",
                   title = "Responses to Elevation by Individual")

check_model(elev_uni_fit)



#Global model








################################ iSSF Method #################################

elev_uni_issf <- fit_issf(steps, case_ ~ elevation + strata(step_id_))
summary(elev_uni_issf)

gpp_uni_issf <- fit_issf(steps, case_ ~ gpp + strata(step_id_))
summary(gpp_uni_issf)

npp_uni_issf <- fit_issf(steps, case_ ~ npp + strata(step_id_))
summary(nppuni_issf)

ndvi_uni_issf <- fit_issf(steps, case_ ~ ndvi + strata(step_id_))
summary(ndviuni_issf)

evi_uni_issf <- fit_issf(steps, case_ ~ evi + strata(step_id_))
summary(eviuni_issf)

#ndvi has lowest AIC, but GPP has highest PCA score
AIC(gpp_uni_issf)
AIC(npp_uni_issf)
AIC(ndvi_uni_issf)
AIC(evi_uni_issf)

tab_model(gpp_uni_issf, npp_uni_issf, ndvi_uni_issf, evi_uni_issf)








