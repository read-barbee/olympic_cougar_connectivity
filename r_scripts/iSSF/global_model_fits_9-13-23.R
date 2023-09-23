#### Global Model Fitting ####

# Author: Read Barbee

# Date:2023-09-13 

# Purpose:


################################ libraries #################################
library(tidyverse)
library(glmmTMB)
library(amt)
#library(MuMIn)
#library(INLA)
library(beepr)

#########################################################################
##
## 1. Import and format step data
##
##########################################################################

#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv") 

#%>% mutate(ndvi = ndvi*0.0001)

#set all negative elevations to 0 and filter dispersal tracks for post dispersal event
steps_scaled <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad)) %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric)) %>% 
  filter(is.na(dispersing) | dispersing==TRUE)

#########################################################################
##
## 2. Fit global model and auto-dredge (not working)
##
##########################################################################

#steps_scaled_no_na <- steps_scaled %>% na.omit()

#Fit global glmmTMB model: takes about ~ 17 min
global <- glmmTMB(case_ ~ -1 +
                    #fixed effects
                    tree_cover_hansen + I(tree_cover_hansen^2) +
                    gpp + I(gpp^2) +
                    perc_tree_cover + I(perc_tree_cover^2) +
                    landuse_hii + I(landuse_hii^2) +
                    ndvi + I(ndvi^2) +
                    popdens_hii + I(popdens_hii^2) +
                    northing +
                    easting + 
                    rails_hii +
                    infra_hii +
                    #stratum-based intercept
                    (1|step_id_) +
                    #random slopes
                    (0 + tree_cover_hansen | animal_id) +
                    (0 + gpp | animal_id) +
                    (0 + perc_tree_cover | animal_id) +
                    (0 + landuse_hii | animal_id) +
                    (0 + ndvi | animal_id) +
                    (0 + popdens_hii | animal_id) +
                    (0 + northing | animal_id) +
                    (0 + easting | animal_id) +
                    (0 + rails_hii | animal_id) +
                    (0 + infra_hii | animal_id),
                  family=poisson,
                  data = steps_scaled,
                  doFit=FALSE); 
global$parameters$theta[1] <- log(1e3)
global$mapArg <- list(theta=factor(c(NA, 1:10)))

system.time(fit_quad <- fitTMB(global))#; beep("fanfare")

saveRDS(fit_quad, "muff_global_fit_9-13-23.rds")

fit_quad <- readRDS("muff_global_fit_9-13-23.rds")


#Fit global glmmTMB model: takes about ~ 17 min
global_no_quad <- glmmTMB(case_ ~ -1 +
                    #fixed effects
                    tree_cover_hansen +
                    gpp + 
                    perc_tree_cover + 
                    landuse_hii + 
                    ndvi + 
                    popdens_hii + 
                    northing +
                    easting + 
                    rails_hii +
                    infra_hii +
                    #stratum-based intercept
                    (1|step_id_) +
                    #random slopes
                    (0 + tree_cover_hansen | animal_id) +
                    (0 + gpp | animal_id) +
                    (0 + perc_tree_cover | animal_id) +
                    (0 + landuse_hii | animal_id) +
                    (0 + ndvi | animal_id) +
                    (0 + popdens_hii | animal_id) +
                    (0 + northing | animal_id) +
                    (0 + easting | animal_id) +
                    (0 + rails_hii | animal_id) +
                    (0 + infra_hii | animal_id),
                  family=poisson,
                  data = steps_scaled,
                  doFit=FALSE); 
global_no_quad$parameters$theta[1] <- log(1e3)
global_no_quad$mapArg <- list(theta=factor(c(NA, 1:10)))

system.time(fit_no_quad <- fitTMB(global_no_quad))#; beep("fanfare")

saveRDS(fit_no_quad, "muff_global_fit_no_quad_9-13-23.rds")

fit_no_quad <- readRDS("muff_global_fit_no_quad_9-13-23.rds")


#Fit global glmmTMB model: takes about ~ 18 min
global_aic_mod <- glmmTMB(case_ ~ -1 +
                            #fixed effects
                            roads_hii +
                            popdens_hii + 
                            landuse_hii +
                            infra_hii + 
                            rails_hii + 
                            power_hii +
                            tpi +
                            perc_nonveg + 
                            perc_tree_cover +
                            npp +
                            ndvi +
                            tree_cover_hansen +
                            northing +
                            easting +
                            slope +
                            precip +
                            #stratum-based intercept
                            (1|step_id_) +
                            #random slopes
                            (0 + roads_hii | animal_id) +
                            (0 + popdens_hii | animal_id) +
                            (0 + landuse_hii | animal_id) +
                            (0 + infra_hii | animal_id) +
                            (0 + rails_hii | animal_id) +
                            (0 + power_hii | animal_id) +
                            (0 + tpi | animal_id) +
                            (0 + perc_nonveg | animal_id) +
                            (0 + perc_tree_cover | animal_id) +
                            (0 + npp | animal_id) +
                            (0 + ndvi | animal_id) +
                            (0 + tree_cover_hansen | animal_id) +
                            (0 + northing | animal_id) +
                            (0 + easting | animal_id) +
                            (0 + slope | animal_id) +
                            (0 + precip | animal_id),
                          family=poisson,
                          data = steps_scaled,
                          doFit=FALSE); 
global_aic_mod$parameters$theta[1] <- log(1e3)
global_aic_mod$mapArg <- list(theta=factor(c(NA, 1:16)))

system.time(fit_aic_mod <- fitTMB(global_aic_mod))#; beep("fanfare")

saveRDS(fit_aic_mod, "muff_global_fit_aic_mod_9-13-23.rds")

fit_aic_mod <- readRDS("muff_global_fit_aic_mod_9-13-23.rds")





#Null Model
null <- glmmTMB(case_ ~ -1 + 
                  #stratum-based intercept
                  (1|step_id_),
                family=poisson,
                data = steps_scaled,
                doFit=FALSE); 
null$parameters$theta[1] <- log(1e3)
null$mapArg <- list(theta=factor(c(NA)))

system.time(null_fit <- fitTMB(null))#; beep("fanfare")


#Fit global issf to each individual 

run_global <- function (dat){
  #library(survival)  survival::clogit
  mod <- amt::fit_issf(formula = case_ ~ 
                         tree_cover_hansen + I(tree_cover_hansen^2) +
                         gpp + I(gpp^2) +
                         perc_tree_cover + I(perc_tree_cover^2) +
                         landuse_hii + I(landuse_hii^2) +
                         ndvi + I(ndvi^2) +
                         popdens_hii + I(popdens_hii^2) +
                         northing +
                         easting + 
                         rails_hii +
                         infra_hii +
                         #steps
                         sl_ + log_sl_ + cos_ta_ +
                         #strata
                         strata(step_id_),
                       data = dat,
                       na.action = "na.omit",
                       model = TRUE)
  return(mod)
}

run_global_no_quad <- function (dat){
  #library(survival)  survival::clogit
  mod <- amt::fit_issf(formula = case_ ~ 
                         tree_cover_hansen + 
                         gpp + 
                         perc_tree_cover +
                         landuse_hii + 
                         ndvi + 
                         popdens_hii +
                         northing +
                         easting + 
                         rails_hii +
                         infra_hii +
                         #steps
                         sl_ + log_sl_ + cos_ta_ +
                         #strata
                         strata(step_id_),
                       data = dat,
                       na.action = "na.omit",
                       model = TRUE)
  return(mod)
}


run_global_aic_mod <- function (dat){
  #library(survival)  survival::clogit
  mod <- amt::fit_issf(formula = case_ ~ 
                         roads_hii +
                         popdens_hii + 
                         landuse_hii +
                         infra_hii + 
                         rails_hii + 
                         power_hii +
                         tpi +
                         perc_nonveg + 
                         perc_tree_cover +
                         npp +
                         ndvi +
                         tree_cover_hansen +
                         northing +
                         easting +
                         slope +
                         precip +
                         #steps
                         sl_ + log_sl_ + cos_ta_ +
                         #strata
                         strata(step_id_),
                       data = dat,
                       na.action = "na.omit",
                       model = TRUE)
  return(mod)
}

#########################################################################
##
## 1. Fit global iSSF model to each individual
##
##########################################################################

#50 convergence issues with imputed set for clogit. 19 warnings for fit_issf
global_fits <- steps %>%  
  nest_by(animal_id, sex, dispersal_status) %>% 
  rename(steps=data) %>% 
  pull(steps) %>% 
  map(run_global)

global_fits_no_quad <- steps %>%  
  nest_by(animal_id, sex, dispersal_status) %>% 
  rename(steps=data) %>% 
  pull(steps) %>% 
  map(run_global_no_quad)

steps$global_fit <- global_fits


test <- global_fits[[1]]
test2 <- global_fits_no_quad[[1]]
update_sl_distr(test2)

indiv_coeffs <- bind_rows(global_fits2 %>% map(function(dat) dat$model$coefficients))

indiv_coeffs %>% pivot_longer(cols = everything(), names_to = "cov", values_to = "value") %>%  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~cov, scales = "free")

start <- make_start(c(steps[1, ]$x1_))
k1 <- redistribution_kernel(m_0, map = env, start = start)
k1 <- redistribution_kernel(m_0, map = env, start = start,
                            landscape = "continuous", tolerance.outside = 0.2, 
                            n.control = 1e4)

hist(indiv_coeffs$tree_cover_hansen)
hist(exp(indiv_coeffs$gpp))
