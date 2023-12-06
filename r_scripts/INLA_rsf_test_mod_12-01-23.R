
library(tidyverse)
library(INLA)

dat <- read_csv("rsf_used_avail_static_covs_five_to_one_12-01-23.csv")

dat2 <- dat %>% 
  mutate(across(elevation:distance_water, scale)) %>% 
  mutate(across(elevation:distance_water, as.numeric)) %>% 
  mutate(case_=as.numeric(case_),
         year = as.numeric(as.factor(year)))

form <- case_ ~ elevation + f(year, model='iid')
test_mod <- inla(form, family = "logistic", data = dat2, control.compute=list(dic=TRUE, 
                                                                  cpo=TRUE, 
                                                                  waic=TRUE, 
                                                                  return.marginals.predictor=TRUE))

#probably don't need random effects in RSF because I assume the probability of use and the relationship between each covariate and probability of use doesn't vary from year to year