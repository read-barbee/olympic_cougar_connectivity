#### Global Resident Muff Model ####

# Author: Read Barbee

# Date:2023-09-18

# Purpose:


################################ libraries #################################
library(tidyverse)
library(glmmTMB)
library(amt)
#library(MuMIn)
#library(INLA)
library(beepr)
library(buildmer)
library(INLA)

#########################################################################
##
##  Specify Model Parameters
##
##########################################################################

params <- c("roads_hii",
            "popdens_hii",
            "landuse_hii",
            "infra_hii", 
            "rails_hii", 
            "power_hii",
            "tpi",
            "npp",
            "perc_tree_cover",
            "ndvi",
            "tree_cover_hansen",
            "northing",
            "easting",
            "tri",
            "precip")

quad_params <- c("roads_hii",
                "tpi",
                "npp",
                "perc_tree_cover",
                "ndvi",
                "tree_cover_hansen",
                "northing",
                "easting",
                "tri",
                "precip")


rand_terms <- vector()

for(i in 1:length(params)){
  rand_terms[i] <- paste0("(0 + ", params[i], " | animal_id)")
}


quad_terms <- vector()

for(i in 1:length(quad_params)){
  quad_terms[i] <- paste0("I(",quad_params[i],"^2)")
}


#rand_terms_quad <- vector()

# for(i in 1:length(quad_terms)){
#   rand_terms_quad[i] <- paste0("(0 + ", quad_terms[i], " | animal_id)")
# }

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
  filter(dispersal_status=="resident") %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad)) %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric))

#########################################################################
##
## 2. Fit global model without quadratics
##
##########################################################################

#Fit global glmmTMB model: takes about ~ 15 min for 15 covs
run_global <- function (dat){
  #library(survival)  survival::clogit
  form <- as.formula(paste0("case_ ~ -1 + ", #remove standard intercept to be replace with stratum-based intercept
                            #fixed effects
                            paste(params, collapse = " + "),
                            #random intercept (strata)
                            " + (1|step_id_) +",
                            #random slopes
                            paste(rand_terms, collapse = " + ")))
  
  
  mod <- glmmTMB(form, family =poisson, data = dat, doFit = FALSE)
  
  mod$parameters$theta[1] <-log(1e3)
  mod$mapArg <-list(theta=factor(c(NA, 1:length(rand_terms))))
  fit <- glmmTMB::fitTMB(mod)
return(fit)
}


system.time(global_fit <- run_global(steps_scaled))


summary(global_fit)

#saveRDS(global_fit, "muff_global_fit_res_no_quad_9-19-23.rds")

#global_fit <- readRDS("muff_global_fit_res_no_quad_9-19-23.rds")


#########################################################################
##
## 2. Fit global model with quadratics
##
##########################################################################

run_global_quad <- function (dat){
  #library(survival)  survival::clogit
  form <- as.formula(paste0("case_ ~ -1 + ", #remove standard intercept to be replace with stratum-based intercept
                            #fixed effects
                            paste(c(params, quad_terms), collapse = " + "),
                            #random intercept (strata)
                            " + (1|step_id_) +",
                            #random slopes
                            paste(rand_terms, collapse = " + ")))
  
  
  mod <- glmmTMB(form, family =poisson, data = dat, doFit = FALSE)
  
  mod$parameters$theta[1] <-log(1e3)
  mod$mapArg <-list(theta=factor(c(NA, 1:length(rand_terms))))
  fit <- glmmTMB::fitTMB(mod)
  return(fit)
}


#Fit takes ~ 22 min with 15 linear terms and 10 quadratic terms; no random slopes on quads
system.time(global_fit_quad <- run_global_quad(steps_scaled))

#saveRDS(global_fit_quad, "muff_global_fit_quad_res_9-19-23.rds")

global_fit_quad <- readRDS("muff_global_fit_quad_res_9-19-23.rds")


#significant quadratic terms: ndvi, perc_tree_cover.

#including quadratic terms improves AIC

#########################################################################
##
## 2. Fit global model with quadratics in INLA
##
##########################################################################

rand_terms_inla <- vector()
n_indiv = steps_scaled %>% distinct(animal_id) %>% count() %>% pull()

for(i in 1:length(params)){
  rand_terms_inla[i] <- paste0("f(", paste0("id", i), ", ", params[i], ", ", "values = ", paste0("1:",n_indiv), ", ", 'model = "iid", ', "hyper = list(theta = list(initial = log(1), fixed = FALSE, ", 'prior = "pc.prec", ', "param = c(1, .05))))")
}



run_global_inla <- function (dat){
  
  #convert response from categorical to numeric
  dat <- dat %>% mutate(case_ = as.numeric(case_))
  
  #create separate animal_id columns for each random effect
  for(i in 1:length(rand_terms_inla)){
    name <- as.name(paste0("id", i))
    dat[[name]] <- as.numeric(factor(dat$animal_id))
  }
  
  #construct the model formula with quadratic terms
  form <- as.formula(paste("case_ ~ ",  
                           #fixed effects
                           paste(c(params, quad_terms), collapse = " + "), "+",
                           #random intercept (strata)
                           "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T))) + ",
                           #random slopes
                           paste(rand_terms_inla, collapse = " + ")))
  
  #fit the model    
  fit <- inla(form, family ="poisson", data=dat)
  
  return(fit)
}

system.time(inla_fit <- run_global_inla(steps_scaled))


#########################################################################
##
## 2. LOO Cross-Validation in INLA
##
##########################################################################

inla_cv <- inla.group.cv(inla_fit)
ULOOCV = mean(inla_cv$cv)

#control.fixed = list(
# mean = mean.beta,
# prec = list(default = prec.beta)), control.compute = list(dic = TRUE)

#########################################################################
##
## 2. Remove non-significant quadratic terms 
##
##########################################################################
quad_terms2 <- c("I(roads_hii^2)",
                 "I(ndvi^2)",
                 "I(perc_tree_cover^2)")

# rand_terms_quad2 <- c("(0 + I(ndvi^2)|animal_id",
#                       "(0 + I(perc_tree_cover^2) | animal_id")

form <- as.formula(paste0("case_ ~ -1 + ", #remove standard intercept to be replace with stratum-based intercept
                          #fixed effects
                          paste(c(params, quad_terms2), collapse = " + "),
                          #random intercept (strata)
                          " + (1|step_id_) +",
                          #random slopes
                          paste(c(rand_terms), collapse = " + ")))


mod <- glmmTMB(form, family =poisson, data = steps_scaled, doFit = FALSE)

mod$parameters$theta[1] <-log(1e3)
mod$mapArg <-list(theta=factor(c(NA, 1:length(rand_terms))))
system.time(global_fit_quad2 <- glmmTMB::fitTMB(mod))

#saveRDS(global_fit_quad2, "muff_global_fit_quad2_res_9-19-23.rds")

#global_fit_quad2 <- readRDS("muff_global_fit_quad2_res_9-19-23.rds")



#########################################################################
##
## 3. Automated stepwise (buildmer)
##
##########################################################################

form <- global_fit_quad2$modelInfo$allForm$formula

test <- buildmer(f = form, data = steps_scaled, buildmerControl = )








#########################################################################
##
## 3. Manual_backward_stepwise
##
##########################################################################

#find term with highest p-value to remove
start_mod <- summary(stepwise4)
tab <- start_mod$coefficients$cond %>% 
  as.data.frame() %>% 
  #filter(`Pr(>|z|)` > 0.05) %>% 
  arrange(desc(`Pr(>|z|)`))

tab2 <- start_mod$coefficients$cond %>% 
  as.data.frame() %>% 
  #filter(`Pr(>|z|)` > 0.05) %>% 
  arrange(abs(Estimate))


form <- as.formula(case_ ~ -1 + 
                     roads_hii +  I(roads_hii^2) +
                     popdens_hii + 
                     #landuse_hii + 
                     infra_hii +  
                     #rails_hii + 
                     #power_hii + 
                     tpi + 
                     npp + 
                     perc_tree_cover + I(perc_tree_cover^2) + 
                     ndvi +  I(ndvi^2) +
                     tree_cover_hansen + 
                     northing + 
                     easting + 
                     tri + 
                     precip  +  
                    (1 | step_id_) + 
                     (0 +  roads_hii | animal_id) + 
                     (0 + popdens_hii | animal_id) +  
                     #(0 + landuse_hii | animal_id) + 
                     (0 + infra_hii | animal_id) +  
                     #(0 + rails_hii | animal_id) + 
                     #(0 + power_hii | animal_id) +  
                     (0 + tpi | animal_id) + 
                     (0 + npp | animal_id) + 
                     (0 + perc_tree_cover | animal_id) + 
                     (0 + ndvi | animal_id) + 
                     (0 + tree_cover_hansen | animal_id) + 
                     (0 + northing | animal_id) + 
                     (0 + easting | animal_id) + 
                     (0 + tri | animal_id) + 
                     (0 + precip | animal_id))


mod <- glmmTMB(form, family =poisson, data = steps_scaled, doFit = FALSE)

mod$parameters$theta[1] <-log(1e3)
mod$mapArg <-list(theta=factor(c(NA, 1:12)))
system.time(stepwise4 <- glmmTMB::fitTMB(mod))

AIC(stepwise3) < AIC(stepwise4)

#top: stepwise4: AIC 274869

#saveRDS(stepwise4, "muff_top_global_fit_res_9-19-23.rds")

#top_mod <- readRDS("muff_top_global_fit_res_9-19-23.rds")

#removing landuse makes it worse

















## Removing nonsignificant linear terms increases AIC by 711 points

#Remove non-significant terms
params_sig <- c("popdens_hii",
                "landuse_hii",
                "infra_hii",
                "npp",
                "ndvi",
                "tree_cover_hansen",
                "northing",
                "tri",
                "precip")

rand_terms_sig <- vector()

for(i in 1:length(params_sig)){
  rand_terms_sig[i] <- paste0(params_sig[i], " + (0 + ", params_sig[i], " | animal_id)")
}


#fit model (~5 min)
system.time(global_fit_sig <- run_global(steps_scaled, params_sig, rand_terms_sig))

summary(global_fit_sig)


#########################################################################
##
## 3. Manual dredge (glmmTMB) -- CRASHES with 15 covs
##
##########################################################################
vars <- c("roads_hii + I(roads_hii^2) + (0 + roads_hii | animal_id)",
          "perc_tree_cover + I(perc_tree_cover^2) + (0 + perc_tree_cover | animal_id)",
          "ndvi + I(ndvi^2) + (0 + ndvi | animal_id)",
          "popdens_hii + (0 + popdens_hii | animal_id)",
          "landuse_hii + (0 + landuse_hii | animal_id)",
          "infra_hii + (0 + infra_hii | animal_id)", 
          "rails_hii + (0 + rails_hii | animal_id)", 
          "power_hii + (0 + power_hii | animal_id)",
          "tpi + (0 + tpi | animal_id)",
          "npp + (0 + npp | animal_id)",
          "tree_cover_hansen + (0 + tree_cover_hansen | animal_id)",
          "northing + (0 + northing | animal_id)",
          "easting + (0 + easting | animal_id)",
          "tri + (0 + tri | animal_id)",
          "precip + (0 + precip | animal_id)")


all_comb <- do.call("c", lapply(seq_along(vars), function(i) combn(vars, i, FUN = list)))


#~ 10 min for 2,097,151 combinations of 21 covariates
forms <- list()

#CRASHES after a couple of hours
for (i in 1:length(all_comb)){
  var_i <- all_comb[[i]]
  forms[[i]] <- as.formula(paste("case_",  paste("-1", "(1|step_id_)", paste(var_i, collapse="+"), sep="+"), sep="~"))
}

mods <- list()
system.time(for (i in 1:length(forms)){
  form <- forms[[i]]
  mod <- glmmTMB(form, family=poisson, data=steps_scaled, doFit=FALSE)
  mod$parameters$theta[1] <- log(1e3)
  map_length <- 1:(length(mod$parameters$theta)-1)
  mod$mapArg <- list(theta=factor(c(NA, map_length)))
  
  mods[[i]] <- fitTMB(mod)
})












#########################################################################
##
## 4. Backward stepwise (glmmTMB)
##
##########################################################################




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
