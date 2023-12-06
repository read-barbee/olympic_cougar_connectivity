#### SSF Global Model Validation INLA####

# Author: Read Barbee

# Date:2023-11-22

# Purpose:


################################ Libraries #################################
library(tidyverse)
library(terra)
library(INLA)
library(INLAutils)
library(beepr)
library(sf)

################################ Helper Functions #################################
source("r_scripts/02_iSSF/08_cross_validation/k_fold_cv_function_INLA_v2_quad_11-27-23.R")

make_form_rand_quad <- function(params, quad_params, n_indiv){
  rand_terms_inla <- vector()
  for(i in 1:length(params)){
    rand_terms_inla[i] <- paste0("f(", paste0("id", i), ", ", params[i], ", ", "values = ", paste0("1:",n_indiv), ", ", 'model = "iid", ', "hyper = list(theta = list(initial = log(1), fixed = FALSE, ", 'prior = "pc.prec", ', "param = c(3, .05))))")
  }
  
  if (is.null(quad_params)){
    quad_terms <- NULL
  } else{
    quad_terms <- vector()
    for(i in 1:length(quad_params)){
      quad_terms[i] <- paste0(" + f(", paste0("id", (i + length(params)), ", ", "I(",quad_params[i],"^2)", ", ", "values = ", paste0("1:",n_indiv), ", ", 'model = "iid", ', "hyper = list(theta = list(initial = log(1), fixed = FALSE, ", 'prior = "pc.prec", ', "param = c(3, .05))))"))
    }
  }
  
  form <- as.formula(paste("case_ ~ -1 + ",  
                           #fixed effects
                           paste(c(params, quad_terms), collapse = " + "), "+",
                           #random intercept (strata)
                           "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T))) + ",
                           #random slopes
                           paste(rand_terms_inla, collapse = " + ")))
  
  return(form)
}

make_form <- function(params, quad_params, n_indiv){
  rand_terms_inla <- vector()
  for(i in 1:length(params)){
    rand_terms_inla[i] <- paste0("f(", paste0("id", i), ", ", params[i], ", ", "values = ", paste0("1:",n_indiv), ", ", 'model = "iid", ', "hyper = list(theta = list(initial = log(1), fixed = FALSE, ", 'prior = "pc.prec", ', "param = c(3, .05))))")
  }
  
  if (is.null(quad_params)){
    quad_terms <- NULL
  } else{
    quad_terms <- vector()
    for(i in 1:length(quad_params)){
      quad_terms[i] <- paste0("I(",quad_params[i],"^2)")
    }
  }
  
  form <- as.formula(paste("case_ ~ -1 + ",  
                           #fixed effects
                           paste(c(params, quad_terms), collapse = " + "), "+",
                           #random intercept (strata)
                           "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T))) + ",
                           #random slopes
                           paste(rand_terms_inla, collapse = " + ")))
  
  return(form)
}

make_sel_table <- function(fit_list, mod_names){
  mod_table <- list()
  for (i in 1:length(fit_list)){
    name <- mod_names[i]
    
    mod_table[[i]] <- tibble(name = name,
                             waic = fit_list[[i]]$waic$waic,
                             p.eff_waic = fit_list[[i]]$waic$p.eff,
                             dic = fit_list[[i]]$dic$dic,
                             p.eff_dic =fit_list[[i]]$dic$p.eff,
                             mean_cpo = mean(fit_list[[i]]$cpo$cpo),
                             mean_pit = mean(fit_list[[i]]$cpo$pit),
                             mlik = fit_list[[i]]$mlik[1])
    
    #print(i)
  }
  
  mod_table <- bind_rows(mod_table)
  
  return(mod_table)
}

################################ User-defined Parameters #################################
#prior params from https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40#inla-1
mean.beta <- 0
prec.beta <- 1e-4 

#significant quadratics
#ndvi2
#perc_nontree_veg2
#perc_tree_cov2
#popdens_hii2
#roads_hii2
#distance_water2
#landuse_hii2

#########################################################################
##
## 1. Import and format step data and cov stack
##
##########################################################################

#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_annual_cov_11-16-2023.csv")


#covariates
cov_stack <- rast("/Users/tb201494/Desktop/1km_buffer/cov_stack_pred_300m_11-21-2023.tif")

#cov_stack_pred <- rast("/Users/tb201494/Desktop/1km_buffer/cov_stack_pred_11-21-2023.tif")



#########################################################################
##
## 2. Scale continuous covariates for model comparison and correlation analysis
##
##########################################################################

#DON'T Scale step length and turn angle

steps_scaled <- steps %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dist_all_roads_annual), scale)) %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dist_all_roads_annual), as.numeric)) %>%
  rename(perc_tree_cov_annual = perc_tree_cov,
         perc_nontree_veg_annual = perc_nontree_veg,
         perc_nonveg_annual = perc_nonveg)

#DataExplorer::plot_histogram(steps_scaled %>% select(elevation:dist_all_roads_annual))

#convert response from categorical to numeric
mod_dat <- steps_scaled %>% mutate(case_ = as.numeric(case_))

#create separate animal_id columns for each random effect
for(i in 1:20){
  name <- as.name(paste0("id", i))
  mod_dat[[name]] <- as.numeric(factor(mod_dat$animal_id))
}

#get number of individuals from dataframe
n_indiv <-  mod_dat %>% distinct(animal_id) %>% count() %>% pull()

#########################################################################
##
## 3. Top model waic HII_combined
##
##########################################################################

params1 <- c("ndvi_annual",
             "tpi",
             "slope",
             "aspect_northness",
             "perc_nonveg_annual",
             "perc_tree_cov_annual",
             "hii_annual",
             "distance_water")

form1 <- make_form(params1, NULL, n_indiv)

top_waic_comb <- inla(form1, family ="Poisson", 
                      data=mod_dat,
                      control.fixed = list(
                      mean = mean.beta,
                      prec = list(default = prec.beta)),
                      control.compute = list(waic=TRUE, dic=TRUE, cpo = TRUE),
                      safe = TRUE)


#########################################################################
##
## 4. Top model waic HII_separate
##
##########################################################################
params2 <- c("ndvi_annual",
             "tpi",
             "slope",
             "aspect_northness",
             "perc_nonveg_annual",
             "perc_tree_cov_annual",
             "roads_hii_annual",
             "popdens_hii_annual",
             "distance_water")

form2 <- make_form(params2, NULL, n_indiv)

top_waic_sep <- inla(form2, family ="Poisson", 
                      data=mod_dat,
                      control.fixed = list(
                        mean = mean.beta,
                        prec = list(default = prec.beta)),
                      control.compute = list(waic=TRUE, dic=TRUE, cpo = TRUE),
                      safe = TRUE)

#########################################################################
##
## 5. Top model waic_rho HII_combined
##
##########################################################################
params3 <- c("ndvi_annual",
             "tpi",
             "slope",
             "perc_nonveg_annual",
             "perc_tree_cov_annual",
             "hii_annual")

form3 <- make_form(params3, NULL, n_indiv)

top_rho_comb <- inla(form3, family ="Poisson", 
                     data=mod_dat,
                     control.fixed = list(
                       mean = mean.beta,
                       prec = list(default = prec.beta)),
                     control.compute = list(waic=TRUE, dic=TRUE, cpo = TRUE),
                     safe = TRUE)

#########################################################################
##
## 6. Top model waic_rho HII_separate
##
##########################################################################
params4 <- c("ndvi_annual",
             "perc_nonveg_annual",
             "perc_tree_cov_annual",
             "roads_hii_annual",
             "popdens_hii_annual")

form4 <- make_form(params4, NULL, n_indiv)

top_rho_sep <- inla(form4, family ="Poisson", 
                     data=mod_dat,
                     control.fixed = list(
                       mean = mean.beta,
                       prec = list(default = prec.beta)),
                     control.compute = list(waic=TRUE, dic=TRUE, cpo = TRUE),
                     safe = TRUE)


#########################################################################
##
## 7. Top model waic HII_separate performs best. try adding quadratic
##
##########################################################################
params5 <- c("ndvi_annual",
             "tpi",
             "slope",
             "aspect_northness",
             "perc_nonveg_annual",
             "perc_tree_cov_annual",
             "roads_hii_annual",
             "popdens_hii_annual",
             "distance_water")

quad_terms5 <- c("ndvi_annual")

form5 <- make_form_rand_quad(params5, quad_terms5, n_indiv)

form5_no_rand <- make_form(params5, quad_terms5, n_indiv)

ndvi_quad <- inla(form5, family ="Poisson", 
                     data=mod_dat,
                     control.fixed = list(
                       mean = mean.beta,
                       prec = list(default = prec.beta)),
                     control.compute = list(waic=TRUE, dic=TRUE, cpo = TRUE),
                     safe = TRUE)

# ndvi_quad_no_rand <- inla(form5_no_rand, family ="Poisson", 
#                   data=mod_dat,
#                   control.fixed = list(
#                     mean = mean.beta,
#                     prec = list(default = prec.beta)),
#                   control.compute = list(waic=TRUE, dic=TRUE, cpo = TRUE),
#                   safe = TRUE)

#########################################################################
##
## 8. Add popdens quadratic
##
##########################################################################
params6 <- c("ndvi_annual",
             "tpi",
             "slope",
             "aspect_northness",
             "perc_nonveg_annual",
             "perc_tree_cov_annual",
             "roads_hii_annual",
             "popdens_hii_annual",
             "distance_water")

quad_terms6 <- c("popdens_hii_annual")

form6 <- make_form_rand_quad(params6, quad_terms6, n_indiv)

popdens_quad <- inla(form6, family ="Poisson", 
                  data=mod_dat,
                  control.fixed = list(
                    mean = mean.beta,
                    prec = list(default = prec.beta)),
                  control.compute = list(waic=TRUE, dic=TRUE, cpo = TRUE),
                  safe = TRUE)


#########################################################################
##
## 9. Add tree_cov quadratic
##
##########################################################################
params7 <- c("ndvi_annual",
             "tpi",
             "slope",
             "aspect_northness",
             "perc_nonveg_annual",
             "perc_tree_cov_annual",
             "roads_hii_annual",
             "popdens_hii_annual",
             "distance_water")

quad_terms7 <- c("perc_tree_cov_annual") #"ndvi_annual", "popdens_hii_annual"

form7 <- make_form_rand_quad(params7, quad_terms7, n_indiv)

tree_cov_quad <- inla(form7, family ="Poisson", 
                          data=mod_dat,
                          control.fixed = list(
                            mean = mean.beta,
                            prec = list(default = prec.beta)),
                          control.compute = list(waic=TRUE, dic=TRUE, cpo = TRUE),
                          safe = TRUE)

#cv7 <- k_fold_inla(mod_dat, form7, cov_stack, 4)

#4-fold mean rho: 0.5333334


#########################################################################
##
## 9. Ndvi and tree_cov quadratics
##
##########################################################################
params8 <- c("ndvi_annual",
             "tpi",
             "slope",
             "aspect_northness",
             "perc_nonveg_annual",
             "perc_tree_cov_annual",
             "roads_hii_annual",
             "popdens_hii_annual",
             "distance_water")

quad_terms8 <- c("ndvi_annual","perc_tree_cov_annual") 

form8 <- make_form_rand_quad(params8, quad_terms8, n_indiv)

ndvi_tree_cov_quad <- inla(form8, family ="Poisson", 
                      data=mod_dat,
                      control.fixed = list(
                        mean = mean.beta,
                        prec = list(default = prec.beta)),
                      control.compute = list(waic=TRUE, dic=TRUE, cpo = TRUE),
                      safe = TRUE)


#########################################################################
##
## 9. Ndvi and tree_cov quadratics
##
##########################################################################
params10 <- c("ndvi_annual",
             "tpi",
             "slope",
             "aspect_northness",
             "perc_nonveg_annual",
             "perc_tree_cov_annual",
             "roads_hii_annual",
             "popdens_hii_annual",
             "distance_water")

quad_terms10 <- c("ndvi_annual","perc_tree_cov_annual", "popdens_hii_annual") 

form10 <- make_form_rand_quad(params10, quad_terms10, n_indiv)

all_quad <- inla(form10, family ="Poisson", 
                      data=mod_dat,
                      control.fixed = list(
                        mean = mean.beta,
                        prec = list(default = prec.beta)),
                      control.compute = list(waic=TRUE, dic=TRUE, cpo = TRUE),
                      safe = TRUE)


#########################################################################
##
## 9. no ndvi
##
##########################################################################
params9 <- c(#"ndvi_annual",
             "tpi",
             "slope",
             "aspect_northness",
             "perc_nonveg_annual",
             "perc_tree_cov_annual",
             "roads_hii_annual",
             "popdens_hii_annual",
             "distance_water")

quad_terms9 <- c("perc_tree_cov_annual") #"ndvi_annual", "popdens_hii_annual"

form9 <- make_form_rand_quad(params9, quad_terms9, n_indiv)

no_ndvi <- inla(form9, family ="Poisson", 
                      data=mod_dat,
                      control.fixed = list(
                        mean = mean.beta,
                        prec = list(default = prec.beta)),
                      control.compute = list(waic=TRUE, dic=TRUE, cpo = TRUE),
                      safe = TRUE)


#########################################################################
##
## 9. Model Diagnostics
##
##########################################################################

# fit_list <- list(top_waic_comb, top_waic_sep, top_rho_comb, top_rho_sep, ndvi_quad, ndvi_popdens_quad, all_quad)
# mod_names <- c("top_waic_comb", "top_waic_sep", "top_rho_comb", "top_rho_sep", "ndvi_quad", "ndvi_popdens_quad", "all_quad")

fit_list <- list(top_waic_sep, ndvi_quad, tree_cov_quad, popdens_quad, ndvi_tree_cov_quad, ndvi_popdens_quad, all_quad)
mod_names <- c("top_linear", "ndvi_quad", "tree_cov_quad", "popdens_quad", "ndvi_tree_cov_quad", "ndvi_popdens_quad", "all_quad")


mod_table <- make_sel_table(fit_list, mod_names)

saveRDS(ndvi_quad, "top_ssf_ndvi_quad_12_01_23.rds")

#~ 6 min 4 fold cv = 0.7; 10-fold cv = 0.68 ~ 12 minutes
system.time(top_linear_cv <- k_fold_inla(mod_dat, form2, cov_stack, n_folds=10))

# 4 fold didn't converge. 10 fold cv = -0.7236364
system.time(ndvi_quad_cv <- k_fold_inla(mod_dat, form5, cov_stack, n_folds=10))

#10 fold cv = -0.7236364 ~15 min
system.time(all_quad_cv <- k_fold_inla(mod_dat, form7, cov_stack, n_folds=10))

sum(unlist(all_quad_cv[1:10]))/10

