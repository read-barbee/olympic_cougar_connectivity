#### Global Resident Muff Model ####

# Author: Read Barbee

# Date:2023-09-18

# Purpose:


################################ libraries #################################
library(tidyverse)
#library(glmmTMB)
#library(amt)
##library(MuMIn)
#library(INLA)
library(beepr)
#library(buildmer)
library(INLA)
library(INLAutils)
library(terra)
library(sf)

source("r_scripts/02_iSSF/08_cross_validation/k_fold_cv_function_INLA_11-20-23.R")
#########################################################################
##
##  Specify Model Parameters
##
##########################################################################

#waic selected params
params <- c("ndvi_annual",
            "tpi" ,
            #"elevation",
            "slope",
            "aspect_northness",
            "perc_nonveg_annual",
            "perc_tree_cov_annual",
            "roads_hii_annual",
            "popdens_hii_annual",
            "distance_water",
            "landuse_hii_annual",
            "infra_hii_annual",
            "aspect_eastness"
            )

# no quads, tpi. waic = 169363.51 rho = 0.9 4_fold_rho: 0.830303
#with ndvi quad, waic = 198116.98, rho = -0.8389097 
#with perc_tree_cov quad: waic = 193730.79, rho = -0.9636364. Also screws things up
#popdens quadratic doesn't converge. Crashes INLA


quad_params <- c(#"ndvi_annual"#,
                 #"perc_tree_cov_annual"
                #"popdens_hii_annual"
                )
                 


#effect size selected params
params <- c("ndvi_annual",
            "aspect_northness",
            "slope",
            "distance_water",
            "hii_annual",
            #mtpi,
            "tpi" ,
            "aspect_eastness",
            "npp_annual",
            "perc_tree_cov_annual",
            "perc_nonveg_annual"
)

#all rho: 0.1272727, waic = 221854.63
#no npp: rho 0.2 , waic = 162647.80
#all with tpi instead of mtpi: rho = 0.5151515, waic = 189518.19

quad_params <- c(#"ndvi_annual"#,
  #"perc_tree_cov_annual"
  #"popdens_hii_annual"
)




#prior params from https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40#inla-1
mean.beta <- 0
prec.beta <- 1e-4 


rand_terms <- vector()

for(i in 1:length(params)){
  rand_terms[i] <- paste0("(0 + ", params[i], " | animal_id)")
}


quad_terms <- vector()

for(i in 1:length(quad_params)){
  quad_terms[i] <- paste0("I(",quad_params[i],"^2)")
}

quad_terms <- NULL

#rand_terms_quad <- vector()

# for(i in 1:length(quad_terms)){
#   rand_terms_quad[i] <- paste0("(0 + ", quad_terms[i], " | animal_id)")
# }

########################################################################
##
## 1. Import and format step data
##
########################################################################

#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_annual_cov_11-16-2023.csv")

#DataExplorer::plot_histogram(steps %>% filter(case_==TRUE) %>% select(elevation:dist_all_roads_annual))

#scatterplots. not super useful
# DataExplorer::plot_scatterplot(steps %>% select(case_, elevation:dist_all_roads_annual), by= "case_")


#set all negative elevations to 0 and filter dispersal tracks for post dispersal event
steps_scaled <- steps %>% 
  #mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  #select(-c(land_use_usfs, land_cover_usfs, dispersing:disp_qual, season:calving_season)) %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dist_all_roads_annual), scale)) %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dist_all_roads_annual), as.numeric)) %>%
  rename(perc_tree_cov_annual = perc_tree_cov,
         perc_nontree_veg_annual = perc_nontree_veg,
         perc_nonveg_annual = perc_nonveg)



#########################################################################
##
## 2. Fit global model with quadratics in INLA
##
##########################################################################


n_indiv = steps_scaled %>% distinct(animal_id) %>% count() %>% pull()

rand_terms_inla <- vector()
for(i in 1:length(params)){
  rand_terms_inla[i] <- paste0("f(", paste0("id", i), ", ", params[i], ", ", "values = ", paste0("1:",n_indiv), ", ", 'model = "iid", ', "hyper = list(theta = list(initial = log(1), fixed = FALSE, ", 'prior = "pc.prec", ', "param = c(3, .05))))")
}


#convert response from categorical to numeric
dat <- steps_scaled %>% mutate(case_ = as.numeric(case_))

#create separate animal_id columns for each random effect
for(i in 1:length(rand_terms_inla)){
  name <- as.name(paste0("id", i))
  dat[[name]] <- as.numeric(factor(dat$animal_id))
}


#construct the model formula with quadratic terms
form <- as.formula(paste("case_ ~ -1 + ",  
                         #fixed effects
                         paste(c(params, quad_terms), collapse = " + "), "+",
                         #random intercept (strata)
                         "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T))) + ",
                         #random slopes
                         paste(rand_terms_inla, collapse = " + ")))


cov_stack <- rast("/Users/tb201494/Desktop/1km_buffer/cov_stack_pred_11-20-2023.tif")

mtpi <- rast("/Users/tb201494/Desktop/1km_buffer/static/mtpi_1km_buffer.tif")

cov_stack$mtpi <- mtpi

cov_stack <- cov_stack[[c(1:5, 32, 6:31 )]]

cov_stack <- terra::aggregate(cov_stack, fact = 10, fun = "mean", cores = 3)

#run_global_inla <- function (dat){
  
  #convert response from categorical to numeric
  dat <- dat %>% mutate(case_ = as.numeric(case_))
  
  #create separate animal_id columns for each random effect
  for(i in 1:length(rand_terms_inla)){
    name <- as.name(paste0("id", i))
    dat[[name]] <- as.numeric(factor(dat$animal_id))
  }
  
  #construct the model formula with quadratic terms
  form <- as.formula(paste("case_ ~ -1 + ",  
                           #fixed effects
                           paste(c(params, quad_terms), collapse = " + "), "+",
                           #random intercept (strata)
                           "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T))) + ",
                           #random slopes
                           paste(rand_terms_inla, collapse = " + ")))
  
  #fit the model    
  fit <-  inla(form, family ="Poisson", data=dat,
               control.fixed = list(
               mean = mean.beta,
               prec = list(default = prec.beta)),
               control.compute = list(waic=TRUE, dic = TRUE, cpo = FALSE)) #
  
  return(fit)
}

#15 min for 4 folds and 12 covariates
system.time(cv2 <- k_fold_inla(dat, form, cov_stack, n_folds = 4))


