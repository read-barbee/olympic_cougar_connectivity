
### NOTE: requires rsf surface for each model rather than incorporating estimated model coefficients?


#1. Select individual set of coefficients from top Muff model
#2. Create raster surface by applying model to landscape
#3. Simulate tracks for that individual
#4. Repeat for each individual 

library(tidyverse)
library(terra)
library(sf)
library(INLA)

#simulation function
source("r_scripts/02_iSSF/07_simulation/sim_path_function_thompson_1-7-23.R")

mod_dat <- read_csv("top_ssf_dat_1-10-24.rds")

animal_info <- mod_dat %>% distinct(animal_id, .keep_all = TRUE) %>% 
  select(id1, animal_id:dispersal_status) %>% 
  rename(ID = id1)

#########################################################################
##
## 1.Import or fit top global model
##
##########################################################################

top_mod <- readRDS("top_ssf_1-10-24.rds")

#extract individual-level coefficients

coeff_names <-  top_mod$names.fixed

mod_summ <- top_mod$summary.random[names(top_mod$summary.random) != "step_id_"]

names(mod_summ) <- coeff_names

betas <- list()
for(i in 1:length(mod_summ)){
  name <- names(mod_summ)[i]
  betas[[i]] <- mod_summ[[i]] %>% as_tibble() %>% select(ID, mean) %>% 
    rename(!!name := mean)
}

betas <- reduce(betas, full_join, by = "ID")

#dataframe of estimated beta coefficients for each individual with demographics
indiv_coeffs <- animal_info %>% left_join(betas, by = join_by(ID))


#### RESUME HERE. NEXT STEP IS TO DERIVE STEP LENGTH AND TURN ANGLE DISTRIBUTIONS FOR EACH INDIVIDUAL ####


### Population Level Posterior Sampling (Optional)
#list variable name tags for posterior sampling
# top_mod$misc$configs$contents
# 
# #sampling from population-level posterior distributions 
# test <- inla.posterior.sample(n=5, result = top_mod) #, selection = "elevation"
# test2 <- inla.posterior.sample.eval("elevation", test)

#########################################################################
##
## 2.Import or fit top global model
##
##########################################################################


move_pars <- c("gamma_shape" = 1, "gamma_rate" = 1, "vm_mu" = 1, "vm_kappa" = 1)

R_rsf <- rast("/Users/tb201494/Desktop/1km_buffer/cov_stack_pred_30m_11-30-2023.tif")

r_mask <- rast("data/Habitat_Covariates/op_water_mask2_01-07-24.tif")

start_zone <- vect("data/Habitat_Covariates/study_area_polys/sim_start_zone_water_masked_01-07-24.shp")


# tmp <- rast(R_rsf$elevation, vals=0)
# 
# water_mask <- mask(tmp, r_mask, updatevalue=1, inverse=TRUE)
# 
# writeRaster(water_mask, "data/Habitat_Covariates/op_water_mask2_01-07-24.tif")

indivs = 1
paths_per_indiv = 3
steps_per_path = 100

test <- sim_paths(
         n_lists = indivs, 
         n_paths_per_list = paths_per_indiv, 
         n_steps_per_path = steps_per_path, 
         n_rand = 10,
         move_pars = move_pars,
         R_rsf = R_rsf, 
         R_step = NULL,
         R_log_step = NULL,
         R_cos_angle = NULL,
         R_daylight = NULL,
         R_mask = r_mask,
         step_cos_angle_par = 0,
         log_step_cos_angle_par = 0,
         step_day_par = 0,
         log_step_day_par = 0,
         cos_angle_day_par = 0,
         x0_bounds = start_zone, 
         daylight_var = numeric(steps_per_path),
         season_var = rep(1, steps_per_path),
         n_check_mask = 1,
         step_factor = 1,
         n_cores = indivs, 
         n_print = steps_per_path + 1,
         ret_list = (indivs > 1)) 


test %>% st_as_sf(coords = c("x", "y"), crs = 5070) %>% mapview::mapview(zcol = "path_id")
