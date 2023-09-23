#### Global Occupancy Model Selection ####

# Author: Read Barbee

# Date:2023-09-19

# Purpose:


################################ libraries #################################
library(tidyverse)
library(unmarked)
library(ubms)
library(camtrapR)
library(beepr)
library(doParallel)
library(DataExplorer)
library(stocc)
library(sf)

#########################################################################
##
##  Specify Model Parameters
##
##########################################################################

params <- c("popdens_hii",
            "easting",
            "npp",
            "roads_hii",
            "precip",
            "dist_water",
            "perc_tree_cover",
            "infra_hii")



quad_params <- c("dist_water",
                 "roads_hii")


quad_terms <- vector()

for(i in 1:length(quad_params)){
  quad_terms[i] <- paste0("I(",quad_params[i],"^2)")
}


#########################################################################
##
## 1. Import and format step data
##
##########################################################################

occ_dat <- read_csv("data/Camera_Data/master/ocp_occ_dat_9-18-23.csv") %>% 
  mutate(across(station_id_year:year, as.factor)) %>% 
  mutate(aspect_rad = (pi*aspect)/180, .after=aspect) %>%
  mutate(northing = cos(aspect_rad),
         easting = sin(aspect_rad), .after=aspect_rad) %>% 
  rename(aspect_deg = aspect) %>% 
  select(-c(aspect_deg, aspect_rad, land_cover_usfs, land_use_usfs))

proj_coords <- occ_dat %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  st_transform(crs=5070) %>% 
  st_coordinates() %>% as.data.frame

occ_dat <- occ_dat %>% 
  mutate(x = proj_coords[,1],
         y = proj_coords[,2], .after = lat)

#scale covariates
occ_dat_scaled <- occ_dat %>% 
  mutate(across(tree_cover_hansen:dist_water, scale)) %>% 
  mutate(across(tree_cover_hansen:dist_water, as.numeric))


#remove station rows with missing covariate values
complete_cases <- occ_dat_scaled %>% select(tree_cover_hansen:dist_water) %>% complete.cases()

occ_dat_complete <- occ_dat_scaled %>% 
  mutate(comp = complete_cases) %>% 
  filter(comp==TRUE)


#convert to unmarked dataframe
umf <- unmarkedFrameOccu(y = occ_dat_complete %>% select(d_1:d_365),
                         siteCovs = occ_dat_complete %>% select(station_id, tree_cover_hansen:dist_water))

umf_cell <- unmarkedFrameOccu(y = occ_dat_complete %>% select(d_1:d_365),
                         siteCovs = occ_dat_complete %>% select(cell_id, tree_cover_hansen:dist_water))

#########################################################################
##
## 2. Fit global model with quadratics
##
##########################################################################

#Fit global glmmTMB model: takes about ~ 15 min for 15 covs
run_global <- function (dat){
  #library(survival)  survival::clogit
  form <- as.formula(paste0("~1 ~", #remove standard intercept to be replace with stratum-based intercept
                            #fixed effects
                            paste(c(params, quad_terms), collapse = " + "), "+",
                            #random intercept (strata)
                            "(1|station_id)"))
  
  
  fit <- stan_occu(form, data=dat, chains=3, iter=1000)
  return(fit)
}

run_global_cell <- function (dat){
  #library(survival)  survival::clogit
  form <- as.formula(paste0("~1 ~", #remove standard intercept to be replace with stratum-based intercept
                            #fixed effects
                            paste(c(params, quad_terms), collapse = " + "), "+",
                            #random intercept (strata)
                            "(1|cell_id)"))
  
  
  fit <- stan_occu(form, data=dat, chains=3, iter=1000)
  return(fit)
}


#system.time(global_fit <- run_global(umf))
summary(global_fit, "state")

#saveRDS(global_fit, "occu_global_fit_9-19-23.rds")

#global_fit <- readRDS("occu_global_fit_9-19-23.rds")

system.time(global_fit_cell <- run_global_cell(umf_cell))

summary(global_fit_cell, "state")

#saveRDS(global_fit_cell, "occu_global_fit_cell_9-22-23.rds")

#global_fit_cell <- readRDS("occu_global_fit_cell_9-22-23.rds")

#########################################################################
##
## 3. Model diagnostics
##
##########################################################################

#traceplots
traceplot(global_fit_cell)

#residual plots
plot_residuals(global_fit_cell, submodel="state")
plot_residuals(global_fit_cell, submodel="state", covariate="popdens_hii")

#caluclate and plot actual values against values simulated from the model
global_fit_gof <- gof(global_fit_cell, draws=100, quiet=TRUE)

plot(global_fit_gof) #doesn't look good


#simulate new datasets to calculate fit statistic
sim_y <- posterior_predict(global_fit, "y", draws=100)
dim(sim_y)

#calculate proportion of 0s in each simulated dataset
prop0 <- apply(sim_y, 1, function(x) mean(x==0, na.rm=TRUE))
actual_prop0 <- mean(getY(global_fit) == 0, na.rm=TRUE)

#Compare proportion of 0s in each simulated dataset to proportion of 0s in actual dataset
hist(prop0, col='gray')
abline(v=actual_prop0, col='red', lwd=2)

#########################################################################
##
## 4. Model inference
##
##########################################################################
#marginal covariate effects
ubms::plot_effects(global_fit_cell, "state")

#predict occupancy values for each site
psi <- predict(global_fit_cell, submodel="state")



#plot predictid occupancy across the study area
plot_spatial(global_fit)

#plot spatial neighbors to individual site. Consider fitting spatial model?
RSR(proj_coords$X, proj_coords$Y, threshold = 2000, plot_site = 10)


#can't have both regular and spatial random effect in model
# #Fit global glmmTMB model: takes about ~ 15 min for 15 covs
# run_global_spatial <- function (dat, threshold){
#   #library(survival)  survival::clogit
#   form <- as.formula(paste0("~1 ~", #remove standard intercept to be replace with stratum-based intercept
#                             #fixed effects
#                             paste(c(params, quad_terms), collapse = " + "), " + ",
#                             "RSR(x, y, threshold = ", threshold, ") + ",
#                             #random intercept (strata)
#                             "(1|station_id)"))
#   
#   
#   fit <- stan_occu(form, data=dat, chains=3, iter=1000)
#   return(fit)
# }
# 
# 
# fit_global_spatial <- run_global_spatial(umf, threshold = 2000)

#########################################################################
##
## 5. Map predictions
##
##########################################################################

cov_stack <- terra::rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2.tif")

names(cov_stack) <- c("tree_cover_hansen",
                      "gpp",
                      "infra_hii",
                      "landuse_hii",
                      "land_cover_usfs",
                      "land_use_usfs",
                      "npp",
                      "popdens_hii",
                      "power_hii",
                      "precip",
                      "rails_hii",
                      "roads_hii",
                      "elevation",
                      "slope",
                      "aspect",
                      "tri",
                      "tpi",
                      "perc_tree_cover",
                      "perc_nontree_veg",
                      "perc_nonveg",
                      "ndvi",
                      "evi",
                      "dist_water")

cov_stack$northing <- cos((pi*cov_stack$aspect)/180)
cov_stack$easting <- sin((pi*cov_stack$aspect)/180)
#cov_stack$station_id = as.factor("LEKT_Station1")

cov_stack_resampled <- terra::aggregate(cov_stack, fact=10, cores=5)

cov_stack_rast <- raster::stack(cov_stack_resampled)

cov_stack_rast$station_id <- as.factor("LEKT_Station1")

mean_station_locs <- occ_dat %>% 
  group_by(cell_id) %>% 
  summarize(x = mean(x), y = mean(y)) %>% 
  st_as_sf(coords = c("x", "y"), crs = 5070)


station_locs <- occ_dat %>%  st_as_sf(coords = c("x", "y"), crs = 5070)


#PROBLEM: LEKT Station 73 and others are in wildly different locations between years

mapview::mapview(mean_station_locs)
mapview::mapview(station_locs %>% filter(station_id == "LEKT_Station73"))


#calculate maximum distances between stations with the same ids between years
occ_dat %>% group_by(station_id)

station_ids <- occ_dat %>% distinct(station_id) %>% pull()

dists <- vector()
for(i in 1:length(station_ids)){
  rows <- occ_dat %>% 
    filter(station_id == station_ids[i]) %>% 
    st_as_sf(coords=c("x", "y"), crs = 5070) 
  
  dist_pairs <- vector()
  for(j in 1:(nrow(rows)-1)){
    row1 <- rows[j,]
    row2 <- rows[j+1,]
    
  dist_pairs[j] <- st_distance(row1, row2)
  }
  
  first_last <- st_distance(rows[1,], rows[nrow(rows),])
  
  dists[i] <- max(c(dist_pairs, first_last))
  
}

max_dists <- tibble(station_id = station_ids,
                    max_distance = dists)

 max_dists %>% filter(max_distance > 30) %>% 
  filter(str_detect(station_id, "LEKT"))



#`sigma_state[sigma [1|station_id]]`

#works with SpatRaster object but doesn't return raster. Doesn't work with RasterStack

#takes ~23min at 300m resolution. Only seems to work when re.form = NA (random effects not included in predictions). ~7 min with cell_id random eff
system.time(surface_rast <- ubms::predict(object = global_fit_cell, 
                   submodel = "state", 
                   newdata = cov_stack_rast,
                   transform = TRUE,
                   re.form = NA,
                   level = 0.95))

mapview::mapview(surface_rast$Predicted)


#works. see code below to project prediction dataframe into raster surface
system.time(surface <- ubms::predict(object = global_fit, 
                                     submodel = "state", 
                                     newdata = cov_stack_resampled,
                                     transform = TRUE,
                                     re.form = NULL,
                                     level = 0.95))


surface <- surface %>% rownames_to_column("cell")

df <- terra::as.data.frame(cov_stack_resampled, xy = TRUE, cells = TRUE, na.rm = FALSE) %>% 
  select(cell:y, station_id) %>% 
  mutate(cell = as.character(cell))

preds <- surface %>% left_join(df, by=join_by(cell)) %>% select(x,y,everything())

#projection working, but weird
test <- raster::rasterFromXYZ(preds, crs = cov_stack_rast)
mapview::mapview(test$Predicted)




system.time(surface2 <- ubms::predict(object = global_fit, 
                                     submodel = "state", 
                                     #newdata = cov_stack_resampled,
                                     transform = TRUE,
                                     re.form = NULL,
                                     level = 0.95))







#global_fit@submodels@submodels$state@formula

df_dat <- raster::as.data.frame(cov_stack_rast, xy=TRUE)


out <- cbind(df_dat[,1:2,drop=FALSE],
             raster::predict(global_fit, "state", newdata=df_dat, transform=TRUE,
                     re.form=NULL, level=0.95))














