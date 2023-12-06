#### amt Population Level Simulation ####

# Author: Read Barbee

# Date:2023-09-19

# Last updated:2023-12-01

# Purpose:


################################ libraries #################################
library(tidyverse)
library(glmmTMB)
library(amt)
#library(MuMIn)
library(INLA)
library(beepr)
library(sf)
library(terra)
library(scico)
library(foreach)
library(doParallel)


#import top model fit
top_mod <- readRDS("top_ssf_ndvi_quad_12_01_23.rds")

########################################################################
##
## 1. Import and format step data
##
##########################################################################

#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_annual_cov_11-30-2023.csv")

#locs <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_7-11-2023.csv")


#select relevant columns and scale
steps_scaled <- steps %>% 
  select(!contains("usfs"), -aspect) %>% 
  mutate(across(elevation:dens_all_roads_annual, scale)) %>% 
  mutate(across(elevation:dens_all_roads_annual, as.numeric)) 


########################################################################
##
## 2. Import covariate data
##
##########################################################################

cov_stack <- terra::rast("/Users/tb201494/Desktop/1km_buffer/cov_stack_pred_300m_11-30-2023.tif")

cov_stack$`I(ndvi_annual^2)` <- cov_stack[["ndvi_annual"]]^2
#resample cov_stack to lower resolution for simulation
#cov_stack2 <- terra::aggregate(cov_stack, fact=10, cores=5)

########################################################################
##
## 2. Set up iSSF model for simulation
##
##########################################################################

#extract coefficient estimates from the fitted muff-style global model
top_estimates <- top_mod$summary.fixed %>% as.data.frame() %>% rownames_to_column("term") %>% select(term, mean)


coefs <- vector()
names <- vector()
for(i in 1:nrow(top_estimates)){
  coefs[i] <- top_estimates$mean[i]
  names[i] <- top_estimates$term[i] %>% 
    str_replace(coll("^2)"), "2") %>% 
    str_remove(coll("I(")) %>% 
    str_c("_end")
}

names(coefs) <- names

#subset linear terms only
# coefs_linear <- coefs[str_detect(names, "2", negate = TRUE)] 
# 
# names_linear <- names[str_detect(names, "2", negate = TRUE)] %>% str_remove(., coll("_end"))
# 
# cov_stack <- cov_stack[[c(names_linear)]]


# OptionaL: artificially construct gamma distribution
# x <- seq(0.1, 10, 0.1)
# plot(x, dgamma(x, shape = 2, scale = 2), type = "l")

#Fit a global issf style model using the muff-estimated coefficients.
mod <- make_issf_model(coefs = c(coefs_linear, sl_ = median(steps$sl_, na.rm=T), ta_ = median(steps$ta_, na.rm=T)),
                       sl = fit_distr(steps$sl_, "gamma", na.rm = TRUE),
                       ta =fit_distr(steps$ta_, "vonmises", na.rm = TRUE))

########################################################################
##
## 3. Sample individual models
##
##########################################################################
#extract individual coefficients
n_indiv <- top_mod$summary.random[[1]] %>% distinct(ID) %>% count() %>% pull()

id_list <- 1:n_indiv

id_sample <- sample(id_list, size = 1, replace=FALSE)

coeff_estimates <- map(top_mod$summary.random, function(x){return(x[id_sample, 'mean'])}) %>% bind_cols() %>% select(c(id1:id9, id10)) %>% pivot_longer(c(id1:id9, id10), names_to = "name", values_to = "value") %>% deframe()

names(coeff_estimates) <- rownames(top_mod$summary.fixed) 
  


top_mod$summary.random


########################################################################
##
## 3. Prepare barrier polygons and start zone for simulation
##
##########################################################################

## Create non-habitat zones where simulation starts shouldnt be generated (i.e. water bodies)

## Study area boundary ##
op_poly <- st_read("data/Habitat_Covariates/study_area_polys/ocp_study_area_poly_wa_only_10-30-23.shp")

## Water Polygons ##
water_polys <- st_read("data/Habitat_Covariates/washington_water_polygons/op_water_polys_with_salt.shp") %>% st_transform(crs = 5070)

water_polys_cropped <-  water_polys %>% 
  filter(WB_PERIOD_ =="PER" ) %>% #only include permanent water bodies
  sf::st_crop(op_poly)
  
water_polys_start <- water_polys_cropped %>% 
  sf::st_union() %>% 
  st_sf()

#create mask layer
water_barriers <- water_polys_cropped %>% 
  filter(SHAPEAREA >= 2000000) %>% 
  sf::st_union() %>% 
  st_sf()

#mask large water bodies out of cov_stack_raster
#cov_stack_masked <- terra::mask(cov_stack, water_barriers, inverse = TRUE)

#mapview::mapview(water_polys_cropped)


#import study area poly to use as generation zone for starting points
start_zone <- st_read("data/Habitat_Covariates/study_area_polys/sim_start_zone_w_I5_9-25-23.shp")

start_zone<- start_zone$geometry

########################################################################
##
## 4. Simulate
##
##########################################################################

#Function to generate n start locations outside of water. way faster than the way Sarah did it
generate_start_nodes <- function (n, start_zone, barriers)
{
  nodes <- st_sample(start_zone, size = n + 1000, type = "random", exact = FALSE)
  
  pts<- st_as_sf(nodes, coords = c("X", "Y"), crs = 5070) %>% 
    mutate(intersects = lengths(st_intersects(., barriers)) > 0)
  
  filtered <- pts %>% 
    filter(intersects==FALSE) %>% 
    st_coordinates() %>% 
    as.data.frame()
  
  samp <- sample_n(filtered, n)
  
  return(samp)
}


#Function to simulate n paths of desired length based on a single global model
sim_paths_general <- function(n_paths, n_steps, mod, covs, barriers){
  
  ext <- start_zone
  start_nodes <- generate_start_nodes(n_paths, ext, barriers)
  #9 minutes for 3 paths of 100 steps at 300m resolution
  
  #doParallel::registerDoParallel(cores = parallel::detectCores()-1)
  
  paths <- list()
  system.time(for(i in 1:n_paths){
    x_coord <- start_nodes$X[i]
    y_coord <- start_nodes$Y[i]
    
    start <- make_start(x= c(x_coord, y_coord),
                        time = ymd_hms("2020-01-01 00:00:00"),
                        ta = 0,
                        dt = hours(2),
                        crs = 5070)
    
    
    k1 <- amt::redistribution_kernel(mod, map = cov_stack, start = start, landscape = "continuous") #cov_stack_masked
    
    paths[[i]] <- amt::simulate_path(k1, n.steps = n_steps)
    
    print(paste0(i, "/", n_paths ))
  })
  
  #stopImplicitCluster()
  
  paths <- bind_rows(paths, .id = "path_id")
  
  return(paths)
  
}

#run simulation ~ 2 hours for one year long track (4,400 steps), 8.4 hours for 4 tracks of 4,400 steps. Eventually want to do 100 tracks per individual lion
system.time(paths <- sim_paths_general(4, 4400, mod, cov_stack, water_barriers)) #cov_stack_masked #water_polys_start

#plot simulated paths
paths_sf <- paths %>% st_as_sf(coords = c("x_", "y_"), crs = 5070)

mapview::mapview(paths_sf)


########################################################################
##
## 4. Tally number of steps taken per raster cell
##
##########################################################################

#template raster for count calculation
tmp_rast <- cov_stack2$tree_cover_hansen
values(tmp_rast) <- 0


#generate raster of the number of times each cell was used
count_rast <- raster::rasterize(paths_sf, tmp_rast, fun = "count")


#bin counts by quantile
pred_vals <- terra::values(count_rast)
breaks <- quantile(pred_vals, probs = 0:10/10, na.rm = T) #obtain 10 quantile values for preds
#breaks_j <- breaks + (seq_along(breaks) * .Machine$double.eps) #add jitter to breaks
binned <- terra::classify(count_rast, rcl=as.vector(breaks_j))


#layercor


########################################################################
##
## 4. Calculate UD from simulated paths (amt) -maybe unnecessary
##
##########################################################################

test2 <- split(test, test$path_id)

n <- 5
#system.time(p1 <- replicate(n, simulate_path(rdk.1a, n = 30), simplify = FALSE))

# UD--note, this is based around homeranging behavior
uds <- lapply(c(5, 2000), function(i) {
  tibble(
    rep = 1:n, 
    path = map(test2, ~ dplyr::slice(.x, i))
  ) |> unnest(cols = path) |> filter(!is.na(x_)) |> 
    make_track(x_, y_) |> hr_kde(trast = cov_stack2_rast[[1]])
})


#Plot UD
xy <- as.data.frame(uds[[2]]$ud, xy = TRUE)
names(xy) <- c("x", "y", "w")

pl3 <- ggplot(xy, aes(x, y, fill = tree_cover_hansen)) + geom_raster()  +
  #geom_path(aes(x, y), data = geom(p, df = TRUE), lty = 2, col = "grey30", 
  #inherit.aes = FALSE) +
  coord_equal() + 
  theme_light() +
  #geom_point(x = 0, y = 0, col = "red") +
  theme(legend.position = "none") +
  scale_fill_scico(palette = "lajolla") +
  #scale_y_continuous(limits = c(-80, 80)) + scale_x_continuous(
  #limits = c(-80, 80)) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  labs(x = "x")





################################ Graveyard #################################


# coefs <- vector()
# for(i in 1:nrow(top_estimates)){
#   tm <- top_estimates$term[i] %>% 
#     str_replace(coll("^2)"), "2") %>% 
#     str_remove(coll("I("))
#   
#   coefs[i] <- paste0(tm, " = ", top_estimates$Estimate[i])
# }


# start <- make_start(x= c(-2114054, 3007263),
#                     #time = ymd_hms("2022-04-05 05:00:35"),
#                     ta = 0,
#                     dt = hours(2),
#                     crs = 5070)
# 
# 
# k1 <- amt::redistribution_kernel(mod, map = cov_stack2, start = start)
# 
# s1 <- simulate_path(k1, n.steps = 1000)
# 
# terra::plot(cov_stack2$tree_cover_hansen)
# lines(s1$x_, s1$y_, col = "red")

# 
# k1 <- redistribution_kernel(mod, map = cov_stack, start = start,
#                             landscape = "continuous", tolerance.outside = 0.2, 
#                             n.control = 1e4)


# Ensure that start nodes are drawn from within start_zone and save matrix of the start nodes to use for each iteration:

# generate_start_nodes <- function (n, start_zone, barriers)
# {
#   x1 <- start_zone$xmin
#   x2 <- start_zone$xmax
#   y1 <- start_zone$ymin
#   y2 <- start_zone$ymax
#   
#   nodes <- list()
#   for(i in 1:(n + 100)){
#     x_rand <- runif(1, x1, x2)
#     y_rand <- runif(1, y1, y2)
#     
#     nodes[[i]] <- tibble(X = x_rand, Y = y_rand)
#     
#   }
#   nodes <- bind_rows(nodes)
#   pts<- st_as_sf(nodes, coords = c("X", "Y"), crs = 5070) %>% 
#     mutate(intersects = lengths(st_intersects(., barriers)) > 0)
#   
#   filtered <- pts %>% 
#     filter(intersects==FALSE) %>% 
#     st_coordinates() %>% 
#     as.data.frame()
#   
#   samp <- sample_n(filtered, n)
#   
#   return(samp)
# }







#cov_stack$roads_hii2 <- (cov_stack$roads_hii^2)


#reordering doesn't seem to make a difference
#na.omit doesn't help
#converting to raster stack doesn't help

#cov_stack3 <- na.omit(cov_stack2)

#check to make sure point is within raster area--it is
#terra::plot(cov_stack2[[1]])
#terra::points(terra::vect("POINT (-2114054 3007263)", crs="+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"))
#reducing cov_stack_resolution doesn't work either

# mod <- make_issf_model(coefs = c(coefs_linear, sl_ = median(steps$sl_), ta_ = median(steps$ta_), strata = steps$step_id_),
#                        sl = make_gamma_distr(shape = 2, scale = 2, vcov = NULL),
#                        ta =make_vonmises_distr(kappa = 4, vcov = NULL))

#cov_stack2_scaled <- terra::scale(cov_stack2)

# cov_stack2$sl <- 1
# cov_stack2$ta <- 0

#cov_stack2_rast <- raster::stack(cov_stack2)




#sells method. works, but very slow. Probably more efficient just to generate a bunch of random points and then filter once based on barriers
# generate_start_nodes <- function (n, start_zone, barriers)
# {
#   x1 <- start_zone$xmin
#   x2 <- start_zone$xmax
#   y1 <- start_zone$ymin
#   y2 <- start_zone$ymax
#   
#   nodes <- list()
#   for(i in 1:n){
#     repeat{
#       x_rand <- runif(1, x1, x2)
#       y_rand <- runif(1, y1, y2)
#       
#       pt <- tibble(X = x_rand, Y = y_rand)
#       
#       rand_point <- st_as_sf(pt, coords = c("X", "Y"), crs = 5070)
#       
#       if(lengths(st_intersects(rand_point, barriers)) == 0)
#         break
#     }
#     nodes[[i]] <- pt
#   }
#   
#   nodes <- bind_rows(nodes)
#   return(nodes)
# }
