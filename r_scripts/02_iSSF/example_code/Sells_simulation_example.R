################################################################################
### Grizzly Bear CRW Step Selection Simulation ###
### S Sells, MTCWRU ###

# Change as needed:
n.clusters <- 25 # number of clusters for parallel computing
n.steps <- 5000 # number of steps
n.runs <- 100 # number of iterations
n.runs.per.fit <- n.runs/10 # will refit model 10 times to capture variation
n.fits <- n.runs / n.runs.per.fit


# Load data and libraries:
library(amt)
library(purrr)
library(raster)
library(tidyverse)
library(gplots)
library(sf)
library(sp)
library(fitdistrplus)
library(circular)
library(doParallel)
library(boot)
library(adehabitatHR)
library(grizmove)
options(scipen=999)


# Load spatial data:
start_zone <- load_crw_start(where = "NCDE") # load spatial polygon of study area
cov_stack <- load_cov_data(sex = "Male") %>% # load the raster layers for males
  raster::crop(NCDE_300k) %>% raster::mask(NCDE_300k) # and crop to area needed to fit models: NCDE_300k buffer (study area plus 300km buffer, same as from hypothesis testing phase)

# Crop rasters to area needed to speed simulations: 100km buffer of the start_zone
crw_area <- raster::buffer(start_zone, 100000) %>% raster::crop(StudyArea) # crop back to the maximal study area boundary for areas extending beyond after buffering
cov_stack_sim <- cov_stack %>% raster::crop(crw_area) %>% raster::mask(crw_area) # this will be the raster stack used for simulations

# Double check things look right by plotting:
raster::plot(cov_stack_sim[[1]])
plot(start_zone, add = T) # same

# Subset data to appropriate sex:
dat_all <- grizmove::bear_df %>% filter(sex == "Male")



# Run Simulations by Individual:
for (i in 1:length(unique(dat_all$sex.id))) { # sex.id is 1:x for # of males
  # Prep bear data for individual i:
  dat_bear <- dat_all %>%
    filter(sex.id == i) %>%
    arrange(t) %>%
    make_bear_track() # --> see below
  
  # Cycle through n.fits (refitting model n.fits times to capture variation):
  for (m in 1:n.fits) {
    my.cluster <- parallel::makeCluster(n.clusters, type = "PSOCK")
    doParallel::registerDoParallel(cl = my.cluster, cores = n.clusters)
    
    # Fit Individual Model:
    dat <- dat_bear %>%
      sample_bear_steps() # --> see below
    if(Sex == "Male"){ mod <- run_male_annual(dat, i) } # --> see below
    if(Sex == "Female"){ mod <- run_female_annual(dat, i) }
    
    # Prepare to run simulation by extracting mod info:
    mov <- dat %>% dplyr::filter(case_ == T) %>% dplyr::select(sl_, log_sl_, cos_ta_) # save the movement parameters from mod
    new_gamma <- update_sl_distr(mod) # amt fx
    new_vm <- update_ta_distr(mod) # ""
    Scale <- as.numeric(new_gamma$params[2])
    Shape <- as.numeric(new_gamma$params[1])
    ta_mu <- as.numeric(new_vm$params$mu)
    ta_kappa <- new_vm$params$kappa
    vals <- make_bear_map(mod) # --> see below
    
    s_node <- matrix(NA, ncol = 2, nrow = n.runs.per.fit * 3) # prepare to define start nodes by creating empty matrix
    for(j in 1:nrow(s_node)){ s_node[j,] <- generate_start_nodes(start_zone = start_zone) } # --> see below
    
    # Run simulation:
    clusterExport(my.cluster, c("Scale", "Shape", "crw_area", "vals", "ta_mu", "ta_kappa"))
    crw_bear <- foreach(k = 1:n.runs.per.fit,
                        .packages = c("tidyverse", "sp", "sf", "circular", "grizmove"),
                        .errorhandling = "pass") %dopar% {
                          step_sim_fun(n_steps = n.steps, # --> see below
                                       x_start = s_node[k,1],
                                       y_start = s_node[k,2])
                        }
    # Save:
    if(m == 1){ crw_bear_all <- crw_bear } else { crw_bear_all <- rbind(crw_bear_all, crw_bear) }
    parallel::stopCluster(cl = my.cluster)
  }
  save(crw_bear_all, file = paste(Sex, i, Where, "CRW.Rdata", sep ="_"))
}


#### Functions:
# Prepare data for individual, make tracks:
make_bear_track <- function (dat)
{
  dat_ <- dat %>% make_track(x, y, t, crs = sp::CRS("+proj=utm +zone=12 +datum=WGS84")) %>%
    amt::transform_coords(raster.crs) %>% track_resample(rate = hours(3),
                                                         tolerance = minutes(45)) %>% steps_by_burst()
  dat_ <- dat_ %>% mutate(bear_name = first(dat$bear_name),
                          sex.id = first(dat$sex.id)) %>% 
    filter(sl_ >= 100) %>%
    filter(sl_ <= 15000)
  return(dat_)
}

# Sample random steps:
sample_bear_steps <- function (dat)
{
  dat <- dat %>% dplyr::arrange(burst_, t1_, .by_group = T) %>%
    amt::random_steps(10) %>% amt::extract_covariates(cov_stack, where = "end") %>% 
    dplyr::mutate(log_sl_ = log(sl_), cos_ta_ = cos(ta_))
  return(dat)
}

# Run model for individual: this is where each individual's model is specified. These were found through the kfolds CV methods discussed previously. 
# After identifying the best model per bear, I saved their model specifications into this function. Their ids are specified as sex_id in the code above, corresponding to i in the loop:
run_male_annual <- function (dat, i) 
{
  if (i == 1) {
    mod <- dat %>% fit_issf(case_ ~ ndvi_e + ruggedness_e +
                              d2forestedge_e + I(d2forestedge_e^2) + densriparian_e +
                              sl_ + log_sl_ + cos_ta_ + strata(step_id_), model = TRUE)
  }
  if (i == 2) {
    mod <- dat %>% fit_issf(case_ ~ ndvi_e + I(ndvi_e^2) +
                              ruggedness_e + I(ruggedness_e^2) + d2forestedge_e +
                              I(d2forestedge_e^2) + densforestedge_e + I(densforestedge_e^2) +
                              densriparian_e + I(densriparian_e^2) + densbuildings_e +
                              I(densbuildings_e^2) + d2core_e + I(d2core_e^2) +
                              sl_ + log_sl_ + cos_ta_ + strata(step_id_), model = TRUE)
  }
  if (i == 3) {
    mod <- dat %>% fit_issf(case_ ~ ndvi_e + ruggedness_e +
                              d2forestedge_e + I(d2forestedge_e^2) + densforestedge_e +
                              I(densforestedge_e^2) + densriparian_e + densbuildings_e +
                              I(densbuildings_e^2) + d2core_e + sl_ + log_sl_ +
                              cos_ta_ + strata(step_id_), model = TRUE)
    ...
  }
}
  # Make a raster by applying the model to the landscape: output is "vals" for the simulation
  make_bear_map <- function (mod)
  {
    cov_stack_values <- values(cov_stack_sim) %>% data.frame %>%
      mutate(sl_ = median(mov$sl_), # these movement parameters are required in the df since they're in the model
             log_sl_ = median(mov$log_sl_),
             cos_ta_ = median(mov$cos_ta_, na.rm = T), 
             step_id_ = mod$model$model$`strata(step_id_)`[1])
    predictions <- raster::predict(mod$model, newdata = cov_stack_values, type = "lp", allow.new.levels = TRUE)
    predictions <- exp(predictions)
    x.min <- quantile(predictions, 0.025, na.rm = T) # omit the extremes
    x.max <- quantile(predictions, 0.975, na.rm = T)
    predictions[predictions > x.max] <- x.max
    predictions[predictions < x.min] <- x.min
    x.range <- x.max - x.min
    predictions. <- (predictions - x.min)/x.range
    vals <- cov_stack_sim[[1]] %>% raster::setValues(predictions.)
    return(vals)
  }
  
  # Ensure that start nodes are drawn from within start_zone and save matrix of the start nodes to use for each iteration:
  generate_start_nodes <- function (start_zone)
  {
    x1 <- start_zone@bbox[1, 1]
    x2 <- start_zone@bbox[1, 2]
    y1 <- start_zone@bbox[2, 1]
    y2 <- start_zone@bbox[2, 2]
    repeat {
      x_rand <- runif(1, x1, x2)
      y_rand <- runif(1, y1, y2)
      if (sp::point.in.polygon(point.x = x_rand, point.y = y_rand,
                               pol.x = start_zone@polygons[[1]]@Polygons[[1]]@coords[,
                                                                                     1], pol.y = start_zone@polygons[[1]]@Polygons[[1]]@coords[,
                                                                                                                                               2]) == 1)
        break
    }
    return(cbind(x_rand, y_rand))
  }
  
  # The main simulation function to have the simulated agent move:
  step_sim_fun <- function (n_steps, x_start, y_start)
  {
    mat <- matrix(NA, nrow = n_steps, ncol = 3)
    colnames(mat) <- c("x", "y", "abs_angle")
    for (j in 1:n_steps) {
      if (j == 1) { proposed_coords <- matrix(replicate(11, first_pt_gen(x = x_start, y = y_start)), ncol = 3, byrow = T) } # --> see below
      else { proposed_coords <- matrix(replicate(11, avail_step_gen(prev_abs_angle = mat[j - 1, "abs_angle"], x = mat[j - 1, "x"], y = mat[j - 1, "y"])), ncol = 3, byrow = T) } # --> see below
      proposed_coords <- as.data.frame(proposed_coords)
      proposed_coords <- sf::st_as_sf(proposed_coords, coords = c(1:2))
      wts <- raster::extract(vals, proposed_coords)
      selected <- proposed_coords[sample(1:11, 1, prob = wts/sum(wts)),
      ]
      mat[j, ] <- do.call(cbind, sf::st_geometry(selected)) %>%
        tibble() %>% t() %>% cbind(selected$V3)
    }
    return(mat)
  }
  
  # For first step, move to a likely starting step from a sample of 11 points near the start_node for true first step:
  first_pt_gen <- function (x, y) # use the movement parameters from the model to help sample steps nearby, and ensure they're in the study area:
  {
    repeat {
      samp_angle <- runif(1, min = 0, max = 2 * pi)
      samp_dist <- rgamma(1, scale = Scale, shape = Shape)
      x_gen <- x + cos(samp_angle) * samp_dist
      y_gen <- y + sin(samp_angle) * samp_dist
      if (sp::point.in.polygon(point.x = x_gen, point.y = y_gen,
                               pol.x = crw_area@polygons[[1]]@Polygons[[1]]@coords[,
                                                                                   1], pol.y = crw_area@polygons[[1]]@Polygons[[1]]@coords[,
                                                                                                                                           2]) == 1)
        break
    }
    return(cbind(x_gen, y_gen, samp_angle))
  }
  
  # For subsequent steps, continue moving to new steps by sampling 11 nearby, moving to one, through rest of simulation:
  avail_step_gen <- function (prev_abs_angle, x, y) # similar to above, but now also accounting for previous turning angle etc.:
  {
    repeat {
      samp_angle <- rvonmises(1, mu = circular(ta_mu), kappa = ta_kappa,
                              control.circular = list(units = "radians"))
      samp_angle <- ifelse(samp_angle > pi, samp_angle - 2 *
                             pi, samp_angle)
      rand_angle <- samp_angle + prev_abs_angle
      rand_angle <- ifelse(rand_angle < 0, rand_angle + 2 *
                             pi, rand_angle)
      samp_length <- rgamma(1, scale = Scale, shape = Shape)
      x_rand <- x + cos(rand_angle) * samp_length
      y_rand <- y + sin(rand_angle) * samp_length
      proposed_coords <- as.data.frame(cbind(x_rand, y_rand))
      proposed_coords <- sf::st_as_sf(proposed_coords, coords = c(1:2))
      wts <- raster::extract(vals, proposed_coords)
      if (sp::point.in.polygon(point.x = x_rand, point.y = y_rand,
                               pol.x = crw_area@polygons[[1]]@Polygons[[1]]@coords[,
                                                                                   1], pol.y = crw_area@polygons[[1]]@Polygons[[1]]@coords[,
                                                                                                                                           2]) == 1 && !is.na(wts))
        break
    }
    return(cbind(x_rand, y_rand, rand_angle))
  }
  