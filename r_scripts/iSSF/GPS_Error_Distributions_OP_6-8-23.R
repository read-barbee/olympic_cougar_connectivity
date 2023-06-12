#### GPS Error Distributions for the Olympic Peninsula Based on Sager 2007 ####

# Author: Read Barbee

# Date:2023-06-08 

# Purpose: Simulate GPS error data for 2D and 3D classes using Sager 2007 analysis from the Olympic Peninsula. Fit t-dsitribution to simulated data to derive tau and nu parameters for each fix class to use as constants in the state-space model


################################ Libraries #################################
library(tidyverse)
library(sf)
library(geosphere)
library(fitdistrplus)
library(extraDistr)
library(sgt)
# library(rmutil)
# library(spatstat)

#candidate distributions for simulation
#Weibull (exponential): k=0.2, lambda = 18.5 (best so far)
#Gamma dist: k = 1.441, theta = 14.47
#Weibull dist: k=1.23, lambda = 22.5
#levy dist: mu = 4.5, c = 2.1 



################################ Parameters from Sager 2007 #################################
#1,075 GPS locations acquired
#561 3D
#1,144 2D

# 3D error values 
max_error_3D <- 73
error_3D_95 <- 17.7
mean_error_3D <-  6.3 # to the SW

#2D error values
max_error_2D <- 2230
error_2D_95 <- 264.6
mean_error_2D <-  11.1 # to the NW

#number of points to simulate
n_3D <- 561 #561 is Sager's number #10000
n_2D <- 1144

n_dist_sample <- 20000 #number of distances to sample from weibull dist. Should be many more than n_3D or n_2D

#Weibull distribution parameters (distance)
shape_weib = 0.19
scale_weib = 18.5

#VonMises distribution parameters (bearing)
m_vm_3D = 3.92699 #mean bearing in radians (SW)
k_vm_3D = 1.5

m_vm_2D = 5.49779 #mean bearing in radians (NW)
k_vm_2D = 1.5


## Additional parameters

#reference coordinates
ref_lat <- 47.801262
ref_lon <- -123.717987

point_crs <- 4326 #WGS 84
dist_crs <- 5070 #Albers Equal Area Projection (meters)

seed <- 777


t_est <- function(){

################################ Step 1: Simulate 3D fixes #################################
# set.seed(seed)

#generate reference point (true location)
ref_point <- st_point(c(x = ref_lon , y = ref_lat)) #x = 47.801262, y = -123.717987 


######## Sample distances from custom Weibull distribution (3D fixes)  ####################

dist_3D <- rweibull(n = n_dist_sample, shape = shape_weib, scale = scale_weib ) #generate more values than needed for later truncation
dist_3D <- dist_3D %>% 
  as_tibble() %>% 
  filter(value <= max_error_3D) %>% #truncate to values <= 73
  pull(value) %>% 
  sample(size = n_3D) %>% #sample  values so length is equal to bearing vector
  as_tibble() %>% 
  rename(distance = value)

#check distribution
# dist_3D %>% pull(distance) %>% summary()
# dist_3D %>% pull(distance) %>% hist()


#Sample bearings from custom Von Mises distribution (3D fixes)  
dir_3D <- Rfast::rvonmises(n = n_3D, m = m_vm_3D, k = k_vm_3D, rads = TRUE) %>% 
  as_tibble() %>% 
  rename(bearing = value) %>% 
  mutate(bearing = as.numeric(bearing)) %>% 
  mutate(bearing = (bearing*180)/pi) #convert radians back to degrees for distance function
  

# dir_3D %>% pull(bearing) %>% summary()
# dir_3D %>% pull(bearing) %>% hist()


#Bind distance and bearing vectors to make matrix for distance simulation 
dist_matrix <- bind_cols(dist_3D, dir_3D)
  

#Simulate coordinates in lat/long based on simulation matrix and reference point 
sim_pts_3D <- destPoint(p = c(ref_lon, ref_lat), b = dist_matrix$bearing, d = dist_matrix$distance) %>% as_tibble() 

######## Plot Simulated Data ###

#Convert to sf object and reproject to Albers Equal Area for plotting
sim_pts_3D_sf <- st_as_sf(sim_pts_3D, coords = c("lon", "lat"), crs=point_crs) %>% st_transform(crs=dist_crs)

#create a new reference point at 0,0 in the Albers coordinate system for plotting
ref_point_utm <- tibble(X = 0, Y = 0) %>% 
  st_as_sf(coords = c("X", "Y"), crs=dist_crs)

#Reproject simulated points to meters and center to 0
sim_pts_3D_scaled <-  st_as_sf(sim_pts_3D, coords = c("lon", "lat"), crs=point_crs) %>%
  st_transform(crs=dist_crs) %>% 
  st_coordinates() %>% 
  as_tibble() %>% 
  mutate(X = as.numeric(scale(X, scale=FALSE)),
         Y = as.numeric(scale(Y, scale = FALSE)))

#convert scaled points to sf object for plotting
sim_pts_3D_scaled_sf <- sim_pts_3D_scaled %>%  st_as_sf(coords = c("X", "Y"), crs=dist_crs)

#calculate the radius of the circle encompassing 95% of the simulated data
radius <- quantile(sqrt(sim_pts_3D_scaled$X^2 + sim_pts_3D_scaled$Y^2), 0.95)

#plot simulated data with circle around 95% of data
# sim_plot_3D <- ggplot() +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   geom_sf(data =sim_pts_3D_scaled_sf, aes(geometry=geometry)) +
#   geom_sf(data=ref_point_utm, aes(geometry=geometry, col="red")) +
#   ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = radius), color = "red", fill = NA) +
#   coord_sf(datum=st_crs(dist_crs)) +
#   labs(x = "Error easting (m)", y="Error northing (m)", title = paste0("Simulated Error of 3D GPS Locations n = ",n_3D)) +
#   theme(legend.position = "none") +
#   xlim(-100, 100) +
#   ylim (-100, 100)



################################ Step 2: Fit t-distribution to simulated 3D fixes to estimate parameters for the state-space model #################################

#doesn't seem to fit any common distributions
# descdist(sim_pts_3D$lon, discrete = FALSE)
# hist(sim_pts_3D$lon, freq = FALSE)
# lines(density(sim_pts_3D$lon), add= TRUE)

#sim_pts_3D_utm <- sim_pts_3D_sf %>% st_coordinates %>% as_tibble()

#ref_pt_utm2 <- tibble(X = -123.717987, Y = 47.801262) %>% st_as_sf(coords = c("X", "Y"), crs=4326) %>% st_transform(crs=5070)

#with mean fixed to true location
t_fit_lon_3D <- fitdist(sim_pts_3D$lon, "lst",
                        method = "mle",
                        start = list(sigma = 1, 
                                     df = 3), 
                        fix.arg = list(mu = ref_lon),
                        lower=c(0.000097, 1)) #0.00005

# plot(t_fit_lon_3D)
# gofstat(t_fit_lon_3D)


# #without fixed mean (performs worse)
# t_fit_lon_3D <- fitdist(sim_pts_3D$lon, "lst",
#                         start = list(mu = mean(sim_pts_3D$lon),
#                                      sigma = sd(sim_pts_3D$lon),
#                                      df = 1),
#                         lower=c(-125, 0.00007, 0))


#try skewed generalized t: fails
# t_fit_lon_3D <- fitdist(sim_pts_3D$lon, "sgt", 
#                         start = list(sigma = sd(sim_pts_3D$lon),
#                                      lambda = -0.5,
#                                      p = 0.5,
#                                      q = 0.5),
#                         fix.arg = list(mu =-123.717987),
#                         lower=c(0, -1, 0, 0))


### Repeat for Longitude ###
t_fit_lat_3D <-  fitdist(sim_pts_3D$lat, "lst", 
                         start = list(sigma = 1, 
                                      df = 1), 
                         fix.arg = list(mu = ref_lat),
                         lower=c(0.000042, 0)) #0.00005


# plot(t_fit_lat_3D)
# gofstat(t_fit_lat_3D)


#likelihood surface plot
# llsurface(sim_pts_3D$lon, distr ="lst", plot.arg = c("df", "sigma"), fix.arg = list(mu = ref_lon), min.arg = c(0, 0), max.arg = c(2, 2))
# points(t_fit_lon_3D$estimate[2], t_fit_lon_3D$estimate[1], pch="x", col="red")



################################ Step 3: Simulate 2D fixes #################################
######## Sample distances from custom Weibull distribution (3D fixes)  ####################

dist_2D <- rweibull(n = n_dist_sample, shape = shape_weib, scale = scale_weib ) #generate more values than needed for later truncation
dist_2D <- dist_2D %>% 
  as_tibble() %>% 
  filter(value <= error_2D_95) %>%
  pull(value) %>% 
  sample(size = n_2D) %>% #sample  values so length is equal to bearing vector
  as_tibble() %>% 
  rename(distance = value)

#check distribution
# dist_2D %>% pull(distance) %>% summary()
# dist_2D %>% pull(distance) %>% hist()


#Sample bearings from custom Von Mises distribution (3D fixes)  
dir_2D <- Rfast::rvonmises(n = n_2D, m = m_vm_2D, k = k_vm_2D, rads = TRUE) %>% 
  as_tibble() %>% 
  rename(bearing = value) %>% 
  mutate(bearing = as.numeric(bearing)) %>% 
  mutate(bearing = (bearing*180)/pi) #convert radians back to degrees for distance function


# dir_2D %>% pull(bearing) %>% summary()
# dir_2D %>% pull(bearing) %>% hist()


#Bind distance and bearing vectors to make matrix for distance simulation 
dist_matrix <- bind_cols(dist_2D, dir_2D)


#Simulate coordinates in lat/long based on simulation matrix and reference point 
sim_pts_2D <- destPoint(p = c(ref_lon, ref_lat), b = dist_matrix$bearing, d = dist_matrix$distance) %>% as_tibble() 

######## Plot Simulated Data ###

#Convert to sf object and reproject to Albers Equal Area for plotting
sim_pts_2D_sf <- st_as_sf(sim_pts_2D, coords = c("lon", "lat"), crs=point_crs) %>% st_transform(crs=dist_crs)


#Reproject simulated points to meters and center to 0
sim_pts_2D_scaled <-  st_as_sf(sim_pts_2D, coords = c("lon", "lat"), crs=point_crs) %>%
  st_transform(crs=dist_crs) %>% 
  st_coordinates() %>% 
  as_tibble() %>% 
  mutate(X = as.numeric(scale(X, scale=FALSE)),
         Y = as.numeric(scale(Y, scale = FALSE)))

#convert scaled points to sf object for plotting
sim_pts_2D_scaled_sf <- sim_pts_2D_scaled %>%  st_as_sf(coords = c("X", "Y"), crs=dist_crs)

#calculate the radius of the circle encompassing 95% of the simulated data
radius <- quantile(sqrt(sim_pts_2D_scaled$X^2 + sim_pts_2D_scaled$Y^2), 0.95)

#plot simulated data with circle around 95% of data
# sim_plot_2D <- ggplot() +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   geom_sf(data =sim_pts_2D_scaled_sf, aes(geometry=geometry)) +
#   geom_sf(data=ref_point_utm, aes(geometry=geometry, col="red")) +
#   ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = radius), color = "red", fill = NA) +
#   coord_sf(datum=st_crs(dist_crs)) +
#   labs(x = "Error easting (m)", y="Error northing (m)", title = paste0("Simulated Error of 2D GPS Locations n = ",n_2D)) +
#   theme(legend.position = "none") +
#   xlim(-1000, 1000) +
#   ylim(-1000, 1000)




################################ Step 4: Fit t-distribution to simulated 2D fixes to estimate parameters for the state-space model #################################
#fit to longitude with mean fixed to true location
t_fit_lon_2D <- fitdist(sim_pts_2D$lon, "lst",
                        method = "mle",
                        start = list(sigma = 1, 
                                     df = 3), 
                        fix.arg = list(mu = ref_lon),
                        lower=c(0.00027, 0)) #0.00005

# plot(t_fit_lon_2D)
# gofstat(t_fit_lon_2D)



### Repeat for Longitude ###
t_fit_lat_2D <-  fitdist(sim_pts_2D$lat, "lst", 
                         start = list(sigma = 1, 
                                      df = 1), 
                         fix.arg = list(mu = ref_lat),
                         lower=c(0.000205, 0)) #0.00005


# plot(t_fit_lat_2D)
# gofstat(t_fit_lat_2D)


#likelihood surface plot
# llsurface(sim_pts_3D$lon, distr ="lst", plot.arg = c("df", "sigma"), fix.arg = list(mu = ref_lon), min.arg = c(0, 0), max.arg = c(2, 2))
# points(t_fit_lon_3D$estimate[2], t_fit_lon_3D$estimate[1], pch="x", col="red")

return(tibble(tau_lon_3D = t_fit_lon_3D$estimate[1],
              nu_lon_3D = t_fit_lon_3D$estimate[2],
              tau_lat_3D = t_fit_lat_3D$estimate[1],
              nu_lat_3D = t_fit_lat_3D$estimate[2],
              tau_lon_2D = t_fit_lon_2D$estimate[1],
              nu_lon_2D = t_fit_lon_2D$estimate[2],
              tau_lat_2D = t_fit_lat_2D$estimate[1],
              nu_lat_2D = t_fit_lat_2D$estimate[2]
              ))

}

### conduct many simulations and average parameter estimates

t_est_dat <- tibble(tau_lon_3D = NA,
                    nu_lon_3D = NA,
                    tau_lat_3D = NA,
                    nu_lat_3D =NA,
                    tau_lon_2D = NA,
                    nu_lon_2D = NA,
                    tau_lat_2D = NA,
                    nu_lat_2D = NA)

#for loop to simulate estimates
for (i in 1:nsim){
  #t_est_dat[i-1,] = t_est()[1,]
  t_est_dat[i,] = t_est()[1,]
  
}

#calculate mean and median parameter estimates across all simulations
means <- t_est_dat %>%
  summarize(across(everything(), mean, .names = "mean_{.col}"))
medians <- t_est_dat %>%
  summarize(across(everything(), median, .names = "median_{.col}"))

###############################################################################  


