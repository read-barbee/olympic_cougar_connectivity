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
n_3D <- 561 #561 is Sager's number
n_2D <- 1144


#Weibull distribution parameters (distance)
shape_weib = 0.19
scale_weib = 18.5

#VonMises distribution parameters (bearing)
m_vm = 3.92699 #mean bearing in radians
k_vm = 1.5


################################ Step 1: Simulate 3D fixes #################################
set.seed(777)

#generate reference point (true location)
ref_point <- st_point(c(x = -123.717987 , y = 47.801262)) #x = 47.801262, y = -123.717987 

######## Sample distances from custom Weibull distribution (3D fixes)  ####################

dist_3D <- rweibull(n = 20000, shape = shape_weib, scale = scale_weib ) #generate more values than needed for later truncation
dist_3D <- dist_3D %>% 
  as_tibble() %>% 
  filter(value <= 73) %>% #truncate to values <= 73
  pull(value) %>% 
  sample(size = n_3D) %>% #sample  values so length is equal to bearing vector
  as_tibble() %>% 
  rename(distance = value)

#check distribution
dist_3D %>% pull(distance) %>% summary()
dist_3D %>% pull(distance) %>% hist()


#Sample bearings from custom Von Mises distribution (3D fixes)  
dir_3D <- Rfast::rvonmises(n = n_3D, m = m_vm, k = k_vm, rads = TRUE) %>% 
  as_tibble() %>% 
  rename(bearing = value) %>% 
  mutate(bearing = as.numeric(bearing)) %>% 
  mutate(bearing = (bearing*180)/pi) #convert radians back to degrees for distance function
  

dir_3D %>% pull(bearing) %>% summary()
dir_3D %>% pull(bearing) %>% hist()


#Bind distance and bearing vectors to make matrix for distance simulation 
dist_matrix <- bind_cols(dist_3D, dir_3D)
  

#Simulate coordinates in lat/long based on simulation matrix and reference point 
sim_pts_3D <- destPoint(p = c(-123.717987, 47.801262), b = dist_matrix$bearing, d = dist_matrix$distance) %>% as_tibble() 


######## Plot Simulated Data ###

#Convert to sf object and reproject to Albers Equal Area for plotting
sim_pts_3D_sf <- st_as_sf(sim_pts_3D, coords = c("lon", "lat"), crs=4326) %>% st_transform(crs=5070)

#create a new reference point at 0,0 in the Albers coordinate system for plotting
ref_point_utm <- tibble(X = 0, Y = 0) %>% 
  st_as_sf(coords = c("X", "Y"), crs=5070)

#Reproject simulated points to meters and center to 0
sim_pts_3D_scaled <-  st_as_sf(sim_pts_3D, coords = c("lon", "lat"), crs=4326) %>%
  st_transform(crs=5070) %>% 
  st_coordinates() %>% 
  as_tibble() %>% 
  mutate(X = as.numeric(scale(X, scale=FALSE)),
         Y = as.numeric(scale(Y, scale = FALSE)))

#convert scaled points to sf object for plotting
sim_pts_3D_scaled_sf <- sim_pts_3D_scaled %>%  st_as_sf(coords = c("X", "Y"), crs=5070)

#calculate the radius of the circle encompassing 95% of the simulated data
radius <- quantile(sqrt(sim_pts_3D_scaled$X^2 + sim_pts_3D_scaled$Y^2), 0.95)

#plot simulated data with circle around 95% of data
ggplot() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_sf(data =sim_pts_3D_scaled_sf, aes(geometry=geometry)) +
  geom_sf(data=ref_point_utm, aes(geometry=geometry, col="red")) +
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = radius), color = "red", fill = NA) +
  coord_sf(datum=st_crs(5070)) +
  labs(x = "Error easting (m)", y="Error northing (m)", title = paste0("Simulated Error of 3D GPS Locations n = ",n_3D)) +
  theme(legend.position = "none")



######## Fit t-distribution to simulated data to estimate parameters for the state-space model ###

t_fit_lon_3D <- fitdist(sim_pts_3D$lon, "lst", 
                        start = list(mu = mean(sim_pts_3D$lon),
                                     sigma = sd(sim_pts_3D$lon), 
                                     df = 1), 
                        lower=c(-125, 0.00007, 0))

t_fit_lat_3D <- fitdist(sim_pts_3D$lat, "t", start = list(df = 1))

fitdist(sim_pts_3D$lon, "t", start = list(df = 1))

llsurface(sim_pts_3D$lon, "t", plot.arg = c())




################################ Step 2: Fit t distribution to simulated fixes #################################



###############################################################################  


