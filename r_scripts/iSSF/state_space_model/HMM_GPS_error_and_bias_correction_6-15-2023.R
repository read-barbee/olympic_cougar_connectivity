#### Addressing location error and habitat bias with momentuHMM ####

# Author: Read Barbee

# Date:2023-06-15 

# Purpose:

#Albers Equal Area (5070) or WGS84/UTM zone 10N? (32610) for projected coords?

#### Libraries ####
library(tidyverse)
library(momentuHMM)
library(sf)


data <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_master_5-16-2023.csv", col_types = list(fix_type = col_character())) %>% 
  mutate(date_time_local = force_tz(date_time_local, tzone="US/Pacific"))

#calculate utm coordinates
data_sf <- data %>% 
  na.omit() %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs=4326) %>% 
  st_transform(crs = 5070)

#append utm coordinates to original dataframe
data2 <- data %>% 
  na.omit() %>% 
  rename(lat_wgs84 = latitude,
         lon_wgs84 = longitude) %>% 
  mutate(lat_utm = st_coordinates(data_sf)[,2],
         lon_utm = st_coordinates(data_sf)[,1], .after=lon_wgs84)



# set up t-distribution parameters and create look-up values for the different location classess
sigma.lon <- c()
sigma.lat <- c()
nu.lon <- c()
nu.lat <- c()
sd.lon <- c()
sd.lat<- c()

# parameters for t-distributions (obtained by MLE); estimated from Sager et al 2007

# sigma (tau) and nu estimates in meters
sigma.lon[1] <- 0.1530714 #3D
sigma.lon[2] <- 0.0002700001 #2D


sigma.lat[1] <- 5.274257e-05 #3D
sigma.lat[2] <- 0.0002059334 #2D


nu.lon[1] <- 1 #3D
nu.lon[2] <- 0.12603 #2D


nu.lat[1] <- 0.1151372 #3D
nu.lat[2] <- 0.1197635 #2D

sd.lon[1] <- 	11.26639 #3D
sd.lon[2] <- 34.4252 #2D

sd.lat[1] <- 	10.22842 #3D
sd.lat[2] <-  37.18543#2D



#Add error values based on t distribution
data2<- data2 %>% mutate(sd.lon = case_when(fix_type =="3D" ~ sd.lon[1],
                                                    fix_type =="2D"~ sd.lon[2],
                                                    is.na(fix_type) ~ max(sd.lon[2])),
                              sd.lat = case_when(fix_type =="3D" ~ sd.lat[1],
                                                    fix_type =="2D"~ sd.lat[2],
                                                    is.na(fix_type) ~ max(sd.lat[2])))


#subset only locations for Al and set error correlation to 0
data_al <- data2 %>% 
  filter(animal_id=="Al") %>% 
  select(animal_id, date_time_utc, lon_utm, lat_utm, sd.lon, sd.lat) %>% 
  rename(ID = animal_id,
         time = date_time_utc,
         x = lon_utm,
         y = lat_utm) %>% 
  mutate(error.corr = 0,
         time = round_date(time, unit = "hour"))


#Define prior distributions for errors of 2D and 3D fixes

# priors <-  function(p) { 
#     dt(p[1], sigma.lon[1], nu.lon[1] , log = FALSE) + #3D lon
#     dt(p[2], sigma.lat[1], nu.lat[1] , log = FALSE) + #3D lat
#     dt(p[3], sigma.lon[2], nu.lon[2], log = FALSE) + #2D lon
#     dt(p[4], sigma.lat[2], nu.lat[2], log = FALSE) + #2D lat
#     dnorm(p[5], -4, 2 , log = FALSE) #beta
# }


#Fit a continuous time corelated random walk (CTCRW) model to estimate locations for each regular time interval (2 hours)
crawl_hmm <- crawlWrap(obsData = data_al,
                       timeStep = "2 hours",
                       err.model = list(x = ~ log(sd.lon)-1,
                                        y = ~ log(sd.lat) -1,
                                        rho = ~ error.corr),
                       initialSANN = list(
                         maxit = 1500,
                         trace = 0),
                       attempts = 15,
                       ncores = 3)

crawl_hmm_naive<- crawlWrap(obsData = data_al,
                       timeStep = "2 hours",
                       initialSANN = list(
                         maxit = 1500,
                         trace = 0),
                       attempts = 15,
                       ncores = 3)

#bsam not working

# al_bsam <- data_al %>% 
#   mutate(lc= "G", .after = date_time_utc) %>% 
#   rename(id = animal_id, date = date_time_utc, lon = lon_wgs84, lat = lat_wgs84) %>% 
#   select(id, date, lc, lon, lat)
# 
# library(bsam)
# bsam_test <- bsam::fit_ssm(data = al_bsam,
#               model = "DCRW",
#               tstep = 0.08333333,
#               adapt = 10000,
#               samples = 5000,
#               thin = 5,
#               span = 0.2)





#compare predicted locations to original
al_sf <- data_al %>% mutate(type = "observed") %>% st_as_sf(coords=c("x", "y"), crs=5070)
smm_sf <-crawl_hmm$crwPredict %>%  as_tibble() %>% mutate(type = "predicted") %>% st_as_sf(coords=c("mu.x", "mu.y"), crs=5070)

sf_combined <- rbind(smm_sf %>% select(time, type, geometry), al_sf %>% select(time, type, geometry))

mapview::mapview(sf_combined, zcol = "type" )









#format simulated data
prep_dat <- prepData(crawl_hmm)

# add cosinor covariate based on hour of day
prep_dat$hour <- as.integer(strftime(prep_dat$time, format = "%H", tz="UTC"))


#ACF plot of step lengths
acf(prep_dat$step[!is.na(prep_dat$step)],lag.max=300)

#all plots (acf, step length, turning angle distributions)
plot(prep_dat)


################################ HMMs for fun #################################
#fit naive three state HMM model: need diff starting values
# label states
stateNames <- c("resting", "moving","feeding")
# distributions for observation processes
dist = list(step = "gamma", angle = "wrpcauchy")
# initial parameters
Par0_m1 <- list(step=c(100,600,100,500,100,200),angle=c(0.1,0.3,0.7))
# fit model
m1 <- fitHMM(data = prep_dat, nbStates = 3, dist = dist, Par0 = Par0_m1,
             estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 3)


#can add covariates to transition probabilities (i.e. time of day)
formula <- ~ cosinor(hour, period = 24)
# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formula)
# fit model
m2 <- fitHMM(data = prep_dat, nbStates = 3, dist = dist, Par0 = Par0_m2$Par,
             beta0=Par0_m2$beta, stateNames = stateNames, formula=formula)

#assess model fit
# compute pseudo-residuals for the steps and the angles
pr <- pseudoRes(m2)
# plot the ACF of step pseudo-residuals
acf(pr$stepRes[!is.na(pr$stepRes)],lag.max = 300)



###############################################################################  




## orthogonalize W based on locations ----
W_ortho <- W
W_path <- extract(x = W, y = matrix(buffalo_proj@coords, ncol = 2)) 
obstimes <- as.numeric(buffalo_proj$POSIX) / 3600 # numeric hours 
W_tilde <- apply(W_path * c(0, diff(obstimes)), 2, cumsum) 
W_tilde_svd <- svd(W_tilde)
W_tilde_proj_mat <- W_tilde_svd$v %*% diag(W_tilde_svd$d^(-1)) 
W_mat <- as.matrix(W)
W_mat_proj <- W_mat %*% W_tilde_proj_mat
for(layer in 1:ncol(W_mat)){
  values(W_ortho[[layer]]) <- W_mat_proj[, layer]
  names(W_ortho[[layer]]) <- paste0("svd", layer) 
  }

buffaloData <- data.frame(ID = 1,
                          time = obstimes,
                          x = buffalo_proj@coords[, 1],
                          y = buffalo_proj@coords[, 2],
                          ln.sd.x = lnError$ln.sd.x,
                          ln.sd.y = lnError$ln.sd.y,
                          error.corr = lnError$error.corr)













