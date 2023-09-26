dat4jags <- function(data, timestep = 1)
{
#
#	Code to convert typical Argos data into proper format for WinBUGS
#
#	Created by Ian Jonsen
#
#		created on: 01/04/2005
#	  last modified on: 05/26/2005
#
#
#		Note data file must include the following named columns in the first 5 positions: date, time, lc, lat, lon
# 		
#			date must be a numeric object of class "dates times", this can be accomplished as follows:
#				data$date <- chron(as.character(data$date), format='d/m/y', out.format='m/d/y')
  #   data$date <- chron(as.character(data$date), format='m/d/y', out.format='m/d/y')
#				where format (eg, ='d/m/y') must be set to whatever format the date is in, and out.format is set to 'm/d/y'
#
#			time must be converted to "minutes after midnight" format, this can be accomplished as follows:
#				data$time <- as.numeric(times(as.character(data$time))) * 1440				
#
#		timestep units are in days
#
#
#		try to convert date and time to correct format if not already so
#

library(chron)
library(tidyverse)
library(lubridate)
  
example_data <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/Reference Literature/Spatial data sampling and management/Correcting for habitat bias/State Space Models/Jonsen Supplements/Supplement1/hseal_tst.csv")

data <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location_Data/Source_Files/locations_master/gps_locs_master_5-16-2023.csv", col_types = list(fix_type = col_character())) %>% 
  filter(animal_id == "Al")

#1853 observations
#1526 NAs

#3379 total lines

#separate date and time columns and select for relevant columns
data <- data %>% 
  mutate(date_time_local = force_tz(date_time_local, tz = "US/Pacific"),
         date = date(date_time_local),
         time = hms::as_hms(date_time_local)) %>% 
  dplyr::select(date, time, fix_type, latitude, longitude, date_time_local, deployment_id)

#convert time column to minutes after midnight
data <- data %>% mutate(time_mam = as.numeric(hour(time) * 60 + minute(time) + second(time) / 60)) %>% filter(!is.na(latitude))

newdata <- data

time <- newdata$time
lon <- newdata$longitude
lat <- newdata$latitude
date <- newdata$date
date_time <- newdata$date_time_local

interval_hours <- 2

timestep = hours(interval_hours)



# group observations into 'windows' of duration = timestep and create index (idx) for WinBUGS
start_time <- round_date(min(date_time), unit = "hour")
end_time <- round_date(max(date_time), unit = "hour")
time_interval <- end_time - start_time
n_intervals <- ceiling(as.period(time_interval) / timestep)
steps <- seq(from = start_time, to =(start_time + n_intervals * timestep), by = "2 hours")
row.out <- sapply(2:n_intervals, function(i){row.names(newdata)[date_time < steps[i] & date_time >= steps[i-1] & is.na(newdata$latitude)==FALSE]}) #assigns observed locations to regular intervals
tmp <- lapply(1:length(row.out), function(i){if(length(row.out[[i]]) == 0) row.out[[i]] <- NA else row.out[[i]]})# set regular intervals with no observations to NA, else use what's already there (the observation)
idx <- cumsum(c(1, sapply(tmp, length))) #cumulative sum of number of observations in each regular interval
newdata <- newdata[unlist(tmp), ] #adding rows for regular intervals (days) with missing observations
newdata <- data.frame(newdata, row.names = NULL)

#data <- data %>% mutate(unique_id = 1:nrow(data), .before = date)

#idx <- data$unique_id


# j isn't actually necessary because there's already a data point for every time step even if there's no actual fix. 

# determine fraction of the regular interval at which each observation is made, j
j <- c()
for(i in 1:(length(idx) - 1)){
	t <- newdata$time[idx[i]:(idx[i+1]-1)]
	d <- as.numeric(newdata$date[idx[i]:(idx[i+1]-1)] - newdata$date[idx[i]])
	j1 <- (t + d * 1440) / (timestep * 1440)
	j <- c(j, j1)
}

# t <- newdata$time_mam[idx[3]:(idx[3+1]-1)]
# d <- as.numeric(newdata$date[idx[1]:(idx[1+1]-1)] - newdata$date[idx[1]])
# j1 <- (t + interval_hours * 60) / (interval_hours * 60)
# j <- c(j, j1)

# set up t-distribution parameters and create look-up values for the different location classess
sigma.lon <- c()
sigma.lat <- c()
nu.lon <- c()
nu.lat <- c()

# parameters for t-distributions (obtained by MLE); estimated from Sager et al 2007

#   these values are in Appendix A, Table A1
sigma.lon[1] <- 9.700003e-05 #3D
sigma.lon[2] <- 0.00027 #2D


sigma.lat[1] <- 4.200001e-05 #3D
sigma.lat[2] <- 0.000205 #2D


nu.lon[1] <- 2.460458 #3D
nu.lon[2] <- 2.138559 #2D


nu.lat[1] <- 1.653926 #3D
nu.lat[2] <- 2.397886 #2D


# add following code to restrict nu to be >= 2; required for WinBUGS
# nu.lon[nu.lon < 2] <- 2
# nu.lat[nu.lat < 2] <- 2

# convert sigma's from km back to degrees-- **should already be in degrees
# sigma.lon <- (sigma.lon/6366.71 * 180)/pi
# sigma.lat <- (sigma.lat/6366.71 * 180)/pi

# convert standard errors to precisions (se^-2); required for WinBUGS

newdata <- newdata %>% mutate(itau2.lon = case_when(fix_type =="3D" ~ sigma.lon[1]^-2,
                                      fix_type =="2D"~ sigma.lon[2]^-2,
                                      is.na(fix_type) ~ max(sigma.lon[2]^-2)),
                itau2.lat = case_when(fix_type =="3D" ~ sigma.lat[1]^-2,
                                      fix_type =="2D"~ sigma.lat[2]^-2,
                                      is.na(fix_type) ~ max(sigma.lat[2]^-2)),
                nu.lon = case_when(fix_type =="3D" ~ nu.lon[1],
                                      fix_type =="2D"~ nu.lon[2],
                                      is.na(fix_type) ~ 2),
                nu.lat = case_when(fix_type =="3D" ~ nu.lat[1],
                                   fix_type =="2D"~ nu.lat[2],
                                   is.na(fix_type) ~ 2))

# Note following code deals with timesteps that have no observations by interpolating a single observation at mid-day
# this could also be addressed by putting priors on the missing data for itau2, nu, and j.

# convert NA's in nu to smallest nu allowed by Winbugs


# set missing observations to be at mid-day
#j[is.na(j)] <- 0.5

return(list(y = cbind(newdata$longitude, newdata$latitude),
            itau2 = cbind(newdata$itau2.lon, newdata$itau2.lat),
            nu = cbind(newdata$nu.lon, newdata$nu.lat),
            idx = idx,
            j = rep(0, nrow(newdata)), 
            RegN = nrow(newdata)-1))

test_dat <- list(y = cbind(newdata$longitude, newdata$latitude),
                 itau2 = cbind(newdata$itau2.lon, newdata$itau2.lat),
                 nu = cbind(newdata$nu.lon, newdata$nu.lat),
                 idx = idx,
                 j = rep(0, nrow(newdata)), 
                 RegN = nrow(newdata)-1)


#original form
# return(list(y = cbind(newdata$lon,newdata$lat), 
#             itau2 = cbind(itau2.lon,itau2.lat), 
#             nu = cbind(nu.lon,nu.lat), 
#             idx = idx, 
#             j = j, 
#             RegN = length(idx) - 1))

}
