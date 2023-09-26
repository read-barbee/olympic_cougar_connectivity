dat4bugs <- function(data, timestep = 1)
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
if(is.factor(data$date) || is.character(data$date)) data$date <- chron(as.character(data$date), out.format='m/d/y')
else if(class(data$date)[1] != "dates") {stop('Error - date is not in correct format, see comments in dat4bugs')}

if(!is.numeric(data$time)) data$time <- as.numeric(times(as.character(data$time))) * 1440
else if(max(data$time) > 1440 || min(data$time) < 0){stop('Error - time is not in correct format, see comments in dat4bugs')}

newdata <- data

time <- newdata$time
lon <- newdata$lon
lat <- newdata$lat
date <- newdata$date

# group observations into 'windows' of duration = timestep and create index (idx) for WinBUGS
start.date <- date[1]
end.date <- date[nrow(newdata)]
numday <- end.date - start.date
numsteps <- ceiling(numday / timestep)
steps <- seq(start.date, dates(as.numeric(start.date) + numsteps * timestep), by = timestep)
row.out <- sapply(2:numsteps, function(i){row.names(newdata)[date < steps[i] & date >= steps[i-1]]})
tmp <- sapply(1:length(row.out), function(i){if(length(row.out[[i]]) == 0) row.out[[i]] <- NA else row.out[[i]]})
idx <- cumsum(c(1, sapply(tmp, length)))
newdata <- newdata[unlist(tmp), ]
newdata <- data.frame(newdata, row.names = NULL)

# determine fraction of day at which each observation is made, j
j <- c()
for(i in 1:(length(idx) - 1)){
	t <- newdata$time[idx[i]:(idx[i+1]-1)]
	d <- as.numeric(newdata$date[idx[i]:(idx[i+1]-1)] - newdata$date[idx[i]])
	j1 <- (t + d * 1440) / (timestep * 1440)
	j <- c(j, j1)
	}

# set up t-distribution parameters and create look-up values for the different location classess
sigma.lon <- c()
sigma.lat <- c()
nu.lon <- c()
nu.lat <- c()

# parameters for t-distributions (obtained by MLE); estimated from Vincent et al. 2002 Argos data
#   these values are in Appendix A, Table A1
sigma.lon[1] <- 0.2898660
sigma.lon[2] <- 0.3119293
sigma.lon[3] <- 0.9020423
sigma.lon[4] <- 2.1625936
sigma.lon[5] <- 0.5072920
sigma.lon[6] <- 4.2050261

sigma.lat[1] <- 0.1220553
sigma.lat[2] <- 0.2605126
sigma.lat[3] <- 0.4603374
sigma.lat[4] <- 1.607056
sigma.lat[5] <- 0.5105468
sigma.lat[6] <- 3.041276

nu.lon[1] <- 3.070609
nu.lon[2] <- 1.220822
nu.lon[3] <- 2.298819
nu.lon[4] <- 0.9136517
nu.lon[5] <- 0.786954
nu.lon[6] <- 1.079216

nu.lat[1] <- 2.075642
nu.lat[2] <- 6.314726
nu.lat[3] <- 3.896554
nu.lat[4] <- 1.010729
nu.lat[5] <- 1.057779
nu.lat[6] <- 1.331283

# add following code to restrict nu to be >= 2; required for WinBUGS
nu.lon[nu.lon < 2] <- 2
nu.lat[nu.lat < 2] <- 2

# convert sigma's from km back to degrees
sigma.lon <- (sigma.lon/6366.71 * 180)/pi
sigma.lat <- (sigma.lat/6366.71 * 180)/pi

# convert standard errors to precisions (se^-2); required for WinBUGS
itau2.lon <- sigma.lon[newdata$lc] ^ -2
itau2.lat <- sigma.lat[newdata$lc] ^ -2

# Note following code deals with timesteps that have no observations by interpolating a single observation at mid-day
# this could also be addressed by putting priors on the missing data for itau2, nu, and j.
#
# convert NA's in itau to very high standard errors
itau2.lon[is.na(itau2.lat)] <- max(itau2.lon, na.rm=T)
itau2.lat[is.na(itau2.lat)] <- max(itau2.lat, na.rm=T)

# converts NA's in nu to smallest nu allowed by Winbugs
nu.lon <- nu.lon[newdata$lc]
nu.lat <- nu.lat[newdata$lc]

nu.lon[is.na(nu.lon)] <- nu.lon[1]
nu.lat[is.na(nu.lat)] <- nu.lat[1]

# set missing observations to be at mid-day
j[is.na(j)] <- 0.5

return(list(y = cbind(newdata$lon,newdata$lat), itau2 = cbind(itau2.lon,itau2.lat), nu = cbind(nu.lon,nu.lat), idx = idx, j = j, RegN = length(idx) - 1))

}
