#### Addressing habitat bias by imputing locations with ctmcmove ####

# Author: Read Barbee

# Date:2023-06-22 

# Purpose:


#############################################################################

#### Libraries ####
library(tidyverse)
library(sf)
library(terra)
library(ctmcmove)
library(fda)
library(doParallel)
library(foreach)

cores <- 10
registerDoParallel(cores=cores)

#########################################################################
##
## 1. Import and format location data
##
##########################################################################

#all location data
data <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_master_5-16-2023.csv", col_types = list(fix_type = col_character())) %>% 
  mutate(date_time_local = force_tz(date_time_local, tzone="US/Pacific"))

#calculate utm coordinates
data_sf <- data %>% 
  st_as_sf(coords = c("longitude", "latitude"), na.fail = FALSE, remove = FALSE, crs=4326) %>% 
  st_transform(crs = 5070)

#append utm coordinates to original dataframe and remove missing loc
data2 <- data %>% 
  rename(lat_wgs84 = latitude,
         lon_wgs84 = longitude) %>% 
  mutate(lat_utm = st_coordinates(data_sf)[,2],
         lon_utm = st_coordinates(data_sf)[,1], .after=lon_wgs84) %>% 
  mutate(across(lat_utm:lon_utm, ~ ifelse(is.nan(.), NA, .)))

al_dat_full <- data2 %>% 
  filter(animal_id =="Al") %>% 
  mutate(date_time_utc = round_date(date_time_utc, unit = "hour"),
         date_time_local = round_date(date_time_local, unit = "hour"))

#Subset Al's data as a test
al_dat <- al_dat_full %>%  
  na.omit()

#define location and time columns
xyt <-  al_dat %>% dplyr::select(lon_utm, lat_utm, date_time_utc)
xy <- al_dat %>% dplyr::select(lon_utm, lat_utm) %>% as.matrix()
x <- al_dat %>% pull(lon_utm)
y <- al_dat %>% pull(lat_utm)
t <- al_dat %>% pull(date_time_utc)

t_full <- al_dat_full %>% pull(date_time_utc)

#plot points
plot(xyt[,1:2],type="b") #b indicates scatter plot with lines


#########################################################################
##
## 2. Fit Quasi-Continuous Path Model
##
##########################################################################


## Define the knots of the spline expansion (how many segments the spline has).
#Problems with fitting the functional movement model can often be fixed by varying the spacing of the knots.

knots = seq(min(t),max(t),by="2 hours")

## create B-spline basis vectors used to approximate the path
b=create.bspline.basis(c(min(t),max(t)),breaks=knots,norder=3)

## define the sequence of times on which to sample the imputed path
tpred=seq(min(t),max(t),by="2 hours")

## Fit latent Gaussian model using MCMC (Quasi-Continuous Path)
#can set fixed observation error (sigma) with the sigma.fixed argument
#increasing the mcmc iterations seems to tighten up the path
out=mcmc.fmove(xy,t,b,tpred,QQ="CAR",n.mcmc=1000,a=1,r=1,num.paths.save=30)

## plot 3 imputed paths
plot(xy,type="b")
points(out$pathlist[[1]]$xy,col="red",type="l")
points(out$pathlist[[2]]$xy,col="blue",type="l")
points(out$pathlist[[3]]$xy,col="green",type="l")


paths <- out$pathlist
path_frames = list()
for (i in 1:length(paths)){
  path_frames[[i]] <- tibble(lon_utm = paths[[i]]$xy[,1],
                             lat_utm = paths[[i]]$xy[,2],
                             date_time_utc = paths[[i]]$t)
}


paths_sf <- list()
for(i in 1:length(paths)){
  paths_sf[[i]] <- path_frames[[i]] %>%
    mutate(path_id = i) %>%
    st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070)
}


#fill missing locations with imputed points from ctmcmove (remove rows before first successful fix)
al_dat_imp<- al_dat_full[-(1:2),] %>%
  mutate(source = case_when(is.na(lat_utm) ~ "imputed",
                            .default = source),
         fix_type = case_when(is.na(lat_utm) ~ "imputed",
                            .default = fix_type),
         imp_status = case_when(is.na(lat_utm) ~ "imputed",
                            .default = "observed"),
         dop = case_when(is.na(lat_utm) ~ NA,
                            .default = dop),
         lat_utm = case_when(is.na(lat_utm) ~ path_frames[[1]]$lat_utm,
                       .default = lat_utm),
         lon_utm = case_when(is.na(lon_utm) ~ path_frames[[1]]$lon_utm,
                       .default = lon_utm))


#visualize imputed points
al_sf2 <- al_dat_imp %>% st_as_sf(coords=c("lon_utm", "lat_utm"), crs=5070)

#paths_sf_combined <- bind_rows(al_sf2, paths_sf[[1]])

mapview::mapview(al_sf2, zcol = "imp_status")

################################ I can probably stop here for my needs #################################

#Might be able to use the remaining portion for help fitting Muff Poisson models



#########################################################################
##
## 2. import rasters
##
##########################################################################

elev <-  rast("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Habitat_Covariates/puma_cov_stack_v1/elevation_albers_30m.tif")

plot(elev)

al_sf <- al_dat %>% 
  dplyr::select(lon_utm, lat_utm) %>% 
  st_as_sf(coords=c("lon_utm", "lat_utm"), crs=5070)

convex_hull <- st_convex_hull(al_sf) %>% st_buffer(dist=10000)

elev_cropped <- crop(elev, convex_hull)

#set all values of sst raster to 1
int=raster::raster(elev_cropped)
values(int) <- 1

##########################################################################
##
## 3 Impute Quasi-Continuous Paths
##
##########################################################################

#number of paths to impute
P=3 #20

plot(elev_cropped,col=grey.colors(100))
for(i in 1:P){
  points(out$pathlist[[i]]$xy,col=i,type="l",lwd=2)
}
points(xyt[,1:2],type="b",pch=20,cex=2,lwd=2)


##########################################################################
##
## 4. Turn continuous space path into a CTMC discrete space path
##
##########################################################################

#can set water mask values with zero.idx parameter

path=out$pathlist[[1]]
ctmc=path2ctmc(path$xy,as.numeric(path$t), int, method="LinearInterp")
## alternate method, useful if you have impassible barriers, but slower
## ctmc=path2ctmc(path$xy,path$t,int,method="ShortestPath")

str(ctmc)

##########################################################################
##
## 5. Turn CTMC discrete paths into latent Poisson GLM data
##
##########################################################################

#stack of rasters of static location-based covariates
loc.stack=raster::stack(int, raster::raster(elev_cropped))
names(loc.stack) <- c("intercept","elevation")

#stack of rasters of directional gradient-based covariates
grad.stack = raster::raster(elev_cropped)

#initialize list of glm objects and create first element
#glm.list=list()
#glm.list[[1]]=ctmc2glm(ctmc,loc.stack,grad.stack)


#Convert each previously generated path to CTMC discrete space path and then to glm format
glm.list <- foreach(i = 1:P) %dopar% {
  cat(i," ")
  path=out$pathlist[[i]]
  ctmc=path2ctmc(path$xy,as.numeric(path$t),int,method="LinearInterp")
  glm.list[[i]]=ctmc2glm(ctmc,loc.stack,grad.stack)
}

## remove transitions that are nearly instantaneous
##  (These are essentially outliers in the following regression analyses)
for(i in 1:P){
  idx.0=which(glm.list[[i]]$tau<10^-5)
  if(length(idx.0)>0){
    glm.list[[i]]=glm.list[[i]][-idx.0,]
  }
  glm.list[[i]]$t=glm.list[[i]]$t-min(glm.list[[i]]$t)
}

## Stack the P imputations together
##

glm.data=glm.list[[1]]
for(i in 2:P){
  glm.data=rbind(glm.data,glm.list[[i]])
}


##########################################################################
##
## 6. Fit Poisson GLM
##    (here we are fitting all "M" paths simultaneously,
##     giving each one a weight of "1/M")
##
##########################################################################

fit.SWL=glm(z~elevation,
            weights=rep(1/P,nrow(glm.data)),family="poisson",offset=log(tau),data=glm.data)
summary(fit.SWL)

beta.hat.SWL=coef(fit.SWL)
beta.se.SWL=summary(fit.SWL)$coef[,2]


##########################################################################
##
## 6. Fit Poisson GLM
##    (here we are fitting using Multiple Imputation)
##
##########################################################################

## Fit each path individually
glm.fits=list()
glm.fits <- foreach(i = 1:P) %dopar% {
  glm.fits[[i]]=glm(z~elevation,
                    family="poisson",offset=log(tau),data=glm.list[[i]])
}

## get point estimates and sd estimates using Rubin's MI combining rules
beta.hat.mat=integer()
beta.se.mat=integer()

for(i in 1:P){
  beta.hat.mat=rbind(beta.hat.mat,coef(glm.fits[[i]]))
  beta.se.mat=rbind(beta.se.mat,summary(glm.fits[[i]])$coef[,2])
}

beta.hat.mat
beta.se.mat

## E(beta) = E_paths(E(beta|path)): mean coefficient estimates
beta.hat.MI=apply(beta.hat.mat,2,mean)
beta.hat.MI

## Var(beta) = E_paths(Var(beta|path))+Var_paths(E(beta|path)): combined variance of coeff estimates
beta.var.MI=apply(beta.se.mat^2,2,mean)+apply(beta.hat.mat,2,var)
beta.se.MI=sqrt(beta.var.MI) #compute standard error from variance

#make table of combined estimates
cbind(beta.hat.MI,beta.se.MI)

##
## compare estimates from MI and Stacked Weighted Likelihood approach
##

## standardize regression coefficients by multiplying by the SE of the X matrix
sds=apply(model.matrix(fit.SWL),2,sd)
sds[1]=1

## plot MI and SWL regression coefficients (Multiple Imputation and Stacked Weighted Likelihood)
par(mfrow=c(1,2))
plot(beta.hat.MI*sds,beta.hat.SWL*sds,main="(a) Coefficient Estimates",
     xlab="Weighted Likelihood Coefficient",
     ylab="Multiple Imputation Coefficient",pch=20,cex=2)
abline(0,1,col="red")
plot(log(beta.se.MI),log(beta.se.SWL),
     main="(b) Estimated log(Standard Errors)",xlab="Weighted Likelihood log(SE)",
     ylab="Multiple Imputation log(SE)",pch=20,cex=2)
abline(0,1,col="red")




###############################################################################  