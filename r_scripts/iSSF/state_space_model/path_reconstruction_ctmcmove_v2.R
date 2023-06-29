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
library(pathroutr)

cores <- 10
registerDoParallel(cores=cores)

n_paths = 3

latlon_crs <- 4326
utm_crs <- 5070

#########################################################################
##
## 1. Import and format location data
##
##########################################################################

#all location data
data <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_master_6-23-2023.csv", col_types = list(fix_type = col_character())) %>% 
  mutate(date_time_local = force_tz(date_time_local, tzone="US/Pacific"))

al_dat_full <- data %>% 
  filter(animal_id =="Al") %>% 
  mutate(date_time_utc = round_date(date_time_utc, unit = "hour"),
         date_time_local = round_date(date_time_local, unit = "hour"))

#Subset Al's data as a test
al_dat <- al_dat_full %>%  
  na.omit()

al_sf <- al_dat %>% st_as_sf(coords= c("lon_utm", "lat_utm"), crs = 5070, remove=FALSE) %>% 
  dplyr::select(animal_id, date_time_utc, lon_utm, lat_utm) %>% 
  rename(id = animal_id, date = date_time_utc, x=lon_utm, y=lat_utm) %>% 
  mutate(type = "original") %>% 
  st_cast("MULTIPOINT")

#define location and time columns
xyt <-  al_dat %>% dplyr::select(lon_utm, lat_utm, date_time_utc)
xy <- al_dat %>% dplyr::select(lon_utm, lat_utm) %>% as.matrix()
x <- al_dat %>% pull(lon_utm)
y <- al_dat %>% pull(lat_utm)
t <- al_dat %>% pull(date_time_utc)

t_full <- al_dat_full %>% pull(date_time_utc)

#plot points
#plot(xyt[,1:2],type="b") #b indicates scatter plot with lines

##################################################################
##
## 2. Import and format barrier polygons
##
##################################################################

#Import water body polygons
water_polys <- st_read("data/Habitat_Covariates/washington_water_polygons/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation.shp") %>% st_transform(crs = 5070)

water_polys_filtered <- water_polys %>% filter(WB_PERIOD_ =="PER" ) #|WB_PERIOD_ =="INT"

#Filter water body polygons to retain those within convex hull of all of Al's points
water_polys_cropped <- sf::st_buffer(al_sf, dist = 10000) %>%
  sf::st_union() %>%
  sf::st_convex_hull() %>%
  sf::st_intersection(water_polys_filtered) %>%
  st_collection_extract('POLYGON') %>%
  st_union() %>% 
  st_sf()

#view the water polygon features retained
mapview::mapview(water_polys_cropped)


#########################################################################
##
## 3. Fit Quasi-Continuous Path Model
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
out=mcmc.fmove(xy,t,b,tpred,QQ="CAR",n.mcmc=1000,a=1,r=1, num.paths.save=n_paths)

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


##################################################################
##
## 4. Re-route path so predicted points don't end up in the middle of water bodies
##
##################################################################

#create buffer around barrier objects as visgraph for rerouting function (connects all verticies of barrier polygon with Delaunay triangle mesh and removes any edges that cross the barrier). Essentially it creates a roadmap of traversible terrain.
visgraph <- pathroutr::prt_visgraph(water_polys_cropped, buffer = 15)

#Create table of all consecutive track points that intersect with a barrier polygon and calculate the shortest path through the visibility network between the non-intersecting points.
# segs_tbl <- pathroutr::get_barrier_segments(predicted_sf2, water_polys_cropped) %>% 
#   prt_shortpath(visgraph, blend=TRUE)
# This is simplified by using the prt_reroute() funciton below:

#Reroute the path based on the visibility network
# rerouted <- pathroutr::prt_reroute(paths_sf[[1]], water_polys_cropped, visgraph, blend = TRUE) %>% 
#   pathroutr::prt_update_points(paths_sf[[1]])



path_sims <- foreach(i = 1:n_paths) %dopar% {
  cat(i," ")
  path_sims <- list()
  path <- paths_sf[[i]]
  path_sims[[i]] <- prt_reroute(path, water_polys_cropped, visgraph) %>% 
    prt_update_points(path) 
}

mapview::mapview(bind_rows(path_sims), zcol="path_id")


##################################################################
##
## 5. Fill in missing locations using rerouted state space model
##
##################################################################

#fill missing locations with imputed points from ctmcmove (remove rows before first successful fix)

# Find the index of the first non-NA value in col1
# first_non_na <- which(!is.na(al_dat_full$lat_utm))[1]
# 
# # Remove rows before the first non-NA value
# al_dat_trimmed <- al_dat_full %>%
#   slice(first_non_na:n())

al_dat_trimmed <- al_dat_full %>% 
  filter(date_time_utc >= min(path_sims[[1]]$date_time_utc) & date_time_utc <= max(path_sims[[1]]$date_time_utc))


combine_imputed2 <- foreach(i = 1:length(path_sims)) %dopar% {
  
  #initialize output list
  dat_imp <- list()
  
  #extract coordinates of imputed points for replacement in al_dat dataframe
  rerouted_coords_utm <- path_sims[[i]] %>% st_coordinates() %>% as_tibble()
  rerouted_coords_latlon <- path_sims[[i]] %>% st_transform(crs=latlon_crs) %>% st_coordinates() %>% as_tibble()
  
  #Combine observed and imputed points
  dat_imp[[i]]<- al_dat_trimmed %>%
    mutate(source = case_when(is.na(lat_utm) ~ "imputed", .default = source),
           fix_type = case_when(is.na(lat_utm) ~ "imputed", .default = fix_type),
           imp_status = case_when(is.na(lat_utm) ~ "imputed", .default = "observed"),
           dop = case_when(is.na(lat_utm) ~ NA, .default = dop),
           #geometry = case_when(is.na(lat_utm) ~ rerouted$geometry,.default = geometry),
           lat_utm = case_when(is.na(lat_utm) ~ rerouted_coords_utm$Y,.default = lat_utm),
           lon_utm = case_when(is.na(lon_utm) ~ rerouted_coords_utm$X, .default = lon_utm),
           lat_wgs84 = case_when(is.na(lat_wgs84) ~ rerouted_coords_latlon$Y,.default = lat_wgs84),
           lon_wgs84 = case_when(is.na(lon_wgs84) ~ rerouted_coords_latlon$X, .default = lon_wgs84),
           path_id = i)
  
  
}


combine_imputed <- function(path_sim, og_dat_trimmed){
 
  #extract coordinates of imputed points for replacement in al_dat dataframe
  rerouted_coords <- path_sim %>% st_coordinates() %>% as_tibble()
  
  #Combine observed and imputed points
  dat_imp<- og_dat_trimmed %>%
    mutate(source = case_when(is.na(lat_utm) ~ "imputed", .default = source),
           fix_type = case_when(is.na(lat_utm) ~ "imputed", .default = fix_type),
           imp_status = case_when(is.na(lat_utm) ~ "imputed", .default = "observed"),
           dop = case_when(is.na(lat_utm) ~ NA, .default = dop),
           #geometry = case_when(is.na(lat_utm) ~ rerouted$geometry,.default = geometry),
           lat_utm = case_when(is.na(lat_utm) ~ rerouted_coords$Y,.default = lat_utm),
           lon_utm = case_when(is.na(lon_utm) ~ rerouted_coords$X, .default = lon_utm))
  
  
}

test <- map(path_sims, combine_imputed, og_dat_trimmed = al_dat_trimmed)



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
  glm.list <- list()
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




################################ GRAVEYARD #################################

#extract coordinates of imputed points for replacement in al_dat dataframe
rerouted_coords <- path_sims[[1]] %>% st_coordinates() %>% as_tibble()



#Combine observed and imputed points
al_dat_imp<- al_dat_trimmed %>%
  mutate(source = case_when(is.na(lat_utm) ~ "imputed", .default = source),
         fix_type = case_when(is.na(lat_utm) ~ "imputed", .default = fix_type),
         imp_status = case_when(is.na(lat_utm) ~ "imputed", .default = "observed"),
         dop = case_when(is.na(lat_utm) ~ NA, .default = dop),
         #geometry = case_when(is.na(lat_utm) ~ rerouted$geometry,.default = geometry),
         lat_utm = case_when(is.na(lat_utm) ~ rerouted_coords$Y,.default = lat_utm),
         lon_utm = case_when(is.na(lon_utm) ~ rerouted_coords$X, .default = lon_utm))



###############################################################################  