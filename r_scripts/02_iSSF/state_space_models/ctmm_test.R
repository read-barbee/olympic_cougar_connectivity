library(ctmm)

#ctmm.boot() could be helpful here, but it takes a long time

#meta() meta analysis of multiple movement models with error estimation
#mean.ctmm() averages multiple movement models
#mean.variogram()

#convert Al's data to movebank format, then ctmm telemetry format

all_locs <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_master_6-29-2023.csv", col_types = list(fix_type = col_character())) %>% 
  mutate(date_time_local = with_tz(date_time_local, tzone = "US/Pacific"))

al_dat_full <- all_locs %>% filter(animal_id =="Al")

al_telem <- al_dat_full %>% 
  dplyr::select(animal_id:date_time_utc, lat_wgs84,lon_wgs84) %>% 
  rename(individual.local.identifier = animal_id,
         tag.local.identifier = collar_id,
         timestamp = date_time_utc,
         location.long = lon_wgs84,
         location.lat = lat_wgs84) %>% 
  as.telemetry(projection = "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs")

#plot of sampling intervals
dt.plot(al_telem)

#create and plot variograms
SVF <- variogram(al_telem)
level <- c(0.5,0.95) # 50% and 95% CIs
xlim <- c(0,12 %#% "hour") # 0-12 hour window
plot(SVF,xlim=xlim,level=level)
title("zoomed in")
plot(SVF,fraction=0.65,level=level)
title("zoomed out")

#fit a variogram with gui. Save as GUESS
variogram.fit(SVF)

#Perform model selection to find which model fits the variogram best
FITZ <- ctmm.select(al_telem,GUESS,verbose=TRUE,cores=8)
summary(FITZ)

fit <- ctmm.fit(al_telem)

boot <- ctmm.boot(al_telem, FITZ[[1]], cores=8)

#none of them fit well because all of them assume home range movement, but Al is dispersing
plot(SVF,CTMM=FITZ,col.CTMM=c("red","purple","blue"),fraction=0.65,level=0.5)
title("zoomed out")
plot(SVF,CTMM=FITZ,col.CTMM=c("red","purple","blue"),xlim=xlim,level=0.5)
title("zoomed in")

#emulate function. Not sure why I would need this
emulate_test <- emulate(FITZ[[1]], data= al_telem, fast = TRUE)


meta_test <- meta(FITZ)

FITZ[[1]]$error <- FALSE

#generate continuous time occurrence distribution. Can add covariate constraints, polygon boundaries, etc.
occ_test <- occurrence(data = al_telem, CTMM = FITZ[[1]])

#high speeds in blue and distant locations in red (designed for home range-based analysis)--not super useful
outlie(al_telem)
plot(periodogram(al_telem))
plot.variogram(al_telem)

#split data with time slider
cleave(al_telem)

sim_test <- simulate(object = fit, data =al_telem, complete=TRUE)

al_telem_sf <- al_telem %>% as_tibble() %>% st_as_sf(coords = c("x", "y"), crs = 5070) %>% 
  mutate(type= "og")

sim_test_sf <- sim_test %>% as_tibble() %>% st_as_sf(coords = c("x", "y"), crs = 5070) %>% 
  mutate(type ="sim")

mapview::mapview(bind_rows(sim_test_sf, al_telem_sf), zcol="type")


