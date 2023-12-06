library(terra)
library(tidyverse)
library(amt)
library(sf)
library(INLA)


################################ User-defined parameters #################################

scales <- c(seq(30, 90, 10), seq(100, 600, 100))


project_crs=5070

#prior params from https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40#inla-1
mean.beta <- 0
prec.beta <- 1e-4 

#function to make inla formula
make_form_rand_quad <- function(params, quad_params, n_indiv){
  rand_terms_inla <- vector()
  for(i in 1:length(params)){
    rand_terms_inla[i] <- paste0("f(", paste0("id", i), ", ", params[i], ", ", "values = ", paste0("1:",n_indiv), ", ", 'model = "iid", ', "hyper = list(theta = list(initial = log(1), fixed = FALSE, ", 'prior = "pc.prec", ', "param = c(1, .05))))")
  }
  
  if (is.null(quad_params)){
    quad_terms <- NULL
  } else{
    quad_terms <- vector()
    for(i in 1:length(quad_params)){
      quad_terms[i] <- paste0(" + f(", paste0("id", (i + length(params)), ", ", "I(",quad_params[i],"^2)", ", ", "values = ", paste0("1:",n_indiv), ", ", 'model = "iid", ', "hyper = list(theta = list(initial = log(1), fixed = FALSE, ", 'prior = "pc.prec", ', "param = c(3, .05))))"))
    }
  }
  
  form <- as.formula(paste("case_ ~ -1 + ",  
                           #fixed effects
                           paste(c(params, quad_terms), collapse = " + "), "+",
                           #random intercept (strata)
                           "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T))) + ",
                           #random slopes
                           paste(rand_terms_inla, collapse = " + ")))
  
  return(form)
}


#########################################################################
##
## 1. Make 2h steps
##
##########################################################################

locs_screened <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_10-02-2023.csv") 

multi_track <- function(d){
  make_track(d, lon_utm, lat_utm, date_time_utc, check_duplicates = TRUE, all_cols = TRUE, crs = project_crs) 
  #%>% transform_coords(crs_to = project_crs)
}

#Nest locations by animal_id 
locs_nested <- locs_screened %>% 
  nest_by(animal_id, sex, dispersal_status) 


#create tracks
locs_nested$tracks <-map(locs_nested$data, multi_track)

steps_calc <- function(x){
    part1 <- x %>% 
      amt::track_resample(rate = hours(2), tolerance = minutes(15)) %>% #resample 
      amt::filter_min_n_burst() %>% #divide into bursts of min 3 pts
      amt::steps_by_burst() %>%  #convert to steps
      dplyr::filter(sl_ >=100) %>% 
      dplyr::filter(sl_<=20000) #remove stationary and unreasonable step lengths
    blank <- list()
    
    if(nrow(part1) < 3){ return(blank)}
    
    else{
      
      #generate 10 times the amount of desired available steps
      part2 <- part1 %>% 
        amt::random_steps(n_control = 10) %>% #generate random steps per used step
        #amt::extract_covariates(cov_stack, where = "end") %>% #extract covariates at end of step
        amt::time_of_day(include.crepuscule = FALSE) %>% #calculate time of day for each step--not working
        mutate(unique_step = paste(burst_,step_id_,sep="_"), .before=step_id_) %>%  #add column for unique step id }
        mutate(log_sl_ = log(sl_),    
               cos_ta_ = cos(ta_))
    }
    
      return(part2)
}

indiv_to_remove <- c("Andre", "Angie", "Anne", "Brian_Miller", "Cato","Cookie", "Elena", "Freddy Cougar", "Hadleigh", "Hank", "Jasper", "Jax", "Julie", "Junior", "Keelie", "Kiki", "Lyla", "Makah-F01", "Makah-F02", "Makah-F03", "Makah-F04", "Makah-F05", "Makah-M01", "Makah-M02", "Makah-M03", "Melodie", "Moses", "Nick", "Richard", "Roxy", "Salmon River", "Sampson", "Sophie", "Tswift")

locs_nested <- locs_nested %>% 
  filter(!(animal_id %in% indiv_to_remove)) %>% 
  filter(!(animal_id %in% c("Butch", "Crash", "Freya", "Hera", "Otook_Tom", "Promise")))

locs_nested$steps <- map(locs_nested$tracks, steps_calc)



#########################################################################
##
## 2. Make covariate layers at different scales
##
##########################################################################

dem <- rast("/Users/tb201494/Desktop/1km_buffer/static/elevation_1km_buffer.tif")

#replace -999 values with NA
dem2 <-  mask(dem, dem, maskvalues = -999, updatevalue=NA)

#create template rasters for resampling
tmp_rasts <- list()
for(i in 1:length(scales)){
  scale <- scales[i]
  tmp_rasts[[i]] <- rast(vals =1, crs = "EPSG:5070", extent = ext(dem2), resolution = scale)
}
names <- vector()
for(i in 1:length(scales)){
  scale <- scales[i]
  if(scale < 1000){
    names[i] <- paste0("tpi_", scale, "m")
  }
  else{
    names[i] <- paste0("tpi_", (scale/1000), "km")
  }
}

names(tmp_rasts) <- names


#resample original raster to each new resolution

resampled_rasts <- map(tmp_rasts, function(x){return(terra::resample(dem2, x, method = "bilinear"))})


#calculate tpi at each resolution
tpi_scales <- map(resampled_rasts, function(x){return(terra::terrain(x, "TPI"))})


#########################################################################
##
## 3. Create step frames with values extracted at each scale
##
##########################################################################

amt_steps <- locs_nested %>% 
  select(animal_id:dispersal_status,
         steps) %>% 
  filter(length(steps)>0)

#Sample covariates at each scale
step_frames <- list()
step_names <- vector()
for(i in 1:length(scales)){
  amt_steps_copy <- amt_steps
  amt_steps_copy$steps <- map(amt_steps_copy$steps, extract_covariates, covariates = tpi_scales[[i]])
  
  amt_steps_copy$steps <- map(amt_steps_copy$steps, function(x){return(x %>% mutate(name = names[i]) %>% select(name, everything()))})
  step_names[i] <- paste0("steps_", names[i])
  
  step_frames[[i]] <- amt_steps_copy
}

names(step_frames) <- step_names


step_frames_unnest <- map(step_frames, unnest, cols=steps)

n_indiv <- step_frames_unnest[[1]] %>% distinct(animal_id) %>% nrow()

step_frames_unnest2 <- map(step_frames_unnest, function(x){return(x %>%
                                                                    ungroup() %>% 
                                                                    mutate(case_=as.numeric(case_),
                                                                           animal_id = as.numeric(as.factor(animal_id))) %>% 
                                                                    rename(id1=animal_id))})


#########################################################################
##
## 4.Fit univariate models at each scale
##
##########################################################################
inla_form <- make_form_rand_quad(params = c("TPI"), quad_params = NULL, n_indiv = n_indiv)

mods <- list()
mod_names <- vector()
for(i in 1:length(step_frames_unnest2)){
  
  mods[[i]] <- inla(inla_form, family ="Poisson", 
                   data=step_frames_unnest2[[i]],
                   control.fixed = list(
                   mean = mean.beta,
                   prec = list(default = prec.beta)),
                   control.compute = list(waic=TRUE, dic=FALSE, cpo = FALSE),
                   safe = TRUE)
 
  mod_names[[i]] <- step_frames_unnest2[[i]]$name[1]
  
  
  print(paste0(step_frames_unnest2[[i]]$name[1],"   ", i, "/", length(step_frames_unnest2)))
  
}

#########################################################################
##
## 5. Visualize trends in model coefficients with scale
##
##########################################################################

posteriors <- list()

for(i in 1:length(mods)){
  posteriors[[i]] <- mods[[i]]$summary.fixed %>% as_tibble() %>% mutate(name = names[i], .before = mean)
}

posteriors <- bind_rows(posteriors)
posteriors$name <- factor(posteriors$name, levels = names)


indiv_posteriors <- list()
for(i in 1:length(mods)){
  indiv_posteriors[[i]] <- mods[[i]]$summary.random$id1 %>% as_tibble() %>% mutate(name = names[i], .before = mean)
}

indiv_posteriors <- bind_rows(indiv_posteriors)
indiv_posteriors$name <- factor(indiv_posteriors$name, levels = names)




#pop level plots
ggplot(posteriors, aes(x = name, y = mean)) +
  geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`), size = 1) +
  labs(x = "Model", y = "Mean") +
  ggtitle("Mean TPI with 95% Credible Intervals for Each Model") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed

ggplot(posteriors, aes(x = name, y = mean)) +
  geom_point(size = 1) +
  labs(x = "Model", y = "Mean") +
  ggtitle("TPI Mean Pop Coefficients by Scale") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed


#indiv level plots

ggplot(indiv_posteriors %>% filter(ID!=1), aes(x = name, y = mean)) +
  geom_boxplot() +
  labs(x = "Model", y = "Mean") +
  ggtitle("TPI Indiv Mean Coefficients by Scale") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# tri_30m <- terra::terrain(dem2, v="TRI")
# tririley_30m <- terra::terrain(dem2, v="TRIriley")
# trirmsd_30m <- terra::terrain(dem2, v="TRIrmsd")
# roughness_30m <- terra::terrain(dem2, v="roughness")
# 


################################ Graveyard #################################



#alternate loop that adds each scale extraction as a nested column to the step frame
# step_frames <- list()
# step_names <- vector()
# 
# amt_steps_copy <- amt_steps
# for(i in 1:length(scales)){
#   
#   step_name <- paste0("steps_", names[i])
#   amt_steps_copy[[step_name]] <- map(amt_steps_copy$steps, extract_covariates, covariates = tpi_scales[[i]])
#}