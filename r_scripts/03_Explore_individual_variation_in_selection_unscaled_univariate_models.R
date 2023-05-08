#### OCP iSSF Module_03: Explore Inidivual Variation in Selection--Unscaled, Quadratic ####

# Author: Read Barbee

# Date:2023-05-02 

# Inputs:
#   •	GPS2: Unscaled amt data frame for fitting iSSFs
# 
# Outputs:
#   •	Figures
# •	Fitted Models for each individual
# 
# Steps
# •	Import scaled and unscaled amt data frame for fitting iSSFs
# •	Fit global iSSF to each individual
# •	Explore summary stats and figures
# •	Calculate log-RSS and CI for range of each covariate for each individual
# •	Compare proportions of individuals overall and in each demographic category with positive, negative or neutral response to each covariate

#note: removing landuse from the global model didn't change the estimates


################################ Libraries #################################
library(tidyverse)
library(terra)
library(lubridate)
library(amt)
library(janitor)

################################ User-Defined Parameters #################################


############################### Import amt-formatted Data #################################
# Import SSF data for 100 individual mountain lions (November 2022)
steps_raw <- read_csv("6h_steps_cov_4-28-2023.csv")


################################ Check covariate distributions #################################
#pivot hsi longer for graphing
steps_long <- steps_raw %>% 
  pivot_longer(elev_end:landuse_hii_end, names_to="cov", values_to="cov_val")

#examine distributions of covariate values for steps-- landuse looks real bad
steps_raw %>% 
  pivot_longer(elev_end:landuse_hii_end, names_to="cov", values_to="cov_val") %>% 
  ggplot() +
  geom_histogram(aes(x=cov_val))+
  facet_wrap(~cov, scales="free") 

#boxplot of used values for all continuous variables faceted by sex and dispersal status
steps_long %>% 
  filter(case_==TRUE) %>% 
  ggplot() +
  geom_boxplot(aes(x=as.factor(dispersal_status), y=cov_val, fill=sex))+
  facet_wrap(~cov, scales="free") 

#boxplot of used and available values for all continuous variables faceted by sex and dispersal status**
used_unused_dem_cats <- steps_long %>%
  filter(cov!="landuse_hii_end") %>% 
  ggplot() +
  geom_boxplot(aes(x=as.factor(dispersal_status), y=cov_val, fill=sex, alpha=case_), notch = T)+
  facet_wrap(~cov, scales="free")  + scale_alpha_manual(values =c(0.5, 1)) +
  ylab("Covariate Values") +
  xlab("Demographic Category")

#ggsave(filename= "used_unused_dem_cats_5-03-2023.png", plot= used_unused_dem_cats)

used_unused_dem_cats_violin <- steps_long %>%
  filter(cov!="landuse_hii_end") %>% 
  ggplot() +
  geom_violin(aes(x=as.factor(dispersal_status), y=cov_val, fill=sex, alpha=case_))+
  facet_wrap(~cov, scales="free")  + scale_alpha_manual(values =c(0.5, 1)) +
  ylab("Covariate Values") +
  xlab("Demographic Category")

#ggsave(filename= "used_unused_dem_cats_violin_5-03-2023.png", plot= used_unused_dem_cats_violin)


#histogram of used values for all continuous variables faceted by sex and dispersal status
steps_long %>% 
  filter(case_==TRUE) %>% 
  ggplot() +
  geom_histogram(aes(x=cov_val, fill=sex))+
  facet_wrap(~cov, scales="free") 

#histogram of used and available values for all continuous variables faceted by sex and dispersal status
steps_long %>% 
  filter(cov!="landuse_hii_end") %>% 
  ggplot() +
  geom_histogram(aes(x=cov_val, fill=sex, alpha=case_))+
  facet_wrap(~cov, scales="free")  + scale_alpha_manual(values =c(0.5, 1))




######Log transforming landuse, elevation, dist water and ndvi doesn't seem to help much...###
# steps_log_cov <- steps_raw %>%
#   mutate(across(c(dist_water_end, elev_end, landuse_hii_end, ndvi_end), log)) %>%
#   na.omit() %>%
#   pivot_longer(elev_end:landuse_hii_end, names_to="cov", values_to="cov_val")
# 
# 
# #used and available locations
# steps_log_cov %>% 
#   ggplot() +
#   geom_boxplot(aes(x=as.factor(dispersal_status), y=cov_val, fill=sex))+
#   facet_wrap(vars(cov, case_), scales="free") 


################################ Prepare data to fit models #################################

#unscaled step data
steps_unscaled_nested <- steps_raw %>%
  na.omit() %>% 
  nest_by(animal_id, sex, dispersal_status) %>% 
  rename(steps=data) %>%  group_by(sex, dispersal_status)

################################ Fit univariate iSSF to each individual #################################


steps_unscaled_nested$fit_elev <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ elev_end  + strata(step_id_), model = TRUE))

steps_unscaled_nested$fit_ndvi <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ ndvi_end  + strata(step_id_), model = TRUE))

steps_unscaled_nested$fit_dist_water <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ dist_water_end  + strata(step_id_), model = TRUE))

steps_unscaled_nested$fit_roads <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ roads_hii_end  + strata(step_id_), model = TRUE))

steps_unscaled_nested$fit_forest <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ forest_end  + strata(step_id_), model = TRUE))

# steps_unscaled_nested$fit_landuse <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ landuse_hii_end  + strata(step_id_), model = TRUE))


#function to calculate log-rss for each univariate model applied to each individual
l_rss_uni <- function(dat, indiv, curr_param){
  indiv_dat <- steps_unscaled_nested %>% 
    na.omit() %>% 
    filter(animal_id == indiv) %>% 
    unnest(cols=c(steps))
  
  if (curr_param == "elev_end"){
    s1 <- data.frame(
      cov <-seq(from = min(indiv_dat$elev_end, na.rm=T), to = max(indiv_dat$elev_end, na.rm=T), length.out = 200)) %>% 
      rename(elev_end = 1)
    
    
    #data frame with means of all covariates encountered by Al
    s2 <- data.frame(
      elev_end <-mean(indiv_dat$elev_end)) %>% 
      rename(elev_end = 1)
    
    
    indiv_dat_nested <- steps_unscaled_nested %>% 
      na.omit() %>% 
      filter(animal_id == indiv)
    
    mod <- indiv_dat_nested$fit_elev[[1]]
  }
  if (curr_param == "ndvi_end"){
    s1 <- data.frame(
      cov <-seq(from = min(indiv_dat$ndvi_end, na.rm=T), to = max(indiv_dat$ndvi_end, na.rm=T), length.out = 200)) %>% 
      rename(ndvi_end = 1)
    
    
    #data frame with means of all covariates encountered by Al
    s2 <- data.frame(
      elev_end <-mean(indiv_dat$ndvi_end)) %>% 
      rename(ndvi_end = 1)
    
    
    indiv_dat_nested <- steps_unscaled_nested %>% 
      na.omit() %>% 
      filter(animal_id == indiv)
    
    mod <- indiv_dat_nested$fit_ndvi[[1]]
  }
  if (curr_param == "dist_water_end"){
    s1 <- data.frame(
      cov <-seq(from = min(indiv_dat$dist_water_end, na.rm=T), to = max(indiv_dat$dist_water_end, na.rm=T), length.out = 200)) %>% 
      rename(dist_water_end = 1)
    
    
    #data frame with means of all covariates encountered by Al
    s2 <- data.frame(
      elev_end <-mean(indiv_dat$dist_water_end)) %>% 
      rename(dist_water_end = 1)
    
    
    indiv_dat_nested <- steps_unscaled_nested %>% 
      na.omit() %>% 
      filter(animal_id == indiv)
    
    mod <- indiv_dat_nested$fit_dist_water[[1]]
  }
  if (curr_param == "roads_hii_end"){
    s1 <- data.frame(
      cov <-seq(from = min(indiv_dat$roads_hii_end, na.rm=T), to = max(indiv_dat$roads_hii_end, na.rm=T), length.out = 200)) %>% 
      rename(roads_hii_end = 1)
    
    
    #data frame with means of all covariates encountered by Al
    s2 <- data.frame(
      elev_end <-mean(indiv_dat$roads_hii_end)) %>% 
      rename(roads_hii_end = 1)
    
    
    indiv_dat_nested <- steps_unscaled_nested %>% 
      na.omit() %>% 
      filter(animal_id == indiv)
    
    mod <- indiv_dat_nested$fit_roads[[1]]
  }
  if (curr_param == "forest_end"){
    s1 <- data.frame(
      cov <-seq(from = min(indiv_dat$forest_end, na.rm=T), to = max(indiv_dat$forest_end, na.rm=T), length.out = 200)) %>% 
      rename(forest_end = 1)
    
    
    #data frame with means of all covariates encountered by Al
    s2 <- data.frame(
      elev_end <-mean(indiv_dat$forest_end)) %>% 
      rename(forest_end = 1)
    
    
    indiv_dat_nested <- steps_unscaled_nested %>% 
      na.omit() %>% 
      filter(animal_id == indiv)
    
    mod <- indiv_dat_nested$fit_forest[[1]]
  }
  
  ### Working. variable names have to be the same across all data frames and model
  l_rss_indiv <- amt::log_rss(mod, s1, s2, ci = "se", ci_level = 0.95)
  
  return(l_rss_indiv$df)
}

indivs <- steps_unscaled_nested %>% pull(animal_id)

# #add rss tables to main data frame
steps_unscaled_nested$elev_rss_uni <- map(indivs, l_rss_uni, dat=steps_unscaled_nested, curr_param="elev_end")

steps_unscaled_nested$ndvi_rss_uni <-map(indivs, l_rss_uni, dat=steps_unscaled_nested, curr_param="ndvi_end")

steps_unscaled_nested$forest_rss_uni <- map(indivs, l_rss_uni, dat=steps_unscaled_nested, curr_param="forest_end")

steps_unscaled_nested$dist_water_rss_uni <-  map(indivs, l_rss_uni, dat=steps_unscaled_nested, curr_param="dist_water_end")

steps_unscaled_nested$roads_rss_uni <- map(indivs, l_rss_uni, dat=steps_unscaled_nested, curr_param="roads_hii_end")


#pivot data for faceting
plot_dat <- steps_unscaled_nested %>% 
  select(-c(steps, fit_elev:fit_forest)) %>% 
  pivot_longer(elev_rss_uni:roads_rss_uni, names_to= "cov", values_to = "rss_val")

#initialize blank lists for loop
ls = list()
ls2=list()



for(i in 1:nrow(plot_dat)){
  
  if(plot_dat$cov[i] == "elev_rss_uni"){
    ls[[i]] <- plot_dat$rss_val[[i]]$elev_end_x1
    ls2[[i]] = plot_dat$rss_val[[i]]$log_rss
  }
  if(plot_dat$cov[i] == "ndvi_rss_uni"){
    ls[[i]] <- plot_dat$rss_val[[i]]$ndvi_end_x1
    ls2[[i]] = plot_dat$rss_val[[i]]$log_rss
  }
  if(plot_dat$cov[i] == "forest_rss_uni"){
    ls[[i]] <- plot_dat$rss_val[[i]]$forest_end_x1
    ls2[[i]] = plot_dat$rss_val[[i]]$log_rss
  }
  if(plot_dat$cov[i] == "dist_water_rss_uni"){
    ls[[i]] <- plot_dat$rss_val[[i]]$dist_water_end_x1
    ls2[[i]] = plot_dat$rss_val[[i]]$log_rss
  }
  if(plot_dat$cov[i] == "roads_rss_uni"){
    ls[[i]] <- plot_dat$rss_val[[i]]$roads_hii_end_x1
    ls2[[i]] = plot_dat$rss_val[[i]]$log_rss
  }
}


#add lists from for loop to data frame and facet plot by sex (Melodie's RSS values don't make any sense)

plot_dat$cov_vals <- ls
plot_dat$rss_vals <- ls2

rss_elev_uni <- plot_dat %>% 
  # filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  # filter(animal_id!="Sampson") %>%
  # filter(animal_id!="Kingsley") %>%
  select(-rss_val) %>% 
  filter(cov=="elev_rss_uni") %>% 
  filter(sex=="Male") %>% 
  unnest(cols=c(cov_vals, rss_vals)) %>% 
  ggplot(., aes(x = cov_vals, y = rss_vals)) +
  geom_smooth(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Covariate Values") +
  ylab("log-RSS vs Mean Covariate Value") +
  theme_bw() +
  facet_wrap(~animal_id, scales = "free")

plotly::ggplotly(rss_elev_uni)



# #for loop to extract just covariate and log rss-values for each covariate and indiviaual
# for(i in 1:nrow(plot_dat)){
#   ls[[i]] = case_when(
#     plot_dat$cov[i] == "elev_rss_uni" ~ plot_dat$rss_val[[i]]$elev_end_x1,
#     plot_dat$cov[i] == "ndvi_rss_uni" ~ plot_dat$rss_val[[i]]$ndvi_end_x1,
#     plot_dat$cov[i] == "forest_rss_uni" ~ plot_dat$rss_val[[i]]$forest_end_x1,
#     plot_dat$cov[i] == "dist_water_rss_uni" ~ plot_dat$rss_val[[i]]$dist_water_end_x1,
#     plot_dat$cov[i] == "roads_rss_uni" ~ plot_dat$rss_val[[i]]$roads_hii_end_x1)
#   
#   ls2[[i]] = plot_dat$rss_val[[i]]$log_rss
# }

################################ Fit iSSF to each individual #################################

#Global iSSF with step length, turn angle, and quadratics for each covariate
steps_unscaled_nested$fit2 <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ elev_end + I(elev_end^2) + ndvi_end + I(ndvi_end^2) + dist_water_end + I(dist_water_end^2) + roads_hii_end + I(roads_hii_end^2) +forest_end + I(forest_end^2) + landuse_hii_end  + I(landuse_hii_end^2) + sl_ + log(sl_) + cos(ta_) + strata(step_id_), model = TRUE))

#inspect fitted object
steps_unscaled_nested$fit2

#inspect model fit for first individual
steps_unscaled_nested$fit2[[1]]$model


################################ Log RSS for all individuals and covariates #################################

#function to calculate log_rss object for elevation for each individual incorporating quadratic terms, sl and ta
l_rss <- function(dat, indiv, curr_param){
  indiv_dat <- dat %>% 
    na.omit() %>% 
    filter(animal_id == indiv) %>% 
    unnest(cols=c(steps))
  
  #data frame varying elevation from min value to max value encountered by Al, holding all other covariates at the mean
  s1 <- data.frame(
    elev_end <-seq(from = -2, to =2, length.out = 200),
    
    elev_end2 <-mean(indiv_dat$elev_end^2),
    
    ndvi_end <- mean(indiv_dat$ndvi_end),
    
    ndvi_end2 <-mean(indiv_dat$ndvi_end^2),
    
    dist_water_end <- mean(indiv_dat$dist_water_end),
    
    dist_water_end2 <-mean(indiv_dat$dist_water_end^2),
    
    roads_hii_end <- mean(indiv_dat$roads_hii_end),
    
    roads_hii_end2 <-mean(indiv_dat$roads_hii_end^2),
    
    forest_end <- mean(indiv_dat$forest_end),
    
    forest_end2 <-mean(indiv_dat$forest_end^2),
    
    landuse_hii_end <- mean(indiv_dat$landuse_hii_end),
    
    landuse_hii_end2 <-mean(indiv_dat$landuse_hii_end^2),
    
    sl <- mean(indiv_dat$sl_),
    
    log_sl <- mean(log(indiv_dat$sl_)),
    
    cos_ta <- mean(cos(indiv_dat$ta_))
  ) %>% 
    rename(elev_end = 1,
           elev_end2 = 2,
           ndvi_end = 3,
           ndvi_end2 = 4,
           dist_water_end = 5,
           dist_water_end2 = 6,
           roads_hii_end = 7,
           roads_hii_end2 = 8,
           forest_end = 9,
           forest_end2 = 10,
           landuse_hii_end = 11,
           landuse_hii_end2 = 12,
           sl_= 13,
           log_sl = 14,
           ta_ = 15
    )
  
  s1$elev_end <- mean(indiv_dat$elev_end)
  
  
  
  #data frame with means of all covariates encountered by Al
  s2 <- data.frame(
    elev_end <-mean(indiv_dat$elev_end),
    
    elev_end2 <-mean(indiv_dat$elev_end^2),
    
    ndvi_end <- mean(indiv_dat$ndvi_end),
    
    ndvi_end2 <-mean(indiv_dat$ndvi_end^2),
    
    dist_water_end <- mean(indiv_dat$dist_water_end),
    
    dist_water_end2 <-mean(indiv_dat$dist_water_end^2),
    
    roads_hii_end <- mean(indiv_dat$roads_hii_end),
    
    roads_hii_end2 <-mean(indiv_dat$roads_hii_end^2),
    
    forest_end <- mean(indiv_dat$forest_end),
    
    forest_end2 <-mean(indiv_dat$forest_end^2),
    
    landuse_hii_end <- mean(indiv_dat$landuse_hii_end),
    
    landuse_hii_end2 <-mean(indiv_dat$landuse_hii_end^2),
    
    sl <- mean(indiv_dat$sl_),
    
    log_sl <- mean(log(indiv_dat$sl_)),
    
    cos_ta <- mean(cos(indiv_dat$ta_))
  ) %>% 
    rename(elev_end = 1,
           elev_end2 = 2,
           ndvi_end = 3,
           ndvi_end2 = 4,
           dist_water_end = 5,
           dist_water_end2 = 6,
           roads_hii_end = 7,
           roads_hii_end2 = 8,
           forest_end = 9,
           forest_end2 = 10,
           landuse_hii_end = 11,
           landuse_hii_end2 = 12,
           sl_= 13,
           log_sl = 14,
           ta_ = 15
    )
  
  if (curr_param == "elev_end"){
    s1$elev_end <- seq(from = min(indiv_dat$elev_end, na.rm=T), to = max(indiv_dat$elev_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "elev_end2"){
    s1$elev_end2 <- seq(from = min((indiv_dat$elev_end)^2, na.rm=T), to = max((indiv_dat$elev_end)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "ndvi_end"){
    s1$ndvi_end <- seq(from = min(indiv_dat$ndvi_end, na.rm=T), to = max(indiv_dat$ndvi_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "ndvi_end2"){
    s1$ndvi_end2 <- seq(from = min((indiv_dat$ndvi_end)^2, na.rm=T), to = max((indiv_dat$ndvi_end)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "forest_end"){
    s1$forest_end <- seq(from = min(indiv_dat$forest_end, na.rm=T), to = max(indiv_dat$forest_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "forest_end2"){
    s1$forest_end2 <- seq(from = min((indiv_dat$forest_end)^2, na.rm=T), to = max((indiv_dat$forest_end)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "dist_water_end"){
    s1$dist_water_end <- seq(from = min(indiv_dat$dist_water_end, na.rm=T), to = max(indiv_dat$dist_water_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "dist_water_end2"){
    s1$dist_water_end2 <- seq(from = min((indiv_dat$dist_water_end)^2, na.rm=T), to = max((indiv_dat$dist_water_end)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "roads_hii_end"){
    s1$roads_hii_end <- seq(from = min(indiv_dat$roads_hii_end, na.rm=T), to = max(indiv_dat$roads_hii_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "roads_hii_end2"){
    s1$roads_hii_end2 <- seq(from = min((indiv_dat$roads_hii_end)^2, na.rm=T), to = max((indiv_dat$roads_hii_end)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "landuse_hii_end"){
    s1$landuse_hii_end <- seq(from = min(indiv_dat$landuse_hii_end, na.rm=T), to = max(indiv_dat$landuse_hii_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "landuse_hii_end2"){
    s1$landuse_hii_end2 <- seq(from = min((indiv_dat$landuse_hii_end)^2, na.rm=T), to = max((indiv_dat$landuse_hii_end)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "sl_"){
    s1$sl <- seq(from = min(indiv_dat$sl_, na.rm=T), to = max(indiv_dat$sl_, na.rm=T), length.out = 200)
  }
  if (curr_param == "log_sl"){
    s1$log_sl <- seq(from = min(indiv_dat$log_sl, na.rm=T), to = max(indiv_dat$log_sl, na.rm=T), length.out = 200)
  }
  if (curr_param == "ta_"){
    s1$ta_ <- seq(from = min(indiv_dat$ta_, na.rm=T), to = max(indiv_dat$ta_, na.rm=T), length.out = 200)
  }
  
  indiv_dat_nested <- dat %>% 
    filter(animal_id == indiv)
  
  mod <- indiv_dat_nested$fit2[[1]]
  
  ### Working. variable names have to be the same across all data frames and model
  l_rss_indiv <- amt::log_rss(mod, s1, s2, ci = "se", ci_level = 0.95)
  
  return(l_rss_indiv$df)
}



indivs <- steps_unscaled_nested %>% pull(animal_id)


# #add rss tables to main data frame
steps_unscaled_nested$elev_rss <- map(indivs, l_rss, dat=steps_unscaled_nested, curr_param="elev_end")

steps_unscaled_nested$ndvi_rss <-map(indivs, l_rss, dat=steps_unscaled_nested, curr_param="ndvi_end")

steps_unscaled_nested$forest_rss <- map(indivs, l_rss, dat=steps_unscaled_nested, curr_param="forest_end")

steps_unscaled_nested$dist_water_rss <-  map(indivs, l_rss, dat=steps_unscaled_nested, curr_param="dist_water_end")

steps_unscaled_nested$roads_rss <- map(indivs, l_rss, dat=steps_unscaled_nested, curr_param="roads_hii_end")

steps_unscaled_nested$landuse_rss <- map(indivs, l_rss, dat=steps_unscaled_nested, curr_param="landuse_hii_end")


#### Quadratic Terms--don't actually need to plot these

# steps_unscaled_nested$elev2_rss <- map(indivs, l_rss, dat=steps_unscaled_nested, curr_param="elev_end2")
# 
# steps_unscaled_nested$ndvi2_rss <-map(indivs, l_rss, dat=steps_unscaled_nested, curr_param="ndvi_end2")
# 
# steps_unscaled_nested$forest2_rss <- map(indivs, l_rss, dat=steps_unscaled_nested, curr_param="forest_end2")
# 
# steps_unscaled_nested$dist_water2_rss <-  map(indivs, l_rss, dat=steps_unscaled_nested, curr_param="dist_water_end2")
# 
# steps_unscaled_nested$roads2_rss <- map(indivs, l_rss, dat=steps_unscaled_nested, curr_param="roads_hii_end2")
# 
# steps_unscaled_nested$landuse2_rss <- map(indivs, l_rss, dat=steps_unscaled_nested, curr_param="landuse_hii_end2")



################################ Univarite RSS plots #################################
#Elevation Log-RSS plot all individuals
# steps_unscaled_nested %>% 
#   select(animal_id:dispersal_status, elev_rss) %>% 
#   unnest(cols=c(elev_rss)) %>% 
#   ggplot(., aes(x = elev_end_x1, y = log_rss)) +
#   #geom_ribbon(aes(ymin = lwr, ymax = upr), 
#              # linetype = "dashed", 
#               #color = "black", fill = "gray80", alpha = 0.5) +
#   geom_smooth(aes(pch=animal_id, color=sex),size = 1) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
#   xlab("Elevation (SD)") +
#   ylab("log-RSS vs Mean Elevation") +
#   theme_bw()
# 
# #Roads Log-RSS plot all individuals
# steps_unscaled_nested %>% 
#   select(animal_id:dispersal_status, roads_rss) %>% 
#   unnest(cols=c(roads_rss)) %>% 
#   ggplot(., aes(x = roads_hii_end_x1, y = log_rss)) +
#   #geom_ribbon(aes(ymin = lwr, ymax = upr), 
#   # linetype = "dashed", 
#   #color = "black", fill = "gray80", alpha = 0.5) +
#   geom_smooth(aes(pch=animal_id, color=sex),size = 1) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
#   xlab("Roads_Impact(SD)") +
#   ylab("log-RSS vs Mean NDVI") +
#   theme_bw()
# 
# #NDVI Log-RSS plot all individuals
# steps_unscaled_nested %>% 
#   select(animal_id:dispersal_status, ndvi_rss) %>% 
#   unnest(cols=c(ndvi_rss)) %>% 
#   ggplot(., aes(x = ndvi_end_x1, y = log_rss)) +
#   #geom_ribbon(aes(ymin = lwr, ymax = upr), 
#   # linetype = "dashed", 
#   #color = "black", fill = "gray80", alpha = 0.5) +
#   geom_smooth(aes(pch=animal_id, color=sex),size = 1) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
#   xlab("NDVI(SD)") +
#   ylab("log-RSS vs Mean NDVI") +
#   theme_bw()
# 


################################ RSS Plot all indiv all covariates #################################

#pivot data for faceting
plot_dat <- steps_unscaled_nested %>% 
  select(-c(steps, fit2)) %>% 
  pivot_longer(elev_rss:landuse_rss, names_to= "cov", values_to = "rss_val")

#initialize blank lists for loop
ls = list()
ls2=list()

#for loop to extract just covariate and log rss-values for each covariate and indiviaual
for(i in 1:nrow(plot_dat)){
  ls[[i]] = case_when(
    plot_dat$cov[i] == "elev_rss" ~ plot_dat$rss_val[[i]]$elev_end_x1,
    plot_dat$cov[i] == "ndvi_rss" ~ plot_dat$rss_val[[i]]$ndvi_end_x1,
    plot_dat$cov[i] == "dist_water_rss" ~ plot_dat$rss_val[[i]]$dist_water_end_x1,
    plot_dat$cov[i] == "forest_rss" ~ plot_dat$rss_val[[i]]$forest_end_x1,
    plot_dat$cov[i] == "roads_rss" ~ plot_dat$rss_val[[i]]$roads_hii_end_x1,
    plot_dat$cov[i] == "landuse_rss" ~ plot_dat$rss_val[[i]]$landuse_hii_end_x1,
    plot_dat$cov[i] == "elev2_rss" ~ plot_dat$rss_val[[i]]$elev_end2_x1,
    plot_dat$cov[i] == "ndvi2_rss" ~ plot_dat$rss_val[[i]]$ndvi_end2_x1,
    plot_dat$cov[i] == "dist_water2_rss" ~ plot_dat$rss_val[[i]]$dist_water_end2_x1,
    plot_dat$cov[i] == "forest2_rss" ~ plot_dat$rss_val[[i]]$forest_end2_x1,
    plot_dat$cov[i] == "roads2_rss" ~ plot_dat$rss_val[[i]]$roads_hii_end2_x1,
    plot_dat$cov[i] == "landuse2_rss" ~ plot_dat$rss_val[[i]]$landuse_hii_end2_x1)
  
  ls2[[i]] = plot_dat$rss_val[[i]]$log_rss
}



#add lists from for loop to data frame and facet plot by sex (Melodie's RSS values don't make any sense)
rss_plot_sex <- plot_dat %>% 
  mutate(cov_vals = ls,
         rss_vals = ls2) %>% 
  filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  filter(animal_id!="Sampson") %>%
  filter(animal_id!="Kingsley") %>%
  select(-rss_val) %>% 
  unnest(cols=c(cov_vals, rss_vals)) %>% 
  ggplot(., aes(x = cov_vals, y = rss_vals)) +
  geom_smooth(aes(pch=animal_id, color=sex),linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Covariate Values") +
  ylab("log-RSS vs Mean Covariate Value") +
  theme_bw() +
  facet_wrap(~cov, scales = "free")

plotly::ggplotly(rss_plot_sex)

#add lists from for loop to data frame and facet plot by sex (Melodie's RSS values don't make any sense)
rss_plot_sex_facet <- plot_dat %>% 
  mutate(cov_vals = ls,
         rss_vals = ls2) %>% 
  filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  filter(animal_id!="Sampson") %>%
  filter(animal_id!="Kingsley") %>%
  select(-rss_val) %>% 
  unnest(cols=c(cov_vals, rss_vals)) %>% 
  filter(cov!="landuse_rss") %>% 
  ggplot(., aes(x = cov_vals, y = rss_vals)) +
  geom_smooth(aes(pch=animal_id, color=sex),linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Covariate Values") +
  ylab("log-RSS vs Mean Covariate Value") +
  theme_bw() +
  facet_wrap(vars(cov, sex), scales = "free")

plotly::ggplotly(rss_plot_sex_facet)

#add lists from for loop to data frame and facet plot by dispersal status
rss_plot_disp <- plot_dat %>% 
  mutate(cov_vals = ls,
         rss_vals = ls2) %>% 
  filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  filter(animal_id!="Sampson") %>%
  filter(animal_id!="Kingsley") %>%
  select(-rss_val) %>% 
  unnest(cols=c(cov_vals, rss_vals)) %>% 
  ggplot(., aes(x = cov_vals, y = rss_vals)) +
  geom_smooth(aes(pch=animal_id, color=dispersal_status),linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Covariate Values") +
  ylab("log-RSS vs Mean Covariate Value") +
  theme_bw() +
  facet_wrap(~cov, scales = "free")

plotly::ggplotly(rss_plot_disp)

#add lists from for loop to data frame and facet plot by dispersal status
rss_plot_disp_facet <- plot_dat %>% 
  mutate(cov_vals = ls,
         rss_vals = ls2) %>% 
  filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  filter(animal_id!="Sampson") %>%
  filter(animal_id!="Kingsley") %>%
  filter(animal_id!="Bunny") %>%
  select(-rss_val) %>% 
  unnest(cols=c(cov_vals, rss_vals)) %>% 
  filter(cov!="landuse_rss") %>% 
  ggplot(., aes(x = cov_vals, y = rss_vals)) +
  geom_smooth(aes(pch=animal_id, color=dispersal_status),linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Covariate Values") +
  ylab("log-RSS vs Mean Covariate Value") +
  theme_bw() +
  facet_wrap(vars(cov, dispersal_status), scales = "free")

plotly::ggplotly(rss_plot_disp_facet)

# ggsave(filename= "male_female_rss_5-02-2023.png", plot= rss_plot_sex_facet)
# ggsave(filename= "res_disp_rss_5-02-2023.png", plot= rss_plot_disp_facet)
# ggsave(filename= "disp_rss_5-02-2023.png", plot= rss_plot_disp)
# ggsave(filename= "mf_rss_5-02-2023.png", plot= rss_plot_sex)



################################ Individual Selection Plots #################################

#nest and scale steps and remove locations with missing data
steps_scaled_nested <- steps_raw %>%
  mutate(across(elev_start:landuse_hii_end, scale)) %>%
  mutate(across(elev_start:landuse_hii_end, as.numeric))%>%
  na.omit() %>%
  nest_by(animal_id, sex, dispersal_status) %>%
  rename(steps=data)

#Refit models with scaled coefficients
steps_scaled_nested$fit <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ elev_end + I(elev_end^2) + ndvi_end + I(ndvi_end^2) + dist_water_end + I(dist_water_end^2) + roads_hii_end + I(roads_hii_end^2) +forest_end + I(forest_end^2) + landuse_hii_end  + I(landuse_hii_end^2) + sl_ + log(sl_) + cos(ta_) + strata(step_id_), model = TRUE))

#inspect fitted object
steps_scaled_nested$fit

#inspect model fit for first individual
steps_scaled_nested$fit[[1]]$model


#summary table of mean coefficient estimates by sex
steps_scaled_nested %>% dplyr::mutate(coef = list(map(fit, ~ broom::tidy(fit$model))))%>%
  dplyr::select(animal_id, sex, dispersal_status, coef) %>% 
  unnest(cols = coef) %>%
  unnest(cols = coef) %>% 
  filter(case_when(!str_detect(term, "\\(") ~ TRUE )) %>% 
  filter(!(term %in% c("sl_", "landuse_hii_end"))) %>% 
  mutate(animal_id = factor(animal_id)) %>%
  group_by(term,sex) %>% 
  summarize(mean=mean(estimate), 
            sd= sd(estimate))

#plot of mean coefficient estimates by sex and dispersal status
ds1 <- steps_scaled_nested %>% dplyr::mutate(coef = list(map(fit, ~ broom::tidy(fit$model))))%>%
  dplyr::select(animal_id, sex, dispersal_status, coef) %>% 
  unnest(cols = coef) %>%
  unnest(cols = coef) %>% 
  filter(case_when(!str_detect(term, "\\(") ~ TRUE )) %>% 
  filter(!(term %in% c("sl_", "landuse_hii_end"))) %>% 
  mutate(animal_id = factor(animal_id),
         term=as.factor(term),
         sex=as.factor(sex),
         dispersal_status=as.factor(dispersal_status)) %>%
  group_by(term, dispersal_status, sex) %>% 
  summarize(mean=mean(estimate), 
            sd= sd(estimate),
            ymin = mean - 1.96* sd,
            ymax = mean + 1.96* sd) %>% 
  ggplot(., aes(x = term, y = mean, col = dispersal_status, pch = sex)) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax), position = position_dodge(width = 0.7), size = 0.8) +
  labs(x = "Covariate", y = "Mean Coefficient Estimate") + 
  theme_light() +
  coord_cartesian(ylim=c(-0.7, 0.7))

plotly::ggplotly(ds1)

#ggsave(filename= "mean_coeffs_sex_disp_5-03-2023.png", plot= ds1)


  

#data for plotting
d2 <- steps_scaled_nested %>% dplyr::mutate(coef = list(map(fit, ~ broom::tidy(fit$model))))%>%
  dplyr::select(animal_id, sex, dispersal_status, coef) %>% 
  unnest(cols = coef) %>%
  unnest(cols = coef) %>% 
  filter(case_when(!str_detect(term, "\\(") ~ TRUE )) %>% 
  filter(term != "sl_") %>% 
  mutate(animal_id = factor(animal_id)) %>%
  group_by(term) %>% 
  summarize(
    mean = mean(estimate),
    ymin = mean - 1.96* sd(estimate), #- 1.96
    ymax = mean + 1.96* sd(estimate)) #1.96

d2$x <- 1:nrow(d2)
d2

#plot of relative selection strength by individual and sex
p1 <- steps_scaled_nested %>% 
  mutate(coef = list(map(fit, ~ broom::tidy(fit$model, conf.int = T)))) %>%
  dplyr::select(animal_id, sex, dispersal_status, coef) %>% 
  unnest(cols = coef) %>% 
  unnest(cols = coef) %>% 
  filter(case_when(!str_detect(term, "\\(") ~ TRUE )) %>% 
  filter(term != "sl_") %>% 
  mutate(animal_id = factor(animal_id)) %>%
  ggplot(., aes(x = term, y = estimate, group = animal_id, col = sex, pch = sex)) +
  geom_rect(mapping = aes(xmin = x - .4, xmax = x + .4, ymin = ymin, ymax = ymax), data = d2, inherit.aes = FALSE,fill = "grey90") +
  geom_segment(mapping = aes(x = x - .4, xend = x + .4,y = mean, yend = mean), data = d2, inherit.aes = FALSE, size = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.7), size = 0.8) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Covariate", y = "Relative selection Strength") + theme_light() +
  scale_x_discrete(labels =d2$term) +
  coord_cartesian(ylim=c(-2,2))


plotly::ggplotly(p1)


#plot of relative selection strength by dispersal status and sex
p2 <- steps_scaled_nested %>% 
  mutate(coef = list(map(fit, ~ broom::tidy(fit$model, conf.int = T)))) %>%
  dplyr::select(animal_id, sex, dispersal_status, coef) %>% 
  unnest(cols = coef) %>% 
  unnest(cols = coef) %>% 
  mutate(animal_id = factor(animal_id)) %>%
  filter(case_when(!str_detect(term, "\\(") ~ TRUE )) %>% 
  filter(term != "sl_") %>% 
  filter(animal_id!="Melodie") %>% 
  ggplot(., aes(x = term, y = estimate, group = animal_id, col = dispersal_status, pch = dispersal_status)) +
  geom_rect(mapping = aes(xmin = x - .4, xmax = x + .4, ymin = ymin, ymax = ymax), data = d2, inherit.aes = FALSE,fill = "grey90") +
  geom_segment(mapping = aes(x = x - .4, xend = x + .4,y = mean, yend = mean), data = d2, inherit.aes = FALSE, size = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.7), size = 0.8) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Covariate", y = "Relative selection Strength") + theme_light() +
  scale_x_discrete(labels =d2$term) +
  coord_cartesian(ylim=c(-2,2))


plotly::ggplotly(p2)





########## TEST ZONE ###################

################################ Fit univariate iSSF to each individual (quadratic) #################################

#unscaled step data
steps_unscaled_nested <- steps_raw %>%
  na.omit() %>% 
  nest_by(animal_id, sex, dispersal_status) %>% 
  rename(steps=data) %>%  group_by(sex, dispersal_status)

steps_unscaled_nested$fit_elev <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ elev_end  + I(elev_end^2) + strata(step_id_), model = TRUE))

steps_unscaled_nested$fit_ndvi <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ ndvi_end + I(ndvi_end^2) + strata(step_id_), model = TRUE))

steps_unscaled_nested$fit_dist_water <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ dist_water_end + I(dist_water_end^2) + strata(step_id_), model = TRUE))

steps_unscaled_nested$fit_roads <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ roads_hii_end  + I(roads_hii_end^2) + strata(step_id_), model = TRUE))

steps_unscaled_nested$fit_forest <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ forest_end + I(forest_end^2) + strata(step_id_), model = TRUE))

# steps_unscaled_nested$fit_landuse <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ landuse_hii_end  + strata(step_id_), model = TRUE))


#function to calculate log-rss for each univariate model applied to each individual
l_rss_uni <- function(dat, indiv, curr_param){
  indiv_dat <- steps_unscaled_nested %>% 
    na.omit() %>% 
    filter(animal_id == indiv) %>% 
    unnest(cols=c(steps))
  
  if (curr_param == "elev_end"){
    s1 <- data.frame(
      cov <-seq(from = min(indiv_dat$elev_end, na.rm=T), to = max(indiv_dat$elev_end, na.rm=T), length.out = 200)) %>% 
      rename(elev_end = 1)
    
    
    #data frame with means of all covariates encountered by Al
    s2 <- data.frame(
      elev_end <-mean(indiv_dat$elev_end)) %>% 
      rename(elev_end = 1)
    
    
    indiv_dat_nested <- steps_unscaled_nested %>% 
      na.omit() %>% 
      filter(animal_id == indiv)
    
    mod <- indiv_dat_nested$fit_elev[[1]]
  }
  if (curr_param == "ndvi_end"){
    s1 <- data.frame(
      cov <-seq(from = min(indiv_dat$ndvi_end, na.rm=T), to = max(indiv_dat$ndvi_end, na.rm=T), length.out = 200)) %>% 
      rename(ndvi_end = 1)
    
    
    #data frame with means of all covariates encountered by Al
    s2 <- data.frame(
      elev_end <-mean(indiv_dat$ndvi_end)) %>% 
      rename(ndvi_end = 1)
    
    
    indiv_dat_nested <- steps_unscaled_nested %>% 
      na.omit() %>% 
      filter(animal_id == indiv)
    
    mod <- indiv_dat_nested$fit_ndvi[[1]]
  }
  if (curr_param == "dist_water_end"){
    s1 <- data.frame(
      cov <-seq(from = min(indiv_dat$dist_water_end, na.rm=T), to = max(indiv_dat$dist_water_end, na.rm=T), length.out = 200)) %>% 
      rename(dist_water_end = 1)
    
    
    #data frame with means of all covariates encountered by Al
    s2 <- data.frame(
      elev_end <-mean(indiv_dat$dist_water_end)) %>% 
      rename(dist_water_end = 1)
    
    
    indiv_dat_nested <- steps_unscaled_nested %>% 
      na.omit() %>% 
      filter(animal_id == indiv)
    
    mod <- indiv_dat_nested$fit_dist_water[[1]]
  }
  if (curr_param == "roads_hii_end"){
    s1 <- data.frame(
      cov <-seq(from = min(indiv_dat$roads_hii_end, na.rm=T), to = max(indiv_dat$roads_hii_end, na.rm=T), length.out = 200)) %>% 
      rename(roads_hii_end = 1)
    
    
    #data frame with means of all covariates encountered by Al
    s2 <- data.frame(
      elev_end <-mean(indiv_dat$roads_hii_end)) %>% 
      rename(roads_hii_end = 1)
    
    
    indiv_dat_nested <- steps_unscaled_nested %>% 
      na.omit() %>% 
      filter(animal_id == indiv)
    
    mod <- indiv_dat_nested$fit_roads[[1]]
  }
  if (curr_param == "forest_end"){
    s1 <- data.frame(
      cov <-seq(from = min(indiv_dat$forest_end, na.rm=T), to = max(indiv_dat$forest_end, na.rm=T), length.out = 200)) %>% 
      rename(forest_end = 1)
    
    
    #data frame with means of all covariates encountered by Al
    s2 <- data.frame(
      elev_end <-mean(indiv_dat$forest_end)) %>% 
      rename(forest_end = 1)
    
    
    indiv_dat_nested <- steps_unscaled_nested %>% 
      na.omit() %>% 
      filter(animal_id == indiv)
    
    mod <- indiv_dat_nested$fit_forest[[1]]
  }
  
  ### Working. variable names have to be the same across all data frames and model
  l_rss_indiv <- amt::log_rss(mod, s1, s2, ci = "se", ci_level = 0.95)
  
  return(l_rss_indiv$df)
}

indivs <- steps_unscaled_nested %>% pull(animal_id)

# #add rss tables to main data frame
steps_unscaled_nested$elev_rss_uni <- map(indivs, l_rss_uni, dat=steps_unscaled_nested, curr_param="elev_end")

steps_unscaled_nested$ndvi_rss_uni <-map(indivs, l_rss_uni, dat=steps_unscaled_nested, curr_param="ndvi_end")

steps_unscaled_nested$forest_rss_uni <- map(indivs, l_rss_uni, dat=steps_unscaled_nested, curr_param="forest_end")

steps_unscaled_nested$dist_water_rss_uni <-  map(indivs, l_rss_uni, dat=steps_unscaled_nested, curr_param="dist_water_end")

steps_unscaled_nested$roads_rss_uni <- map(indivs, l_rss_uni, dat=steps_unscaled_nested, curr_param="roads_hii_end")


#pivot data for faceting
plot_dat <- steps_unscaled_nested %>% 
  select(-c(steps, fit_elev:fit_forest)) %>% 
  pivot_longer(elev_rss_uni:roads_rss_uni, names_to= "cov", values_to = "rss_val")

#initialize blank lists for loop
ls = list()
ls2=list()



for(i in 1:nrow(plot_dat)){
  
  if(plot_dat$cov[i] == "elev_rss_uni"){
    ls[[i]] <- plot_dat$rss_val[[i]]$elev_end_x1
    ls2[[i]] = plot_dat$rss_val[[i]]$log_rss
  }
  if(plot_dat$cov[i] == "ndvi_rss_uni"){
    ls[[i]] <- plot_dat$rss_val[[i]]$ndvi_end_x1
    ls2[[i]] = plot_dat$rss_val[[i]]$log_rss
  }
  if(plot_dat$cov[i] == "forest_rss_uni"){
    ls[[i]] <- plot_dat$rss_val[[i]]$forest_end_x1
    ls2[[i]] = plot_dat$rss_val[[i]]$log_rss
  }
  if(plot_dat$cov[i] == "dist_water_rss_uni"){
    ls[[i]] <- plot_dat$rss_val[[i]]$dist_water_end_x1
    ls2[[i]] = plot_dat$rss_val[[i]]$log_rss
  }
  if(plot_dat$cov[i] == "roads_rss_uni"){
    ls[[i]] <- plot_dat$rss_val[[i]]$roads_hii_end_x1
    ls2[[i]] = plot_dat$rss_val[[i]]$log_rss
  }
}


#add lists from for loop to data frame and facet plot by sex (Melodie's RSS values don't make any sense)

plot_dat$cov_vals <- ls
plot_dat$rss_vals <- ls2

rss_elev_uni <- plot_dat %>% 
  # filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  # filter(animal_id!="Sampson") %>%
  # filter(animal_id!="Kingsley") %>%
  select(-rss_val) %>% 
  filter(cov=="elev_rss_uni") %>% 
  filter(sex=="Male") %>% 
  unnest(cols=c(cov_vals, rss_vals)) %>% 
  ggplot(., aes(x = cov_vals, y = rss_vals)) +
  geom_smooth(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Covariate Values") +
  ylab("log-RSS vs Mean Covariate Value") +
  theme_bw() +
  facet_wrap(~animal_id, scales = "free")

plotly::ggplotly(rss_elev_uni)


