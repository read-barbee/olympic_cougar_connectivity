#### OCP iSSF Module_03: Explore Inidivual Variation in Selection ####

# Author: Read Barbee

# Date:2023-04-28 

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


#boxplot of used values for all continuous variables faceted by sex and dispersal status
steps_long %>% 
  filter(case_==TRUE) %>% 
  ggplot() +
  geom_boxplot(aes(x=as.factor(dispersal_status), y=cov_val, fill=sex))+
  facet_wrap(~cov, scales="free") 

#histogram of used values for all continuous variables faceted by sex and dispersal status
steps_long %>% 
  filter(case_==TRUE) %>% 
  ggplot() +
  geom_histogram(aes(x=cov_val, fill=sex))+
  facet_wrap(~cov, scales="free") 

################################ Prepare data to fit models #################################

#nest and scale steps and remove locations with missing data
# steps_unscaled_nested <- steps_raw %>% 
#   mutate(across(elev_start:landuse_hii_end, scale)) %>% 
#   mutate(across(elev_start:landuse_hii_end, as.numeric))%>% 
#   na.omit() %>% 
#   nest_by(animal_id, sex, dispersal_status) %>% 
#   rename(steps=data)

#unscaled model data
steps_unscaled_nested <- steps_raw %>% 
  na.omit() %>% 
  nest_by(animal_id, sex, dispersal_status) %>% 
  rename(steps=data)

################################ Fit iSSF to each individual #################################

#basic model with just covariates. No step length, turn angle or quadratics
steps_unscaled_nested$fit <-  map(steps_unscaled_nested$steps, ~ amt::fit_issf(., case_ ~ elev_end + ndvi_end + dist_water_end + roads_hii_end + forest_end + landuse_hii_end  + strata(step_id_), model = TRUE))


#inspect fitted object
steps_unscaled_nested$fit

#inspect model fit for first individual
steps_unscaled_nested$fit[[1]]$model


################################ TEST: Calculate log_RSS for one covariate #################################

# #function to calculate log_rss object for elevation for each individual. Working.
# elev_rss <- function(dat, indiv){
#   indiv_dat <- dat %>% 
#     na.omit() %>% 
#     filter(animal_id == indiv) %>% 
#     unnest(cols=c(steps))
#   
#   #data frame varying elevation from min value to max value encountered by Al, holding all other covariates at the mean
#   s1 <- data.frame(
#     elev_end <- seq(from = min(indiv_dat$elev_end), to = max(indiv_dat$elev_end), length.out = 200),
#     
#     ndvi_end <- mean(indiv_dat$ndvi_end),
#     
#     dist_water_end <- mean(indiv_dat$dist_water_end),
#     
#     roads_hii_end <- mean(indiv_dat$roads_hii_end),
#     
#     forest_end <- mean(indiv_dat$forest_end),
#     
#     landuse_hii_end <- mean(indiv_dat$landuse_hii_end)
#   ) %>% 
#     rename(elev_end = elev_end....seq.from...min.indiv_dat.elev_end...to...max.indiv_dat.elev_end...,
#            ndvi_end = ndvi_end....mean.indiv_dat.ndvi_end.,
#            dist_water_end= dist_water_end....mean.indiv_dat.dist_water_end.,
#            roads_hii_end = roads_hii_end....mean.indiv_dat.roads_hii_end.,
#            forest_end = forest_end....mean.indiv_dat.forest_end.,
#            landuse_hii_end = landuse_hii_end....mean.indiv_dat.landuse_hii_end.)
#   
#   #data frame with means of all covariates encountered by Al
#   s2 <- data.frame(
#     elev_end <- mean(indiv_dat$elev_end),
#     
#     ndvi_end <- mean(indiv_dat$ndvi_end),
#     
#     dist_water_end <- mean(indiv_dat$dist_water_end),
#     
#     roads_hii_end <- mean(indiv_dat$roads_hii_end),
#     
#     forest_end <- mean(indiv_dat$forest_end),
#     
#     landuse_hii_end <- mean(indiv_dat$landuse_hii_end)
#   ) %>% 
#     rename(elev_end = elev_end....mean.indiv_dat.elev_end.,
#            ndvi_end = ndvi_end....mean.indiv_dat.ndvi_end.,
#            dist_water_end = dist_water_end....mean.indiv_dat.dist_water_end.,
#            roads_hii_end = roads_hii_end....mean.indiv_dat.roads_hii_end.,
#            forest_end = forest_end....mean.indiv_dat.forest_end.,
#            landuse_hii_end = landuse_hii_end....mean.indiv_dat.landuse_hii_end.
#     )
#   
#   
#   ### Working. variable names have to be the same across all data frames and model
#   l_rss_al <- amt::log_rss(dat$fit[[1]], s1, s2, ci = "se", ci_level = 0.95)
#   
#   return(l_rss_al)
# }
# 
# elev_rss(dat=steps_unscaled_nested, indiv = "Zebra")
# 
# 
# # Plot using ggplot2
# ggplot(l_rss_al$df, aes(x = elev_end_x1, y = log_rss)) +
#   geom_line(size = 1) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
#   xlab("Elevation (SD)") +
#   ylab("log-RSS vs Mean Elevation") +
#   theme_bw()
# 
# 
# #plot with 95% large-sample confidence intervals
# ggplot(l_rss_al$df, aes(x = elev_end_x1, y = log_rss)) +
#   geom_ribbon(aes(ymin = lwr, ymax = upr), 
#               linetype = "dashed", 
#               color = "black", fill = "gray80", alpha = 0.5) +
#   geom_line(size = 1) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
#   xlab("Elevation (SD)") +
#   ylab("log-RSS vs Mean Elevation") +
#   theme_bw()

################################ Log RSS for all individuals and covariates #################################

##function to calculate log_rss object for elevation for each individual. This version only returns the df component of the log_rss object. only for model without quadratic terms
l_rss2 <- function(dat, indiv, curr_param){
  indiv_dat <- dat %>% 
    na.omit() %>% 
    filter(animal_id == indiv) %>% 
    unnest(cols=c(steps))
  
  #data frame varying elevation from min value to max value encountered by Al, holding all other covariates at the mean
  s1 <- data.frame(
    elev_end <-seq(from = -2, to =2, length.out = 200),
    
    ndvi_end <- mean(indiv_dat$ndvi_end),
    
    dist_water_end <- mean(indiv_dat$dist_water_end),
    
    roads_hii_end <- mean(indiv_dat$roads_hii_end),
    
    forest_end <- mean(indiv_dat$forest_end),
    
    landuse_hii_end <- mean(indiv_dat$landuse_hii_end)
  ) %>% 
    rename(elev_end = 1,
           ndvi_end = 2,
           dist_water_end = 3,
           roads_hii_end = 4,
           forest_end = 5,
           landuse_hii_end = 6)
  
  s1$elev_end <- mean(indiv_dat$elev_end)
  
  
  
  #data frame with means of all covariates encountered by Al
  s2 <- data.frame(
    elev_end <-mean(indiv_dat$elev_end),
    
    ndvi_end <- mean(indiv_dat$ndvi_end),
    
    dist_water_end <- mean(indiv_dat$dist_water_end),
    
    roads_hii_end <- mean(indiv_dat$roads_hii_end),
    
    forest_end <- mean(indiv_dat$forest_end),
    
    landuse_hii_end <- mean(indiv_dat$landuse_hii_end)
  ) %>% 
    rename(elev_end = 1,
           ndvi_end = 2,
           dist_water_end = 3,
           roads_hii_end = 4,
           forest_end = 5,
           landuse_hii_end = 6)
  
  if (curr_param == "elev_end"){
    s1$elev_end <- seq(from = min(indiv_dat$elev_end, na.rm=T), to = max(indiv_dat$elev_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "ndvi_end"){
    s1$ndvi_end <- seq(from = min(indiv_dat$ndvi_end, na.rm=T), to = max(indiv_dat$ndvi_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "forest_end"){
    s1$forest_end <- seq(from = min(indiv_dat$forest_end, na.rm=T), to = max(indiv_dat$forest_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "dist_water_end"){
    s1$dist_water_end <- seq(from = min(indiv_dat$dist_water_end, na.rm=T), to = max(indiv_dat$dist_water_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "roads_hii_end"){
    s1$roads_hii_end <- seq(from = min(indiv_dat$roads_hii_end, na.rm=T), to = max(indiv_dat$roads_hii_end, na.rm=T), length.out = 200)
  }
  if (curr_param == "landuse_hii_end"){
    s1$landuse_hii_end <- seq(from = min(indiv_dat$landuse_hii_end, na.rm=T), to = max(indiv_dat$landuse_hii_end, na.rm=T), length.out = 200)
  }
  
  indiv_dat_nested <- dat %>% 
    filter(animal_id == indiv)
  
  mod <- indiv_dat_nested$fit[[1]]
  
  ### Working. variable names have to be the same across all data frames and model
  l_rss_indiv <- amt::log_rss(mod, s1, s2, ci = "se", ci_level = 0.95)
  
  return(l_rss_indiv$df)
}


indivs <- steps_unscaled_nested %>% pull(animal_id)

# #add rss tables to main data frame
steps_unscaled_nested$elev_rss <- map(indivs, l_rss2, dat=steps_unscaled_nested, curr_param="elev_end")

steps_unscaled_nested$ndvi_rss <-map(indivs, l_rss2, dat=steps_unscaled_nested, curr_param="ndvi_end")

steps_unscaled_nested$forest_rss <- map(indivs, l_rss2, dat=steps_unscaled_nested, curr_param="forest_end")

steps_unscaled_nested$dist_water_rss <-  map(indivs, l_rss2, dat=steps_unscaled_nested, curr_param="dist_water_end")

steps_unscaled_nested$roads_rss <- map(indivs, l_rss2, dat=steps_unscaled_nested, curr_param="roads_hii_end")

steps_unscaled_nested$landuse_rss <- map(indivs, l_rss2, dat=steps_unscaled_nested, curr_param="landuse_hii_end")


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

#code for fit without quadratics
#pivot data frame to get a row for every individual and every covariate
plot_dat <- steps_unscaled_nested %>%
  select(-c(steps, fit)) %>%
  pivot_longer(elev_rss:landuse_rss, names_to= "cov", values_to = "rss_val")

#initialize blank lists for for loop
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
    plot_dat$cov[i] == "landuse_rss" ~ plot_dat$rss_val[[i]]$landuse_hii_end_x1)

  ls2[[i]] = plot_dat$rss_val[[i]]$log_rss
}


#add lists from for loop to data frame and facet plot by sex 
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
  # xlab("Elevation (SD)") +
  # ylab("log-RSS vs Mean Elevation") +
  theme_bw() +
  facet_wrap(~cov, scales = "free")

plotly::ggplotly(rss_plot_sex)

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
  # xlab("Elevation (SD)") +
  # ylab("log-RSS vs Mean Elevation") +
  theme_bw() +
  facet_wrap(~cov, scales = "free")

plotly::ggplotly(rss_plot_disp)




################################ Individual Selection Plots #################################

#These will probably work better with scaled values

d2 <- steps_unscaled_nested %>% dplyr::mutate(coef = list(map(fit, ~ broom::tidy(fit$model))))%>%
  dplyr::select(animal_id, sex, dispersal_status, coef) %>% 
  unnest(cols = coef) %>%
  unnest(cols = coef) %>% 
  filter(term!="landuse_hii_end") %>% 
  mutate(animal_id = factor(animal_id)) %>% group_by(term) %>%
  summarize(
    mean = mean(estimate),
    ymin = mean - 1.96 * sd(estimate),
    ymax = mean + 1.96 * sd(estimate))

d2$x <- 1:nrow(d2)
d2

#plot of relative selection strength by individual and sex
p1 <- steps_unscaled_nested %>% 
  mutate(coef = list(map(fit, ~ broom::tidy(fit$model, conf.int = T)))) %>%
  dplyr::select(animal_id, sex, dispersal_status, coef) %>% 
  unnest(cols = coef) %>% 
  unnest(cols = coef) %>%
  filter(term!="landuse_hii_end") %>% 
  mutate(animal_id = factor(animal_id)) %>%
  ggplot(., aes(x = term, y = estimate, group = animal_id, col = sex, pch = sex)) +
  geom_rect(mapping = aes(xmin = x - .4, xmax = x + .4, ymin = ymin, ymax = ymax), data = d2, inherit.aes = FALSE,fill = "grey90") +
  geom_segment(mapping = aes(x = x - .4, xend = x + .4,y = mean, yend = mean), data = d2, inherit.aes = FALSE, size = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.7), size = 0.8) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Covariate", y = "Relative selection Strength") + theme_light() +
  scale_x_discrete(labels =d2$term) 
#+coord_cartesian(ylim=c(-2,2))


plotly::ggplotly(p1)


#plot of relative selection strength by dispersal status and sex
p2 <- steps_unscaled_nested %>% 
  mutate(coef = list(map(fit, ~ broom::tidy(fit$model, conf.int = T)))) %>%
  dplyr::select(animal_id, sex, dispersal_status, coef) %>% 
  unnest(cols = coef) %>% 
  unnest(cols = coef) %>% 
  filter(term!="landuse_hii_end") %>% 
  mutate(animal_id = factor(animal_id)) %>%
  ggplot(., aes(x = term, y = estimate, group = animal_id, col = dispersal_status, pch = sex)) +
  geom_rect(mapping = aes(xmin = x - .4, xmax = x + .4, ymin = ymin, ymax = ymax), data = d2, inherit.aes = FALSE,fill = "grey90") +
  geom_segment(mapping = aes(x = x - .4, xend = x + .4,y = mean, yend = mean), data = d2, inherit.aes = FALSE, size = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.7), size = 0.8) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Covariate", y = "Relative selection Strength") + theme_light() +
  scale_x_discrete(labels =d2$term) +
  coord_cartesian(ylim=c(-0.2,0.2))


plotly::ggplotly(p2)









