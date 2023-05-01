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



################################ Prepare data to fit models #################################

#nest and scale steps and remove locations with missing data
steps_scaled_nested <- steps_raw %>% 
  mutate(across(elev_start:landuse_hii_end, scale)) %>% 
  mutate(across(elev_start:landuse_hii_end, as.numeric))%>% 
  na.omit() %>% 
  nest_by(animal_id, sex, dispersal_status) %>% 
  rename(steps=data)

#unscaled model data
steps_unscaled_nested <- steps_raw %>% 
  na.omit() %>% 
  nest_by(animal_id, sex, dispersal_status) %>% 
  rename(steps=data)

################################ Fit iSSF to each individual #################################

#the piping workflow isn't working for some reason
# indiv_issfs_global <- steps_scaled_nested %>% dplyr::mutate(fit = map(steps, ~ amt::fit_issf(., case_ ~ elev_end + ndvi_end + dist_water_end + roads_hii_end + forest_end + landuse_hii_end +  strata(unique_step))))

steps_scaled_nested$fit <-  map(steps_scaled_nested$steps, ~ amt::fit_issf(., case_ ~ elev_end + ndvi_end + dist_water_end + roads_hii_end + forest_end + landuse_hii_end  + strata(step_id_), model = TRUE))

#+ sl_ + log(sl_) --not converging for some reason

#inspect fitted object
steps_scaled_nested$fit

#inspect model fit for first individual
steps_scaled_nested$fit[[1]]$model



################################ Calculate log_RSS for one individual #################################

al <- steps_raw %>% 
  na.omit() %>% 
  filter(animal_id =="Al")

al2 <- steps_scaled_nested %>% 
  na.omit() %>% 
  filter(animal_id =="Al") %>% 
  unnest(cols=c(steps))

#data frame varying elevation from min value to max value encountered by Al, holding all other covariates at the mean
s1 <- data.frame(
elev_end <- seq(from = min(al2$elev_end), to = max(al2$elev_end), length.out = 200),

ndvi_end <- mean(al2$ndvi_end),

dist_water_end <- mean(al2$dist_water_end),

roads_hii_end <- mean(al2$roads_hii_end),

forest_end <- mean(al2$forest_end),

landuse_hii_end <- mean(al2$landuse_hii_end)
) %>% 
  rename(elev_end = elev_end....seq.from...min.al2.elev_end...to...max.al2.elev_end...,
         ndvi_end = ndvi_end....mean.al2.ndvi_end.,
         dist_water_end= dist_water_end....mean.al2.dist_water_end.,
         roads_hii_end = roads_hii_end....mean.al2.roads_hii_end.,
         forest_end = forest_end....mean.al2.forest_end.,
         landuse_hii_end = landuse_hii_end....mean.al2.landuse_hii_end.)

#data frame with means of all covariates encountered by Al
s2 <- data.frame(
  elev_end <- mean(al2$elev_end),
  
  ndvi_end <- mean(al2$ndvi_end),
  
  dist_water_end <- mean(al2$dist_water_end),
  
  roads_hii_end <- mean(al2$roads_hii_end),
  
  forest_end <- mean(al2$forest_end),
  
  landuse_hii_end <- mean(al2$landuse_hii_end)
) %>% 
  rename(elev_end = elev_end....mean.al2.elev_end.,
         ndvi_end = ndvi_end....mean.al2.ndvi_end.,
         dist_water_end = dist_water_end....mean.al2.dist_water_end.,
         roads_hii_end = roads_hii_end....mean.al2.roads_hii_end.,
         forest_end = forest_end....mean.al2.forest_end.,
         landuse_hii_end = landuse_hii_end....mean.al2.landuse_hii_end.
         )


### Working. variable names have to be the same across all data frames and model
l_rss_al <- amt::log_rss(steps_scaled_nested$fit[[1]], s1, s2, ci = "se", ci_level = 0.95)

# Plot using ggplot2
ggplot(l_rss_al$df, aes(x = elev_end_x1, y = log_rss)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Elevation (SD)") +
  ylab("log-RSS vs Mean Elevation") +
  theme_bw()


#plot with 95% large-sample confidence intervals
ggplot(l_rss_al$df, aes(x = elev_end_x1, y = log_rss)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              linetype = "dashed", 
              color = "black", fill = "gray80", alpha = 0.5) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Elevation (SD)") +
  ylab("log-RSS vs Mean Elevation") +
  theme_bw()


d2 <- steps_scaled_nested %>% dplyr::mutate(coef = list(map(fit, ~ broom::tidy(fit$model))))%>%
  dplyr::select(animal_id, sex, dispersal_status, coef) %>% 
  unnest(cols = coef) %>%
  unnest(cols = coef) %>% 
  mutate(animal_id = factor(animal_id)) %>% group_by(term) %>%
  summarize(
    mean = mean(estimate),
    ymin = mean - 1.96 * sd(estimate),
    ymax = mean + 1.96 * sd(estimate))

d2$x <- 1:nrow(d2)
d2

#plot of relative selection strength by individual and sex
p1 <- steps_scaled_nested %>% 
  mutate(coef = list(map(fit, ~ broom::tidy(fit$model, conf.int = T)))) %>%
  dplyr::select(animal_id, sex, dispersal_status, coef) %>% 
  unnest(cols = coef) %>% 
  unnest(cols = coef) %>% 
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
  ggplot(., aes(x = term, y = estimate, group = animal_id, col = dispersal_status, pch = sex)) +
  geom_rect(mapping = aes(xmin = x - .4, xmax = x + .4, ymin = ymin, ymax = ymax), data = d2, inherit.aes = FALSE,fill = "grey90") +
  geom_segment(mapping = aes(x = x - .4, xend = x + .4,y = mean, yend = mean), data = d2, inherit.aes = FALSE, size = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.7), size = 0.8) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Covariate", y = "Relative selection Strength") + theme_light() +
  scale_x_discrete(labels =d2$term) +
  coord_cartesian(ylim=c(-2,2))


plotly::ggplotly(p2)

