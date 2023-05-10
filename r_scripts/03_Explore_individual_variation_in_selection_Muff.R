#### OCP iSSF Module_03: Explore Inidivual Variation in Selection--Muff ####

# Author: Read Barbee

# Date:2023-05-08 

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
library(GGally)
library(glmmTMB)
library(broom.mixed)
library(sjPlot)

################################ User-Defined Parameters #################################


############################### Import amt-formatted Data #################################
# Import SSF data for 100 individual mountain lions (November 2022)
steps_raw <- read_csv("6h_steps_unscaled_cov_5-09-2023.csv")



################################ Prepare data to fit models #################################

date_filter <-  interval(start=min(steps_raw$t1_), end=max(steps_raw$t2_))

#individuals with <100 steps or <100 days of steps
# removal_list <- steps_raw %>%
#   group_by(animal_id) %>%
#   distinct(step_id_) %>%
#   count() %>%
#   filter(n<100) %>%
#   pull(animal_id)


#make list of individuals with < 100 steps and <100 days of steps
removal_list <- steps_raw %>%
  group_by(animal_id) %>%
  mutate(date_range=interval(start=min(t1_), end=max(t2_)),
         step_days = as.duration(date_range)/ddays(1)) %>% 
  distinct(step_id_, .keep_all = TRUE) %>% 
  summarize(n=n(), step_days = round(first(step_days), 0)) %>%
  filter(n<100 & step_days <= 100) %>%
  pull(animal_id)


#unscaled step data
steps_unscaled <- steps_raw %>%
  na.omit() %>% 
  filter(!(animal_id %in% removal_list)) %>%
  filter(sl_ >= 100) #%>% remove steps of less than 100m
  # nest_by(animal_id, sex, dispersal_status) %>% 
  # rename(steps=data) %>%  group_by(sex, dispersal_status)

#scaled step data
steps_scaled<- steps_raw %>%
  mutate(across(elev_start:landuse_hii_end, scale)) %>%
  mutate(across(elev_start:landuse_hii_end, as.numeric))%>%
  na.omit() %>%  
  filter(!(animal_id %in% removal_list)) #%>% 
  #filter(sl_ >= 100)


males_scaled <- steps_scaled %>% 
  filter(sex=="Male")

females_scaled <- steps_scaled %>% 
  filter(sex=="Female")

dispersers_scaled <- steps_scaled %>% 
  filter(dispersal_status=="disperser")

residents_scaled <- steps_scaled %>% 
  filter(dispersal_status=="resident")


####################### Poisson Point Process/glmTMB (Muff) #################################

################################ Set model formulas #################################

#Removed step length and log step length due to convergence issue for males and residents

#All individuals
TMB_full <- glmmTMB(case_ ~ elev_end + ndvi_end + dist_water_end + roads_hii_end + forest_end  + (1|animal_id) + (0 + elev_end|animal_id) + (0 + ndvi_end|animal_id) + (0 + dist_water_end|animal_id) + (0 + roads_hii_end|animal_id) + (0 + forest_end|animal_id), family = poisson, data = steps_scaled, doFit = FALSE)

#Males
TMB_males <- glmmTMB(case_ ~ elev_end + ndvi_end + dist_water_end + roads_hii_end + forest_end +  (1|animal_id) + (0 + elev_end|animal_id) + (0 + ndvi_end|animal_id) + (0 + dist_water_end|animal_id) + (0 + roads_hii_end|animal_id) + (0 + forest_end|animal_id), family = poisson, data = males_scaled, doFit = FALSE)

#Females
TMB_females <- glmmTMB(case_ ~ elev_end + ndvi_end + dist_water_end + roads_hii_end + forest_end  + (1|animal_id) + (0 + elev_end|animal_id) + (0 + ndvi_end|animal_id) + (0 + dist_water_end|animal_id) + (0 + roads_hii_end|animal_id) + (0 + forest_end|animal_id), family = poisson, data = females_scaled, doFit = FALSE)

#Residents
TMB_residents <- glmmTMB(case_ ~ elev_end + ndvi_end + dist_water_end + roads_hii_end + forest_end +  (1|animal_id) + (0 + elev_end|animal_id) + (0 + ndvi_end|animal_id) + (0 + dist_water_end|animal_id) + (0 + roads_hii_end|animal_id) + (0 + forest_end|animal_id), family = poisson, data = residents_scaled, doFit = FALSE)

#Dispersers
TMB_dispersers <- glmmTMB(case_ ~ elev_end + ndvi_end + dist_water_end + roads_hii_end + forest_end  + (1|animal_id) + (0 + elev_end|animal_id) + (0 + ndvi_end|animal_id) + (0 + dist_water_end|animal_id) + (0 + roads_hii_end|animal_id) + (0 + forest_end|animal_id), family = poisson, data = dispersers_scaled, doFit = FALSE)


#set variance parameter theta to extremely large value
TMB_full$parameters$theta[1] <-log(1e3)
TMB_males$parameters$theta[1] <-log(1e3)
TMB_females$parameters$theta[1] <-log(1e3)
TMB_residents$parameters$theta[1] <-log(1e3)
TMB_dispersers$parameters$theta[1] <-log(1e3)

#match mapArg to number of theta parameters (5 in this case)
TMB_full$mapArg <-list(theta=factor(c(1:6)))
TMB_males$mapArg <-list(theta=factor(c(1:6)))
TMB_females$mapArg <-list(theta=factor(c(1:6)))
TMB_residents$mapArg <-list(theta=factor(c(1:6)))
TMB_dispersers$mapArg <-list(theta=factor(c(1:6)))

################################ Fit Models #################################
#All individuals
full_fit <- glmmTMB::fitTMB(TMB_full)
summary(full_fit)

#Males
males_fit <- glmmTMB::fitTMB(TMB_males)
summary(males_fit)

#Females
females_fit <- glmmTMB::fitTMB(TMB_females)
summary(females_fit)

#Residents
residents_fit <- glmmTMB::fitTMB(TMB_residents)
summary(residents_fit)

#Dispersers
dispersers_fit <- glmmTMB::fitTMB(TMB_dispersers)
summary(dispersers_fit)

################################ Plots #################################

#intercept and coefficient estimates for all males
sjPlot::plot_model(males_fit, type="re",
                   transform=NULL,
                   axis.title = "Coefficient Estimates (Untransformed)",
                   title = "Individual Coefficient Estimates for Males (Untransformed)")

#intercept and coefficient estimates for all females
sjPlot::plot_model(females_fit, type="re",
                   transform=NULL,
                   axis.title = "Coefficient Estimates (Untransformed)",
                   title = "Individual Coefficient Estimates for Females (Untransformed)")

#intercept and coefficient estimates for all residents
sjPlot::plot_model(residents_fit, type="re",
                   transform=NULL,
                   axis.title = "Coefficient Estimates (Untransformed)",
                   title = "Individual Coefficient Estimates for Residents (Untransformed)")

#intercept and coefficient estimates for all dispersers
sjPlot::plot_model(dispersers_fit, type="re",
                   transform=NULL,
                   axis.title = "Coefficient Estimates (Untransformed)",
                   title = "Individual Coefficient Estimates for Dispersers (Untransformed)")

#coefficient estimates for all individuals
plot_model(full_fit, 
           rm.terms=c("sl_", "log(sl_)"), 
           p.shape=TRUE, 
           transform=NULL, 
           show.values = TRUE, 
           grid=FALSE,
           show.p=FALSE,
           m.labels = c("Males", "Females"), 
           legend.title = "Sex",
           title = "Coefficient Estimates by Sex (Untransformed)",
           axis.title = "Coefficient Estimates (Untransformed)",
           axis.labels = c("elev_end" = "Elevation (m)", 
                           "ndvi_end" ="NDVI", 
                           "dist_water_end" = "Distance to Water (m)", 
                           "roads_hii_end"= "Road Impact", 
                           "forest_end"= "Forest Cover"))


#multi-model plot of unscaled coefficient estimates by sex
plot_models(males_fit, females_fit, 
            rm.terms=c("sl_", "log(sl_)"), 
            p.shape=TRUE, 
            transform=NULL, 
            show.values = TRUE, 
            grid=FALSE,
            show.p=FALSE,
            m.labels = c("Males", "Females"), 
            legend.title = "Sex",
            title = "Coefficient Estimates by Sex (Untransformed)",
            axis.title = "Coefficient Estimates (Untransformed)",
            axis.labels = c("elev_end" = "Elevation (m)", 
                            "ndvi_end" ="NDVI", 
                            "dist_water_end" = "Distance to Water (m)", 
                            "roads_hii_end"= "Road Impact", 
                            "forest_end"= "Forest Cover"))

#multi-model plot of unscaled coefficient estimates by sex (faceted)
plot_models(males_fit, females_fit, 
            rm.terms=c("sl_", "log(sl_)"), 
            p.shape=TRUE, 
            transform=NULL, 
            show.values = TRUE, 
            grid=TRUE,
            show.p=FALSE,
            m.labels = c("Males", "Females"), 
            legend.title = "Sex",
            title = "Coefficient Estimates by Sex (Untransformed)",
            axis.title = "Coefficient Estimates (Untransformed)",
            axis.labels = c("elev_end" = "Elevation (m)", 
                            "ndvi_end" ="NDVI", 
                            "dist_water_end" = "Distance to Water (m)", 
                            "roads_hii_end"= "Road Impact", 
                            "forest_end"= "Forest Cover"))


#multi-model plot of unscaled coefficient estimates by dispersal status
plot_models(residents_fit, dispersers_fit, 
            rm.terms=c("sl_", "log(sl_)"), 
            p.shape=TRUE, 
            transform=NULL, 
            show.values = TRUE, 
            grid=FALSE,
            show.p=FALSE,
            m.labels = c("Residents", "Dispersers"), 
            legend.title = "Dispersal Status",
            title = "Coefficient Estimates by Dispersal Status (Untransformed)",
            axis.title = "Coefficient Estimates (Untransformed)",
            axis.labels = c("elev_end" = "Elevation (m)", 
                            "ndvi_end" ="NDVI", 
                            "dist_water_end" = "Distance to Water (m)", 
                            "roads_hii_end"= "Road Impact", 
                            "forest_end"= "Forest Cover"))




#try out model table

tab_model(males_fit, females_fit, residents_fit, dispersers_fit,
          transform=NULL,
          dv.labels = c("Males", "Females", "Residents", "Dispersers"),
          string.est = "Estimate (95% CI)",
          collapse.ci=TRUE,
          show.ngroups = TRUE,
          pred.labels = c("elev_end" = "Elevation (m)", 
                          "ndvi_end" ="NDVI", 
                          "dist_water_end" = "Distance to Water (m)", 
                          "roads_hii_end"= "Road Impact", 
                          "forest_end"= "Forest Cover"),
          file="muff_demographic_table1_full_screening_5-9-2023.png")




