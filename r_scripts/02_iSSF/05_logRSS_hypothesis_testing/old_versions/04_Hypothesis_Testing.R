#### iSSF Hypothesis Testing ####

# Author: Read Barbee

# Date:2023-07-07 

# Purpose: Examine log-RSS responses of each individual to each covariate and look for differences between sexes and dispersal groups


################################ Libraries #################################
library(tidyverse)
library(terra)
library(amt)
library(janitor)
library(GGally)

#########################################################################
##
## 0. Helper Functions
##
##########################################################################
run_global <- function (dat){
  #library(survival)  survival::clogit
  mod <- amt::fit_issf(formula = case_ ~ 
                                        #prey 
                                        evi + I(evi^2) +
                                        dist_water + I(dist_water^2) +
                                        #cover
                                        perc_tree_cover + I(perc_tree_cover^2) +
                                        #terrain
                                        tpi + I(tpi^2) +
                                        #humans
                                        popdens_hii + I(popdens_hii^2) +
                                        roads_hii + I(roads_hii^2) +
                                        infra_hii + I(infra_hii^2) +
                                        #steps
                                        sl_ + log_sl_ + cos_ta_ +
                                        #strata
                                        strata(step_id_),
                                      data = dat,
                                      na.action = "na.omit",
                                      model = TRUE)
  return(mod)
}

set_df_s1 <- function (indiv_dat)
{
  #data frame varying elevation from min value to max value encountered by individual, holding all other covariates at the mean
  s1 <- data.frame(
    evi <-seq(from = -2, to =2, length.out = 200),
    
    evi2 <-mean(indiv_dat$evi^2),
    
    dist_water <- mean(indiv_dat$dist_water),
    
    dist_water2 <-mean(indiv_dat$dist_water^2),
    
    perc_tree_cover <- mean(indiv_dat$perc_tree_cover),
    
    perc_tree_cover2 <-mean(indiv_dat$perc_tree_cover^2),
    
    tpi <- mean(indiv_dat$tpi),
    
    tpi2 <-mean(indiv_dat$tpi^2),
    
    popdens_hii <- mean(indiv_dat$popdens_hii),
    
    popdens_hii2 <-mean(indiv_dat$popdens_hii^2),
    
    roads_hii <- mean(indiv_dat$roads_hii),
    
    roads_hii2 <-mean(indiv_dat$roads_hii^2),
    
    infra_hii <-mean(indiv_dat$infra_hii),
    
    infra_hii2 <-mean(indiv_dat$infra_hii^2),
    
    sl_ <- mean(indiv_dat$sl_),
    
    log_sl_ <- mean(indiv_dat$log_sl_),
    
    cos_ta_ <- mean(indiv_dat$cos_ta_)
  ) %>% 
    rename(evi = 1,
           evi2 = 2,
           dist_water = 3,
           dist_water2 = 4,
           perc_tree_cover = 5,
           perc_tree_cover2 = 6,
           tpi = 7,
           tpi2 = 8,
           popdens_hii = 9,
           popdens_hii2 = 10,
           roads_hii = 11,
           roads_hii2 = 12,
           infra_hii = 13,
           infra_hii2 = 14,
           sl_= 15,
           log_sl_ = 16,
           cos_ta_ = 17
    )
  
  s1$evi <- mean(indiv_dat$evi)
  return(s1)
}

set_df_s2 <- function (indiv_dat)
{
  #data frame with means of all covariates encountered by individual
  s2 <-data.frame(
    evi <-mean(indiv_dat$evi),
    
    evi2 <-mean(indiv_dat$evi^2),
    
    dist_water <- mean(indiv_dat$dist_water),
    
    dist_water2 <-mean(indiv_dat$dist_water^2),
    
    perc_tree_cover <- mean(indiv_dat$perc_tree_cover),
    
    perc_tree_cover2 <-mean(indiv_dat$perc_tree_cover^2),
    
    tpi <- mean(indiv_dat$tpi),
    
    tpi2 <-mean(indiv_dat$tpi^2),
    
    popdens_hii <- mean(indiv_dat$popdens_hii),
    
    popdens_hii2 <-mean(indiv_dat$popdens_hii^2),
    
    roads_hii <- mean(indiv_dat$roads_hii),
    
    roads_hii2 <-mean(indiv_dat$roads_hii^2),
    
    infra_hii <-mean(indiv_dat$infra_hii),
    
    infra_hii2 <-mean(indiv_dat$infra_hii^2),
    
    sl_ <- mean(indiv_dat$sl_),
    
    log_sl_ <- mean(indiv_dat$log_sl_),
    
    cos_ta_ <- mean(indiv_dat$cos_ta_)
  ) %>% 
    rename(evi = 1,
           evi2 = 2,
           dist_water = 3,
           dist_water2 = 4,
           perc_tree_cover = 5,
           perc_tree_cover2 = 6,
           tpi = 7,
           tpi2 = 8,
           popdens_hii = 9,
           popdens_hii2 = 10,
           roads_hii = 11,
           roads_hii2 = 12,
           infra_hii = 13,
           infra_hii2 = 14,
           sl_= 15,
           log_sl_ = 16,
           cos_ta_ = 17
    )
  return(s2)
}

classify_results <- function (tracks, curr_param)
{
  classifications <- data.frame(bear_name = NA, uncertain = NA,
                                positive = NA)
  classifications_ <- classifications
  for (i in 1:nrow(tracks)) {
    mod <- tracks$mod[[i]]
    data <- tracks$data[[i]]
    s1 <- set_df_s1()
    if (curr_param == "NDVI") {
      s1$ndvi_e <- seq(from = min(data$ndvi_e, na.rm = T),
                       to = max(data$ndvi_e, na.rm = T), length.out = 200)
    }
    if (curr_param == "Ruggedness") {
      s1$ruggedness_e <- seq(from = min(data$ruggedness_e,
                                        na.rm = T), to = max(data$ruggedness_e, na.rm = T),
                             length.out = 200)
    }
    if (curr_param == "D2forestedge") {
      s1$d2forestedge_e <- seq(from = min(data$d2forestedge_e,
                                          na.rm = T), to = max(data$d2forestedge_e, na.rm = T),
                               length.out = 200)
    }
    if (curr_param == "Densforestedge") {
      s1$densforestedge_e <- seq(from = min(data$densforestedge_e,
                                            na.rm = T), to = max(data$densforestedge_e,
                                                                 na.rm = T), length.out = 200)
    }
    if (curr_param == "Densriparian") {
      s1$densriparian_e <- seq(from = min(data$densriparian_e,
                                          na.rm = T), to = max(data$densriparian_e, na.rm = T),
                               length.out = 200)
    }
    if (curr_param == "Densbuildings") {
      s1$densbuildings_e <- seq(from = min(data$densbuildings_e,
                                           na.rm = T), to = max(data$densbuildings_e, na.rm = T),
                                length.out = 200)
    }
    if (curr_param == "D2core") {
      s1$d2core_e <- seq(from = min(data$d2core_e, na.rm = T),
                         to = max(data$d2core_e, na.rm = T), length.out = 200)
    }
    lr._ci_se <- amt::log_rss(mod, s1, s2, ci = "se", ci_level = 0.95)
    classifications_$bear_name <- tracks$bear_name[[i]]
    count_pos <- lr._ci_se$df %>% summarize(count_pos = sum(log_rss >
                                                              0))
    classifications_$positive <- ifelse(count_pos > 100,
                                        T, F)
    count_overlap <- lr._ci_se$df %>% rowwise %>% dplyr::summarize(count_overlap = sum(lwr <=
                                                                                         0 && upr >= 0)) %>% dplyr::filter(count_overlap ==
                                                                                                                             1)
    count_overlap <- nrow(count_overlap)
    classifications_$uncertain <- ifelse(count_overlap >
                                           100, T, F)
    classifications <- rbind(classifications, classifications_)
  }
  classifications <- classifications[-1, ]
  classifications <- classifications %>% dplyr::group_by(bear_name) %>%
    dplyr::arrange(bear_name)
  classifications$plot.order <- 1:nrow(classifications)
  return(classifications)
}

#function to calculate log_rss object for elevation for each individual incorporating quadratic terms, sl and ta
l_rss <- function(indiv, dat, curr_param){
  indiv_dat <- dat %>% 
    filter(animal_id == indiv) %>% 
    unnest(cols=c(steps)) %>% 
    na.omit()
  
  s1 <- set_df_s1(indiv_dat)
  s2 <- set_df_s2(indiv_dat)
  
  if (curr_param == "evi"){
    s1$evi <- seq(from = min(indiv_dat$evi, na.rm=T), to = max(indiv_dat$evi, na.rm=T), length.out = 200)
  }
  if (curr_param == "evi2"){
    s1$evi2 <- seq(from = min((indiv_dat$evi)^2, na.rm=T), to = max((indiv_dat$evi)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "dist_water"){
    s1$dist_water <- seq(from = min(indiv_dat$dist_water, na.rm=T), to = max(indiv_dat$dist_water, na.rm=T), length.out = 200)
  }
  if (curr_param == "dist_water2"){
    s1$dist_water2 <- seq(from = min((indiv_dat$dist_water)^2, na.rm=T), to = max((indiv_dat$dist_water)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "perc_tree_cover"){
    s1$perc_tree_cover <- seq(from = min(indiv_dat$perc_tree_cover, na.rm=T), to = max(indiv_dat$perc_tree_cover, na.rm=T), length.out = 200)
  }
  if (curr_param == "perc_tree_cover2"){
    s1$perc_tree_cover2 <- seq(from = min((indiv_dat$perc_tree_cover)^2, na.rm=T), to = max((indiv_dat$perc_tree_cover)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "tpi"){
    s1$tpi <- seq(from = min(indiv_dat$tpi, na.rm=T), to = max(indiv_dat$tpi, na.rm=T), length.out = 200)
  }
  if (curr_param == "tpi2"){
    s1$tpi2 <- seq(from = min((indiv_dat$tpi)^2, na.rm=T), to = max((indiv_dat$tpi)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "popdens_hii"){
    s1$popdens_hii <- seq(from = min(indiv_dat$popdens_hii, na.rm=T), to = max(indiv_dat$popdens_hii, na.rm=T), length.out = 200)
  }
  if (curr_param == "popdens_hii2"){
    s1$popdens_hii2 <- seq(from = min((indiv_dat$popdens_hii)^2, na.rm=T), to = max((indiv_dat$popdens_hii)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "roads_hii"){
    s1$roads_hii <- seq(from = min(indiv_dat$roads_hii, na.rm=T), to = max(indiv_dat$roads_hii, na.rm=T), length.out = 200)
  }
  if (curr_param == "roads_hii2"){
    s1$roads_hii2 <- seq(from = min((indiv_dat$roads_hii)^2, na.rm=T), to = max((indiv_dat$roads_hii)^2, na.rm=T), length.out = 200)
  }
  if (curr_param == "infra_hii"){
    s1$infra_hii <- seq(from = min(indiv_dat$infra_hii, na.rm=T), to = max(indiv_dat$infra_hii, na.rm=T), length.out = 200)
  }
  if (curr_param == "infra_hii2"){
    s1$infra_hii2 <- seq(from = min((indiv_dat$infra_hii)^2, na.rm=T), to = max((indiv_dat$infra_hii)^2, na.rm=T), length.out = 200)
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
  
  model <- indiv_dat_nested$global_fit[[1]]
  
  ### Working. variable names have to be the same across all data frames and model
  l_rss_indiv <- amt::log_rss(model, s1, s2, ci = "se", ci_level = 0.95)
  
  return(l_rss_indiv$df)
}


#########################################################################
##
## 1. Import and format step data
##
##########################################################################

steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_imputed_7-12-2023.csv")

#set all negative elevations to 0
steps <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad))


#########################################################################
##
## 2. Reduce land cover and land use categories
##
##########################################################################
steps <- steps %>% mutate(land_cover_usfs_lumped = fct_collapse(land_cover_usfs, trees = c("trees", "tall_trees_shrubs", "gfh_tree_mix", "tree_shrub_mix",  "barren_tree_mix"),
                                                                shrubs = c("tall_shrubs", "shrubs", "gfh_shrub_mix", "barren_shrub_mix"),
                                                                gfh = c("gfh", "barren_gfh_mix"),
                                                                barren = c("barren_impervious", "snow_ice")),
                          land_use_usfs_lumped = fct_collapse(land_use_usfs, agriculture = c("agriculture", "rangeland_pasture")), .after = land_cover_usfs) %>% 
  select(-c(dispersing, disp_date_nsd, disp_qual)) %>% 
  # dummify(select=c("land_cover_usfs_lumped",
  #                  "land_use_usfs_lumped",
  #                  "season",
  #                  "hunting_season",
  #                  "calving_season")) %>% 
  nest_by(animal_id, sex, dispersal_status) %>% 
  rename(steps=data)

#plot_bar(steps %>% select(gpp:calving_season))

#########################################################################
##
## 3. Scale continuous covariates for model comparison and analysis
##
##########################################################################
# 
# steps_scaled <- steps %>% 
#   mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
#   mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric)) %>% 
#   nest_by(animal_id) %>% rename(steps=data)

#plot_histogram(steps_scaled %>% select(sl_, ta_, gpp:power_hii))

#plot_boxplot(steps_scaled %>% select(case_, sl_, ta_, gpp:power_hii), by = "case_")

#########################################################################
##
## 1. Fit global iSSF model to each individual
##
##########################################################################

#9 convergence issues with imputed set for clogit. 19 warnings for fit_issf
global_fits <- steps %>%  
  pull(steps) %>% 
  map(run_global)

steps$global_fit <- global_fits


#########################################################################
##
## 4. Calculate log-RSS and classify individual responses to each covariate
##
##########################################################################

indivs <- steps %>% pull(animal_id)

steps$evi_rss <- map(indivs, l_rss, dat=steps, curr_param="evi")
steps$dist_water_rss <- map(indivs, l_rss, dat=steps, curr_param="dist_water")
steps$perc_tree_cover_rss <- map(indivs, l_rss, dat=steps, curr_param="perc_tree_cover")
steps$tpi_rss <- map(indivs, l_rss, dat=steps, curr_param="tpi")
steps$popdens_hii_rss <- map(indivs, l_rss, dat=steps, curr_param="pop_dens_hii")
steps$roads_hii_rss <- map(indivs, l_rss, dat=steps, curr_param="roads_hii")
steps$infra_hii_rss <- map(indivs, l_rss, dat=steps, curr_param="infra_hii")


#########################################################################
##
##5. Make RSS plots
##
##########################################################################


#pivot data for faceting
rss_plot_dat <- steps %>% 
  pivot_longer(evi_rss:infra_hii_rss, names_to= "cov", values_to = "rss_val")


#initialize blank lists for loop
ls = list()
ls2=list()


#extract covariate ranges and log_rss values for each covariate and individual
for(i in 1:nrow(rss_plot_dat)){
  
  if(rss_plot_dat$cov[i] == "evi_rss"){
    ls[[i]] <- rss_plot_dat$rss_val[[i]]$evi_x1
    ls2[[i]] = rss_plot_dat$rss_val[[i]]$log_rss
  }
  if(rss_plot_dat$cov[i] == "dist_water_rss"){
    ls[[i]] <- rss_plot_dat$rss_val[[i]]$dist_water_x1
    ls2[[i]] = rss_plot_dat$rss_val[[i]]$log_rss
  }
  if(rss_plot_dat$cov[i] == "perc_tree_cover_rss"){
    ls[[i]] <- rss_plot_dat$rss_val[[i]]$perc_tree_cover_x1
    ls2[[i]] = rss_plot_dat$rss_val[[i]]$log_rss
  }
  if(rss_plot_dat$cov[i] == "tpi_rss"){
    ls[[i]] <- rss_plot_dat$rss_val[[i]]$tpi_x1
    ls2[[i]] = rss_plot_dat$rss_val[[i]]$log_rss
  }
  if(rss_plot_dat$cov[i] == "popdens_hii_rss"){
    ls[[i]] <- rss_plot_dat$rss_val[[i]]$popdens_hii_x1
    ls2[[i]] = rss_plot_dat$rss_val[[i]]$log_rss
  }
  if(rss_plot_dat$cov[i] == "roads_hii_rss"){
    ls[[i]] <- rss_plot_dat$rss_val[[i]]$roads_hii_x1
    ls2[[i]] = rss_plot_dat$rss_val[[i]]$log_rss
  }
  if(rss_plot_dat$cov[i] == "infra_hii_rss"){
    ls[[i]] <- rss_plot_dat$rss_val[[i]]$infra_hii_x1
    ls2[[i]] = rss_plot_dat$rss_val[[i]]$log_rss
  }
}


#add lists from for loop to data frame and facet plot by sex (Melodie's RSS values don't make any sense)

rss_plot_dat$cov_vals <- ls
rss_plot_dat$rss_vals <- ls2


#pretty weird looking

evi_rss_uni <- rss_plot_dat %>% 
  # filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  # filter(animal_id!="Sampson") %>%
  # filter(animal_id!="Kingsley") %>%
  select(-rss_val) %>% 
  filter(cov=="evi_rss") %>% 
  filter(sex=="Male") %>% 
  unnest(cols=c(cov_vals, rss_vals)) %>% 
  ggplot(., aes(x = cov_vals, y = rss_vals)) +
  geom_smooth(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Covariate Values") +
  ylab("log-RSS vs Mean Covariate Value") +
  theme_bw() +
  facet_wrap(~animal_id, scales = "free")

plotly::ggplotly(evi_rss_uni)








