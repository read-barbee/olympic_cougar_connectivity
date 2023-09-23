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
    
    dist_water <- 0,
    
    perc_tree_cover <- 0,
    
    tpi <- 0,
    
    popdens_hii <- 0,
    
    roads_hii <- 0,
    
    infra_hii <-0,
    
    sl_ <- 100,
    
    log_sl_ <- log(100),
    
    cos_ta_ <- 1
  ) %>% 
    rename(evi = 1,
           dist_water = 2,
           perc_tree_cover = 3,
           tpi = 4,
           popdens_hii = 5,
           roads_hii = 6,
           infra_hii = 7,
           sl_= 8,
           log_sl_ = 9,
           cos_ta_ = 10
    )
  
  s1$evi <- 0
  return(s1)
}

set_df_s2 <- function (indiv_dat)
{
  #data frame with means of all covariates encountered by individual
  s2 <-data.frame(
    evi <-0,
    
    dist_water <- 0,
    
    perc_tree_cover <- 0,
    
    tpi <- 0,
    
    popdens_hii <- 0,
    
    roads_hii <- 0,
    
    infra_hii <-0,
    
    sl_ <- 100,
    
    log_sl_ <- log(100),
    
    cos_ta_ <- 1
  ) %>% 
    rename(evi = 1,
           dist_water = 2,
           perc_tree_cover = 3,
           tpi = 4,
           popdens_hii = 5,
           roads_hii = 6,
           infra_hii = 7,
           sl_= 8,
           log_sl_ = 9,
           cos_ta_ = 10
    )
  return(s2)
}

classify_results <- function (curr_param, steps)
{
  classifications <- data.frame(animal_id = NA, uncertain = NA,
                                positive = NA)
  classifications_ <- classifications
  for (i in 1:nrow(steps)) {
    mod <- steps$global_fit[[i]]
    data <- steps$steps[[i]]
    
    rss_name <- as.name(paste0(curr_param, "_rss"))
    
    l_rss <- steps[i, rss_name][[1]] %>% as.data.frame()
    classifications_$animal_id <- steps$animal_id[[i]]
    
    count_pos <- l_rss %>% summarize(count_pos = sum(log_rss >0))
    
    classifications_$positive <- ifelse(count_pos > 100,T, F)
    
    count_overlap <- l_rss %>% 
      rowwise %>% 
      dplyr::summarize(count_overlap = sum(lwr <=0 && upr >= 0)) %>% 
      dplyr::filter(count_overlap == 1)
    
    count_overlap <- nrow(count_overlap)
    
    classifications_$uncertain <- ifelse(count_overlap > 100, T, F)
    
    classifications <- rbind(classifications, classifications_)
  }
  classifications <- classifications[-1, ]
  classifications <- classifications %>% dplyr::group_by(animal_id) %>%
    dplyr::arrange(animal_id)
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
  
  cp_name <- as.name(curr_param)
  
  s1[,curr_param] <- seq(from = min(indiv_dat[,cp_name], na.rm=T), to = max(indiv_dat[,cp_name], na.rm=T), length.out = 200)
  
  
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

#without imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv")

#with imputation
#steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_imputed_7-12-2023.csv")

#set all negative elevations to 0
steps <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad)) %>% 
  filter(is.na(dispersing) | dispersing==TRUE)


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
steps$popdens_hii_rss <- map(indivs, l_rss, dat=steps, curr_param="popdens_hii")
steps$roads_hii_rss <- map(indivs, l_rss, dat=steps, curr_param="roads_hii")
steps$infra_hii_rss <- map(indivs, l_rss, dat=steps, curr_param="infra_hii")

#convert steps to datatable for ease of processing
steps2 <- data.table::as.data.table(steps)


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


#plot of individual log-rss curves for each univariate model by sex
rss_uni_sex <- rss_plot_dat %>% 
  # filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  # filter(animal_id!="Sampson") %>%
  # filter(animal_id!="Kingsley") %>%
  # filter(animal_id!="Bunny") %>%
  select(-rss_val) %>% 
  unnest(cols=c(cov_vals, rss_vals)) %>% 
  ggplot(., aes(x = cov_vals, y = rss_vals)) +
  geom_smooth(aes(col=sex, pch=animal_id), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Covariate Values") +
  ylab("log-RSS vs Mean Covariate Value") +
  theme_bw() +
  facet_wrap(~cov, scales = "free")

plotly::ggplotly(rss_uni_sex)

#ggsave(filename= "mf_rss_uni__5-06-2023.png", plot= rss_uni_sex)

#plot of individual log-rss curves for each univariate model by dispersal status
rss_uni_disp <- rss_plot_dat %>% 
  # filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  # filter(animal_id!="Sampson") %>%
  # filter(animal_id!="Kingsley") %>%
  # filter(animal_id!="Bunny") %>%
  select(-rss_val) %>% 
  unnest(cols=c(cov_vals, rss_vals)) %>% 
  ggplot(., aes(x = cov_vals, y = rss_vals)) +
  geom_smooth(aes(col=dispersal_status, pch=animal_id), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  xlab("Covariate Values") +
  ylab("log-RSS vs Mean Covariate Value") +
  theme_bw() +
  facet_wrap(~cov, scales = "free")

plotly::ggplotly(rss_uni_disp)

#ggsave(filename= "disp_rss_uni_5-06-2023.png", plot= rss_uni_disp)


#faceted by dispersal status
rss_disp_facet_uni <- rss_plot_dat %>% 
  # filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  # filter(animal_id!="Sampson") %>%
  # filter(animal_id!="Kingsley") %>%
  # filter(animal_id!="Bunny") %>%
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

plotly::ggplotly(rss_disp_facet_uni)

#ggsave(filename= "disp_rss_facet_imp_7-14-2023.png", plot= rss_disp_facet_uni)

#faceted by sex
rss_sex_facet_uni <- rss_plot_dat %>% 
  #filter(animal_id!="Melodie") %>% #remove individuals with outlying values skewing plots
  #filter(animal_id!="Sampson") %>%
  #filter(animal_id!="Kingsley") %>%
  #filter(animal_id!="Bunny") %>%
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

plotly::ggplotly(rss_sex_facet_uni)

#ggsave(filename= "sex_rss_facet_imp_7-14-2023.png", plot= rss_sex_facet_uni)











#########################################################################
##
## 6. Classify RSS Values
##
##########################################################################

#make list of covaraites
param_list <- c("evi", "dist_water", "perc_tree_cover", "tpi", "popdens_hii", "roads_hii", "infra_hii")

#classify individual log_rss trends as per Sells et al 2022
classified <- map(param_list, classify_results, steps=steps)
names(classified) <- param_list

#extract sex and dispersal infromation from steps
metadat <- steps %>% select(animal_id, sex, dispersal_status)

#join sex and dispersal status to log_rss classification tables
for(i in 1:length(classified)){
  classified[[i]] <- classified[[i]] %>% left_join(metadat, by=join_by(animal_id)) %>% 
    relocate(c(sex, dispersal_status), .after = animal_id)
}

classified

#########################################################################
##
## 7. Summarize RSS trends
##
##########################################################################

rows <- list()
names <- names(classified)

for (i in 1:length(classified)){
  classes <- classified[[i]] %>% mutate(class =case_when(uncertain==TRUE ~ "Uncertain",
                                                         uncertain==FALSE & positive==TRUE ~ "Positive",
                                                         uncertain==FALSE & positive==FALSE ~ "Negative"))
  
  class_perc <- classes %>% group_by(sex, dispersal_status, class) %>% summarize(count =n()) %>% mutate(freq = (count / sum(count))) %>% 
    mutate(perc = round((freq * 100), digits=2))
  
  rows[[i]] <- class_perc %>% 
    mutate(perc_count = paste0(perc, " (",count,")")) %>% 
    unite("class", sex:class, sep="_") %>% 
    select(-c(count:perc)) %>% 
    pivot_wider(names_from = class, values_from = perc_count) %>% 
    mutate(name = names[i]) %>% 
    select(name, everything())
}

summary_table <- bind_rows(rows) %>% mutate(across(everything(), ~replace_na(., "0.0")))
  
#write_csv(summary_table, "feature_selection/rss_summary_no_imp_7-14-23.csv")
################################ GRAVEYARD #################################

#old l_rss function
#l_rss_old <- function(indiv, dat, curr_param){
#   indiv_dat <- dat %>% 
#     filter(animal_id == indiv) %>% 
#     unnest(cols=c(steps)) %>% 
#     na.omit()
#   
#   s1 <- set_df_s1(indiv_dat)
#   s2 <- set_df_s2(indiv_dat)
#   
#   if (curr_param == "evi"){
#     s1$evi <- seq(from = min(indiv_dat$evi, na.rm=T), to = max(indiv_dat$evi, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "evi2"){
#     s1$evi2 <- seq(from = min((indiv_dat$evi)^2, na.rm=T), to = max((indiv_dat$evi)^2, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "dist_water"){
#     s1$dist_water <- seq(from = min(indiv_dat$dist_water, na.rm=T), to = max(indiv_dat$dist_water, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "dist_water2"){
#     s1$dist_water2 <- seq(from = min((indiv_dat$dist_water)^2, na.rm=T), to = max((indiv_dat$dist_water)^2, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "perc_tree_cover"){
#     s1$perc_tree_cover <- seq(from = min(indiv_dat$perc_tree_cover, na.rm=T), to = max(indiv_dat$perc_tree_cover, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "perc_tree_cover2"){
#     s1$perc_tree_cover2 <- seq(from = min((indiv_dat$perc_tree_cover)^2, na.rm=T), to = max((indiv_dat$perc_tree_cover)^2, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "tpi"){
#     s1$tpi <- seq(from = min(indiv_dat$tpi, na.rm=T), to = max(indiv_dat$tpi, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "tpi2"){
#     s1$tpi2 <- seq(from = min((indiv_dat$tpi)^2, na.rm=T), to = max((indiv_dat$tpi)^2, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "popdens_hii"){
#     s1$popdens_hii <- seq(from = min(indiv_dat$popdens_hii, na.rm=T), to = max(indiv_dat$popdens_hii, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "popdens_hii2"){
#     s1$popdens_hii2 <- seq(from = min((indiv_dat$popdens_hii)^2, na.rm=T), to = max((indiv_dat$popdens_hii)^2, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "roads_hii"){
#     s1$roads_hii <- seq(from = min(indiv_dat$roads_hii, na.rm=T), to = max(indiv_dat$roads_hii, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "roads_hii2"){
#     s1$roads_hii2 <- seq(from = min((indiv_dat$roads_hii)^2, na.rm=T), to = max((indiv_dat$roads_hii)^2, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "infra_hii"){
#     s1$infra_hii <- seq(from = min(indiv_dat$infra_hii, na.rm=T), to = max(indiv_dat$infra_hii, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "infra_hii2"){
#     s1$infra_hii2 <- seq(from = min((indiv_dat$infra_hii)^2, na.rm=T), to = max((indiv_dat$infra_hii)^2, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "sl_"){
#     s1$sl <- seq(from = min(indiv_dat$sl_, na.rm=T), to = max(indiv_dat$sl_, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "log_sl"){
#     s1$log_sl <- seq(from = min(indiv_dat$log_sl, na.rm=T), to = max(indiv_dat$log_sl, na.rm=T), length.out = 200)
#   }
#   if (curr_param == "ta_"){
#     s1$ta_ <- seq(from = min(indiv_dat$ta_, na.rm=T), to = max(indiv_dat$ta_, na.rm=T), length.out = 200)
#   }
#   
#   indiv_dat_nested <- dat %>%
#     filter(animal_id == indiv)
#   
#   model <- indiv_dat_nested$global_fit[[1]]
#   
#   ### Working. variable names have to be the same across all data frames and model
#   l_rss_indiv <- amt::log_rss(model, s1, s2, ci = "se", ci_level = 0.95)
#   
#   return(l_rss_indiv$df)
# }
# 
