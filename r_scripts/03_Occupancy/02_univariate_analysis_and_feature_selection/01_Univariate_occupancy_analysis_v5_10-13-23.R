#### Univariate occupancy analysis ####

# Author: Read Barbee

# Date:2023-10-12 

#Last updated: 2023-10-14

# Purpose:



library(tidyverse)
library(sf)
library(terra)
library(ubms)
library(beepr)
#library(doParallel)

################################ User-defined Parameters #################################

#  Sampling occasion specs
survey_period <- "day" 
number_of_periods <- 14 # how many days you want each survey to be

#if there are 362 distinct julian days represented in data, why are there only 331 in w_obs?
#there are the same number of unique intervals generated as unique dates 
#but any given cell_year doesn't have more than 331 days in its matrix

##test <- dat %>% group_by(cell_id, year) %>% distinct(date) %>% count()
#range(test$n)

#########################################################################
##
##  1. Import stacked occupancy data
##
##########################################################################
#  Source functions to make occupancy surveys of different lenghts
source("r_scripts/03_Occupancy/01_data_prep/Utility/sampling_interval_aggregation/original/make_int_fun.R")

#example data
#load("r_scripts/03_Occupancy/01_data_prep/Utility/sampling_interval_aggregation/original/dat_example.Rdata") # object called dat

make_interval2 <- function(x, dt_col, time_int = "month", increment = 1) {
  stopifnot(tibble::is_tibble(x))
  tmp <- dplyr::pull(x, {{dt_col}})
  if (time_int == "hour") {
    stopifnot(inherits(tmp, c("POSIXct", "POSIXlt")))
  }
  stopifnot(all(!is.na(tmp)))
  stopifnot(time_int %in% c("hour", "day", "week", "month"))
  if (dplyr::is_grouped_df(x)) {
    args_passed <- as.list(match.call())
    stop(
      "Data is grouped, try something like",
      "\n\t",
      paste0(
        "x %>% ",
        "dplyr::do(gbn_make_interval(",
        args_passed$dt_col, ", ",
        "'", args_passed$time_int, "', ",
        args_passed$increment, "))"
      )
    )
  }
  
  int <- lubridate::interval(min(tmp), tmp)
  
  new <- x %>%
    dplyr::mutate(
      rnd_dt = lubridate::round_date(tmp, unit = paste(increment, time_int)),
      interval = switch(
        time_int,
        "month" = int %/% lubridate::period(increment, units = "month"),
        "week" = int %/% lubridate::weeks(increment),
        "day" = int %/% lubridate::days(increment),
        "hour" = int %/% lubridate::hours(increment)
      ) + 1
    ) 
  return(new)
}

#activity/detection data
dat <- read_csv("data/Camera_Data/master/ocp_onp_occ_dat_long_10-13-23.csv")

# Load raster data
cov_stack <- rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2_asp.tif")


#########################################################################
##
##  2. Format data for new survey interval
##
##########################################################################

# Note date column should be date format:
# $ date   : Date, format: "2009-05-02" "2009-08-10" "2009-05-14" ...

#check for missing days
missing_days <- which(!(1:365 %in% dat$j_day ))

missing_dates <- list()
for(i in 1:length(missing_days)){
missing_dates[[i]] <- ymd("2020-01-01") + days(missing_days[i])
}

#missing 4 days march 23-26

# dat %>% group_by(grid_id, year) %>% 
#   summarize(start = min(date),
#             end = max(date)) %>% View()


#Assign each row to an observation interval
obs_tmp1 <- make_interval2(dat, date, survey_period, number_of_periods) 


#test_dat <- obs_tmp1 %>% filter(cell_id == 49402630 & year == 2019)

#add column for survey interval
obs_tmp2 <- obs_tmp1 %>% 
	dplyr::group_by(cell_id, year) %>% 
	dplyr::mutate(survey_interval = interval - min(interval) + 1) %>%
	dplyr::ungroup()

#test_dat <- obs_tmp2 %>% filter(cell_id == 49402630 & year == 2019)

#calculate the maximum number of observations for each survey interval
obs_tmp3 <- obs_tmp2 %>%
	group_by(survey_interval, cell_id, year) %>% #station_id, 
  summarize(station_id = first(station_id),
            grid_id = first(grid_id),
            lon = first(lon),
            lat = first(lat),
            utm_e = first(utm_e),
            utm_n = first(utm_n),
            # annual_effort_correction_factor = first(effort_correction_factor),
            # annual_effort_correction_factor_bait = first(effort_correction_factor_bait),
            # annual_effort_correction_factor_snare = first(effort_correction_factor_snare),
            cam_days = sum(cam_status, na.rm = TRUE),
            bait_days = sum(bait_status, na.rm = TRUE),
            snare_days = sum(snare_status, na.rm = TRUE),
            cougar_detections_binary = case_when(all(is.na(cougar_detections_binary))==T ~NA,
                                                 sum(cougar_detections_binary, na.rm = TRUE) == 0 ~ 0,
                                                 sum(cougar_detections_binary, na.rm = TRUE) > 0 ~ 1))



w_obs <- obs_tmp3 %>%
	arrange(survey_interval, cell_id, year) %>%
	pivot_wider(values_from = c(cam_days, bait_days, snare_days, cougar_detections_binary), names_from = survey_interval) %>% 
 mutate(across(!contains("detections") , \(x) replace_na(x, 0))) %>% 
	as.data.frame()
	
#w_obs[w_obs == "-Inf"] <- NA

#id_cols = c(cell_id, year),


#########################################################################
##
##  2. Make objects for ubms
##
##########################################################################


sf_pts <- w_obs %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(crs = 5070)

#mapview::mapview(sf_pts)

#Extract cell numbers for each station location
occ_dat <-terra::extract(cov_stack, sf_pts, bind = TRUE) %>% 
  as.data.frame()

occ_dat_scaled <- occ_dat %>% 
  select(-land_cover_usfs, -land_use_usfs) %>% 
  mutate(across(tree_cover_hansen:easting, scale)) %>% 
  mutate(across(tree_cover_hansen:easting, as.numeric))

#remove station rows with missing covariate values
complete_cases <- occ_dat_scaled %>% select(tree_cover_hansen:easting) %>% 
  complete.cases()

occ_dat_complete <- occ_dat_scaled %>% 
  mutate(comp = complete_cases) %>% 
  unite("cell_id_year", cell_id, year, remove = FALSE) %>% 
  mutate(cell_id_year = as.factor(cell_id_year)) %>% 
  filter(comp==TRUE)
  

#construct umf stack

nsite <- nrow(occ_dat_complete)
y <- occ_dat_complete %>% 
	dplyr::select(contains("detections")) %>% 
	as.matrix()

# Number of surveys detected per site
summary(rowSums(y, na.rm = TRUE))	

eff <- occ_dat_complete %>%
  dplyr::select(contains("cam")) %>%
	as.matrix()

bait <- occ_dat_complete %>%
  dplyr::select(contains("bait")) %>%
  as.matrix()

snare <- occ_dat_complete %>%
  dplyr::select(contains("snare")) %>%
  as.matrix()

covs_scaled <- occ_dat_complete %>% select(cell_id_year, tree_cover_hansen:easting)

umf_stack <- unmarkedFrameOccu(y = y, 
	siteCovs = covs_scaled,
	obsCovs = list(eff = eff,
	               bait = bait,
	               snare = snare)) #, survey = surveyID))
head(umf_stack)


#########################################################################
##
##  2. Inspect correlation of covariates
##
##########################################################################
corr_dat <- covs %>% 
	dplyr::select(-ID)
cor_matrix <- cor(corr_dat, use= "complete.obs")

#corrplot::corrplot(corr_dat, method = "number", type = "upper")
#ggsave(file = "covariate_correlation_plot.png")

threshold <- 0.5

high_cor_pairs <- which(abs(cor_matrix) >= threshold, arr.ind = TRUE)
high_cor_pairs <- high_cor_pairs[high_cor_pairs[, 1] != high_cor_pairs[, 2], ]

var1 <- vector()
var2 <- vector()
correlation <- vector()

for (i in 1:nrow(high_cor_pairs)) {
  var1[i] <- rownames(cor_matrix)[high_cor_pairs[i, 1]]
  var2[i] <- colnames(cor_matrix)[high_cor_pairs[i, 2]]
  correlation[i] <- cor_matrix[high_cor_pairs[i, 1], high_cor_pairs[i, 2]]
  #print(paste("Variables:", var1, "and", var2, "- Correlation:", correlation))
}

#ordered list of correlated covariates
cor_dat <- tibble(variables = paste0(var1, "_", var2), corr = correlation) %>% 
  #filter(corr < 0.977) %>% 
  arrange(desc(corr))

#the GGAlly package corrplot is more informative for numeric variables
GGally::ggcorr(corr_dat, label = TRUE)

#write_csv(cor_dat, "feature_selection/pairwise_cov_corr_occ_12month_10-14-23.csv")

#########################################################################
##
##  2. Determine best model for detection
##
##########################################################################

#null model
null <- stan_occu(~1 ~1, 
   data = umf_stack, chains = 3, iter = 1500, cores = 3)
null


cam_days <- stan_occu(~scale(eff)  ~1, 
   data = umf_stack, chains = 3, iter = 1500, cores = 3)
cam_days

bait_days <- stan_occu(~scale(bait) ~1, 
   data = umf_stack, chains = 3, iter = 1500, cores = 3)
bait_days

snare_days <- stan_occu(~scale(snare) ~1, 
   data = umf_stack, chains = 3, iter = 1500, cores = 3)
snare_days

cam_bait_snare <- stan_occu(~scale(eff) + scale(bait) + scale(snare) ~1, 
                   data = umf_stack, chains = 3, iter = 1500, cores = 3)
cam_bait_snare

cam_bait <- stan_occu(~scale(eff) + scale(bait) ~1, 
                   data = umf_stack, chains = 3, iter = 1500, cores = 3)
cam_bait

cam_snare <- stan_occu(~scale(eff) + scale(snare) ~1, 
                   data = umf_stack, chains = 3, iter = 1500, cores = 3)
cam_snare

det_fits <- c(null, cam_days, bait_days, snare_days, cam_bait_snare, cam_bait, cam_snare)

names(det_fits) <- c("null", "cam_days", "bait_days", "snare_days", "cam_bait_snare", "cam_bait", "cam_snare")

fit_list <- fitList(det_fits)

mod_sel <- modSel(fit_list, nullmod="null") %>% 
  rownames_to_column(var = "model")



estimates <- list()
waic <- list()
looic <- list ()
for(i in 1:length(det_fits)){
  estimates[[i]] <- summary(det_fits[[i]], "det") %>% 
    rownames_to_column(var = "term") %>% 
    filter(str_detect(term, "Intercept")== FALSE) %>% 
    mutate(model = names(det_fits)[i], .before = term)
  
  waic[[i]] <- waic(det_fits[[i]])$estimates %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "metric") %>% 
    pivot_wider(names_from = "metric", values_from = c("Estimate", "SE")) %>% 
    mutate(model = names(det_fits)[i], .before = Estimate_elpd_waic)
  
  looic[[i]] <- tibble(
    model = names(det_fits)[i],
    looic = det_fits[[i]]@loo$looic 
  )
}

estimates <- bind_rows(estimates)
waic <- bind_rows(waic)
looic <- bind_rows(looic)


occu_det_summ <- mod_sel %>% 
  left_join(estimates, by=join_by(model)) %>% 
  left_join(waic, by = join_by(model)) %>% 
  left_join(looic, by = join_by(model)) %>%
  relocate(elpd:weight, .after=Rhat)

#write_csv(occu_det_summ, "/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/feature_selection/Occupancy/occu_det_model_fits_12month_10-14-23.csv")


#  Significant effects of: ___ so take this detection
#   model forward
# save(fitd3, file = "Single_covariate_models/fitd3.Rdata")


#########################################################################
##
##  2. Run univariate models with best detection model
##
##########################################################################


#### LINEAR MODELS ###
#register doParallel backend. Note:the cores argument implies forking with doMC backend, but you can check with the getDoParName() function
#doParallel::registerDoParallel(cores = 9) 

#### NULL MODEL ###

null_fit <- stan_occu(~1 ~ (1|cell_id_year), data=umf_stack, chains=3, iter=1500)

traceplot(null_fit)

traceplot(null_fit, pars = "sigma_state")

null_fit_list <- list(null_fit)

names(null_fit_list) <- "null"



cov_names <- covs_scaled %>% select(-cell_id_year) %>% names()

#parallel loop doesn't seem to work

#takes 2 hours without parallel loop. Try parallel next time.
uni_fits <- list()
system.time(for(i in 1:length(cov_names)){
  cov <- cov_names[i]
  
  #construct the model formula with linear terms only
  form <- as.formula(paste0("~ scale(eff) + scale(bait) + scale(snare) ~ ", cov, " + (1|cell_id_year)"))
  
  
  uni_fits[[i]] <- stan_occu(form, data = umf_stack, chains = 3, iter = 1500, cores = 3)
  
  print(paste0(i, "/", length(cov_names)))
  
}) ; beep("fanfare")


names(uni_fits) <- cov_names

uni_fits <- c(uni_fits, null_fit_list)


fit_list <- fitList(uni_fits)

mod_sel <- modSel(fit_list,  nullmod="null") %>%   #, nullmod="null"
  rownames_to_column(var = "model")



estimates <- list()
waic <- list()
looic <- list ()
for(i in 1:length(uni_fits)){
  estimates[[i]] <- summary(uni_fits[[i]], "state") %>% 
    rownames_to_column(var = "term") %>% 
    filter(str_detect(term, "Intercept")== FALSE) %>% 
    mutate(model = names(uni_fits)[i], .before = term)
  
  waic[[i]] <- waic(uni_fits[[i]])$estimates %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "metric") %>% 
    pivot_wider(names_from = "metric", values_from = c("Estimate", "SE")) %>% 
    mutate(model = names(uni_fits)[i], .before = Estimate_elpd_waic)
  
  looic[[i]] <- tibble(
    model = names(uni_fits)[i],
    looic = uni_fits[[i]]@loo$looic 
  )
}

estimates <- bind_rows(estimates)
waic <- bind_rows(waic)
looic <- bind_rows(looic)


occu_uni_summ <- mod_sel %>% 
  left_join(estimates, by=join_by(model)) %>% 
  left_join(waic, by = join_by(model)) %>% 
  left_join(looic, by = join_by(model)) %>%
  relocate(elpd:weight, .after=Rhat) %>% 
  filter(!str_detect(term, "sigma") | model=="null")

#write_csv(occu_uni_summ, "feature_selection/Occupancy/occu_uni_fits_14day_cam_bait_snare_10-15-23.csv")
#RESUME HERE. SAVE MODEL STATS, THEN DO THE SAME WITH QUADRATICS
#CONSIDER DETECTION COVARIATE FOR PARK VS OCP INSTEAD OF BAIT AND SNARE


#########################################################################
##
##  4. Fit quadratic univariate stacked occupancy models
##
##########################################################################

#parallelized for loop: fits in ~ 10 minutes with 3 chains and 100 iterations; 35 min for 1000 iterations
uni_fits_quad <- list()
system.time(uni_fits_quad <- for (i in 1:length(cov_names)){
  
  cov <- cov_names[i]
  
  form <- as.formula(paste0("~ scale(eff) + scale(bait) + scale(snare) ~", #remove standard intercept to be replace with stratum-based intercept
                            #fixed effects
                            cov, "+", "I(", cov,  "^2) + ",
                            #random intercept (strata)
                            "(1|cell_id_year)"))
  
  
  uni_fits_quad[[i]] <- stan_occu(form, data=umf_stack, chains=3, iter=1500)
  
  print(paste0(i, "/", length(cov_names)))
  
}); beep("fanfare")

#rename qudaratic covariates with "2" at the end
quad_names <- vector()
for(i in 1:length(cov_names)){
  old_name <- cov_names[i]
  quad_names[i] <- paste0(old_name, "2")
}

names(uni_fits_quad) <- quad_names



#########################################################################
##
##  4. Fit top model from SSA
##
##########################################################################

top_ssa <- stan_occu(~ scale(eff) + scale(bait) + scale(snare) ~ tree_cover_hansen + ndvi + northing + tri + perc_tree_cover + perc_nonveg + roads_hii + elevation + popdens_hii + landuse_hii + I(elevation^2) + I(tree_cover_hansen^2) + I(ndvi^2) + I(perc_tree_cover^2) + I(roads_hii^2) + (1|cell_id_year), data=umf_stack, chains=3, iter=1500) ; beep("fanfare")


#extract fixed effect beta values from fitted_model
betas <- inla_fit_quad3$summary.fixed$mean

#reff <- ranef(global_fit)




#manual predictions
preds <- exp((betas[1] * cov_stack_pred$tree_cover_hansen) + 
               (betas[2] * cov_stack_pred$ndvi) + 
               #(betas[3] * cov_stack_pred$tpi) + 
               (betas[3] * cov_stack_pred$northing) + 
               #(betas[4] * cov_stack_pred$tri) + 
               (betas[4] * cov_stack_pred$perc_tree_cover) +
               (betas[5] * cov_stack_pred$perc_nonveg) +
               (betas[6] * cov_stack_pred$roads_hii) + 
               (betas[7] * cov_stack_pred$elevation) +
               (betas[8] * cov_stack_pred$popdens_hii) +
               (betas[9] * cov_stack_pred$landuse_hii) + 
               #(betas[11] * cov_stack_pred$rails_hii) +
               #(betas[12] * cov_stack_pred$infra_hii))/
               #(betas[13] * cov_stack_pred$easting) +
               (betas[10] * cov_stack_pred$elevation^2) +
               (betas[11] * cov_stack_pred$tree_cover_hansen^2) +
               (betas[12] * cov_stack_pred$ndvi^2) + 
               #(betas[16] * cov_stack_pred$tpi^2) +
               #(betas[17] * cov_stack_pred$northing^2) +
               #(betas[14] * cov_stack_pred$tri^2) +
               (betas[13] * cov_stack_pred$perc_tree_cover^2) +
               (betas[14] * cov_stack_pred$roads_hii^2))/
              # (betas[17] * cov_stack_pred$popdens_hii^2))/
              #(betas[22] * cov_stack_pred$landuse_hii^2))/
  (1 + exp((betas[1] * cov_stack_pred$tree_cover_hansen) + 
             (betas[2] * cov_stack_pred$ndvi) + 
             #(betas[3] * cov_stack_pred$tpi) + 
             (betas[3] * cov_stack_pred$northing) + 
             #(betas[4] * cov_stack_pred$tri) + 
             (betas[4] * cov_stack_pred$perc_tree_cover) +
             (betas[5] * cov_stack_pred$perc_nonveg) +
             (betas[6] * cov_stack_pred$roads_hii) + 
             (betas[7] * cov_stack_pred$elevation) +
             (betas[8] * cov_stack_pred$popdens_hii) +
             (betas[9] * cov_stack_pred$landuse_hii) + 
             #(betas[11] * cov_stack_pred$rails_hii) +
             #(betas[12] * cov_stack_pred$infra_hii))/
             #(betas[13] * cov_stack_pred$easting) +
             (betas[10] * cov_stack_pred$elevation^2) +
             (betas[11] * cov_stack_pred$tree_cover_hansen^2) +
             (betas[12] * cov_stack_pred$ndvi^2) + 
             #(betas[16] * cov_stack_pred$tpi^2) +
             #(betas[17] * cov_stack_pred$northing^2) +
             #(betas[14] * cov_stack_pred$tri^2) +
             (betas[13] * cov_stack_pred$perc_tree_cover^2) +
             (betas[14] * cov_stack_pred$roads_hii^2)))
     # (betas[17] * cov_stack_pred$popdens_hii^2))/
     #(betas[22] * cov_stack_pred$landuse_hii^2))/
   
terra::plot(preds)

#######################################################################
##
## 3. Bin predictions
##
##########################################################################

pred_vals <- terra::values(preds)

breaks <- quantile(pred_vals, probs = 0:10/10, na.rm = T)

bins <- cut(pred_vals, breaks, include.lowest=TRUE)

binned <- terra::classify(preds, rcl=as.vector(breaks))

terra::plot(binned)

levels(binned) <- 1:10



########################################################################
##
## 3. Internal validation (residents)
##
##########################################################################
locs <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_10-02-2023.csv")

residents <- locs %>% 
  filter(dispersal_status =="resident") %>% 
  select(-c(disp_date_nsd:dispersing))

eval_res <- residents %>% st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070)

#extract binned predicted probability of use values for each used location
loc_preds<- terra::extract(binned, eval_res)[,2] %>% as.numeric() %>% na.omit() 

#plot proportion of test points in each bin. 
ggplot() +
  geom_bar(aes(x=loc_preds, y = after_stat(prop))) 

n_bins <- length(unique(loc_preds))

#calculate pearson correlation between the bin number and the proportion of used locations in each bin
cor_cv_res <- cor.test(x = 1:n_bins, y = as.numeric(table(loc_preds))/(length(loc_preds)/n_bins), method = "spearman", exact = FALSE)$estimate

cor_cv_res 


########################################################################
##
## 3. External validation (dispersers)
##
##########################################################################

dispersers <- locs %>% 
  filter(dispersal_status =="disperser") %>% 
  filter(dispersing == TRUE)

eval_disp <- dispersers %>% st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070)

#extract binned predicted probability of use values for each used location
loc_preds_disp<- terra::extract(binned, eval_disp)[,2] %>% as.numeric() %>% na.omit() 

#plot proportion of test points in each bin. 
ggplot() +
  geom_bar(aes(x=loc_preds_disp, y = after_stat(prop))) 

n_bins <- length(unique(loc_preds_disp))

#calculate pearson correlation between the bin number and the proportion of used locations in each bin
cor_cv_res_disp <- cor.test(x = 1:n_bins, y = as.numeric(table(loc_preds_disp))/(length(loc_preds_disp)/n_bins), method = "spearman", exact = FALSE)$estimate

cor_cv_res_disp

########################################################################
##
## 3. Mapping
##
##########################################################################

#viridis magma
ggplot()+
  tidyterra::geom_spatraster(data=binned, mapping=aes()) +
  #geom_sf(data=red_deer_used, aes(geometry=geometry))+
  #coord_sf(datum = st_crs(5070))+
  scale_fill_manual(values = viridis::magma(10), na.value = NA, guide = guide_legend(reverse = TRUE), na.translate=FALSE)  #rev() to reverse pallete

mapview::mapview(raster::raster(binned))


########################################################################
# Get summary stats for each single covariate model
s <- rstan::summary(fit1@stanfit, probs = c(0.05, 0.95))
s$summary[1:2,]
s <- rstan::summary(fit2@stanfit, probs = c(0.05, 0.95))
s$summary[1:2,]
s <- rstan::summary(fit3@stanfit, probs = c(0.05, 0.95))
s$summary[1:2,]
s <- rstan::summary(fit4@stanfit, probs = c(0.05, 0.95))
s$summary[1:2,]
s <- rstan::summary(fit5@stanfit, probs = c(0.05, 0.95))
s$summary[1:2,]



########################################################################
#  Test quadratics for certain covariates 

#  Fit each covaraite that is reasonable to have a quadratic relationship
sq_fit1 <- stan_occu(~ scale(eff) + log(n_station) 
		~1 + scale(asp) + I(scale(asp)^2) + (1|id), 
   data = umf_stack, chains = 3, iter = 15000, cores = 3)
sq_fit2 <- stan_occu(~ scale(eff) + log(n_station) 
		~1 + scale(elv) + I(scale(elv)^2) + (1|id), 
   data = umf_stack, chains = 3, iter = 15000, cores = 3)
sq_fit3 <- stan_occu(~ scale(eff) + log(n_station) 
		~1 + scale(slp) + I(scale(slp)^2) + (1|id), 
   data = umf_stack, chains = 3, iter = 15000, cores = 3)
sq_fit4 <- stan_occu(~ scale(eff) + log(n_station) 
		~1 + scale(bc_prc) + I(scale(bc_prc)^2) + (1|id), 
   data = umf_stack, chains = 3, iter = 15000, cores = 3)
sq_fit5 <- stan_occu(~ scale(eff) + log(n_station) 
		~1 +scale(bc_tmp_mu) + I(scale(bc_tmp_mu)^2) + (1|id), 
   data = umf_stack, chains = 3, iter = 15000, cores = 3)
   
# save(sq_fit1, file = "sq_fit1.Rdata")
# save(sq_fit2, file = "sq_fit2.Rdata")
# save(sq_fit3, file = "sq_fit3.Rdata")
# save(sq_fit3b, file = "sq_fit3b.Rdata")
# save(sq_fit4, file = "sq_fit4.Rdata")
# save(sq_fit5, file = "sq_fit5.Rdata")
# save(sq_fit10, file = "sq_fit10.Rdata")
# save(sq_fit17, file = "sq_fit17.Rdata")
# save(sq_fit18, file = "sq_fit18.Rdata")
# save(sq_fit19, file = "sq_fit19.Rdata")


s <- rstan::summary(sq_fit1@stanfit, probs = c(0.05, 0.95))
s$summary[1:2,]
s <- rstan::summary(sq_fit2@stanfit, probs = c(0.05, 0.95))
s$summary[1:2,]
s <- rstan::summary(sq_fit3@stanfit, probs = c(0.05, 0.95))
s$summary[1:2,]
s <- rstan::summary(sq_fit4@stanfit, probs = c(0.05, 0.95))
s$summary[1:2,]
s <- rstan::summary(sq_fit5@stanfit, probs = c(0.05, 0.95))
s$summary[1:2,]


# Compare linear and quadratic for these variables
#  Check for evidence that quadratic is needed
#  Models where ELPD is improved by > 4 with quadratic term


f <- fitList(fit1, sq_fit1)
modSel(f)
f <- fitList(fit2, sq_fit2)
modSel(f)
f <- fitList(fit3, sq_fit3)
modSel(f)
f <- fitList(fit4, sq_fit4)
modSel(f)
f <- fitList(fit5, sq_fit5)
modSel(f)


################################ Graveyard #################################


# #for loop to convert julian day into date. Working but VERY slow
# date <- list()
# for(i in 1:nrow(dat_form)){
#   
#   origin_date = ymd(paste0(dat_form$year[i], "-01-01"))
#   
#   day(origin_date) <- as.integer(dat_form$j_day[i])
#   
#   date[[i]] <- origin_date
#   
#   print(paste0(i, "/", nrow(dat_form)))
# }