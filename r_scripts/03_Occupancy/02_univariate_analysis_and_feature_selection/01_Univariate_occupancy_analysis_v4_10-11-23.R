
library(tidyverse)
library(sf)
library(terra)
library(ubms)

################################ User-defined Parameters #################################

#  Sampling occasion specs
survey_period <- "day" 
number_of_periods <- "7"  # how many days you want each survey to be


#########################################################################
##
##  1. Import stacked occupancy data
##
##########################################################################
#  Source functions to make occupancy surveys of different lenghts
source("r_scripts/03_Occupancy/01_data_prep/Utility/sampling_interval_aggregation/original/make_int_fun.R")

#example data
#load("r_scripts/03_Occupancy/01_data_prep/Utility/sampling_interval_aggregation/original/dat_example.Rdata") # object called dat

#activity/detection data
dat <- read_csv("data/Camera_Data/master/ocp_onp_occ_dat_long_10-11-23.csv")

# Load raster data
cov_stack <- rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2_asp.tif")


#########################################################################
##
##  2. Format data for new survey interval
##
##########################################################################

# Note date column should be date format:
# $ date   : Date, format: "2009-05-02" "2009-08-10" "2009-05-14" ...

#Assign each row to an observation interval
obs_tmp1 <- make_interval(dat, date, survey_period, number_of_periods) 

#add column for survey interval
obs_tmp2 <- obs_tmp1 %>% 
	dplyr::group_by(cell_id, year) %>% 
	dplyr::mutate(survey_interval = interval - min(interval) + 1) %>%
	dplyr::ungroup()

#calculate the maximum number of observations for each survey interval
obs_tmp3 <- obs_tmp2 %>%
	group_by(survey_interval, cell_id, year) %>% #station_id, 
  summarize(station_id = first(station_id),
            grid_id = first(grid_id),
            lon = first(lon),
            lat = first(lat),
            utm_e = first(utm_e),
            utm_n = first(utm_n),
            annual_effort_correction_factor = first(effort_correction_factor),
            annual_effort_correction_factor_bait = first(effort_correction_factor_bait),
            annual_effort_correction_factor_snare = first(effort_correction_factor_snare),
            cam_days = sum(cam_status, na.rm = TRUE),
            bait_days = sum(bait_status, na.rm = TRUE),
            snare_days = sum(snare_status, na.rm = TRUE),
            cougar_detections_binary = case_when(all(is.na(cougar_detections_binary))==T ~NA,
                                                 sum(cougar_detections_binary, na.rm = TRUE) == 0 ~ 0,
                                                 sum(cougar_detections_binary, na.rm = TRUE) > 0 ~ 1))



w_obs <- obs_tmp3 %>%
	arrange(survey_interval, cell_id, year) %>%
	pivot_wider(values_from = c(cam_days, bait_days, snare_days, cougar_detections_binary), names_from = survey_interval) %>% 
 mutate(across(!contains("detections") , replace_na, 0)) %>% 
	as.data.frame()
	
#w_obs[w_obs == "-Inf"] <- NA

#id_cols = c(cell_id, year),


#########################################################################
##
##  2. Make objects for ubms
##
##########################################################################


sf_pts <- w_obs %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  st_transform(crs = 5070)

#mapview::mapview(sf_pts)

#Extract cell numbers for each station location
covs <-terra::extract(cov_stack, sf_pts)


#construct umf stack

nsite <- nrow(w_obs)
y <- w_obs %>% 
	dplyr::select(contains("detections")) %>% 
	as.matrix()

# Number of surveys detected per site
summary(rowSums(y, na.rm = TRUE))	

eff <- w_obs %>%
  dplyr::select(contains("cam")& !contains("correction")) %>%
	as.matrix()

bait <- w_obs %>%
  dplyr::select(contains("bait")& !contains("correction")) %>%
  as.matrix()

snare <- w_obs %>%
  dplyr::select(contains("snare") & !contains("correction")) %>%
  as.matrix()

umf_stack <- unmarkedFrameOccu(y = y, 
	siteCovs = covs,
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

#write_csv(cor_dat, "feature_selection/pairwise_cov_corr_occ_10-12-23.csv")

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

#write_csv(occu_det_summ, "/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/feature_selection/Occupancy/occu_det_model_fits_daily_10-12-23.csv")


#  Significant effects of: ___ so take this detection
#   model forward
# save(fitd3, file = "Single_covariate_models/fitd3.Rdata")


########################################################################
#  Run single covariate models
 
fit1 <- stan_occu(~ scale(eff) + log(n_station) 
		~ scale(asp) + (1|id), 
   data = umf_stack, chains = 3, iter = 15000, cores = 3)
fit2 <- stan_occu(~ scale(eff)+ log(n_station) 
		~ scale(elv) + (1|id), 
   data = umf_stack, chains = 3, iter = 15000, cores = 3)
fit3 <- stan_occu(~ scale(eff) + log(n_station) 
		~ scale(slp) + (1|id), 
   data = umf_stack, chains = 3, iter = 15000, cores = 3)
fit4 <- stan_occu(~ scale(eff) + log(n_station) 
		~ scale(bc_prc) + (1|id), 
   data = umf_stack, chains = 3, iter = 15000, cores = 3)
fit5 <- stan_occu(~ scale(eff) + log(n_station) 
		~ scale(bc_tmp_mu) + (1|id), 
   data = umf_stack, chains = 3, iter = 15000, cores = 3)  

# save(fit1, file = "fit1.Rdata")
# save(fit2, file = "fit2.Rdata")
# save(fit3, file = "fit3.Rdata") 
# save(fit4, file = "fit4.Rdata")
# save(fit5, file = "fit5.Rdata")


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