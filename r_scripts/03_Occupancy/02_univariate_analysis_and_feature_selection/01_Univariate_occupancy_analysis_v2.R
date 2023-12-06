
library(tidyverse)
library(sf)
library(ubms)

################################ User-defined Parameters #################################

#  Sampling occasion specs
survey_period <- "day" 
number_of_periods <- "5"  # how many days you want each survey to be


#########################################################################
##
##  1. Import stacked occupancy data
##
##########################################################################
#  Source functions to make occupancy surveys of different lenghts
source("r_scripts/03_Occupancy/01_data_prep/Utility/sampling_interval_aggregation/original/make_int_fun.R")

# Load data
load("r_scripts/03_Occupancy/01_data_prep/Utility/sampling_interval_aggregation/original/dat_example.Rdata") # object called dat

#activity/detection data
dat2 <- read_csv("data/Camera_Data/master/ocp_onp_occ_dat_10-06-23.csv")

#bait matrix
bait_mat <- read_csv("data/Camera_Data/Olympic_National_Park/cam_act/onp_fisher_2013_2016_bait_mat_10-11-23.csv")

snare_mat <- read_csv("data/Camera_Data/Olympic_National_Park/cam_act/onp_fisher_2013_2016_snare_mat_10-11-23.csv")


#########################################################################
##
##  2. Format data for interval function
##
##########################################################################

#select necessary columns, calculate effort column, and pivot longer
days_active_by_cell_year <- dat2 %>% 
  select(station_id_year:effort_correction_factor, d_1:d_366) %>% 
  pivot_longer(cols = d_1:d_366, names_to = "obs_p", values_to = "obs") %>% 
  group_by(cell_id, year) %>% 
  summarise(cell_year_days_active = sum(is.na(obs)==F))

dat_form <- dat2 %>% 
  left_join(days_active_by_cell_year, by = join_by(cell_id, year)) %>% 
  mutate(cell_year_effort = cell_year_days_active + effort_correction_factor, .after = effort_correction_factor) %>% 
  select(station_id_year:snare_days_good, cell_year_days_active, effort_correction_factor, cell_year_effort, d_1:d_366) %>% 
  pivot_longer(cols = d_1:d_366, names_to = "obs_p", values_to = "obs") %>% 
  mutate(j_day = str_remove(obs_p, coll("d_")), .after = obs_p,
       daily_eff = case_when(is.na(obs)==TRUE ~ 0,
                       obs >= 0 ~ 1)) 



# Create a Julian day vector and year vector
julian_days <- as.numeric(dat_form$j_day)
years <- dat_form$year

# Create origin dates for all years
origin_dates <- ymd(paste0(years, "-01-01"))

# Calculate the date by adding days to origin dates
dates <- origin_dates + days(julian_days - 1)  # Subtract 1 because Julian days start from 

#add new date column to original dataframe
dat_form <- dat_form %>% 
  mutate(date = dates)


#########################################################################
##
##  2. Format observation matrix for new survey interval
##
##########################################################################

# Note date column should be date format:
# $ date   : Date, format: "2009-05-02" "2009-08-10" "2009-05-14" ...


#Select necessary columns
obs_int <- dat_form %>%
	dplyr::select(station_id_year, cell_id, year, date, obs) %>% 
	as_tibble()

#Assign each row to an observation interval
obs_tmp1 <- make_interval(obs_int, date, survey_period, number_of_periods) 

#add column for survey interval
obs_tmp2 <- obs_tmp1 %>% 
	dplyr::group_by(cell_id, year) %>% 
	dplyr::mutate(survey_interval = interval - min(interval) + 1) %>%
	dplyr::ungroup()

#calculate the maximum number of observations for each survey interval
obs_tmp3 <- obs_tmp2 %>%
	dplyr::select(station_id_year, cell_id, year, obs, survey_interval) %>%  
	group_by(survey_interval, cell_id, year) %>% 
	mutate(obs = max(obs, na.rm = TRUE)) %>%
	slice(1) %>% #get the first row of each group
	as.data.frame() 	

#pivot from long format to wide format
w_obs <- obs_tmp3 %>%
	#dplyr::select(station_id_year, cell_id, year, obs, survey_interval) %>%  
	arrange(survey_interval, cell_id, year) %>%
	pivot_wider(values_from = obs, names_from = survey_interval) %>% 
	as.data.frame()
w_obs[w_obs == "-Inf"] <- NA

#id_cols = c(cell_id, year),

#########################################################################
##
##  2. Format effort matrix for new survey interval
##
##########################################################################

#  Effort/Detection covariate: total number of days all cameras within a cell was operational during each interval

#Select necessary columns
eff_int <- dat_form %>%
	dplyr::select(station_id_year, cell_id, year, date, daily_eff) %>%
	as_tibble()

#Assign each row to an observation interval
eff_tmp1 <- make_interval(eff_int, date, survey_period, number_of_periods) 

#add column for survey interval
eff_tmp2 <- eff_tmp1 %>% 
	dplyr::group_by(cell_id, year) %>% 
	dplyr::mutate(survey_interval = interval - min(interval) + 1) %>%
	dplyr::ungroup()

#calculate the total number of camera_days for each survey interval
eff_tmp3 <- eff_tmp2 %>%
	dplyr::select(station_id_year, cell_id, year, daily_eff, survey_interval) %>%  
	group_by(survey_interval, cell_id, year) %>% 
	mutate(camdays = sum(daily_eff, na.rm = TRUE)) %>%
	slice(1) %>%
	as.data.frame() 

#pivot from long format to wide format
w_eff <- eff_tmp3 %>%
	dplyr::select(station_id_year, cell_id, year, camdays, survey_interval) %>%  
	group_by(cell_id, survey_interval, year) %>% 
	mutate(camdays = sum(camdays, na.rm = TRUE)) %>%
	slice(1) %>%
	as.data.frame() %>%
	arrange(survey_interval, cell_id, year) %>%
	pivot_wider(values_from = camdays, names_from = survey_interval) %>%
	as.data.frame()
w_eff[is.na(w_eff)] <- 0	

#id_cols = c(cell_id, year),


#########################################################################
##
##  2. Format bait matrix for new survey interval
##
##########################################################################

#  Effort/Detection covariate: total number of days all cameras within a cell was operational during each interval

#Select necessary columns
eff_int <- dat_form %>%
  dplyr::select(station_id_year, cell_id, year, date, daily_eff) %>%
  as_tibble()

#Assign each row to an observation interval
eff_tmp1 <- make_interval(eff_int, date, survey_period, number_of_periods) 

#add column for survey interval
eff_tmp2 <- eff_tmp1 %>% 
  dplyr::group_by(cell_id, year) %>% 
  dplyr::mutate(survey_interval = interval - min(interval) + 1) %>%
  dplyr::ungroup()

#calculate the total number of camera_days for each survey interval
eff_tmp3 <- eff_tmp2 %>%
  dplyr::select(station_id_year, cell_id, year, daily_eff, survey_interval) %>%  
  group_by(survey_interval, cell_id, year) %>% 
  mutate(camdays = sum(daily_eff, na.rm = TRUE)) %>%
  slice(1) %>%
  as.data.frame() 

#pivot from long format to wide format
w_eff <- eff_tmp3 %>%
  dplyr::select(station_id_year, cell_id, year, camdays, survey_interval) %>%  
  group_by(cell_id, survey_interval, year) %>% 
  mutate(camdays = sum(camdays, na.rm = TRUE)) %>%
  slice(1) %>%
  as.data.frame() %>%
  arrange(survey_interval, cell_id, year) %>%
  pivot_wider(values_from = camdays, names_from = survey_interval) %>%
  as.data.frame()
w_eff[is.na(w_eff)] <- 0	

#id_cols = c(cell_id, year),


#########################################################################
##
##  2. Format snare matrix for new survey interval
##
##########################################################################

#  Effort/Detection covariate: total number of days all cameras within a cell was operational during each interval

#Select necessary columns
eff_int <- dat_form %>%
  dplyr::select(station_id_year, cell_id, year, date, daily_eff) %>%
  as_tibble()

#Assign each row to an observation interval
eff_tmp1 <- make_interval(eff_int, date, survey_period, number_of_periods) 

#add column for survey interval
eff_tmp2 <- eff_tmp1 %>% 
  dplyr::group_by(cell_id, year) %>% 
  dplyr::mutate(survey_interval = interval - min(interval) + 1) %>%
  dplyr::ungroup()

#calculate the total number of camera_days for each survey interval
eff_tmp3 <- eff_tmp2 %>%
  dplyr::select(station_id_year, cell_id, year, daily_eff, survey_interval) %>%  
  group_by(survey_interval, cell_id, year) %>% 
  mutate(camdays = sum(daily_eff, na.rm = TRUE)) %>%
  slice(1) %>%
  as.data.frame() 

#pivot from long format to wide format
w_eff <- eff_tmp3 %>%
  dplyr::select(station_id_year, cell_id, year, camdays, survey_interval) %>%  
  group_by(cell_id, survey_interval, year) %>% 
  mutate(camdays = sum(camdays, na.rm = TRUE)) %>%
  slice(1) %>%
  as.data.frame() %>%
  arrange(survey_interval, cell_id, year) %>%
  pivot_wider(values_from = camdays, names_from = survey_interval) %>%
  as.data.frame()
w_eff[is.na(w_eff)] <- 0	

#id_cols = c(cell_id, year),
#########################################################################
##
##  2. Make objects for ubms
##
##########################################################################

nsite <- nrow(w_obs)
y <- w_obs %>% 
	dplyr::select(-c(station_id_year:year)) %>% 
	as.matrix()

# Number of surveys detected per site
summary(rowSums(y, na.rm = TRUE))	

eff <- w_eff %>%
  dplyr::select(-c(station_id_year:year)) %>%
	as.matrix()

umf_stack <- unmarkedFrameOccu(y = y, 
	siteCovs = covs,
	obsCovs = list(eff = eff)) #, survey = surveyID))
head(umf_stack)


########################################################################
#  Corrletation plot of covaraites 
corr_dat <- covs %>% 
	dplyr::select(-id, -year, -country, -n_station)
corr_dat2 <- cor(corr_dat, use= "complete.obs")

corrplot::corrplot(corr_dat2, method = "number", type = "upper")
#ggsave(file = "covariate_correlation_plot.png")


########################################################################
# Determine best model for detection
fitd0 <- stan_occu(~1 ~1, 
   data = umf_stack, chains = 3, iter = 15000, cores = 3)
fitd0

fitd1 <- stan_occu(~scale(eff)  ~1, 
   data = umf_stack, chains = 3, iter = 15000, cores = 3)
fitd1

fitd2 <- stan_occu(~log(n_station) ~1, 
   data = umf_stack, chains = 3, iter = 15000, cores = 3)
fitd2

fitd3 <- stan_occu(~scale(eff) + log(n_station) ~1, 
   data = umf_stack, chains = 3, iter = 15000, cores = 3)
fitd3

fitd0@loo$looic
fitd1@loo$looic
fitd2@loo$looic
fitd3@loo$looic


#  Significant effects of effort and log number of stations - so take this detection
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