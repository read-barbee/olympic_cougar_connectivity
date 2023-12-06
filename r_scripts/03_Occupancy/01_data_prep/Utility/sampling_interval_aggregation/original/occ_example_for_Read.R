
library(dplyr)
library(lubridate)
library(tidyr)
library(sf)
library(ubms)



#  Source functions to make occupancy surveys of different lenghts
source("make_int_fun.R")

# Load data
load("dat_example.Rdata") # object called dat


########################################################################
#  Sampling occasion specs
survey_period <- "day" 
number_of_periods <- "5"  # how many days you want each survey to be

# Note date column should be date format:
# $ date   : Date, format: "2009-05-02" "2009-08-10" "2009-05-14" ...


#  Format encounter history 
#   Observations
obs_int <- dat %>%
	dplyr::select(id, year, date, obs) %>% 
	as_tibble()
obs_tmp1 <- make_interval(obs_int, date, survey_period, number_of_periods) 

obs_tmp2 <- obs_tmp1 %>% 
	dplyr::group_by(id, year) %>% 
	dplyr::mutate(survey_interval = interval - min(interval) + 1) %>%
	dplyr::ungroup()
obs_tmp3 <- obs_tmp2 %>%
	dplyr::select(id, year, obs, survey_interval) %>%  
	group_by(survey_interval, id, year) %>% 
	mutate(obs = max(obs, na.rm = TRUE)) %>%
	slice(1) %>%
	as.data.frame() 	
w_obs <- obs_tmp3 %>%
	dplyr::select(id, year, obs, survey_interval) %>%  
	arrange(survey_interval, id, year) %>%
	pivot_wider(id_cols = c(id, year),
		values_from = obs, names_from = survey_interval) %>%
	as.data.frame()
w_obs[w_obs == "-Inf"] <- NA


#  Effort/Detection covariate: total number of days all cameras within a cell was operational during each interval
eff_int <- dat %>%
	dplyr::select(id, year, date, eff) %>%
	as_tibble()
eff_tmp1 <- make_interval(eff_int, date, survey_period, number_of_periods) 
eff_tmp2 <- eff_tmp1 %>% 
	dplyr::group_by(id, year) %>% 
	dplyr::mutate(survey_interval = interval - min(interval) + 1) %>%
	dplyr::ungroup()
eff_tmp3 <- eff_tmp2 %>%
	dplyr::select(id, year, eff, survey_interval) %>%  
	group_by(survey_interval, id, year) %>% 
	mutate(camdays = sum(eff, na.rm = TRUE)) %>%
	slice(1) %>%
	as.data.frame() 
w_eff <- eff_tmp3 %>%
	dplyr::select(id, year, camdays, survey_interval) %>%  
	group_by(id, survey_interval, year) %>% 
	mutate(camdays = sum(camdays, na.rm = TRUE)) %>%
	slice(1) %>%
	as.data.frame() %>%
	arrange(survey_interval, id, year) %>%
	pivot_wider(id_cols = c(id, year),
		values_from = camdays, names_from = survey_interval) %>%
	as.data.frame()
w_eff[is.na(w_eff)] <- 0	



########################################################################
#  Make objects for ubms
nsite <- nrow(w_obs)
y <- w_obs %>% 
	dplyr::select(-id, -year) %>% 
	as.matrix()

# Number of surveys detected per site
summary(rowSums(y, na.rm = TRUE))	

eff <- w_eff %>%
	dplyr::select(-id, -year) %>% 
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


