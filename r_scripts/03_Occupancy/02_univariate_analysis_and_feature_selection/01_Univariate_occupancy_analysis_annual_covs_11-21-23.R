#### Univariate occupancy analysis ####

# Author: Read Barbee

# Date:2023-10-12 

#Last updated: 2023-12-07

# Purpose:



library(tidyverse)
library(sf)
library(terra)
library(ubms)
library(beepr)
library(DataExplorer)
library(doParallel)

################################ User-defined Parameters #################################


#########################################################################
##
##  1. Import stacked occupancy data
##
######################################################################

#activity/detection data
occ_dat <- read_csv("data/Camera_Data/master/ocp_onp_occ_dat_annual_covs_12-07-23.csv")


#########################################################################
##
## 2. Check covariate distributions
##
##########################################################################

cov_check <- occ_dat %>% 
  select(elevation:dens_all_roads_annual) %>% 
  mutate(across(land_cover_usfs_annual:land_use_change_usfs_annual, as.factor))

#data overview
introduce(cov_check)
plot_intro(cov_check)

#check for missing values
plot_missing(cov_check)


#continuous histograms all
plot_histogram(cov_check)


#categorical bar plots
plot_bar(cov_check %>% select(land_cover_usfs_annual, land_use_usfs_annual, land_use_change_usfs_annual))


#full report with all DataExplorer metrics
#create_report(steps, y="case_")


#########################################################################
##
##  2. Make objects for ubms
##
##########################################################################
occ_dat_scaled <- occ_dat %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dens_all_roads_annual), scale)) %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dens_all_roads_annual), as.numeric))


#remove station rows with missing covariate values
complete_cases <- occ_dat_scaled %>% select(c(elevation:precip_annual, popdens_hii_annual:dens_all_roads_annual)) %>% 
  complete.cases()

occ_dat_complete <- occ_dat_scaled %>% 
  mutate(comp = complete_cases) %>% 
  #unite("cell_id_year", cell_id, year, remove = FALSE) %>% 
  #mutate(cell_id_year = as.factor(cell_id_year)) %>% 
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

covs_scaled <- occ_dat_complete %>% select(cell_id, year, elevation:dens_all_roads_annual)

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
corr_dat <- covs_scaled %>% 
	dplyr::select(-c(cell_id, year))


cont <- covs_scaled %>% 
  select(elevation:dens_all_roads_annual, -aspect) %>% 
  mutate(across(land_cover_usfs_annual:land_use_change_usfs_annual, as.factor)) %>% 
  dummify(select = c("land_cover_usfs_annual", "land_use_usfs_annual", "land_use_change_usfs_annual")) %>%
  #mutate(across(elevation:dist_all_roads_annual, as.numeric)) %>% 
  select(!contains("NA"))

cor_matrix <- cor(cont, use= "complete.obs")

#corrplot::corrplot(corr_dat, method = "number", type = "upper")
#ggsave(file = "covariate_correlation_plot.png")

threshold <- 0.0

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

#write_csv(cor_dat, "feature_selection/pairwise_cov_corr_occ_2_week_12-07-23.csv")

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

#write_csv(occu_det_summ, "/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/feature_selection/occu_det_model_fits_14_day_12-7-23.csv")


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

null_fit <- stan_occu(~ scale(eff) + scale(bait) + scale(snare) ~ (1|cell_id), data=umf_stack, chains=3, iter=1500)

traceplot(null_fit)

traceplot(null_fit, pars = "sigma_state", inc_warmup=TRUE)

null_fit_list <- list(null_fit)

names(null_fit_list) <- "null"



cov_names <- covs_scaled %>% select(-c(cell_id, year)) %>% names()

#fit univariate models
uni_fits <- list()
system.time(for(i in 1:length(cov_names)){
  cov <- cov_names[i]
  
  #construct the model formula with linear terms only
  form <- as.formula(paste0("~ scale(eff) + scale(bait) + scale(snare) ~ ", cov, " + (1|cell_id)"))
  
  
  uni_fits[[i]] <- stan_occu(form, data = umf_stack, chains = 3, iter = 1500, cores = 3)
  
  print(paste0(i, "/", length(cov_names)))
  
}) ; beep("fanfare")


names(uni_fits) <- cov_names

uni_fits <- c(uni_fits, null_fit_list)

#uni_fits <- uni_fits[1:35]

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
  print(paste0(i, "/", length(uni_fits)))
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

#write_csv(occu_uni_summ, "feature_selection/occu_uni_fits_14day_cam_bait_snare_12-7-23.csv")
#RESUME HERE. SAVE MODEL STATS, THEN DO THE SAME WITH QUADRATICS
#CONSIDER DETECTION COVARIATE FOR PARK VS OCP INSTEAD OF BAIT AND SNARE


#########################################################################
##
##  4. Fit quadratic univariate stacked occupancy models
##
##########################################################################

#parallelized for loop: fits in ~ 10 minutes with 3 chains and 100 iterations; 35 min for 1000 iterations

#~ 9 hours for 35 models
uni_fits_quad <- list()
system.time(uni_fits_quad <- for (i in 1:length(cov_names)){
  
  cov <- cov_names[i]
  
  form <- as.formula(paste0("~ scale(eff) + scale(bait) + scale(snare) ~", #remove standard intercept to be replace with stratum-based intercept
                            #fixed effects
                            cov, "+", "I(", cov,  "^2) + ",
                            #random intercept (strata)
                            "(1|cell_id)"))
  
  
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