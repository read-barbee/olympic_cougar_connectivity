#### Global model all subsets selection occupancy annual covariates####

# Author: Read Barbee

# Date:2023-11-26
# Last Updated:2023-12-07

# Purpose:


################################ Libraries #################################
library(tidyverse)
library(sf)
library(terra)
library(ubms)
library(beepr)
library(DataExplorer)


#########################################################################
##
##  1. Import stacked occupancy data
##
######################################################################

#activity/detection data
occ_dat <- read_csv("data/Camera_Data/master/ocp_onp_occ_dat_annual_covs_14_day_period_11-21-23.csv")

#########################################################################
##
##  2. Make objects for ubms
##
##########################################################################
occ_dat_scaled <- occ_dat %>% 
  select(!contains("usfs")) %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dist_all_roads_annual), scale)) %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dist_all_roads_annual), as.numeric))


#remove station rows with missing covariate values
complete_cases <- occ_dat_scaled %>% select(elevation:dist_all_roads_annual) %>% 
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

covs_scaled <- occ_dat_complete %>% select(cell_id, year, elevation:dist_all_roads_annual)

umf_stack <- unmarkedFrameOccu(y = y, 
                               siteCovs = covs_scaled,
                               obsCovs = list(eff = eff,
                                              bait = bait,
                                              snare = snare)) #, survey = surveyID))
head(umf_stack)


#########################################################################
##
## 3. List parameters to run all subsets for
##
##########################################################################

#List of parameters to compute all subset models for
#waic selected params
params <- c("aspect_northness",
            #"aspect_eastness" ,
            "evi_annual",
            "dist_minor_roads_annual",
            "perc_nonveg_annual",
            "tri",
            "hii_annual",
            "mtpi",
            #"gpp_annual",
            "distance_water"
            #"precip_annual"
)


#List of quadratic parameters if necessary 

#quad_params <- c("ndvi_annual"#,
  #"perc_tree_cov_annual"
  #"popdens_hii_annual"
#)


#########################################################################
##
## 6. Make list of all term combinations
##
##########################################################################
#create list of all combinations of variables
all_comb <- do.call("c", lapply(seq_along(params), function(i) combn(params, i, FUN = list)))

# all_comb <- do.call("c", lapply(seq_along(vars), function(i) combn(vars, i, FUN = list)))
# response <- "case_"

#########################################################################
##
## 7. Import cov stack for CV
##
##########################################################################
cov_stack <- rast("/Users/tb201494/Desktop/1km_buffer/cov_stack_pred_300m_11-21-2023.tif")

steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_annual_cov_11-16-2023.csv")

#########################################################################
##
## 8. Function to genrate two_thirds CV score for each model
##
##########################################################################
occu_cv <- function (fitted_model, cov_stack) {

  set.seed(777)

  
  #extract coefficient names and values from the fitted model
  mod_sum <- summary(fitted_model, "state")
  coeffs <- mod_sum %>% as.data.frame() 
  names <- coeffs %>% 
    rownames_to_column("param") %>% 
    filter(str_detect(.$param, coll("(Intercept)"), negate = T)) %>% 
    filter(str_detect(.$param, coll("sigma [1|cell_id_year]"), negate = T)) %>% 
    pull(param)
  #coeffs$param
  
  cov_stack_sel <- cov_stack[[names]]
  
  #Generate predictions by multiplying each covariate layer by its coefficient estimate from the model
  pred <- list()
  for(j in 1:length(names)){
    pred[[j]] <- cov_stack_sel[[j]] * coeffs[names[j], "mean"]
    #print(paste0(j, "/", length(names)))
  }
  predictions <- Reduce("+", pred) #sum layers together
  pred_vals <- terra::values(predictions) #convert raster to values for calculations
  preds_exp <- plogis(pred_vals) #backtransform predictions from logit to prob scale
  breaks <- quantile(preds_exp, probs = 0:10/10, na.rm = T) #obtain 10 quantile values for preds
  breaks_j <- breaks + (seq_along(breaks) * .Machine$double.eps) #add jitter to breaks
  vals <- cov_stack_sel[[1]] %>% terra::setValues(preds_exp) #project those values onto the map
  binned <- terra::classify(vals, rcl=as.vector(breaks_j)) #classify the map based on breaks
  
  #Subset used steps from the test data
  test_steps <- dplyr::filter(steps, case_ == T)
  
  #convert end points of used steps to sf points and extract values from the binned raster predictions
  test_sf <- st_as_sf(test_steps, coords = c("x2_", "y2_"), crs=5070)
  test_vals <- as.numeric(terra::extract(binned, test_sf)[,2]) %>% na.omit()
  
  #plot proportion of test points in each bin. Currently storing all of the same plot for some reason.
  plot <-  ggplot() +
    geom_bar(aes(x=test_vals, y = after_stat(prop)))
  
  n_bins <- length(unique(test_vals))
  
  #calculate pearson correlation between the bin number and the proportion of used locations in each bin
  cor_cv <- cor.test(x = 1:n_bins, y = as.numeric(table(test_vals))/(length(test_vals)/n_bins), method = "spearman", exact = FALSE)$estimate #y = num points in each bin / total num points / num bins
  #cv_mat[i] <- cor_cv #save the rho value from each iteration to a matrix
  #cor_cv
  #print(paste0(i, "/", n_folds))
  
  
  return(list(cor_cv, plot))
}


#########################################################################
##
## 9. Run all subsets (~ 24 hours for 255 combinations)
##
##########################################################################

#initialize lists to store model formulas and waic scores
mod_names <- list()
waic_list <- list()
looic_list <- list()
cv_list <- list()
## Took 32273 seconds (~9hrs) on Legion for 11 covariates (2047 combinations)
system.time(
  for(i in 1:length(all_comb)){ #length(all_comb)
    #start_time <- Sys.time() # Record start time
    var_i <- all_comb[[i]]
    
    form <- as.formula(paste0("~ scale(eff) + scale(bait) + scale(snare) ~ ", paste(var_i, collapse="+"), " + (1|cell_id_year)"))
    
    # form <- as.formula(paste(response, paste("-1", "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T)))", paste(var_i, collapse="+"), sep="+"), sep="~"))
    
    mod <- stan_occu(form, data = umf_stack, chains = 3, iter = 1500, cores = 3)
    
    mod_names[[i]] <- paste(var_i, collapse="+")
    
    waic_list[[i]] <-waic(mod)$estimates %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "metric") %>% 
      pivot_wider(names_from = "metric", values_from = c("Estimate", "SE")) 
    
    looic_list[[i]] <-loo(mod)$estimates %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "metric") %>% 
      pivot_wider(names_from = "metric", values_from = c("Estimate", "SE")) 
    
    cv_list[[i]] <- occu_cv(mod, cov_stack)
    
    print(paste0(i, "/", length(all_comb)))
  }
)



 #########################################################################
##
## 10. Make table of output summaries for export
##
##########################################################################

# mod_names2 <- mod_names[485:511]
# waic_list2 <- waic_list[485:511]
# cv_list2 <- cv_list[485:511]

waic_tab <- bind_rows(waic_list)
looic_tab <- bind_rows(looic_list)

aic_tab <- bind_cols(waic_tab, looic_tab)

rho_scores <- vector()

for(i in 1:length(cv_list)){
  rho_scores[i] <- cv_list[[i]][[1]]
}

#dredge selection table
dredge_table_occu <- tibble(model = unlist(as.character(mod_names)),
                           rho = rho_scores) %>% 
                      bind_cols(aic_tab)
  

write_csv(dredge_table_occu, "feature_selection/all_subsets_model_selection_annual_covs_occu_11-28-23.csv")

# saveRDS(dredge_table_inla, "feature_selection/all_subsets_model_selection_annual_covs_inla_hii_subs_485_511_11-22-23.rds")

################################ Graveyard #################################

#time out stuff for INLA calls not working
# tryCatch({
#   mod <- withTimeout(inla(form, family ="Poisson", data=dat,
#                           control.fixed = list(
#                             mean = mean.beta,
#                             prec = list(default = prec.beta)),
#                           control.compute = list(waic=TRUE)),
#                      timeout = 100)}, #dic = TRUE,
#   TimeoutException = function(ex){
#     cat("Operation took more than 5 minutes. Exiting loop.\n")
#     terminated_indices <- c(terminated_indices, i)
#     next 
#   })

# # Check elapsed time
# elapsed_time <- timing["elapsed"]
# 
# # If elapsed time exceeds 5 minutes, break out of the loop
# if (elapsed_time > 300) {
#   cat("Operation took more than 5 minutes. Exiting loop.\n")
#   terminated_indices <- c(terminated_indices, i)
#   next  # Skip to the next iteration
# }

