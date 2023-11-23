#### Global model all subsets selection INLA. All individuals, annual covariates####

# Author: Read Barbee

# Date:2023-11-21

# Purpose:


################################ Libraries #################################
library(tidyverse)
library(terra)
library(INLA)
library(beepr)
library(sf)
#library(R.utils)


#########################################################################
##
## 1. Import and format step data
##
##########################################################################

#choose whether to use steps without imputation, imputation but no rerouting, or iimputation and rerouting

#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_annual_cov_11-16-2023.csv")

#########################################################################
##
## 2. Scale continuous covariates for model comparison and correlation analysis
##
##########################################################################

#DON'T Scale step length and turn angle

steps_scaled <- steps %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dist_all_roads_annual), scale)) %>% 
  mutate(across(c(elevation:precip_annual, popdens_hii_annual:dist_all_roads_annual), as.numeric)) %>%
  rename(perc_tree_cov_annual = perc_tree_cov,
         perc_nontree_veg_annual = perc_nontree_veg,
         perc_nonveg_annual = perc_nonveg)

#DataExplorer::plot_histogram(steps_scaled %>% select(elevation:dist_all_roads_annual))

#########################################################################
##
## 3. List parameters to run all subsets for
##
##########################################################################

#List of parameters to compute all subset models for
#waic selected params
params <- c("ndvi_annual",
            "tpi" ,
            #"elevation",
            "slope",
            "aspect_northness",
            "perc_nonveg_annual",
            "perc_tree_cov_annual",
            "roads_hii_annual",
            "popdens_hii_annual",
            "distance_water" #,
            #"landuse_hii_annual",
            #"infra_hii_annual",
            #"aspect_eastness"
)

#full hii only to reduce combinations ~ 2.32 hours
# params <- c("ndvi_annual",
#             "tpi" ,
#             #"elevation",
#             "slope",
#             "aspect_northness",
#             "perc_nonveg_annual",
#             "perc_tree_cov_annual",
#             "hii_annual",
#             "distance_water"
#             #"aspect_eastness"
# )

#List of quadratic parameters if necessary 

#quad_params <- c("ndvi_annual"#,
  #"perc_tree_cov_annual"
  #"popdens_hii_annual"
#)

#########################################################################
##
## 4. Set prior parameters
##
##########################################################################
#prior params from https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40#inla-1
mean.beta <- 0
prec.beta <- 1e-4 


#########################################################################
##
## 5. Generate INLA model terms for all combination list
##
##########################################################################

#get number of individuals from dataframe
n_indiv = steps_scaled %>% distinct(animal_id) %>% count() %>% pull()


# quad_terms <- vector()
# 
# for(i in 1:length(quad_params)){
#   quad_terms[i] <- paste0("I(",quad_params[i],"^2)")
# }

quad_terms <- NULL

vars <- vector()
for(i in 1:length(params)){
  vars[i] <- paste0(params[i]," + f(", paste0("id", i), ", ", params[i], ", ", "values = ", paste0("1:",n_indiv), ", ", 'model = "iid", ', "hyper = list(theta = list(initial = log(1), fixed = FALSE, ", 'prior = "pc.prec", ', "param = c(3, .05))))")
}

#convert response from categorical to numeric
dat <- steps_scaled %>% mutate(case_ = as.numeric(case_))

#create separate animal_id columns for each random effect
for(i in 1:length(vars)){
  name <- as.name(paste0("id", i))
  dat[[name]] <- as.numeric(factor(dat$animal_id))
}

#########################################################################
##
## 6. Make list of all term combinations
##
##########################################################################
#create list of all combinations of variables
all_comb <- do.call("c", lapply(seq_along(vars), function(i) combn(vars, i, FUN = list)))
response <- "case_"


#########################################################################
##
## 7. Import cov stack for CV
##
##########################################################################
cov_stack <- rast("/Users/tb201494/Desktop/1km_buffer/cov_stack_pred_300m_11-21-2023.tif")


#########################################################################
##
## 8. Function to genrate two_thirds CV score for each model
##
##########################################################################
two_thirds_cv <- function (dat, form, cov_stack, n_folds = 3, mean_beta = 0, prec_beta = 1e-4 ) {
  #set starting parameters
  mean.beta <- mean_beta
  prec.beta <- prec_beta
  
  set.seed(777)
  #Randomly assign each individual to one of n groups
  nested <- dat %>% nest_by(animal_id)
  n <- ceiling(nrow(nested)/n_folds)
  group <- c(replicate(n, sample(1:n_folds, n_folds, replace = F)))
  group <- group[c(1:nrow(nested))]
  nested$group <- group
  trn.tst <- nested %>% unnest(cols=c(data)) %>% ungroup()
  
  # Initialize data storage objects
  cv_mat <- matrix(NA, ncol = n_folds, nrow = 1)
  plots <- list()
  
  #Run the k-fold CV
  
  #partition training and test sets by group
  train <- dplyr::filter(trn.tst, group != 1)
  test <- dplyr::filter(trn.tst, group == 1)
  
  # Fit the model to the subsetted training dataset:
  mod <- inla(form, family ="Poisson", data= train,
              control.fixed = list(
                mean = mean.beta,
                prec = list(default = prec.beta)),
              control.compute = list(waic=TRUE, dic = FALSE, cpo = FALSE)) #; beepr::beep('fanfare')
  
  #extract coefficient names and values from the fitted model
  mod_sum <- summary(mod)
  coeffs <- mod_sum$fixed %>% as.data.frame() #%>% rownames_to_column("param")
  names <- rownames(coeffs)
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
  test_steps <- dplyr::filter(test, case_ == T)
  
  #convert end points of used steps to sf points and extract values from the binned raster predictions
  test_sf <- st_as_sf(test_steps, coords = c("x2_", "y2_"), crs=5070)
  test_vals <- as.numeric(terra::extract(binned, test_sf)[,2]) %>% na.omit()
  
  #plot proportion of test points in each bin. Currently storing all of the same plot for some reason.
  # plots[[i]] <-  ggplot() +
  #   geom_bar(aes(x=test_vals, y = after_stat(prop))) 
  
  n_bins <- length(unique(test_vals))
  
  #calculate pearson correlation between the bin number and the proportion of used locations in each bin
  cor_cv <- cor.test(x = 1:n_bins, y = as.numeric(table(test_vals))/(length(test_vals)/n_bins), method = "spearman", exact = FALSE)$estimate #y = num points in each bin / total num points / num bins
  #cv_mat[i] <- cor_cv #save the rho value from each iteration to a matrix
  #cor_cv
  #print(paste0(i, "/", n_folds))
  
  
  return(cor_cv)
}


#########################################################################
##
## 9. Run all subsets
##
##########################################################################

#initialize lists to store model formulas and waic scores
mod_names <- list()
waic_list <- list()
cv_list <- list()
## Took 32273 seconds (~9hrs) on Legion for 11 covariates (2047 combinations)
system.time(
  for(i in 1:length(all_comb)){ #length(all_comb)
    #start_time <- Sys.time() # Record start time
    var_i <- all_comb[[i]]
    form <- as.formula(paste(response, paste("-1", "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T)))", paste(var_i, collapse="+"), sep="+"), sep="~"))
    
    mod <- inla(form, family ="Poisson", 
                      data=dat,
                      #control.fixed = list(
                      #mean = mean.beta,
                      #prec = list(default = prec.beta)),
                      control.compute = list(waic=TRUE),
                      safe = TRUE) #dic = TRUE,
    
    mod_names[[i]] <- paste(mod$names.fixed, collapse="+")
    
    waic_list[[i]] <- mod$waic$waic
    
    cv_list[[i]] <- two_thirds_cv(dat, form, cov_stack)
    
    print(paste0(i, "/", length(all_comb)))
  }
)



#########################################################################
##
## 10. Make table of output summaries for export
##
##########################################################################

mod_names2 <- mod_names[485:511]
waic_list2 <- waic_list[485:511]
cv_list2 <- cv_list[485:511]

#dredge selection table
dredge_table_inla <- tibble(model = unlist(as.character(mod_names2)),
                           waic = unlist(waic_list2),
                           rho = unlist(cv_list2))
  

write_csv(dredge_table_inla, "feature_selection/all_subsets_model_selection_annual_covs_inla_hii_subs_485_511_11-22-23.csv")

saveRDS(dredge_table_inla, "feature_selection/all_subsets_model_selection_annual_covs_inla_hii_subs_485_511_11-22-23.rds")

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

