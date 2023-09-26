### Sells Individual Predictive iSSF Example Code ###

# I worked with one bear's data subset at a time. The basic idea is to "turn off" various covariates iteratively and run a k-folds cross validation to identify which variation of the model has the highest Spearman rank score for that bear. I selected one of those two models (global or catered) as the final version for that bear depending on which performed best. 
# 
# Firstly, get R set up:

### Grizzly Bear CV ###
### S Sells, MTCWRU ###

library(amt)
library(purrr)
library(raster)
library(tidyverse)
library(broom)
library(doParallel)
library(grizmove)

# Change as needed: while identifying potential variations of the catered model, run fewer iterations for speed.
n.cv <- 5 # number of cross validation iterations (try using same value as # of clusters when you are trialing different catered models)
n.clusters <- 5 # number of clusters for parallel run (much faster to run in parallel)

# Load data and initialize results tables:
cov_stack <- load_cov_data(sex = "Male", extent = "NCDE_300km") # load covariate rasters
dat_all <- grizmove::bear_df %>% filter(sex == "Male") # load bear location dataset
bear_mat <- matrix(NA, ncol = 4, nrow = n.cv) # initialize results table

#Subset and prepare the first individual's dataset:
  
  #### Bear 1: Adam
  # This takes the full dataset, subsets to Adam, and preps the basic amt-style track df including sampled steps:
  track <- dat_all %>%
    filter(bear == "Adam") %>% arrange(t) %>%
    make_track(x, y, t, crs = sp::CRS("+proj=utm +zone=12 +datum=WGS84")) %>%
    amt::transform_coords(raster.crs) %>% # transform to the same coord sys as the raster data, which I saved to easily call as "raster.crs"
    track_resample(rate = hours(3), tolerance = minutes(45)) %>% steps_by_burst() %>% 
    mutate(bear_name = "Adam") %>% 
    filter(sl_ >= 100) %>% # filter steps to at least 100m apart (we wanted movements, not stationary bears)
    filter(sl_ <= 15000) %>% # filter steps that are unreasonably long, as these are probably inaccurate
    dplyr::arrange(burst_, t1_, .by_group = T) %>%
    amt::random_steps(10) %>% 
    amt::extract_covariates(cov_stack, where = "end") %>% 
    dplyr::mutate(log_sl_ = log(sl_), cos_ta_ = cos(ta_))
  
  #Next, prepare to run the k-folds CV, starting with the global model:
    
    # Start with global model: set up the function will use:
    run_global <- function(dat) { mod <- dat %>% amt::fit_issf(case_ ~ 
                                                                 ndvi_e + I(ndvi_e^2) +
                                                                 ruggedness_e + I(ruggedness_e^2) + 
                                                                 d2forestedge_e + I(d2forestedge_e^2) + 
                                                                 densforestedge_e + I(densforestedge_e^2) +
                                                                 densriparian_e + I(densriparian_e^2) + 
                                                                 densbuildings_e + I(densbuildings_e^2) + 
                                                                 d2core_e + I(d2core_e^2) + 
                                                                 sl_ + log_sl_ + cos_ta_ + 
                                                                 strata(step_id_), model = TRUE) }
  
  #Running everything in parallel is super helpful for speed:
    
    # Run global cross-validation ncde_4fold in parallel:
    my.cluster <- parallel::makeCluster(n.clusters, type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster, cores = ncores)
  bear_mat <- foreach(s = 1:n.cv, .combine = rbind, .packages = c("tidyverse", "amt", "raster", "grizmove")) %dopar%
    ncde_4fold(dat = track, cov_stack, mod = "global")
  parallel::stopCluster(cl = my.cluster) # once done, manually stop the cores from continuing to act busy
  
#Here's the guts of the whole operation:  ncde_4fold function:

ncde_4fold <- function (dat, mod, cov_stack) {
  # Split out the used GPS points and the available points:
  case_t <- dplyr::filter(dat, case_ == T)
  case_f <- dplyr::filter(dat, case_ == F) %>% dplyr::mutate(group = NA)
  # Divide the datasets into 4 folds:
  n <- ceiling(nrow(case_t)/4)
  group <- c(replicate(n, sample(1:4, 4, replace = F)))
  group <- group[c(1:nrow(case_t))]
  trn.tst <- cbind(case_t, group)
  trn.tst <- rbind(trn.tst, case_f) %>% dplyr::arrange(t1_) %>% fill(group)
  
  # Run the 4-fold CV:
  cv_mat <- matrix(NA, ncol = 4, nrow = 1)
  for (i in 1:4) {
    train <- dplyr::filter(trn.tst, group != i)
    test <- dplyr::filter(trn.tst, group == i)
  # Fit the model to the subsetted training dataset:
    if (mod == "global") { mod.fit <- run_global(dat = train) }
    if (mod == "catered") { mod.fit <- run_catered(dat = train) }
 
    # Apply model to the landscape as a RSF-style raster map:
    mov <- train %>% dplyr::filter(case_ == T) %>% dplyr::select(sl_, log_sl_, cos_ta_)
    cov_stack_values <- values(cov_stack) %>% data.frame
    cov_stack_values$sl_ <- median(mov$sl_) # we need these items to predict from the model, so take medians.
    cov_stack_values$log_sl_ <- median(mov$log_sl_)
    cov_stack_values$cos_ta_ <- median(mov$cos_ta_, na.rm = T)
    cov_stack_values$step_id_ = mod.fit$model$model$`strata(step_id_)`[1]
    predictions <- raster::predict(mod.fit$model, newdata = cov_stack_values, type = "lp", allow.new.levels = TRUE)
    
    # Now test how well this map predicts the test locations:
    train_ssf_vals_df <- as.data.frame(as.vector(predictions))
    train_ssf_vals_df <- within(train_ssf_vals_df, bins <- cut(predictions, quantile(predictions, probs = 0:10/10, na.rm = T), include.lowest = TRUE))
    test_steps <- dplyr::filter(test, case_ == T)
    test_steps$step_id_ = mod.fit$model$model$`strata(step_id_)`[1]
    test_ssf_vals <- predict(mod.fit$model, newdata = test_steps,  type = "lp", allow.new.levels = TRUE)
    test_ssf_vals_df <- as.data.frame(test_ssf_vals)
    test_ssf_vals_df <- within(test_ssf_vals_df, bins <- cut(test_ssf_vals, quantile(predictions, probs = 0:10/10, na.rm = T), include.lowest = TRUE))
    cor_cv <- cor.test(x = 1:10, y = as.numeric(table(test_ssf_vals_df$bins))/(nrow(test_steps)/10), method = "spearman", exact = FALSE)$estimate
    cv_mat[i] <- cor_cv
    cor_cv
  }
  return(cv_mat)
}

#Summarize the results of CV for the global model:

mean(bear_mat)
# This is the overall mean of however many CV iterations you ran (as specified by n.cv). Save this to a useful location and make 
# sure to label it as the global model to compare with the catered next.

#Repeat all this, but for a catered model:

# Build catered: 
run_catered <- function(dat){
  mod <- dat %>% fit_issf(case_ ~
                            ndvi_e + #I(ndvi_e^2)+
                            ruggedness_e + #I(ruggedness_e^2) +
                            d2forestedge_e + I(d2forestedge_e^2) +
                            # densforestedge_e + #I(densforestedge_e^2) +
                            densriparian_e + #I(densriparian_e^2) +
                            # densbuildings_e + #I(densbuildings_e^2) +
                            # d2core_e + I(d2core_e^2) +
                            sl_ + log_sl_ + cos_ta_ +
                            strata(step_id_), model = TRUE)
}
mod <- run_catered(track) # I found it helpful to look at the parameters overlapping 0; often, these could be dropped to improve CV scores. 

# Run CV Catered:
my.cluster <- parallel::makeCluster(n.clusters, type = "PSOCK")
clusterExport(my.cluster, c("run_catered"))
doParallel::registerDoParallel(cl = my.cluster, cores = ncores)
bear_mat <- foreach(s = 1:n.cv, .combine = rbind,.packages = c("tidyverse", "amt", "raster", "grizmove")) %dopar%
  ncde_4fold(dat = track, cov_stack, mod = "catered")
parallel::stopCluster(cl = my.cluster) # once done, manually stop the cores from continuing to act busy

# Summarize results:
mean(bear_mat) # again save this to a useful location, label as "catered".

#Note that this above section for the catered model is where I would manually turn on/off the parameters using #, and then review model fit and try a quick CV run, using just 5 iterations of the entire CV process (n.cv in the very top end of this code) so that it took a few minutes to run each time. Each iteration of the CV will use ncde_4fold to split up the dataset 4 ways and test fit, so running 5 iterations would give you 5 Spearman rank scores that you can average to see if this is overall likely better than another model iteration. 

#Once I identified a best catered model for that individual, I ran through 100 iterations of n.cv for each the global and catered model and got a final Spearman rank score to differentiate which of the 2 models was best, global or catered.

#Hope that helps? It's a lot to work through.
