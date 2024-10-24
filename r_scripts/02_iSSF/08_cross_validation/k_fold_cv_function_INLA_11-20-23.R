#### K-fold cross-validation function ####

# Author: Read Barbee

# Date:2023-11-20

# Purpose:

#1. Partition the dataset into training and testing sets
#2. Fit model to training data
#3. Project predictions from the trained model to a raster surface (newdata=covariate layers)
#4. Bin predictions into 10 quantiles
#5. Make predictions usting the test data. (predicting whether each location in the test set will be used or unused based on the underlying covariate values)

# system.time(test <- k_fold_tmb(dat=steps_scaled, form =form, cov_stack=cov_stack_scaled, n_folds=5))

k_fold_inla <- function (dat, form, cov_stack, n_folds, mean_beta = 0, prec_beta = 1e-4 ) {
  
  #libraries
  library(tidyverse)
  library(INLA)
  library(terra)
  library(sf)
  
  #set starting parameters
  mean.beta <- mean_beta
  prec.beta <- prec_beta

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
  for (i in 1:n_folds) {
    #partition training and test sets by group
    train <- dplyr::filter(trn.tst, group != i)
    test <- dplyr::filter(trn.tst, group == i)
    
    # Fit the model to the subsetted training dataset:
    mod <- inla(form, family ="Poisson", data= train,
                         control.fixed = list(
                           mean = mean.beta,
                           prec = list(default = prec.beta)),
                         control.compute = list(waic=TRUE, dic = FALSE, cpo = FALSE)) ; beepr::beep('fanfare')
    
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
    plots[[i]] <-  ggplot() +
      geom_bar(aes(x=test_vals, y = after_stat(prop))) 
    
    n_bins <- length(unique(test_vals))
  
  #calculate pearson correlation between the bin number and the proportion of used locations in each bin
    cor_cv <- cor.test(x = 1:n_bins, y = as.numeric(table(test_vals))/(length(test_vals)/n_bins), method = "spearman", exact = FALSE)$estimate #y = num points in each bin / total num points / num bins
    cv_mat[i] <- cor_cv #save the rho value from each iteration to a matrix
    cor_cv
    print(paste0(i, "/", n_folds))
  } # end k-fold cv loop

  return(c(cv_mat, plots))
}

