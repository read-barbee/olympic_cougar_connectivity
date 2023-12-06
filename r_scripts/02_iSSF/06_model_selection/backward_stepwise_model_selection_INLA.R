#### Backward stepwise model selecion in INLA ####

# Author: Read Barbee

# Date:2023-10-13

# Purpose:


################################ libraries #################################
library(tidyverse)
library(INLA)
library(beepr)
library(INLAutils)
library(terra)
library(sf)

#########################################################################
##
##  Specify Model Parameters
##
##########################################################################

params <- c("tree_cover_hansen",
            "ndvi",
            #"tpi",
            #"northing",
            "tri",
            "perc_tree_cover",
            "perc_nonveg",
            "roads_hii",
            "elevation")
            # "popdens_hii",
            # "landuse_hii",
            # "rails_hii",
            # "infra_hii",
            # "easting")

quad_params <- c("tree_cover_hansen",
                 "ndvi",
                 "tri",
                 "perc_tree_cover",
                 "roads_hii",
                 "elevation")
                 # "popdens_hii",#not sure that this makes biological sense
                 # "landuse_hii") #not sure that this makes biological sense
                 
                 
                 #"tpi",
                 #"northing", #not sure that this makes biological sense

             
                 
               



#prior params from https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40#inla-1
mean.beta <- 0
prec.beta <- 1e-4 


quad_terms <- vector()

for(i in 1:length(quad_params)){
  quad_terms[i] <- paste0("I(",quad_params[i],"^2)")
}


#########################################################################
##
## 1. Import and format step data
##
##########################################################################

#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_10-02-2023.csv") %>% 
  mutate(ndvi = ndvi*0.0001)



#set all negative elevations to 0 and filter dispersal tracks for post dispersal event
steps_scaled <- steps %>% 
  filter(dispersal_status=="resident") %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(land_use_usfs, land_cover_usfs, dispersing:disp_qual, season:calving_season)) %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric))



#########################################################################
##
## 2. Fit global model with quadratics in INLA
##
##########################################################################


n_indiv = steps_scaled %>% distinct(animal_id) %>% count() %>% pull()

rand_terms_inla <- vector()
for(i in 1:length(params)){
  rand_terms_inla[i] <- paste0("f(", paste0("id", i), ", ", params[i], ", ", "values = ", paste0("1:",n_indiv), ", ", 'model = "iid", ', "hyper = list(theta = list(initial = log(1), fixed = FALSE, ", 'prior = "pc.prec", ', "param = c(1, .05))))")
}



run_global_inla <- function (dat){
  
  #convert response from categorical to numeric
  dat <- dat %>% mutate(case_ = as.numeric(case_))
  
  #create separate animal_id columns for each random effect
  for(i in 1:length(rand_terms_inla)){
    name <- as.name(paste0("id", i))
    dat[[name]] <- as.numeric(factor(dat$animal_id))
  }
  
  #construct the model formula with quadratic terms
  form <- as.formula(paste("case_ ~ -1 + ",  
                           #fixed effects
                           paste(c(params, quad_terms), collapse = " + "), "+",
                           #random intercept (strata)
                           "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T))) + ",
                           #random slopes
                           paste(rand_terms_inla, collapse = " + ")))
  
  #fit the model    
  fit <-  inla(form, family ="Poisson", data=dat,
               control.fixed = list(
                 mean = mean.beta,
                 prec = list(default = prec.beta)),
               control.compute = list(waic=TRUE, dic = TRUE, cpo = FALSE)) #
  
  return(fit)
}


get_inla_form <- function(dat){
  #convert response from categorical to numeric
  dat <- dat %>% mutate(case_ = as.numeric(case_))
  
  #create separate animal_id columns for each random effect
  for(i in 1:length(rand_terms_inla)){
    name <- as.name(paste0("id", i))
    dat[[name]] <- as.numeric(factor(dat$animal_id))
  }
  
  #construct the model formula with quadratic terms
  form <- as.formula(paste("case_ ~ -1 + ",  
                           #fixed effects
                           paste(c(params, quad_terms), collapse = " + "), "+",
                           #random intercept (strata)
                           "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T))) + ",
                           #random slopes
                           paste(rand_terms_inla, collapse = " + ")))
  return(form)
}


###################################################################
##
## 2. Backward stepwise model selection
##
###################################################################

formula_full <- get_inla_form(steps_scaled)

# Function for backward stepwise model selection based on WAIC
backward_selection_inla <- function(formula_full, dat, quad, threshold) {
  
  #convert response from categorical to numeric
  dat <- dat %>% mutate(case_ = as.numeric(case_))
  
  #create separate animal_id columns for each random effect
  for(i in 1:length(rand_terms_inla)){
    name <- as.name(paste0("id", i))
    dat[[name]] <- as.numeric(factor(dat$animal_id))
  }
  
  # Fit the full model
  model_full <- inla(formula_full, family ="Poisson", data=dat,
       control.fixed = list(
         mean = mean.beta,
         prec = list(default = prec.beta)),
       control.compute = list(waic=TRUE, dic = FALSE, cpo = FALSE))
  
  # Store the initial WAIC
  waic_full <- model_full$waic$waic
  waic_best <- waic_full
  
  #determine if function will perform selection on quadratic or linear terms
  if(quad == TRUE){
    fixed_effects <- quad_terms
  } else if(quad == FALSE){
    fixed_effects <- params
  }
  
  # Initialize the list of variables to be removed
  variables_to_remove <- character(0)
  
  while (length(variables_to_remove) < length(fixed_effects)) {
    
    terms <- terms(formula_full)
    variable_to_remove_best <- NULL
    
    fixed_effects <- setdiff(fixed_effects, variables_to_remove)
    
    # Try removing each fixed effect one by one
    waics <- numeric(length(fixed_effects))
  
    for (i in 1:length(fixed_effects)) {
      variable_to_remove <- fixed_effects[i]
    
      removal_index <- which(attr(terms, "term.labels") == variable_to_remove)
      
      new_terms <- drop.terms(terms, removal_index, keep.response = TRUE)
      
      reduced_formula <- reformulate(attr(new_terms, "term.labels"), response = "case_")
      
      # Fit the model without the selected fixed effect
      reduced_model <- inla(reduced_formula, family ="Poisson", data=dat,
                            control.fixed = list(
                              mean = mean.beta,
                              prec = list(default = prec.beta)),
                            control.compute = list(waic=TRUE, dic = FALSE, cpo = FALSE))
      
      # Calculate the WAIC for the reduced model
      waics[i] <- reduced_model$waic$waic
      
      # If the reduced model has a lower WAIC, update the best model and keep track of the variable to be removed
      if ((waics[i] - waic_best) < -threshold) {
        waic_best <- waics[i]
        variable_to_remove_best <- variable_to_remove
      }
    }
    
    # If none of the models result in a lower WAIC, stop the loop
    if ((min(waics) - waic_best) < -threshold) {
      #variable_to_remove_best <- NULL
   
    
    #if(is.null(variable_to_remove_best) == FALSE){
    # Remove the best variable and update the formula
    variables_to_remove <- c(variables_to_remove, variable_to_remove_best)
    
    removal_index_best <- which(attr(terms, "term.labels") == variable_to_remove_best)
    
    new_terms <- drop.terms(terms, removal_index_best, keep.response = TRUE)
    
    formula_full <- reformulate(attr(new_terms, "term.labels"), response = "case_")
    print(paste0("waic_best: ", min(waics)))
    print(paste0("variable_to_remove_best: ", variable_to_remove_best))
    print(paste0("variables_to_remove: ", variables_to_remove))
    } else{
      variable_to_remove_best <- "none"
      print(paste0("waic_best: ", min(waics)))
      print(paste0("variable_to_remove_best: ", variable_to_remove_best))
      print(paste0("variables_to_remove: ", variables_to_remove))
      break
    }
  } 
  
  # Fit the final selected model
  final_model <- inla(formula_full, family ="Poisson", data=dat,
                      control.fixed = list(
                        mean = mean.beta,
                        prec = list(default = prec.beta)),
                      control.compute = list(waic=TRUE, dic = FALSE, cpo = FALSE))
  
  # Return the final selected model
  return(final_model)
}

# Example usage:
# Define the full formula, data, family, and random/fixed effects
full_formula <- Loc ~ -1 + STAU1 + REST1 + Sohlenbrei + Breaks_Dis + 
  f(str_ID, model = "iid") + 
  f(ANIMAL_ID1, Sohlenbrei, values = 1:9, model = "iid", 
    hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(3, 0.05))))

data <- YourData  # Replace with your data
family <- "poisson"  # Specify the family
fixed_effects <- c("STAU1", "REST1", "Sohlenbrei", "Breaks_Dis")

# Perform backward model selection
final_selected_model <- backward_selection_inla(full_formula, data, family, fixed_effects)

# Display the final selected model summary
summary(final_selected_model)


# Function for backward stepwise model selection based on WAIC
backward_selection_inla <- function(formula_full, data, family, random_effects, fixed_effects, control = list()) {
  # Fit the full model
  model_full <- inla(formula = formula_full, data = data, family = family, control.predictor = control)
  
  # Store the initial WAIC
  waic_full <- inla.waic(model_full)
  waic_best <- waic_full
  
  # Initialize the list of variables to be removed
  variables_to_remove <- character(0)
  
  while (length(variables_to_remove) < length(fixed_effects)) {
    # Try removing each fixed effect one by one
    waics <- numeric(length(fixed_effects))
    
    for (i in 1:length(fixed_effects)) {
      variable_to_remove <- fixed_effects[i]
      reduced_formula <- as.formula(paste(formula_full, paste("- ", variable_to_remove), sep = ""))
      
      # Fit the model without the selected fixed effect
      reduced_model <- inla(reduced_formula, data = data, family = family, control.predictor = control)
      
      # Calculate the WAIC for the reduced model
      waics[i] <- inla.waic(reduced_model)$waic
      
      # If the reduced model has a lower WAIC, update the best model and keep track of the variable to be removed
      if (waics[i] < waic_best) {
        waic_best <- waics[i]
        variable_to_remove_best <- variable_to_remove
      }
    }
    
    # If none of the models result in a lower WAIC, stop the loop
    if (min(waics) >= waic_best) {
      break
    }
    
    # Remove the best variable and update the formula
    variables_to_remove <- c(variables_to_remove, variable_to_remove_best)
    formula_full <- as.formula(paste(formula_full, paste("- ", variable_to_remove_best), sep = ""))
  }
  
  # Fit the final selected model
  final_model <- inla(formula = formula_full, data = data, family = family, control.predictor = control)
  
  # Return the final selected model
  return(final_model)
}




###################################################################
##
## 2. Backward stepwise model selection fixed effects only
##
###################################################################

stepwise_dat <- steps_scaled %>% 
  select(case_, animal_id, step_id_, matches(params)) %>% 
  mutate(case_ = as.numeric(case_),
         animal_id = as.numeric(as.factor(animal_id)))
  

INLAstep(
  fam1 = "poisson",
  dataf = stepwise_dat,
  spde = NULL,
  in_stack = NULL,
  invariant = "0 + Intercept",
  direction = "backwards",
  y = "case_",
  y2 = "case_",
  include = NULL,
  powerl = 2,
  inter = 1,
  thresh = 2,
  num.threads = 3
)





















