#### All Subsets Model Selection ####

# Author: Read Barbee

# Date:2023-07-17 

# Purpose:

################################ Libraries #################################
library(tidyverse)
library(amt)
library(doParallel)
library(terra)
library(sf)

#########################################################################
##
## 1. Import and format step data
##
##########################################################################
#choose whether to use steps without imputation, imputation but no rerouting, or iimputation and rerouting

#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv")

#imputed and rerouted
#steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_imputed_7-12-2023.csv")


#set all negative elevations to 0 and filter dispersal tracks for post dispersal event
steps <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad)) %>% 
  filter(is.na(dispersing) | dispersing==TRUE)


#########################################################################
##
## 3. Reduce land cover and land use categories
##
##########################################################################
steps <- steps %>% mutate(land_cover_usfs_lumped = fct_collapse(land_cover_usfs, trees = c("trees", "tall_trees_shrubs", "gfh_tree_mix", "tree_shrub_mix",  "barren_tree_mix"),
                                                                shrubs = c("tall_shrubs", "shrubs", "gfh_shrub_mix", "barren_shrub_mix"),
                                                                gfh = c("gfh", "barren_gfh_mix"),
                                                                barren = c("barren_impervious", "snow_ice")),
                          land_use_usfs_lumped = fct_collapse(land_use_usfs, agriculture = c("agriculture", "rangeland_pasture")), .after = land_cover_usfs)

#plot_bar(steps %>% select(gpp:calving_season))


#########################################################################
##
## 4. Scale continuous covariates for model comparison and correlation analysis
##
##########################################################################

steps_scaled <- steps %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric)) %>% 
  select(-c(land_cover_usfs, land_use_usfs)) %>% 
  DataExplorer::dummify(select=c("land_use_usfs_lumped", "land_cover_usfs_lumped", "season", "hunting_season", "calving_season")) #%>% 
  #select(case_, animal_id, step_id_, gpp:calving_season_yes) %>% 
  #select(-c(sex:disp_qual))

#plot_histogram(steps_scaled %>% select(sl_, ta_, gpp:power_hii))


#########################################################################
##
## 4. Fit all model subsets (Dredge)
##
##########################################################################


vars <- c("evi",
          "dist_water",
          "perc_tree_cover",
          "tpi",
          "popdens_hii",
          "roads_hii",
          "infra_hii")



all_comb <- do.call("c", lapply(seq_along(vars), function(i) combn(vars, i, FUN = list)))


#n.clusters <- parallel::detectCores()-1

#my_cluster <- parallel::makeCluster(n.clusters, type = "FORK")

doParallel::registerDoParallel(cores = parallel::detectCores()-1) #cl=my_cluster, 

response <- "case_"
mods <- list()
#takes about 6 minutes with 19 cores. ~ 30 min for imputed dataset
system.time(
mods <- foreach(i = 1:length(all_comb), .packages = c("amt"), .verbose=TRUE)%dopar%{
  # aic_list <- list()
  # bic_list <- list()
  
    var_i <- all_comb[[i]]
    form <- as.formula(paste(response, paste(paste(var_i, collapse="+"), " strata(step_id_)", sep="+"), sep="~"))
    
    fit_issf(form, data=steps_scaled)
  }
)

#stopCluster(my_cluster)

names(mods) <- all_comb

tidy_list <- list()
mod_summ_issf <- list()

for (i in 1:length(mods)){
  Name <- names(mods)[[i]]
  tidy_list[[i]] <- broom::tidy(mods[[i]]$model) %>%
    mutate(name = Name, .before=term)
  
  mod_summ_issf[[i]] <- broom::glance(mods[[i]]$model) %>%
    select(r.squared:BIC, statistic.wald, p.value.wald) %>%
    mutate(name = Name, .before=r.squared)
}

tidy_list <- bind_rows(tidy_list)
mod_summ_issf <- bind_rows(mod_summ_issf)

issf_summary_tab_all_terms <- tidy_list %>% left_join(mod_summ_issf, by= join_by(name)) 
#%>% mutate(name= case_when(term!="cov" ~ paste0(name, "_", term), .default = name)) %>% select(-term)

#remove redundant linear terms from quadratic models
issf_summary_tab_mods_only <- issf_summary_tab_all_terms %>% distinct(name, .keep_all = TRUE) %>% select(name, r.squared:p.value.wald)

#write_csv(issf_summary_tab_all_terms , "issf_summary_all_subsets_all_terms_imp_7-18-23.csv")
#write_csv(issf_summary_tab_mods_only, "issf_summary_all_subsets_mods_only_imp_7-18-23.csv")


#########################################################################
##
## 4. Cross Validation
##
##########################################################################

fit_mod <- function(form, dat){
  fit <- fit_issf(form, data=dat)
}


forms <- list()

for (i in 1:length(all_comb)){
  var_i <- all_comb[[i]]
  forms[[i]] <- as.formula(paste(response, paste(paste(var_i, collapse="+"), " strata(step_id_)", sep="+"), sep="~"))
}


ncde_4fold <- function (form, dat, cov_stack) {
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
     mod_fit <- fit_issf(formula = form, data = train)
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




########### TEST ZONE ##############

al <- steps %>% filter(animal_id=="Al") %>% select(-c(land_cover_usfs, land_use_usfs)) %>% 
  DataExplorer::dummify(select=c("land_use_usfs_lumped", "land_cover_usfs_lumped", "season", "hunting_season", "calving_season"))



cov_stack <- terra::rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2.tif")

names(cov_stack) <- c("tree_cover_hansen",
                      "gpp",
                      "infra_hii",
                      "landuse_hii",
                      "land_cover_usfs",
                      "land_use_usfs",
                      "npp",
                      "popdens_hii",
                      "power_hii",
                      "precip",
                      "rails_hii",
                      "roads_hii",
                      "elevation",
                      "slope",
                      "aspect",
                      "tri",
                      "tpi",
                      "perc_tree_cover",
                      "perc_nontree_veg",
                      "perc_nonveg",
                      "ndvi",
                      "evi",
                      "dist_water")

# #remove categorical layers. not working yet.
# cov_stack <- cov_stack[[c(1:4, 7:23)]]
# 
# cov_stack_values <- values(cov_stack[[1:2]]) %>% data.frame() %>% 
#  cov_stack_values %>% na.omit()

#crop cov_stack to extent of individual data to save memory when projecting
indiv_sf <- al %>% st_as_sf(coords=c("x1_", "y1_"), crs = 5070)
indiv_poly <- indiv_sf %>% sf::st_buffer(dist = 10000) %>%
  sf::st_union() %>%
  sf::st_convex_hull()
cov_stack_cropped <- terra::crop(cov_stack, indiv_poly)
  
system.time(tmp <- ncde_4fold(forms[[2]], al, cov_stack_cropped))

ncde_4fold <- function (form, dat, cov_stack) {
  # Split out the used GPS points and the available points:
  case_t <- dplyr::filter(al, case_ == T)
  case_f <- dplyr::filter(al, case_ == F) %>% dplyr::mutate(group = NA)
  # Divide the datasets into 4 folds:
  n <- ceiling(nrow(case_t)/4)
  group <- c(replicate(n, sample(1:4, 4, replace = F))) #assign random group numbers to each row of true observations
  group <- group[c(1:nrow(case_t))]
  trn.tst <- cbind(case_t, group)
  trn.tst <- rbind(trn.tst, case_f) %>% dplyr::arrange(t1_) %>% fill(group)
  
  # Run the 4-fold CV:
  cv_mat <- matrix(NA, ncol = 4, nrow = 1)
  for (i in 1:4) {
    train <- dplyr::filter(trn.tst, group != i)
    test <- dplyr::filter(trn.tst, group == i)
    # Fit the model to the subsetted training dataset:
    mod.fit <- fit_issf(formula = form, data = train, model=TRUE)
    # Apply model to the landscape as a RSF-style raster map:
    mov <- train %>% dplyr::filter(case_ == T) %>% dplyr::select(sl_, log_sl_, cos_ta_)
    cov_stack_values <- terra::values(cov_stack, dataframe=TRUE)
    cov_stack_values$sl_ <- median(mov$sl_) # we need these items to predict from the model, so take medians.
    cov_stack_values$log_sl_ <- median(mov$log_sl_)
    cov_stack_values$cos_ta_ <- median(mov$cos_ta_, na.rm = T)
    cov_stack_values$step_id_ = mod.fit$model$model$`strata(step_id_)`[1]
    predictions <- terra::predict(mod.fit$model, newdata = cov_stack_values, type = "lp", allow.new.levels = TRUE)
    
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

predict(object=mod.fit$model, newdata=cov_stack_values, type="lp")


al_ranks <- list()
mod_names <- list()
system.time(for (i in 1:3){
  al_ranks[[i]] <- ncde_4fold(forms[[i]], al, cov_stack_cropped) %>% as.data.frame()
  mod_names[[i]] <- paste0(forms[[i]][3])
  
  #paste0(i, "/3")
})

df <- bind_rows(al_ranks) %>% mutate(name = mod_names, .before=V1)
