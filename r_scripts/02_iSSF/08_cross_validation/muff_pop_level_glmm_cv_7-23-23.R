steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv")

#set all negative elevations to 0 and filter dispersal tracks for post dispersal event
steps <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad)) %>% 
  filter(is.na(dispersing) | dispersing==TRUE)

steps <- steps %>% mutate(land_cover_usfs_lumped = fct_collapse(land_cover_usfs, trees = c("trees", "tall_trees_shrubs", "gfh_tree_mix", "tree_shrub_mix",  "barren_tree_mix"),
                                                                shrubs = c("tall_shrubs", "shrubs", "gfh_shrub_mix", "barren_shrub_mix"),
                                                                gfh = c("gfh", "barren_gfh_mix"),
                                                                barren = c("barren_impervious", "snow_ice")),
                          land_use_usfs_lumped = fct_collapse(land_use_usfs, agriculture = c("agriculture", "rangeland_pasture")), .after = land_cover_usfs)

steps_scaled <- steps %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric)) %>% 
  select(-c(land_cover_usfs, land_use_usfs)) %>% 
  DataExplorer::dummify(select=c("land_use_usfs_lumped", "land_cover_usfs_lumped", "season", "hunting_season", "calving_season")) 




vars <- c("land_use_usfs_lumped_forest + (1|step_id_) + (0 + land_use_usfs_lumped_forest | animal_id)",
          "ndvi + I(ndvi^2) + (1|step_id_) + (0 + ndvi | animal_id)",
          "tree_cover_hansen + (1|step_id_) + (0 + tree_cover_hansen | animal_id)",
          "gpp + (1|step_id_) + (0 + gpp | animal_id)",
          "npp + (1|step_id_) + (0 + npp | animal_id)",
          "popdens_hii + I(popdens_hii^2) + (1|step_id_) + (0 + popdens_hii | animal_id)",
          "land_cover_usfs_lumped_water + (1|step_id_) + (0 + land_cover_usfs_lumped_water | animal_id)",
          "slope + I(slope^2) + (1|step_id_) + (0 + slope | animal_id)",
          "perc_tree_cover + I(perc_tree_cover^2) + (1|step_id_) + (0 + perc_tree_cover | animal_id)",
          "northing + (1|step_id_) + (0 + northing | animal_id)",
          "land_cover_usfs_lumped_gfh + (1|step_id_) + (0 + land_cover_usfs_lumped_gfh | animal_id)",
          "perc_nonveg + I(perc_nonveg^2) + (1|step_id_) + (0 + perc_nonveg | animal_id)",
          "perc_nontree_veg + I(perc_nontree_veg^2) + (1|step_id_) + (0 + perc_nontree_veg | animal_id)",
          "tri + I(tri^2) + (1|step_id_) + (0 + tri | animal_id)",
          "land_cover_usfs_lumped_shrubs + (1|step_id_) + (0 + land_cover_usfs_lumped_shrubs | animal_id)",
          "landuse_hii + I(landuse_hii^2) + (1|step_id_) + (0 + landuse_hii| animal_id)",
          "easting + I(easting^2) + (1|step_id_) + (0 + easting | animal_id)",
          "roads_hii + I(roads_hii^2) + (1|step_id_) + (0 + roads_hii | animal_id)",
          "dist_water + I(dist_water^2) + (1|step_id_) + (0 + dist_water | animal_id)",
          "rails_hii + (1|step_id_) + (0 + rails_hii | animal_id)",
          "infra_hii + (1|step_id_) + (0 + infra_hii | animal_id)")


all_comb <- do.call("c", lapply(seq_along(vars), function(i) combn(vars, i, FUN = list)))


forms <- list()

for (i in 1:length(all_comb)){
  var_i <- all_comb[[i]]
  forms[[i]] <- as.formula(paste("case_",  paste("-1", paste(var_i, collapse="+"), sep="+"), sep="~"))
}


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

library(sf)

all_sf <- steps_scaled %>% st_as_sf(coords=c("x1_", "y1_"), crs = 5070)
all_poly <- all_sf %>% sf::st_buffer(dist = 10000) %>%
  sf::st_union() %>%
  sf::st_convex_hull()
cov_stack_cropped <- terra::crop(cov_stack, all_poly)

ncde_4fold <- function (form, dat, cov_stack) {
  # Split out the used GPS points and the available points:
  case_t <- dplyr::filter(steps_scaled, case_ == T)
  case_f <- dplyr::filter(steps_scaled, case_ == F) %>% dplyr::mutate(group = NA)
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
    mod.fit <- glmmTMB(formula = forms[[i]], family = poisson, data = train, doFit=FALSE)
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
    breaks <- quantile(predictions, probs = 0:10/10, na.rm = T)
    breaks_j <- breaks + (seq_along(breaks) * .Machine$double.eps)
    train_ssf_vals_df <- within(train_ssf_vals_df, bins <- cut(predictions, breaks_j, include.lowest=TRUE))
    test_steps <- dplyr::filter(test, case_ == T)
    test_steps$step_id_ = mod.fit$model$model$`strata(step_id_)`[1]
    test_ssf_vals <- predict(mod.fit$model, newdata = test_steps,  type = "lp", allow.new.levels = TRUE)
    test_ssf_vals_df <- as.data.frame(test_ssf_vals)
    test_ssf_vals_df <- within(test_ssf_vals_df, bins <- cut(test_ssf_vals, breaks_j, include.lowest = TRUE))
    cor_cv <- cor.test(x = 1:10, y = as.numeric(table(test_ssf_vals_df$bins))/(nrow(test_steps)/10), method = "spearman", exact = FALSE)$estimate
    cv_mat[i] <- cor_cv
    cor_cv
  }
  return(cv_mat)
}



