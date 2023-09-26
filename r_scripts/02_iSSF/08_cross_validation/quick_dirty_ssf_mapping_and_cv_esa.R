library(tidyverse)
library(terra)
library(sf)
library(survival)

#########################################################################
##
## 1. Import and format step data
##
##########################################################################
#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv")


#set all negative elevations to 0 and filter dispersal tracks for post dispersal event
steps <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad)) %>% 
  filter(is.na(dispersing) | dispersing==TRUE)

#collapse landcover categories
steps <- steps %>% mutate(land_cover_usfs_lumped = fct_collapse(land_cover_usfs, trees = c("trees", "tall_trees_shrubs", "gfh_tree_mix", "tree_shrub_mix",  "barren_tree_mix"),
                                                                shrubs = c("tall_shrubs", "shrubs", "gfh_shrub_mix", "barren_shrub_mix"),
                                                                gfh = c("gfh", "barren_gfh_mix"),
                                                                barren = c("barren_impervious", "snow_ice")),
                          land_use_usfs_lumped = fct_collapse(land_use_usfs, agriculture = c("agriculture", "rangeland_pasture")), .after = land_cover_usfs)

#scale and dummify covariates
steps_scaled <- steps %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric)) %>% 
  select(-c(land_cover_usfs, land_use_usfs)) %>% 
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

cov_stack<- terra::aggregate(cov_stack, fact=10, fun="mean", cores=5)

cov_stack$aspect_rad <- (pi*cov_stack$aspect)/180
cov_stack$northing <- cos(cov_stack$aspect_rad)
cov_stack$easting <- sin(cov_stack$aspect_rad)

cov_stack_values <- terra::values(cov_stack2, dataframe=TRUE)
cov_stack$sl_ <- median(steps_scaled$sl_) # we need these items to predict from the model, so take medians.
cov_stack$log_sl_ <- median(steps_scaled$log_sl_)
cov_stack$cos_ta_ <- median(steps_scaled$cos_ta_, na.rm = T)
cov_stack$step_id_ = mod$model$`strata(step_id_)`[1]

#takes about 3 min at 60m resolution
system.time(predictions <- terra::predict(mod$model, newdata = cov_stack2, type = "lp", allow.new.levels = TRUE))

predictions <- exp(predictions)
x.min <- quantile(predictions, 0.025, na.rm = T) # omit the extremes
x.max <- quantile(predictions, 0.975, na.rm = T)
predictions[predictions > x.max] <- x.max
predictions[predictions < x.min] <- x.min
x.range <- x.max - x.min
predictions. <- (predictions - x.min)/x.range
vals <- cov_stack[[1]] %>% terra::setValues(predictions.)


data <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_7-11-2023.csv")



#crop cov_stack to extent of individual data to save memory when projecting
# all_sf <- data %>% st_as_sf(coords=c("lon_utm", "lat_utm"), crs = 5070)
# all_poly <- all_sf %>% sf::st_buffer(dist = 100000) %>% #300km buffer
#   sf::st_union() %>%
#   sf::st_convex_hull()
# cov_stack_cropped <- terra::crop(cov_stack, all_poly)



mod <- amt::fit_issf(case_ ~
                   tree_cover_hansen +
                   gpp +
                   northing +
                   easting +
                   perc_nontree_veg +
                   tpi +
                   popdens_hii + #I(popdens_hii^2) +
                   landuse_hii +
                   ndvi + #I(ndvi^2) +
                   infra_hii +
                   sl_ +
                   log_sl_ +
                   cos_ta_ +
                   strata(step_id_), data=steps_scaled, model=TRUE)

mod2 <- amt::clogit(case_ ~
                       tree_cover_hansen +
                       #gpp +
                       northing +
                       #easting +
                       perc_nontree_veg +
                       #tpi +
                       #popdens_hii + #I(popdens_hii^2) +
                       landuse_hii +
                       ndvi + #I(ndvi^2) +
                       infra_hii +
                       sl_ +
                       log_sl_ +
                       cos_ta_ +
                       strata(step_id_), data=steps_scaled, cluster=animal_id, method="approximate",model=TRUE)


cov_stack$case_ <- 1
cov_stack$step_id_ = mod2$model$`strata(step_id_)`[1]



cov_stack2 <- terra::scale(cov_stack) #scaling before multiplying actually makes the values make sense
#does the same thing as predict. Just wanted to double check
test <- cov_stack2$tree_cover_hansen *mod$model$coefficients["tree_cover_hansen"] +
   cov_stack2$gpp *mod$model$coefficients["gpp"] +
  cov_stack2$northing *mod$model$coefficients["northing"] +
  cov_stack2$easting *mod$model$coefficients["easting"] +
  cov_stack2$perc_nontree_veg *mod$model$coefficients["perc_nontree_veg"] +
  cov_stack2$tpi *mod$model$coefficients["tpi"] +
  cov_stack2$popdens_hii *mod$model$coefficients["popdens_hii"] +
  cov_stack2$landuse_hii *mod$model$coefficients["landuse_hii"] +
  (cov_stack2$ndvi) *mod$model$coefficients["ndvi"] + #*0.0001
  cov_stack2$infra_hii *mod$model$coefficients["infra_hii"]


test_exp <- exp(test) / (1+ exp(test))

exp_preds <- terra::values(test_exp)

breaks <- quantile(exp_preds, probs = 0:10/10, na.rm = T)

bins <- cut(exp_preds, breaks, include.lowest=TRUE)

binned <- terra::classify(test_exp, rcl=as.vector(breaks))


#viridis magma
ggplot()+
  tidyterra::geom_spatraster(data=binned, mapping=aes()) +
  #geom_sf(data=red_deer_used, aes(geometry=geometry))+
  #coord_sf(crs=4269, xlim = c(-116.4, -115.3), ylim=c(51.1, 51.9))+
  scale_fill_manual(values = rev(viridis::magma(10)), na.value = NA, guide = guide_legend(reverse = TRUE), na.translate=FALSE)  #rev() to reverse pallete
  
mapview::mapview(raster::raster(binned))




#### mod2
test2 <- cov_stack2$tree_cover_hansen *mod2$coefficients["tree_cover_hansen"] +
  #cov_stack2$gpp *mod2$model$coefficients["gpp"] +
  cov_stack2$northing *mod2$coefficients["northing"] +
  #cov_stack2$easting *mod2$model$coefficients["easting"] +
  cov_stack2$perc_nontree_veg *mod2$coefficients["perc_nontree_veg"] +
  #cov_stack2$tpi *mod2$model$coefficients["tpi"] +
  #cov_stack2$popdens_hii *mod2$model$coefficients["popdens_hii"] +
  cov_stack2$landuse_hii *mod2$coefficients["landuse_hii"] +
  (cov_stack2$ndvi) *mod2$coefficients["ndvi"] + #*0.0001
  cov_stack2$infra_hii *mod2$coefficients["infra_hii"]


test_exp2 <- exp(test2) / (1+ exp(test2))

exp_preds2 <- terra::values(test_exp2)

breaks2 <- quantile(exp_preds2, probs = 0:10/10, na.rm = T)

binned2 <- terra::classify(test_exp2, rcl=as.vector(breaks2))


#viridis magma
ggplot()+
  tidyterra::geom_spatraster(data=binned2, mapping=aes()) +
  #geom_sf(data=red_deer_used, aes(geometry=geometry))+
  #coord_sf(crs=4269, xlim = c(-116.4, -115.3), ylim=c(51.1, 51.9))+
  scale_fill_manual(values = rev(viridis::magma(10)), na.value = NA, guide = guide_legend(reverse = TRUE), na.translate=FALSE)  #rev() to reverse pallete

mapview::mapview(raster::raster(binned2))


cv <- ncde_4fold(form = as.formula(case_ ~
                                     tree_cover_hansen +
                                     gpp +
                                     northing +
                                     easting +
                                     perc_nontree_veg +
                                     tpi +
                                     popdens_hii + #I(popdens_hii^2) +
                                     landuse_hii +
                                     ndvi + #I(ndvi^2) +
                                     infra_hii +
                                     sl_ +
                                     log_sl_ +
                                     cos_ta_ +
                                     strata(step_id_)),
                 dat=steps_scaled,
                 cov_stack = cov_stack2)


#### CV
ncde_4fold <- function (form, dat, cov_stack) {
  # Split out the used GPS points and the available points:
  case_t <- dplyr::filter(dat, case_ == T)
  case_f <- dplyr::filter(dat, case_ == F) %>% dplyr::mutate(group = NA)
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

