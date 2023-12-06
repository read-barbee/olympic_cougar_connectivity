#### Occupancy global model validation ####

# Author: Read Barbee

# Date:2023-10-15

#Last updated: 2023-10-15

# Purpose:



#########################################################################
##
##  4. Fit top model from SSA
##
##########################################################################

top_ssa <- stan_occu(~ scale(eff) + scale(bait) + scale(snare) ~ tree_cover_hansen + ndvi + northing + tri + perc_tree_cover + perc_nonveg + roads_hii + elevation + popdens_hii + landuse_hii + I(elevation^2) + I(tree_cover_hansen^2) + I(ndvi^2) + I(perc_tree_cover^2) + I(roads_hii^2) + (1|cell_id_year), data=umf_stack, chains=3, iter=1500) ; beep("fanfare")


#extract fixed effect beta values from fitted_model
betas <- inla_fit_quad3$summary.fixed$mean

#reff <- ranef(global_fit)


#########################################################################
##
##  4. with ocp data only
##
##########################################################################
occ_dat_complete_ocp <- occ_dat_scaled %>% 
  mutate(comp = complete_cases) %>% 
  unite("cell_id_year", cell_id, year, remove = FALSE) %>% 
  mutate(cell_id_year = as.factor(cell_id_year)) %>% 
  filter(comp==TRUE) %>% 
  filter(grid_id != "ONP")


#construct umf stack

nsite <- nrow(occ_dat_complete_ocp)
y <- occ_dat_complete_ocp %>% 
  dplyr::select(contains("detections")) %>% 
  as.matrix()

# Number of surveys detected per site
summary(rowSums(y, na.rm = TRUE))	

eff <- occ_dat_complete_ocp %>%
  dplyr::select(contains("cam")) %>%
  as.matrix()

bait <- occ_dat_complete_ocp %>%
  dplyr::select(contains("bait")) %>%
  as.matrix()

snare <- occ_dat_complete_ocp %>%
  dplyr::select(contains("snare")) %>%
  as.matrix()

covs_scaled <- occ_dat_complete_ocp %>% select(cell_id_year, tree_cover_hansen:easting)

umf_stack_ocp <- unmarkedFrameOccu(y = y, 
                               siteCovs = covs_scaled,
                               obsCovs = list(eff = eff,
                                              bait = bait,
                                              snare = snare)) #, survey = surveyID))


top_ssa_ocp <- stan_occu(~ scale(eff) ~ tree_cover_hansen + ndvi + northing + tri + perc_tree_cover + perc_nonveg + roads_hii + elevation + popdens_hii + landuse_hii + I(elevation^2) + I(tree_cover_hansen^2) + I(ndvi^2) + I(perc_tree_cover^2) + I(roads_hii^2) + (1|cell_id_year), data=umf_stack_ocp, chains=3, iter=1500) ; beep("fanfare")

#extract fixed effect beta values from fitted_model
betas <- summary(top_ssa_ocp, submodel = "state")$mean

cov_stack_capped <- terra::rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2_capped_30m.tif")

cov_stack <- terra::rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2_300m.tif")

#remove categorical layers from cov stack
cov_stack2 <- subset(cov_stack, c("land_cover_usfs", "land_use_usfs"), negate = TRUE)

cov_stack_capped <- list()
for (i in 1:length(names(cov_stack2))){
  name <- names(cov_stack2)[[i]]
  
  upper <- max(occ_dat[[name]], na.rm = T)
  lower <- min(occ_dat[[name]], na.rm = T)
  
  cov_stack_capped[[i]] <- clamp(cov_stack2[[i]], lower = lower, upper = upper)
  
  print(paste0(i, "/", length(names(cov_stack2))))
  
}

#add names to clamped layers
names(cov_stack_capped) <- names(cov_stack2)

cov_stack_capped <- rast(cov_stack_capped)

#cov_stack_capped <- terra::aggregate(cov_stack_capped, fact = 10)

cov_stack_pred <- terra::scale(cov_stack_capped)

#ext <- terra::draw()

ext <- ext(-2153256.05769231, -1938727.21153846, 2837586.92307692, 3152175.86538462)

cov_stack_pred <- terra::crop(cov_stack_pred, ext)

#manual predictions
preds <- exp((betas[1] +
                betas[2] * cov_stack_pred$tree_cover_hansen) + 
               (betas[3] * cov_stack_pred$ndvi) + 
               #(betas[3] * cov_stack_pred$tpi) + 
               (betas[4] * cov_stack_pred$northing) + 
               (betas[5] * cov_stack_pred$tri) + 
               (betas[6] * cov_stack_pred$perc_tree_cover) +
               (betas[7] * cov_stack_pred$perc_nonveg) +
               (betas[8] * cov_stack_pred$roads_hii) + 
               (betas[9] * cov_stack_pred$elevation) +
               (betas[10] * cov_stack_pred$popdens_hii) +
               (betas[11] * cov_stack_pred$landuse_hii) + 
               #(betas[11] * cov_stack_pred$rails_hii) +
               #(betas[12] * cov_stack_pred$infra_hii))/
               #(betas[13] * cov_stack_pred$easting) +
               (betas[12] * cov_stack_pred$elevation^2) +
               (betas[13] * cov_stack_pred$tree_cover_hansen^2) +
               (betas[14] * cov_stack_pred$ndvi^2) + 
               #(betas[16] * cov_stack_pred$tpi^2) +
               #(betas[17] * cov_stack_pred$northing^2) +
               #(betas[14] * cov_stack_pred$tri^2) +
               (betas[15] * cov_stack_pred$perc_tree_cover^2) +
               (betas[16] * cov_stack_pred$roads_hii^2))/
  # (betas[17] * cov_stack_pred$popdens_hii^2))/
  #(betas[22] * cov_stack_pred$landuse_hii^2))/
  (1 + exp((betas[1] +
              betas[2] * cov_stack_pred$tree_cover_hansen) + 
             (betas[3] * cov_stack_pred$ndvi) + 
             #(betas[3] * cov_stack_pred$tpi) + 
             (betas[4] * cov_stack_pred$northing) + 
             (betas[5] * cov_stack_pred$tri) + 
             (betas[6] * cov_stack_pred$perc_tree_cover) +
             (betas[7] * cov_stack_pred$perc_nonveg) +
             (betas[8] * cov_stack_pred$roads_hii) + 
             (betas[9] * cov_stack_pred$elevation) +
             (betas[10] * cov_stack_pred$popdens_hii) +
             (betas[11] * cov_stack_pred$landuse_hii) + 
             #(betas[11] * cov_stack_pred$rails_hii) +
             #(betas[12] * cov_stack_pred$infra_hii))/
             #(betas[13] * cov_stack_pred$easting) +
             (betas[12] * cov_stack_pred$elevation^2) +
             (betas[13] * cov_stack_pred$tree_cover_hansen^2) +
             (betas[14] * cov_stack_pred$ndvi^2) + 
             #(betas[16] * cov_stack_pred$tpi^2) +
             #(betas[17] * cov_stack_pred$northing^2) +
             #(betas[14] * cov_stack_pred$tri^2) +
             (betas[15] * cov_stack_pred$perc_tree_cover^2) +
             (betas[16] * cov_stack_pred$roads_hii^2)))
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

# residents <- locs %>% 
#   filter(dispersal_status =="resident") %>% 
#   select(-c(disp_date_nsd:dispersing))

eval_res <- locs %>% st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070)

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

# dispersers <- locs %>% 
#   filter(dispersal_status =="disperser") %>% 
#   filter(dispersing == TRUE)
# 
# eval_disp <- dispersers %>% st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070)
# 
# #extract binned predicted probability of use values for each used location
# loc_preds_disp<- terra::extract(binned, eval_disp)[,2] %>% as.numeric() %>% na.omit() 
# 
# #plot proportion of test points in each bin. 
# ggplot() +
#   geom_bar(aes(x=loc_preds_disp, y = after_stat(prop))) 
# 
# n_bins <- length(unique(loc_preds_disp))
# 
# #calculate pearson correlation between the bin number and the proportion of used locations in each bin
# cor_cv_res_disp <- cor.test(x = 1:n_bins, y = as.numeric(table(loc_preds_disp))/(length(loc_preds_disp)/n_bins), method = "spearman", exact = FALSE)$estimate
# 
# cor_cv_res_disp

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



plot_dat <- summary(top_ssa, submodel = "state") %>% rownames_to_column("variable") %>% 
  filter(!(variable %in% c("sigma [1|cell_id_year]","(Intercept)","elevation", "ndvi", "perc_tree_cover", "roads_hii", "tree_cover_hansen")))

ggplot(plot_dat) +
  aes(x = mean, y = variable, colour = as.factor(variable)) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_hue(direction = 1) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`, height = 0.25)) + # Adjust the width of the error bars
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add a reference line at 0
  labs(title = "Coefficient Posterior Means with 95% Credible Intervals")
