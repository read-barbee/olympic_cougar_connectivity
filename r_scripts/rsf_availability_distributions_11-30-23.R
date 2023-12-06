################################ RSF Availability Distributions #################################

#Date: 11/30/2023


library(tidyverse)
library(terra)
library(sf)

#import steps

project_crs <- 5070 #NAD83 Albers Equal Area Projection


#import covariate layers
cov_stack <- rast("/Users/tb201494/Desktop/1km_buffer/cov_stack_pred_raw_unscaled_30m_11-30-2023.tif")

## Study area boundary ##
op_poly <- st_read("data/Habitat_Covariates/study_area_polys/ocp_study_area_poly_wa_only_10-30-23.shp")
#op_poly_buffer <- st_buffer(op_poly, study_area_buffer_dist)

## Water Polygons ##
water_polys <- st_read("data/Habitat_Covariates/washington_water_polygons/op_water_polys_with_salt.shp") %>% st_transform(crs = project_crs)

## Freshwater only ##
water_polys_cropped <-  water_polys %>%
  filter(WB_PERIOD_ =="PER" | OBJECTID != 1) %>% #only include permanent water bodies with area > 100 m2
  filter(SHAPEAREA >= 100) %>%
  sf::st_crop(op_poly) #crop to study area
## Freshwater only (dissolved) ##
water_polys_mask <- water_polys_cropped %>%
  sf::st_union() %>% #dissolve polygons into single vector mask layer
  st_sf()
## Saltwater only ##
ocean <- water_polys %>% filter(OBJECTID == 1)

#make land mask
land_mask <- sf::st_difference(op_poly, water_polys_mask)


#import used locations
used_locs <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_10-02-2023.csv")

#format usde locations
used_locs2 <- used_locs %>% select(lon_utm, lat_utm) %>% mutate(case_=TRUE, .before= lon_utm)


#sample random points from land mask
rand_pts <- st_sample(land_mask, size = (nrow(used_locs2)*10))

#view subset of random points
mapview::mapview(rand_pts2[1:500])

#format random points
rand_pts2 <- rand_pts %>% 
  st_coordinates() %>% 
  as_tibble() %>% 
  rename(lon_utm = X,
         lat_utm = Y) %>% 
  mutate(case_=FALSE, .before = lon_utm)

#combine used and random locations
pts_combined <- bind_rows(used_locs2, rand_pts2) %>% st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070, remove = FALSE)

#extract covariates at all locations
vals <- terra::extract(cov_stack, pts_combined, bind = TRUE) %>% as_tibble()


vals_pivot <- vals %>% 
  select(!contains("usfs")) %>% 
  pivot_longer(c(elevation:dens_all_roads_annual), names_to = "cov", values_to= "value")


ggplot(vals_pivot, aes(x = value, y = ..density.., fill = as.factor(case_)))+ geom_histogram(position="identity", alpha=0.7) + 
  facet_wrap(vars(cov), scales = "free") +
  xlab("Used vs Available Distributions") + theme(axis.title.x=element_text(size=16))  
  
#write_csv(vals, "rsf_frame_v2_ten_to_one_12-01-23.csv")
  

rsf_dat <- read_csv("rsf_frame_v2_ten_to_one_12-01-23.csv")

scaled <- rsf_dat %>% mutate(across(elevation:dens_all_roads_annual, scale)) %>% 
  mutate(across(elevation:dens_all_roads_annual, as.numeric)) #%>% 
  #na.omit()
  


test_mod_1 <- glm(case_ ~ ndvi_annual + slope + aspect_northness + hii_annual + distance_water, family = binomial(logit), data = scaled)

#########################################################################
##
## 8. Correlation analysis
##
##########################################################################
#Correlation plot
plot_correlation(scaled %>% 
                   select(-c(case_, lat_utm, lon_utm)) %>% 
                   na.omit())

#select covariates of interest and dummy code factors
cont <-scaled %>% 
          select(-c(case_, lat_utm, lon_utm)) %>% 
          na.omit()


#create a correlation matrix
cor_matrix <- cor(cont, use="pairwise.complete.obs", method = "pearson")

# Set correlation threshold
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

rows_to_remove <- seq(from = 2, to = nrow(cor_dat), by = 2)
cor_dat_no_dupes <- cor_dat[-rows_to_remove, ]

#write_csv(cor_dat_no_dupes, "feature_selection/pairwise_cov_corr_rsf_2023_only_12-01-23.csv")

#the GGAlly package corrplot is more informative for numeric variables
ggcorr(steps_scaled %>% 
         select(elevation:dist_all_roads_annual), label = TRUE)





#########################################################################
##
## 8. Null model
##
##########################################################################
run_null <- function(dat){
  #convert response from categorical to numeric
  dat <- dat %>% mutate(case_ = as.numeric(case_))
  
  form <- case_ ~  1
  fit <- glm(form, family = binomial(logit), data = dat)
  
  return(fit)
}


null_fit <- run_null(scaled)

null_fit <- list(null_fit)

names(null_fit) <- "null"


#########################################################################
##
## 9. Linear Models
##
##########################################################################

cov_names <- scaled %>% select(-c(case_, lon_utm, lat_utm)) %>% names()


# ~ 11 min
uni_fits <- list()
system.time(for (i  in 1:length(cov_names)) { 
  
  cov <- cov_names[i]
  
  #construct the model formula with linear terms only
  form <- as.formula(paste0("case_ ~", cov))
  
  
  uni_fits[[i]] <- glm(form, family = binomial(logit), data = scaled)
  print(paste0(i, "/", length(cov_names)))
  
}) ; beep("fanfare")


names(uni_fits) <- cov_names


#########################################################################
##
## 11. Quadratic Models
##
##########################################################################

cov_names_quad <- cov_names

uni_fits_quad <- list()
system.time(for (i  in 1:length(cov_names_quad)){
  
  cov <- cov_names_quad[i]
  
  form <- as.formula(paste0("case_ ~", cov, " + I(",cov,"^2)"))
  
  uni_fits_quad[[i]] <-glm(form, family = binomial(logit), data = scaled)
  
  print(paste0(i, "/", length(cov_names_quad)))
  
}); beep("fanfare")


#rename qudaratic covariates with "2" at the end
quad_names <- vector()
for(i in 1:length(cov_names_quad)){
  old_name <- cov_names_quad[i]
  quad_names[i] <- paste0(old_name, "2")
}

names(uni_fits_quad) <- quad_names



#combine linear and quadratic fits
uni_fits_all <- c(uni_fits, uni_fits_quad, null_fit)

#########################################################################
##
## 12. Create model summary table
##
##########################################################################

mod_names_all <- names(uni_fits_all)



mod_table <- list()
for (i in 1:length(uni_fits_all)){
  name <- mod_names_all[i]
  
  if(str_detect(name, "2")){
    tmp1 <-  broom::tidy(uni_fits_all[[i]]) %>% 
      filter(str_detect(term, "2")==T)
  } else if (name=="null"){
    tmp1 <- broom::tidy(uni_fits_all[[i]])
  } else{
    tmp1 <-   broom::tidy(uni_fits_all[[i]]) %>% 
      filter(str_detect(term, "(Intercept)", negate = T))
  }
  
  tmp2 <- broom::glance(uni_fits_all[[i]])
  
  r2 <- pscl::pR2(uni_fits_all[[i]]) %>% enframe() %>% pivot_wider(names_from = name, values_from = value)
  
  mod_table[[i]] <- bind_cols(tmp1, tmp2, r2)
  print(i)
}

mod_table <- bind_rows(mod_table)


mod_table_renamed <- mod_table %>% 
  mutate(term = str_remove(.$term, coll("I("))) %>% 
  mutate(term = str_remove(.$term, coll(")"))) %>% 
  mutate(term = str_remove(.$term, coll("^"))) %>% 
  mutate(term = str_replace(.$term, coll("(Intercept"), coll("null")))

#write_csv(mod_table_renamed, "feature_selection/uni_rsf_12-01-23.csv")


