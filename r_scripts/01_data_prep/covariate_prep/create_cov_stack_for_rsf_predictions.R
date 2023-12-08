#### Create cov stack for RSF prediction ####

# Author: Read Barbee

# Date:2023-11-20

# Purpose:

################################ libraries #################################
library(tidyverse)
library(terra)


########################################################################
##
## 1. Import and format covariate stacks and step data
##
########################################################################

#most recent annual stack
cov_stack_2023 <- rast("/Users/tb201494/Desktop/annual_cov_stacks_1km_buffer/cov_stack_2023.tif")

#static stack
cov_stack_static <- rast("/Users/tb201494/Desktop/1km_buffer/static_stack_1km_buffer_11-29-23.tif")

#steps (no imputation)
steps <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location_Data/Steps/2h_steps_unscaled_no_imp_annual_cov_11-30-2023.csv")

########################################################################
##
## 2. Rename layers to match location data
##
########################################################################

cov_stack_annual_names_new <- substr(names(cov_stack_2023), 1, nchar(names(cov_stack_2023)) - 5)

names(cov_stack_2023) <- cov_stack_annual_names_new


########################################################################
##
## 3. Combine static and annual stacks
##
########################################################################
cov_stack <- c(cov_stack_static, cov_stack_2023)

cov_stack <- terra::subset(cov_stack, "aspect", negate = TRUE)

########################################################################
##
## 4. Cap raster values based on ranges found in step data
##
########################################################################

cov_stack_capped <- list()
for (i in 1:length(names(cov_stack))){
  name <- names(cov_stack)[[i]]
  
  upper <- max(steps[[name]], na.rm = T)
  lower <- min(steps[[name]], na.rm = T)
  
  cov_stack_capped[[i]] <- clamp(cov_stack[[i]], lower = lower, upper = upper)
  
  print(paste0(i, "/", length(names(cov_stack))))
  
}

cov_stack_capped <- rast(cov_stack_capped)


########################################################################
##
## 5. Calculate quadratic layers
##
########################################################################

# cov_stack_capped$`I(ndvi_annual^2)` <- cov_stack_capped[["ndvi_annual"]]^2
# cov_stack_capped$`I(perc_tree_cov_annual^2)` <- cov_stack_capped[["perc_tree_cov_annual"]]^2
# cov_stack_capped$`I(popdens_hii_annuall^2)` <- cov_stack_capped[["popdens_hii_annual"]]^2


########################################################################
##
## 6. Scale values and export
##
########################################################################

#remove raw aspect layer before scaling

cov_stack_scaled <- scale(cov_stack_capped)


writeRaster(cov_stack_scaled, "/Users/tb201494/Desktop/annual_cov_stacks_1km_buffer/cov_stack_pred_11-30-2023.tif", overwrite=TRUE)


