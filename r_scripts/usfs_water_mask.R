library(terra)

cov_stack <- rast("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Habitat_Covariates/puma_cov_stack_v1/cov_stack1.tif")

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


land_cover <- cov_stack[[5]]

m <- c(0, 14, 1,
       14, 15, 0,
       15,15, NA)

rclmat <- matrix(m, ncol=3, byrow=TRUE)

#lc_reclass <- ifelse(land_cover==15, 0, 1)

lc_reclass <- classify(land_cover, rclmat, include.lowest=TRUE)

#plot(lc_reclass)

writeRaster(lc_reclass, "data/Habitat_Covariates/puma_cov_stack_v1/usfs_water_mask_6-23-23.tif" )
