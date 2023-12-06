library(tidyverse)

cov_stack <- rast("/Users/tb201494/Desktop/1km_buffer/cov_stack_pred_raw_unscaled_30m_11-30-2023.tif")

sims_rerouted <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location_Data/imputed_paths/ctmm_sim_imputed_paths_rerouted_12-04-2023.csv")


imp_rerouted_sf <- sims_rerouted_full_df %>% st_as_sf(coords=c("x_rerouted", "y_rerouted"), crs = 5070)


#extract covariates at all locations
vals <- terra::extract(cov_stack, imp_rerouted_sf, bind = TRUE) %>% as_tibble()

vals_pivot <- vals %>% 
  select(!contains("usfs")) %>% 
  pivot_longer(c(elevation:dens_all_roads_annual), names_to = "cov", values_to= "value")


ggplot(vals_pivot, aes(x = value, y = ..density.., fill = as.factor(imp_status)))+ geom_histogram(position="identity", alpha=0.7) + 
  facet_wrap(vars(cov), scales = "free") +
  xlab("Observed vs Imputed Distributions") + theme(axis.title.x=element_text(size=16)) 


 ggplot(vals_pivot, aes(x = value, fill = as.factor(imp_status)))+ geom_histogram(position="identity", alpha=0.7) + 
  facet_wrap(vars(cov), scales = "free") +
  xlab("Observed vs Imputed Distributions") + theme(axis.title.x=element_text(size=16)) 
 
 
 
 
 #monthly ndvi values
 
 
 monthly_ndvi <- read_csv("monthly_mean_ndvi_values_2010_2022.csv")
 
 monthly_ndvi2 <- read_csv("/Users/tb201494/Desktop/annual_monthly_ndvi_12-04-23.csv")

 monthly_ndvi <- monthly_ndvi %>% mutate(month = as.factor(month))
 
 ggplot(monthly_ndvi) +
   aes(x = as.factor(month), y = ndvi) +
   geom_bar(colour = "#112446") +
   theme_minimal() +
   facet_wrap(vars(year))
 
 ggplot(monthly_ndvi2) +
  aes(x = month, y = mean) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
   geom_errorbar(aes(ymin = (monthly_ndvi2$mean - monthly_ndvi2$sd), ymax = (monthly_ndvi2$mean + monthly_ndvi2$sd)))+
  theme_minimal() +
  facet_wrap(vars(year))
 
 
