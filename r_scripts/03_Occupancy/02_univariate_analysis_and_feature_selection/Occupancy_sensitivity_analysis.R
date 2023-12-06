#### Occupancy sensitivity analysis ####

# Author: Read Barbee

# Date:2023-10-14

#Last updated: 2023-10-14

# Purpose:



################################ Libraries #################################
library(tidyverse)

#########################################################################
##
##  1. Import model selection tables
##
##########################################################################

daily <- read_csv("feature_selection/Occupancy/occu_det_model_fits_daily_10-12-23.csv")
weekly <- read_csv("feature_selection/Occupancy/occu_det_model_fits_weekly_10-14-23.csv")
two_weeks <- read_csv("feature_selection/Occupancy/occu_det_model_fits_2weeks_10-14-23.csv")
monthly <- read_csv("feature_selection/Occupancy/occu_det_model_fits_1month_10-14-23.csv")
two_months <- read_csv("feature_selection/Occupancy/occu_det_model_fits_2month_10-14-23.csv")
three_months <- read_csv("feature_selection/Occupancy/occu_det_model_fits_3month_10-14-23.csv")
four_months <- read_csv("feature_selection/Occupancy/occu_det_model_fits_4month_10-14-23.csv")
twelve_months <- read_csv("feature_selection/Occupancy/occu_det_model_fits_12month_10-14-23.csv")



data_list <- list(daily,
                  weekly,
                  two_weeks,
                  monthly,
                  two_months,
                  three_months,
                  four_months,
                  twelve_months)

names(data_list) <- c("daily",
                      "weekly",
                      "two_weeks",
                      "monthly",
                      "two_months",
                      "three_months",
                      "four_months",
                      "twelve_months")


plot_dat <- list()
for(i in 1:length(data_list)){
  name <- names(data_list)[i]
  
  plot_dat[[i]] <- data_list[[i]] %>% 
    group_by(model) %>% 
    slice(1) %>% 
    select(model, elpd, weight, Estimate_elpd_waic, SE_elpd_waic, Estimate_waic, SE_waic, looic) %>%
    mutate(interval = name, .before = model) %>% 
    pivot_longer(cols = -c(interval, model), names_to = "measure", values_to = "value")
  
}

plot_dat <- bind_rows(plot_dat)

plot_dat2 <- plot_dat %>% 
  mutate(days = case_when(interval == "daily" ~ 1,
                          interval =="weekly" ~ 7,
                          interval == "two_weeks" ~ 14,
                          interval == "monthly" ~ 30,
                          interval == "two_months" ~ 60,
                          interval == "three_months" ~ 90,
                          interval == "four_months" ~ 120,
                          interval == "twelve_months" ~ 365), .after = interval) %>% 
  mutate(model = case_when(model == "bait_days" ~ "bait",
                          model =="cam_days" ~ "cam",
                          model == "snare_days" ~ "snare",
                          .default = model)) %>% 
  pivot_wider(names_from = measure, values_from = value)

plot_dat2$interval <- factor(plot_dat2$interval, levels = names(data_list))

plot_dat2$model <- factor(plot_dat2$model, levels = c("null", "cam", "bait", "snare", "cam_bait", "cam_snare", "cam_bait_snare"))


#make plots

esquisse::esquisser()



#waic
ggplot(plot_dat2) +
  aes(x = model, y = Estimate_waic, colour = model) +
  geom_point(shape = "circle", size = 1.5) +
  geom_errorbar(aes(ymin = Estimate_waic - 2, ymax = Estimate_waic + 2),
                width = 0.25) + # Adjust the width of the error bars
  scale_color_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(interval), scales = "free") +
  scale_y_reverse()+
  theme(axis.text.x = element_text(size = 6)) +
  ggtitle("Detection model selection by survey interval ")

#model weight
ggplot(plot_dat2) +
  aes(x = model, y = weight, colour = model) +
  geom_point(shape = "circle", size = 1.5) +
  #geom_errorbar(aes(ymin = Estimate_waic - 2, ymax = Estimate_waic + 2),
                #width = 0.25) + # Adjust the width of the error bars
  scale_color_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(interval), scales = "free") +
  theme(axis.text.x = element_text(size = 6)) +
  ggtitle("Detection model weight by survey interval ")



# ggplot(plot_dat2) +
#   aes(x = interval, y = weight, fill = model) +
#   geom_col(position = position_dodge(width = 1)) +
#   scale_fill_hue(direction = 1) +
#   theme_minimal() + 
#   scale_x_discrete(limits = names(data_list)) 
# 
# 
# ggplot(plot_dat2) +
#   aes(x = interval, y = Estimate_waic, fill = model) +
#   geom_col(position = position_dodge(width = 1)) +
#   scale_fill_hue(direction = 1) +
#   theme_minimal() + 
#   scale_x_discrete(limits = names(data_list)) 
# 
# 
# 
# ggplot(plot_dat2) +
#   aes(x = interval, y = Estimate_waic, colour = model) +
#   geom_point(shape = "circle", size = 2, position = position_dodge(width = 1)) +
#   geom_errorbar(aes(ymin = Estimate_waic - SE_waic, ymax = Estimate_waic + SE_waic),
#                 width = 0.25,  # Adjust the width of the error bars
#                 position = position_dodge(width = 1)) +
#   scale_color_hue(direction = 1) +
#   theme_minimal() +
#   scale_x_discrete(limits = names(data_list)) 






#function form doesn't allow calling names of list elements:
# plot_form <- function(dat){
#   
#   name <- deparse(substitute(dat))
#   
#   new <- dat %>% 
#     group_by(model) %>% 
#     slice(1) %>% 
#     select(model, elpd, weight, Estimate_elpd_waic, SE_elpd_waic, Estimate_waic, SE_waic, looic) %>%
#     mutate(interval = name, .before = model) %>% 
#     pivot_longer(cols = -c(interval, model), names_to = "measure", values_to = "value") #%>% 
#   #pivot_wider(names_from = c(interval, measure), values_from = value)
#   
#   return(new)
#   
# }

#plot_dat <- map(data_list, plot_form)

