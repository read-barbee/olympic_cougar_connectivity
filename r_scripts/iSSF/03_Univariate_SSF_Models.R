#### OCP iSSF Module_03: Univariate SSF Models ####

# Author: Read Barbee

# Date:2023-06-02 

# Purpose: Fit univariate SSF models to full dataset with random effect for individual to explore individual covariate effects

# Inputs:
#   Unscaled amt step data frame for fitting iSSFs
#
# Outputs:
#   •	univariate regression plots
#   •	univariate regression tables

#Steps
# •	Import amt step data




################################ Libraries #################################
library(tidyverse)
library(terra)
library(lubridate)
library(amt)
library(janitor)
library(DataExplorer)
library(GGally)


############################### Import Steps #################################

steps <- read_csv("data/Location_Data/Steps/6h_steps_unscaled_6-02-2023.csv")


################################ Check covariate distributions #################################

#data overview
introduce(steps)
plot_intro(steps)

#check for missing values
plot_missing(steps)


#continuous histograms
plot_histogram(steps %>% select(sl_, ta_, gpp:power_hii, -c(aspect_deg, aspect_rad)))

#categorical bar plots
plot_bar(steps %>% select(gpp:calving_season))












