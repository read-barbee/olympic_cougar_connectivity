#### Count individual cougars from deployment list before 2018 ####

# Author: Read Barbee

# Date:2023-01-27 

# Purpose:


###########################################################################

#### Library ####
library(tidyverse)
library(lubridate)
library(janitor)


#### Data ####

deployments <- read_csv("data/Location Data/OCP_collar_deployments_master_list_2023-04-16.csv")


## Operations ###

unique <- deployments %>% 
  clean_names() %>% 
  mutate(end_date = ifelse(end_date =="Active", NA, end_date)) %>% 
  mutate(start_date = mdy(start_date),
         end_date = mdy(end_date)) %>% 
  distinct(name, .keep_all = TRUE) %>% 
  filter(start_date <= ymd("2018-01-01"))
