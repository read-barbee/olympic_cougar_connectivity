### Cougar Location Data Formatting Master

##Read Barbee

#Creation Date: 2022-10-04

library(tidyverse)
library(lubridate)
library(janitor)
library(collar)

######################### Step 1: Import raw data from all sources  ##################

lotek_web_raw <- read_csv("data/Location Data/Raw Data/Lotek/Lotek_complete_download_raw_2022-09-27.csv") %>% clean_names()

vec_web_raw <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Vectronics/Vectronic_complete_download_raw_2022-09-27.csv") %>% clean_names()

retrieved_collar_downloads <-
  list.files(
    path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Collar Downloads/Lotek/formatted", 
    pattern = "*.csv",
    full.names = TRUE) %>% 
  map_df(~read_csv(.)) %>% 
  clean_names() %>% 
  mutate(gmt_time = mdy_hm(gmt_time))
