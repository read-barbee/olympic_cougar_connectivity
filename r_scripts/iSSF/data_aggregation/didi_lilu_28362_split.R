library(tidyverse)
library(janitor)
library(lubridate)

dl <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Kim_5-15-2023/historic_with_DOP/Didi_Lilu_28362_20210818132514.csv") %>% 
  clean_names()

deployments <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Metadata/From Teams/Formatted for R/collar_deployments_master_5-11-2023.csv") %>% 
  clean_names() %>% 
  filter(collar_id==28362)

dl_labeled <- dl %>%
  mutate(lmt_date = mdy(lmt_date)) %>%
  filter(lmt_date>= mdy("10/21/18") & lmt_date<= mdy("05/30/20")) %>% 
  mutate(animal_id= case_when(lmt_date <= mdy("02/13/19") ~ "Didi",
                              lmt_date >= mdy("12/02/19") ~ "Lilu",
                              .default = NA_character_)) %>% 
  select(animal_id, collar_id:lmt_time, latitude:fix_type)

dl2 <- dl_labeled %>% 
  mutate(collar_id = as.character(collar_id),
         utc_date_time = mdy_hms(paste0(utc_date," ",utc_time), tz="UTC"),
         lmt_date_time = with_tz(ymd_hms(paste0(lmt_date," ",lmt_time)), tzone="US/Pacific"),
         fix_type = case_when(fix_type %in% c("val. GPS-3D", "GPS-3D") ~ "3D",
                              fix_type == "GPS-2D" ~ "2D",
                              fix_type == "No Fix" ~ NA_character_)) %>% 
  select(animal_id, collar_id, utc_date_time, lmt_date_time, latitude:fix_type) %>% 
  rename(altitude= height_m)


didi <- dl2 %>% filter(animal_id =="Didi")

lilu <- dl2 %>% filter(animal_id =="Lilu")

write_csv(didi, "data/Location Data/Raw Data/Historic Downloads/pre-formatted/Didi_28362_historic_download_raw.csv")

write_csv(lilu, "data/Location Data/Raw Data/Historic Downloads/pre-formatted/Lilu_28362_historic_download_raw.csv")




