##### Collar Download Data Formatting #####


#Author: Read Barbee

#Date: 2022-10-03

#Purpose: Filter for collared dates and combine and format files from retrieved collar downloads. Get everything into Movebank format.

library(tidyverse)
library(lubridate)
library(janitor)
library(collar)
library(naniar)


#Import historical downloads grouped in the 4 formats
hist_form1 <- list.files(
    path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Historic Downloads/format1", 
    pattern = "*.csv",
    full.names = TRUE) %>% 
  map_df(~read_csv(.))

hist_form2 <- list.files(
  path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Historic Downloads/format2", 
  pattern = "*.csv",
  full.names = TRUE) %>% 
  map_df(~read_csv(.))

hist_form3 <- list.files(
  path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Historic Downloads/format3", 
  pattern = "*.csv",
  full.names = TRUE) %>% 
  map_df(~read_csv(.))

hist_form4 <- list.files(
  path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Historic Downloads/format4", 
  pattern = "*.csv",
  full.names = TRUE) %>% 
  map_df(~read_csv(., col_types = list(NAV = col_character())))

hist_form5 <- list.files(
  path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Historic Downloads/format5/", 
  pattern = "*.csv",
  full.names = TRUE) %>% 
  map_df(~read_csv(.))



#######Standardize formatting of historical downloads#########
#Standard fields:
#deployment_id
#animal_id
#collar_id
#date_time_gmt
#latitude
#longitude
#altitude_m
#dop
#fix_type
#temp_c
#main_v
#back_v

f1_formatted <- hist_form1 %>% 
  clean_names() %>% 
  select(animal_id_2, collar_id, utc_date_time, latitude:fix_type, temp_c, main_v) %>% 
  rename(animal_id = animal_id_2,
         date_time_gmt = utc_date_time,
         altitude_m = height_m) %>% 
  mutate(date_time_gmt = mdy_hms(date_time_gmt),
         validated = ifelse(str_detect(fix_type, "val"), "Yes", "No"),
         fix_type = case_when(fix_type == "val. GPS-3D" | fix_type == "GPS-3D" ~ "3D", fix_type=="GPS-2D" ~ "2D", fix_type=="No Fix" ~ NA_character_)
         )



f2_formatted <- hist_form2 %>%  ## Continue formating fix_type collumns for remaining formats
  clean_names() %>%
  select(animal_id:date_time_gmt, latitude:back_v) %>% 
  rename(altitude_m = altitude,
         fix_type = fix_status) %>% 
  mutate(date_time_gmt = mdy(date_time_gmt))%>% 
  relocate(dop, .after = altitude_m) %>% 
  replace_with_na(replace = list(fix_type ="0")) 


f3_formatted <- hist_form3 %>% 
  clean_names() %>% 
  select(-date_time_local) %>% 
  rename(altitude_m=altitude,
         fix_type = fix_status) %>% 
  mutate(date_time_gmt= ymd_hms(date_time_gmt),
         validated = ifelse(str_detect(fix_type, "V"), "Yes", "No"),
         fix_type = case_when(fix_type == "3D-V Fix" | fix_type == "3D Fix" ~ "3D", fix_type=="2D Fix" ~ "2D", fix_type=="No Sats" ~ NA_character_)) %>% 
  relocate(dop, .after = altitude_m)


f4_formatted <- hist_form4 %>% 
  clean_names() %>% 
  mutate(date_time_gmt= mdy_hms(paste(gmt_date, gmt_time))) %>% 
  select(animal_id, collar_id, date_time_gmt, latitude:sats_used, temp, main_vol, bu_vol ) %>% 
  rename(altitude_m = height,
         fix_type = nav,
         temp_c=temp,
         main_v=main_vol,
         back_v= bu_vol) %>% 
  replace_with_na(replace = list(fix_type ="No"))
  relocate(c(sats_used, validated), .after = back_v)


f5_formatted <- hist_form5 %>% 
  clean_names() %>% 
  mutate(date_time_gmt= dmy_hms(paste(gmt_date_dmy, gmt_time))) %>% 
  select(animal_id, collar_id, date_time_gmt, latitude:sats_used, temp, main_vol, bu_vol ) %>% 
  rename(altitude_m = height,
         fix_type = nav,
         temp_c=temp,
         main_v=main_vol,
         back_v= bu_vol) %>% 
  replace_with_na(replace = list(fix_type ="No")) %>% 
  relocate(c(sats_used, validated), .after = back_v)


########## combine all formatted dfs into single df #############
hist_combined <- bind_rows(f1_formatted, f2_formatted, f3_formatted, f4_formatted, f5_formatted)

#add deplooyment_id column

hist_combined <- hist_combined %>% 
  mutate(deployment_id = paste0(animal_id,"_",collar_id)) %>% 
  select(deployment_id, everything())


#check deployments
unique(hist_combined$deployment_id)

write_csv(hist_combined, "/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Source Files/Formatted/hist_downloads_final_2022-10-05.csv")




