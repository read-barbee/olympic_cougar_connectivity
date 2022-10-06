##### Create Master Location Data Frame from Source Files #####


#Author: Read Barbee

#Creation Date: 2022-10-05
#Updated: 2022-10-05

#Purpose: Combine all source files into single master location dataframe and check for completeness

library(tidyverse)
library(lubridate)
library(janitor)


################## Load Source Files ####################

lotek_web <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Source Files/Formatted/lotek_web_final_2022-10-04.csv")

vec_web <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Source Files/Formatted/vec_web_final_2022-10-04.csv")

retrieved <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Source Files/Formatted/retrieved_collars_final_2022-10-05.csv")

historic <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Source Files/Formatted/hist_downloads_final_2022-10-05.csv")

all_locations <- bind_rows(lotek_web, vec_web, retrieved, historic)

all_locations <-all_locations %>% 
  mutate(deployment_id= str_replace_all(deployment_id, "Tswift_38019", "TSwift_38019")) %>% mutate(deployment_id = str_replace_all(deployment_id, " ", "")) %>% 
  mutate(animal_id = str_replace_all(animal_id, " ", ""))

unique(all_locations$deployment_id)

deployments <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/OCP_Cougar_Deployments_9-30-22.csv") %>% 
  clean_names() %>% 
  mutate(deployment_id = paste0(name,"_",collar_id))

setdiff(all_locations$deployment_id, deployments$deployment_id)

### Missing 17 deployments
