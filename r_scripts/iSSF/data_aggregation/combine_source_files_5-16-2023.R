#### Combine all location source files ####

# Author: Read Barbee

# Date:2023-05-16 

# Purpose: Combine all location source files:
#1. web source
#2. hist_source
#3. collar_source

#total deployments: 171


###############################################################################

#### Library ####
library(tidyverse)
library(lubridate)
library(janitor)


################################ Import Source Files #################################
web_source <- read_csv("data/Location Data/Source Files/web_source_5-16-2023.csv", col_types = list(fix_type = col_character(),
                                                                                                    collar_id = col_character())) %>% 
  mutate(source = "web")

hist_source <- read_csv("data/Location Data/Source Files/hist_source_5-16-2023.csv",col_types = list(fix_type = col_character(),
                                                                                                     collar_id = col_character())) %>% 
  mutate(source = "historic")

collar_source <- read_csv("data/Location Data/Source Files/collar_source_5-16-2023.csv",col_types = list(fix_type = col_character(),
                                                                                                         collar_id = col_character())) %>% 
  mutate(source = "collar_download")

deployments_master <- read_csv("data/Location Data/Metadata/From Teams/Formatted for R/collar_deployments_master_5-11-2023.csv") %>%
  mutate(deployment_id = paste0(name,"_",collar_id), .before=name) %>% distinct(deployment_id) %>% pull()

################################ Combine into single data frame #################################
locs_all <- bind_rows(web_source,
                      hist_source,
                      collar_source)



#check for duplicate points
get_dupes(locs_all, deployment_id, date_time_local, latitude, longitude)
get_dupes(locs_all, deployment_id, date_time_local)


#make sure the location file contains all the deployments in the master deployment list
locs_all_deps <- locs_all %>% distinct(deployment_id) %>% pull()

setdiff(deployments_master, locs_all_deps)
setdiff(locs_all_deps, deployments_master)

#total deployments: 171

locs_all <- locs_all %>% 
  mutate(date_time_local = with_tz(date_time_utc, tzone="US/Pacific"))

#with_tz doesn't print to csv unless coerced to a character.
locs_all <- locs_all %>%
  mutate(date_time_utc = as.character(date_time_utc),
         date_time_local = as.character(date_time_local))

write_csv(locs_all, "data/Location Data/Source Files/locations_master/gps_locs_master_5-16-2023.csv")


#create a file detailing the source file for each deployment
deployment_sources <- locs_all %>% distinct(deployment_id, source)

#write_csv(deployment_sources, "data/Location Data/Metadata/deployment_sources_5-16-2023.csv")


###############################################################################  