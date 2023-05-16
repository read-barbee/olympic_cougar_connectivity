#### Combine all location source files ####

# Author: Read Barbee

# Date:2023-05-16 

# Purpose: Combine all location source files:
#1. web source
#2. hist_source
#3. collar_source


###############################################################################

#### Library ####
library(tidyverse)
library(lubridate)
library(janitor)


################################ Import Source Files #################################
web_source <- read_csv("data/Location Data/Source Files/web_source_5-16-2023.csv", col_types = list(fix_type = col_character(),
                                                                                                    collar_id = col_character()))

hist_source <- read_csv("data/Location Data/Source Files/hist_source_5-16-2023.csv",col_types = list(fix_type = col_character(),
                                                                                                     collar_id = col_character()))

collar_source <- read_csv("data/Location Data/Source Files/collar_source_5-16-2023.csv",col_types = list(fix_type = col_character(),
                                                                                                         collar_id = col_character()))

deployments_master <- read_csv("data/Location Data/Metadata/From Teams/Formatted for R/collar_deployments_master_5-11-2023.csv") %>%
  mutate(deployment_id = paste0(name,"_",collar_id), .before=name) %>% distinct(deployment_id) %>% pull()

################################ Combine into single data frame #################################
locs_all <- bind_rows(web_source,
                      hist_source,
                      collar_source)



#check for duplicate points
get_dupes(locs_all, deployment_id, date_time_local, latitude, longitude)

#duplicates in time come from loss of seconds portion of time stamp in lilu 85724 and two missing local timestamps in apollo_28360
get_dupes(locs_all, deployment_id, date_time_local)

#seven locations are missing timestamps
which(is.na(locs_all$date_time_local))

#check deployment list
locs_all_deps <- locs_all %>% distinct(deployment_id) %>% pull()

setdiff(deployments_master, locs_all_deps)
setdiff(locs_all_deps, deployments_master)


###############################################################################  