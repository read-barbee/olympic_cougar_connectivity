#### Comparing ER locations to Deployment Dates####

# Author: Read Barbee

# Date:2023-01-23 

# Purpose: Compare date ranges of EarthRanger locations to collar deployment dates.


###############################################################################
#### Library / Functions / Data ####

#### Library ####
library(tidyverse)
library(lubridate)
library(janitor)
library(naniar)
library(svMisc)
library(foreach)
library(doParallel)


#### Functions ####

extract_deployments <- function(data_p, animal_id_p, collar_id_p, start_date_p, end_date_p) {
  if (is.na(end_date_p) == TRUE){ #don't filter by end date if it's not included
    trimmed_track <- data_p %>% 
      filter(collar_id == collar_id_p) %>% 
      filter(date_time_gmt>= as.POSIXct(start_date_p))
  } else{
    trimmed_track <- data_p %>% 
      filter(collar_id == collar_id_p) %>% 
      filter(date_time_gmt>= as.POSIXct(start_date_p) & date_time_gmt<=as.POSIXct(end_date_p)) 
  }
#extract deployments based on collar ID and deployment dates and append columns of animal_id and deployment_id (Working)
  
  return(trimmed_track)
}

#### Data ####

## Raw ##
er_locs_raw <- read_csv("data/Location Data/All Location Data from EarthRanger/OCP_full_loc_export_2023-01-13.csv")

deployments_raw <- read_csv("data/Location Data/OCP_Cougar_Deployments_1-13-2023.csv")


## Clean ##

er_locs_clean <- er_locs_raw %>% 
  clean_names() %>% 
  select(subject_1:lon) %>% 
  mutate(recorded_at = ymd_hms(recorded_at)) %>% 
  rename(subject = subject_1)

deployments_clean <- deployments_raw %>% 
  clean_names() %>% 
  mutate(deployment_id = paste0(name,"_",collar_id),
         start_date = mdy(start_date),
         end_date = mdy(end_date)) %>% 
  select(deployment_id, everything()) %>% 
  mutate(end_date = replace_na(end_date, mdy("12/31/9999")))



#### Actions ####

#### Step 1: Set up parallel computing cluster ####

#define the number of cores to be used for parallel processing. Leave one free for other tasks
n_cores <- detectCores() - 1

#create parallel computing cluster
cluster <- makeCluster(n_cores,
                       type = "FORK"
                       )
#check cluster definition
print(cluster)

#register cluster with DoParallel
registerDoParallel(cl = cluster)

#check if registered successfully
getDoParRegistered()

#check how many workers are available
getDoParWorkers()


#### Step 2: Assign deployment to every ER location by date ####

#Parallelized loop to popoulate each row of the deployment_id column of er_locs_clean (takes about 13 min to run with 9 cores) 


n <-nrow(er_locs_clean)

loop <- foreach(i=1:n, 
                .combine = 'c', 
                .packages = c('dplyr', 'lubridate')) %dopar% {
dep_select <- deployments_clean %>% 
  filter(name == er_locs_clean$subject[[i]] &
         er_locs_clean$recorded_at[[i]] %within% interval(start_date, end_date))
  
  if (nrow(dep_select) == 1){
    return(dep_select[[1,1]])
  }else if (nrow(dep_select) == 0){
    return(NA)
  }
}

stopCluster(cl=cluster)

#note--loop created fewer rows than contained in dataframe

#these lines aren't working
# er_locs_clean <- er_locs_clean %>% 
#   mutate(deployment_id = loop)
# 
# left_join(er_locs_clean, loop)

#New commands:
# command shift O shows table of contents for script
#command shift 0 restarts R.

###############################################################################  