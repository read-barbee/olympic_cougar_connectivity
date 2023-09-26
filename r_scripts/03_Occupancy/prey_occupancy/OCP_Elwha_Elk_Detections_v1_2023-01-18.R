#### OCP Elwha Elk Detections v1####

# Author: Read Barbee

# Date:2023-01-18 

# Purpose:Collate all elk detections from Elwha grid for Vinny's occupancy model


###############################################################################
#### Library / Functions / Data ####

#### Library ####
library(tidyverse)
library(lubridate)
library(janitor)
library(zoo)
#### Functions ####

#### Data ####
load("data/Camera Data/2019/dtbs.2019.Rdata")

metadata <- read_csv("data/Camera Data/2019/SNA_20190703_20200214_2022.09.19_21.23_metadata_tbl.csv")

#activity sheet with pull dates
camact2 <-  read_csv("data/Camera Data/2019/2019.Activity.Sheet.csv")


### Part I: Camera deployment table ##############################

#Create Deployment ID column formatted like the other datatables
camact3 <- camact2 %>% 
  mutate(CameraID = paste0("CAM",CameraID)) %>%
  mutate(CameraID_2 = case_when(
    Cameras == "Camera1" ~ paste0(CameraID),
    Cameras == "Camera2" ~ paste0(CameraID,"_2"),
    Cameras == "Camera3" ~ paste0(CameraID,"_3")
    )
    ) %>% 
  select(CameraID_2,CameraID, everything(), -Study) %>% 
  unite("DeploymentID", c(Station, Cameras, CameraID_2), remove = FALSE) %>% 
  select(DeploymentID, 
         Station, 
         everything(),
         -c(CameraID, Cameras, Cam.Num)) %>% 
  rename(CameraID = CameraID_2)


#pivot activity history into activitydate and activitystatus columns
camact_pivot <- camact3 %>% 
  pivot_longer(cols = 9:257,
               names_to = "ActivityDate",
               values_to = "ActivityStatus") 

#calculate last active date for each deployment from activity history
last_active <- camact_pivot %>% 
  group_by(DeploymentID, ActivityStatus) %>% 
  summarize(DeploymentDate=min(ActivityDate), 
            DateLastActive= max(ActivityDate)) %>% 
  filter(ActivityStatus == 1) %>% 
  ungroup() %>% 
  select(DeploymentID, DateLastActive)

#camera missing last active date was stolen
View(anti_join(camact3, last_active, by= "DeploymentID"))

#join last active dates to camact3
deployments <- left_join(camact3, last_active, by = "DeploymentID") %>% relocate(DateLastActive, .after = Set_up_Date) %>% 
  rename(DeploymentDate = Set_up_Date,
         PullDate = Pull_Date) %>% 
  select(DeploymentID:PullDate)

#replace NA values with 999, format activity date to POSIXct
camact_pivot2 <- camact_pivot %>% 
  mutate(ActivityStatus = replace_na(camact_pivot$ActivityStatus, 999),
         ActivityDate = mdy(ActivityDate)) 


#calculate dates of activity intervals for each camera based on 0, 1, and NAs in ActivityStatus column

activity_intervals <- camact_pivot2 %>% 
    group_by(CameraID) %>% 
    mutate(cluster_id= na.locf(ifelse(
      sequence(rle(as.character(ActivityStatus))$lengths)==1, row_number(), NA))) %>% 
    group_by(CameraID,cluster_id) %>% 
    summarise(start = min(ActivityDate), 
              end = max(ActivityDate),
              state = min(ActivityStatus)) %>% 
    select(CameraID, state, everything(), -cluster_id)

problem_intervals <- activity_intervals %>% 
  filter(state==0)

#join activity intervals to deployment table
test <- deployments %>% 
  left_join(problem_intervals, by= "CameraID") %>% 
  pivot_wider()


######### Part II: Elk Detection Table ###########################

#### Explore Data Frames ####

#camact: one row for each camera deployed for each date in the season.
# Unique Field: CameraID + Date
#Cameras column redundant to CameraID
get_dupes(camact, c(CameraID, Date))

#indvdl_obsrvtn_tbl: Puma detections with some individual IDs. One row for each attribute of each trigger
# Unique Field: TriggerID + Attribute +ImageNUmber + UserInput + Group + DateTimeInput

get_dupes(indvdl_obsrvtn_tbl, c(TriggerID, ImageNumber, Attribute, UserInput, Group, DateTimeInput))

#spcs_obsrvtn_tbl: Classified species observations. One row for each attribute of each trigger, each mode of entry. Species attribute duplicated for human and machine classifer
# Unique Field: TriggerID + Attribute + UserInput + DateTimeInput + ImageNumber + Group

get_dupes(spcs_obsrvtn_tbl, c(TriggerID, Attribute, UserInput, DateTimeInput, ImageNumber, Group))

#trggr_tbl: one row for each camera file with metadata. May be multiple files for single trigger.
#Unique Field: primary key OR FileName_New OR TriggerID + Image Number
get_dupes(trggr_tbl, primarykey)
get_dupes(trggr_tbl, FileName_New)
get_dupes(trggr_tbl, c(TriggerID, ImageNumber))


#metadata: all metadata for each file. One row for each file. May be multiple files for single trigger. Identical in structure to trggr_tbl but with more metadata fields. 4 pairs of duplicated rows cause discrepancy in rowcount with trggr_tbl

#Unique Field: primary key OR FileName_New 
get_dupes(metadata, c(TriggerID, ImageNumber)) #indicates duplicate primary keys
#remove 4 rows with duplicate primary keys
metadata_no_dupes <- metadata %>% distinct(primarykey, .keep_all = TRUE)
#check for differences in the dataframes
setdiff(metadata$primarykey, trggr_tbl$primarykey)

#verify there are no more duplicate primary keys in metadata
get_dupes(metadata_no_dupes, c(TriggerID, ImageNumber))


#Reduce metadata table to variables of interest. One row for each individual image. May be multiple rows for indiv triggers. Unique by primary key or TriggerID + ImageNumber

images_all <- metadata_no_dupes %>% 
  select(primarykey,
         DeploymentID,
         TriggerID,
         ImageNumber,
         Station,
         CameraID,
         X,
         Y,
         CameraDeploymentBeginDateTime,
         DateTimeOriginal,
         SurveyCheck,
         Species
         ) %>% 
  rename(TriggerDateTime = DateTimeOriginal)

get_dupes(images_all, c(TriggerID, ImageNumber))

#All uniuqe triggers (duplicate images per trigger removed). Unique by primary key or TriggerID
triggers_unique <- images_all %>% distinct(TriggerID, .keep_all = TRUE) %>% select(-ImageNumber)

get_dupes(triggers_unique, TriggerID)


#filter triggers to elk detections. Unique by primarykey or TriggerID
elk_detections <- triggers_unique %>% 
  filter(Species =="Elk_Roosevelt") 

get_dupes(elk_detections, TriggerID)



  

#load("data/Camera Data/2020/dtbs.2020.Rdata")

#load("data/Camera Data/2021/dtbs.2021.Rdata")


### Data Wrangling ###

#filter species observation table for elk observations

#there don't seem to be any duplicate triggers for elk, but I'll use distinct() anyway just to be safe
elk_detections1 <- spcs_obsrvtn_tbl %>% 
  filter(Attribute=="Species") %>% 
  distinct(TriggerID, .keep_all = TRUE) %>% 
  filter(Value=="Elk_Roosevelt")


elk_detections2 <- triggers_all %>% 
  distinct(TriggerID, .keep_all = TRUE) %>% 
  filter(Species=="Elk_Roosevelt")
  

dupes <- get_dupes(triggers_all, TriggerID)

diff <- as_tibble(setdiff(elk_detections1$TriggerID, elk_detections2$TriggerID))

diff2 <- diff <- as_tibble(setdiff(elk_detections2$TriggerID, elk_detections1$TriggerID))


elk_det_all <- left_join(elk_detections1, 
                         elk_detections2, 
                         by= "TriggerID") %>% 
  select(primarykey,
         DeploymentID,
         TriggerID,
         Station:Number)

anti_join()


#Join elk detections to trigger table by TriggerID and remove unnecessary columns

elk_det_join1 <- elk_detections %>% 
  left_join(triggers, by = "TriggerID") 

###############################################################################  