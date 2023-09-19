#### Extract Image Metadata ####

# Author: Read Barbee

# Date:2023-09-15 

# Purpose: Extract metadata from Browning DCL and all other camera images


#INFO:
#The entire camera data workflow, including organizing, timeshifting, and renaming raw images, is available in vignettes at https://cran.r-project.org/web/packages/camtrapR/index.html.

#I find vignette "3. Data Extraction from Images and videos, creating occupancy and secr input" particularly helpful: https://cran.r-project.org/web/packages/camtrapR/vignettes/camtrapr3.pdf

#IMPORTANT: ExifTool must be installed on your computer for this script to work. Installation instructions can be found in vignette "1. Organising raw camera trap images in camtrapR" found here: https://cran.r-project.org/web/packages/camtrapR/vignettes/camtrapr1.pdf

#MacOS requires installation but no additional setup. Windows does not require installation but requires additional setup


################################ Libraries #################################
library(camtrapR)
library(tidyverse)


#########################################################################
##
## 1. Set working directory
##
##########################################################################
#If species are identified, the image folder structure should be Station/Camera/Species

#If you want to keep the "Check" folders, and your structure is Station/Camera/Check, just rename the 'Species' column from the recordTable function as "Check" (as below). 

#If the structure is just Station/Camera, omit the "cameraID" and "camerasIndependent" arguments from the recordTable function below and rename the Species column from the output as "Camera". 


#set your working directory with all images and folder substructure here
wd <- "/your_directory/main_image_folder"

#check number of files to make sure the directory can be accessed
length(list.files(wd, pattern = "JPG", recursive = TRUE))




#########################################################################
##
## 2. Extract metadata from images and save
##
##########################################################################

#Note: This function automatically removes image records from the same camera and station that are less than minDeltaTime apart (in minutes). Default should be 1 minute. Customize by uncommenting the minDeltaTime argument below.

#many more customization options are available, including writing to a csv within the function. see ?recordTable.

rec_tab <- recordTable(inDir = wd, 
                       cameraID = "directory", 
                       camerasIndependent = TRUE,
                       IDfrom = "directory",
                       #minDeltaTime = 1,
                       timeZone = "US/Pacific") %>% 
  rename(Check = Species)


write_csv(rec_tab, "your_output_directory")


