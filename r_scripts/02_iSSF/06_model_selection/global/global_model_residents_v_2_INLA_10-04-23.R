#### Global Resident Muff Model ####

# Author: Read Barbee

# Date:2023-09-18

# Purpose:


################################ libraries #################################
library(tidyverse)
#library(glmmTMB)
#library(amt)
##library(MuMIn)
#library(INLA)
library(beepr)
#library(buildmer)
library(INLA)
library(INLAutils)
library(terra)
library(sf)

#########################################################################
##
##  Specify Model Parameters
##
##########################################################################

params <- c("tree_cover_hansen",
            "ndvi",
            "tpi",
            "northing", 
            "tri", 
            "perc_tree_cover",
            "perc_nonveg",
            "roads_hii",
            "popdens_hii",
            "landuse_hii",
            "rails_hii",
            "infra_hii",
            "easting")

quad_params <- c(#"tree_cover_hansen",
                 "ndvi"
                # "tpi",
                # "northing", #not sure that this makes biological sense
                # "tri",
                #"perc_tree_cover",
                #"roads_hii",
                #"popdens_hii"#not sure that this makes biological sense
                #"landuse_hii" #not sure that this makes biological sense
                ) 


#elevation model
# params <- c("tree_cover_hansen",
#             "ndvi",
#             #"tpi",
#             "northing", 
#             #"tri", 
#             "perc_tree_cover",
#             "roads_hii",
#             "popdens_hii",
#             "landuse_hii",
#             "rails_hii",
#             "infra_hii",
#             "elevation")
# 
# 
# quad_params <- c("tree_cover_hansen",
#                  #"ndvi",
#                  "northing",
#                  "perc_tree_cover",
#                  "roads_hii",
#                  "popdens_hii",
#                  "landuse_hii",
#                  "elevation"
                 # ) 

#prior params from https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40#inla-1
mean.beta <- 0
prec.beta <- 1e-4 


rand_terms <- vector()

for(i in 1:length(params)){
  rand_terms[i] <- paste0("(0 + ", params[i], " | animal_id)")
}


quad_terms <- vector()

for(i in 1:length(quad_params)){
  quad_terms[i] <- paste0("I(",quad_params[i],"^2)")
}

#quad_terms <- NULL

#rand_terms_quad <- vector()

# for(i in 1:length(quad_terms)){
#   rand_terms_quad[i] <- paste0("(0 + ", quad_terms[i], " | animal_id)")
# }

#########################################################################
##
## 1. Import and format step data
##
##########################################################################

#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_10-02-2023.csv") %>% 
  mutate(ndvi = ndvi*0.0001)



#set all negative elevations to 0 and filter dispersal tracks for post dispersal event
steps_scaled <- steps %>% 
  filter(dispersal_status=="resident") %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(land_use_usfs, land_cover_usfs, dispersing:disp_qual, season:calving_season)) %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric))



#########################################################################
##
## 2. Fit global model with quadratics in INLA
##
##########################################################################


n_indiv = steps_scaled %>% distinct(animal_id) %>% count() %>% pull()

rand_terms_inla <- vector()
for(i in 1:length(params)){
  rand_terms_inla[i] <- paste0("f(", paste0("id", i), ", ", params[i], ", ", "values = ", paste0("1:",n_indiv), ", ", 'model = "iid", ', "hyper = list(theta = list(initial = log(1), fixed = FALSE, ", 'prior = "pc.prec", ', "param = c(1, .05))))")
}



run_global_inla <- function (dat){
  
  #convert response from categorical to numeric
  dat <- dat %>% mutate(case_ = as.numeric(case_))
  
  #create separate animal_id columns for each random effect
  for(i in 1:length(rand_terms_inla)){
    name <- as.name(paste0("id", i))
    dat[[name]] <- as.numeric(factor(dat$animal_id))
  }
  
  #construct the model formula with quadratic terms
  form <- as.formula(paste("case_ ~ -1 + ",  
                           #fixed effects
                           paste(c(params, quad_terms), collapse = " + "), "+",
                           #random intercept (strata)
                           "f(step_id_, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T))) + ",
                           #random slopes
                           paste(rand_terms_inla, collapse = " + ")))
  
  #fit the model    
  fit <-  inla(form, family ="Poisson", data=dat,
               control.fixed = list(
               mean = mean.beta,
               prec = list(default = prec.beta)),
               control.compute = list(waic=TRUE, dic = TRUE, cpo = FALSE)) #
  
  return(fit)
}

system.time(inla_fit_linear <- run_global_inla(steps_scaled)) ;beep("fanfare")

p <- autoplot(inla_fit_linear)

inla_fit_linear2 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_linear2)

inla_fit_linear3 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_linear3)

inla_fit_linear4 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_linear4)

inla_fit_linear5 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_linear5)

inla_fit_linear6 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_linear6)

inla_fit_linear7 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_linear7)

inla_fit_linear8 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_linear8)

inla_fit_quad <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_quad)

inla_fit_quad2 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_quad2)

inla_fit_quad3 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_quad2)

# res_plot <- ggplot_inla_residuals(inla_fit_global, observed = steps[1:100, "case_"] , binwidth = NULL)

#inla_fit_linear: all linear terms
#terms: "tree_cover_hansen + ndvi +tpi + northing + tri + perc_tree_cover + perc_nonveg + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii"

#resident rho: 0.9878788 
#disperser rho: 0.9636364 
#waic: 234642.23
#dic: 234958.53

#cross zero: infra_hii, landuse_hii, ndvi, perc_nonveg, popdens_hii, rails_hii, roads_hii, tpi
#remove: rails



#inla_fit_linear2: all linear - rails_hii
#terms: "tree_cover_hansen + ndvi +tpi + northing + tri + perc_tree_cover + perc_nonveg + roads_hii + popdens_hii + landuse_hii + infra_hii"

#resident rho: 0.9636364
#disperser rho: 0.9636364
#waic: 234642.74
#dic: 234959.01

#cross 0: infra_hii, landuse_hii, ndvi, perc_nonveg, popdens_hii, roads_hii, tpi

#remove: roads



#inla_fit_linear3: - rails_hii, roads_hii
#terms: "tree_cover_hansen + ndvi +tpi + northing + tri + perc_tree_cover + perc_nonveg + popdens_hii + landuse_hii + infra_hii"

#resident rho: 0.9636364 
#disperser rho:0.9636364 
#waic: 234718.27
#dic: 235032.08

#cross 0: infra_hii, landuse_hii, ndvi, perc_nonveg, popdens_hii, tpi

#remove: tpi


#inla_fit_linear4: - rails_hii, roads_hii. tpi
#terms: "tree_cover_hansen + ndvi + northing + tri + perc_tree_cover + perc_nonveg + popdens_hii + landuse_hii + infra_hii"

#resident rho: 0.9636364 
#disperser rho: 0.9636364 
#waic: 235023.85
#dic: 235334.49

#cross 0: infra_hii, landuse_hii, ndvi, perc_nonveg, popdens_hii

#remove: perc_nonveg


#inla_fit_linear5: - rails_hii, roads_hii. tpi, perc_nonveg
#terms: "tree_cover_hansen + ndvi + northing + tri + perc_tree_cover + popdens_hii + landuse_hii + infra_hii"

#resident rho: 0.9636364 
#disperser rho: 0.9636364 
#waic: 235063.03
#dic: 235370.56

#cross 0: infra_hii, landuse_hii, ndvi, popdens_hii

#remove: infra_hii


#inla_fit_linear6: - rails_hii, roads_hii. tpi, perc_nonveg, infra_hii
#terms: "tree_cover_hansen + ndvi + northing + tri + perc_tree_cover + popdens_hii + landuse_hii"

#resident rho: 0.9636364 
#disperser rho: 0.9636364 
#waic: 235072.11
#dic: 235379.48

#cross 0: landuse_hii, ndvi, popdens_hii

#remove: landuse_hii

#inla_fit_linear7: - rails_hii, roads_hii. tpi, perc_nonveg, infra_hii, landuse_hii
#terms: "tree_cover_hansen + ndvi + northing + tri + perc_tree_cover + popdens_hii "

#resident rho: 0.9515152
#disperser rho: 0.9636364 
#waic: 235093.52
#dic: 235399.77

#cross 0: ndvi, popdens_hii

#remove: ndvi

#inla_fit_linear8: - rails_hii, roads_hii. tpi, perc_nonveg, infra_hii, landuse_hii, ndvi
#terms: "tree_cover_hansen + northing + tri + perc_tree_cover + popdens_hii "

#resident rho: 
#disperser rho: 
#waic: 235355.79
#dic: 235648.88

#cross 0: popdens_hii

#remove: 


#inla_fit_quad: "tree_cover_hansen + ndvi +tpi + northing + tri + perc_tree_cover + perc_nonveg + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii + tree_cover_hansen2 + ndvi2 + tpi2 + northing2  + tri2 + perc_tree_cover2 + roads_hii2 + popdens_hii2 + landuse_hii2"

#resident rho: 
#disperser rho: 
#waic: 233635.58
#dic: 233938.79

#quads that cross 0: landuse_hii2, northing2, roads_hii2, tpi2, tri2
#linear that cross 0: easting, infra_hii, landuse_hii, perc_nonveg, perc_tree_cover, rails_hii, roads_hii, tpi, tree_cover_hansen, easting

#remove: 

#inla_fit_quad2: "tree_cover_hansen + ndvi +tpi + northing + tri + perc_tree_cover + perc_nonveg + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii + tree_cover_hansen2 + ndvi2 + perc_tree_cover2 + popdens_hii2"

#resident rho: 
#disperser rho: 
#waic: 233635.58
#dic: 233938.79

#quads that cross 0: none
#linear that cross 0: easting, infra_hii, landuse_hii, perc_nonveg, perc_tree_cover, rails_hii, roads_hii, tpi, tree_cover_hansen, easting

#remove: all quads but ndvi


#inla_fit_quad3: "tree_cover_hansen + ndvi +tpi + northing + tri + perc_tree_cover + perc_nonveg + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii + tree_cover_hansen2 + ndvi2 + perc_tree_cover2 + popdens_hii2"

#resident rho: 
#disperser rho: 
#waic: 233635.58
#dic: 233938.79

#quads that cross 0: none
#linear that cross 0: easting, infra_hii, landuse_hii, perc_nonveg, perc_tree_cover, rails_hii, roads_hii, tpi, tree_cover_hansen, easting

#remove: all quads but ndvi


#########################################################################
##
## 2. Old models
##
##########################################################################

#all quad:
#terms: "tree_cover_hansen + ndvi +tpi + northing + tri + perc_tree_cover + perc_nonveg + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii + tree_cover_hansen2 + ndvi2 + tpi2 + northing2  + tri2 + perc_tree_cover2 + roads_hii2 + popdens_hii2 + landuse_hii2"

#rho: -0.9272727  
#waic: 233661.43 
#dic: 233964.37


#non significant quads removed:
#terms: "tree_cover_hansen + ndvi + tpi + northing + tri + perc_tree_cover + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii + tree_cover_hansen2 + ndvi2 + perc_tree_cover2 + popdens_hii2 + landuse_hii2"

#rho: -0.9515152 
#waic: 233665.89
#dic: 233967.75


#Elevation model:
#terms: "tree_cover_hansen + ndvi  + northing  + perc_tree_cover + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii + elevation"

#resident rho: 1
#disperser rho: 1
#waic: 235036.21
#dic: 235345.02

#Elevation2 model:
#terms: "tree_cover_hansen + ndvi + northing  + perc_tree_cover + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii + elevation + elevation2"

#resident rho: 1
#disperser rho: 1
#waic: 235027.76
#dic:  235336.60

#Elevation2 all quads
#terms: "tree_cover_hansen + ndvi + northing + perc_tree_cover + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii + elevation + tree_cover_hansen2 + ndvi2 + northing2  + perc_tree_cover2 + roads_hii2 + popdens_hii2 + landuse_hii2 + elevation2"
#resident rho: -0.9272727 
#disperser rho: -0.8545455 
#waic: 233971.99
#dic: 234270.16

#Elevation2 all quads that don't cross 0
#terms: "tree_cover_hansen + ndvi  + northing  + perc_tree_cover + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii + elevation  + tree_cover_hansen2 + ndvi2 + perc_tree_cover2 + popdens_hii2 + elevation2"

#resident rho: -0.9272727 
#disperser rho: -0.8666667
#waic: 233975.92
#dic: 234273.82

#Elevation and all quads that don't cross 0 - popdens_hii2 (convergence issues)
#terms: "tree_cover_hansen + ndvi  + northing  + perc_tree_cover + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii + elevation  + tree_cover_hansen2 + ndvi2 + perc_tree_cover2 + elevation2"

#resident rho: -0.9151515 
#disperser rho: -0.8424242 
#waic: 234079.94
#dic: 234380.12

#Elevation and all quads that don't cross 0 - popdens_hii2 - tree_cover_hansen2: DOESN'T CONVERGE
#terms: "tree_cover_hansen + ndvi  + northing  + perc_tree_cover + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii + elevation + ndvi2 + perc_tree_cover2 + elevation2"

#resident rho: 
#disperser rho: 
#waic: 
#dic: 

#Elevation and all quads that don't cross 0 - popdens_hii2 - tree_cover_hansen2 - perc_tree_cover2
#terms: "tree_cover_hansen + ndvi  + northing  + perc_tree_cover + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii + elevation   + ndvi2 + elevation2"

#resident rho: -0.9272727 
#disperser rho: -0.8666667 
#waic:234190.48
#dic: 234491.68

#Elevation and all quads except for ndvi2
#terms: "tree_cover_hansen + ndvi + northing + perc_tree_cover + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii + elevation + tree_cover_hansen2 +  northing2  + perc_tree_cover2 + roads_hii2 + popdens_hii2 + landuse_hii2 + elevation2"

#resident rho:-0.7939394
#disperser rho: -0.5878788
#waic:234818.58
#dic: 235124.28

#### Try moving forward from best elevation model-- START HERE TOMORROW

#Elevation2 model:
#terms: "tree_cover_hansen + ndvi + northing  + perc_tree_cover + roads_hii + popdens_hii + landuse_hii + rails_hii + infra_hii + elevation + elevation2"

#resident rho: 1
#disperser rho: 1
#waic: 235027.76
#dic:  235336.60


#########################################################################
##
## 2. Gradually eliminate quadratics
##
##########################################################################

####### all quadratic terms included--takes forever to fit ###
inla_fit2 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit2)


#quadratics that cross 0: landuse_hii, northing, roads_hii, tpi, tri
#Eliminate the weakest one: tri


#all quadratics that cross 0 removed
inla_fit3 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit3)


#Elevation instead of tpi and tri
inla_fit_elev <- run_global_inla(steps_scaled) ;beep("fanfare")

#Elevation with quadratic
inla_fit_elev2 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_elev2)

#Elevation with quadratic all quads
inla_fit_elev_all2 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_elev_all2)

#quads that cross 0: landuse_hii, northing, roads, 
#linear that cross 0: infra_hii, landuse_hii, rails, roads, tree_cover_hansen

#Elevation with sig quads only
inla_fit_elev_sig2 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_elev_sig2)

#quads that cross 0: none
#linear that cross 0: infra_hii, landuse_hii, rails_hii

#remove quadratics one by one to find the problem. start with weakest: popdens2 (actually elevation, but don't want to remove that yet)

inla_fit_elev_prob_search1 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_elev_prob_search1)

#quads that cross 0: none
#linear that cross 0: infra_hii, landuse_hii, popdens_hii, rails_hii, roads_hii, tree_cover_hansen

#next to remove: tree_cover_hansen

inla_fit_elev_prob_search2 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_elev_prob_search2)

#quads that cross 0: none
#linear that cross 0: infra_hii, landuse_hii, popdens_hii, rails_hii, roads_hii, tree_cover_hansen

#next to remove: ndvi

inla_fit_elev_prob_search3 <- run_global_inla(steps_scaled) ;beep("fanfare")
p <- autoplot(inla_fit_elev_prob_search3)

#quads that cross 0: none
#linear that cross 0: infra_hii, landuse_hii, popdens_hii, rails_hii, roads_hii, tree_cover_hansen

#next to remove: ndvi

#########################################################################
##
## 2. Gradually eliminate quadratics
##
##########################################################################


mods <- list()
waic <- list()

inla_fit$waic$waic
inla_fit2$waic$waic

INLAstep(
  faml = "poisson",
  dataf = steps_scaled,
  
)



#########################################################################
##
## 2. LOO Cross-Validation in INLA
##
##########################################################################

inla_cv <- inla.group.cv(inla_fit)
ULOOCV = mean(inla_cv$cv)


#inlasloo() for spatial leave one out cross-validation

#control.fixed = list(
# mean = mean.beta,
# prec = list(default = prec.beta)), control.compute = list(dic = TRUE)


#########################################################################
##
## 3. Plot RSF-style predictive surface
##
##########################################################################

cov_stack <- terra::rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2_300m.tif")

#remove categorical layers from cov stack
cov_stack2 <- subset(cov_stack, c("land_cover_usfs", "land_use_usfs"), negate = TRUE)


#For loop to cap each layer of the covariate stack to the max and min values observed in the step dataset
cov_stack_capped <- list()
for (i in 1:length(names(cov_stack2))){
  name <- names(cov_stack2)[[i]]
  
  upper <- max(steps[[name]], na.rm = T)
  lower <- min(steps[[name]], na.rm = T)
  
  cov_stack_capped[[i]] <- clamp(cov_stack2[[i]], lower = lower, upper = upper)
  
  print(paste0(i, "/", length(names(cov_stack2))))
  
}

#add names to clamped layers
names(cov_stack_capped) <- names(cov_stack2)

#convert from list to SpatRaster
cov_stack_capped <- rast(cov_stack_capped)


#scale capped covariate stack for prediction
cov_stack_pred <- terra::scale(cov_stack_capped)


#extract fixed effect beta values from fitted_model
betas <- inla_fit_quad2$summary.fixed$mean

#reff <- ranef(global_fit)




#manual predictions
preds <- exp((betas[1] * cov_stack_pred$tree_cover_hansen) + 
               (betas[2] * cov_stack_pred$ndvi) + 
               (betas[3] * cov_stack_pred$tpi) + 
               (betas[4] * cov_stack_pred$northing) + 
               (betas[5] * cov_stack_pred$tri) + 
               (betas[6] * cov_stack_pred$perc_tree_cover) +
               (betas[7] * cov_stack_pred$perc_nonveg) +
               (betas[8] * cov_stack_pred$roads_hii) + 
               (betas[9] * cov_stack_pred$popdens_hii) +
               (betas[10] * cov_stack_pred$landuse_hii) +
               (betas[11] * cov_stack_pred$rails_hii) +
               (betas[12] * cov_stack_pred$infra_hii) +
               (betas[13] * cov_stack_pred$easting) +
               #(betas[10] * cov_stack_pred$elevation) +
               (betas[14] * cov_stack_pred$tree_cover_hansen^2) +
               (betas[15] * cov_stack_pred$ndvi^2) +
               #(betas[16] * cov_stack_pred$tpi^2) +
               #(betas[17] * cov_stack_pred$northing^2) +
               #(betas[18] * cov_stack_pred$tri^2) +
               (betas[16] * cov_stack_pred$perc_tree_cover^2) +
               #(betas[20] * cov_stack_pred$roads_hii^2) +
               (betas[17] * cov_stack_pred$popdens_hii^2))/
              #(betas[22] * cov_stack_pred$landuse_hii^2))/
              #(betas[17] * cov_stack_pred$elevation^2) +
  (1 + exp((betas[1] * cov_stack_pred$tree_cover_hansen) + 
             (betas[2] * cov_stack_pred$ndvi) + 
             (betas[3] * cov_stack_pred$tpi) + 
             (betas[4] * cov_stack_pred$northing) + 
             (betas[5] * cov_stack_pred$tri) + 
             (betas[6] * cov_stack_pred$perc_tree_cover) +
             (betas[7] * cov_stack_pred$perc_nonveg) +
             (betas[8] * cov_stack_pred$roads_hii) + 
             (betas[9] * cov_stack_pred$popdens_hii) +
             (betas[10] * cov_stack_pred$landuse_hii) +
             (betas[11] * cov_stack_pred$rails_hii) +
             (betas[12] * cov_stack_pred$infra_hii) +
             (betas[13] * cov_stack_pred$easting) +
             #(betas[10] * cov_stack_pred$elevation) +
             (betas[14] * cov_stack_pred$tree_cover_hansen^2) +
             (betas[15] * cov_stack_pred$ndvi^2) +
             #(betas[16] * cov_stack_pred$tpi^2) +
             #(betas[17] * cov_stack_pred$northing^2) +
             #(betas[18] * cov_stack_pred$tri^2) +
             (betas[16] * cov_stack_pred$perc_tree_cover^2) +
             #(betas[20] * cov_stack_pred$roads_hii^2) +
             (betas[17] * cov_stack_pred$popdens_hii^2)))
     #(betas[22] * cov_stack_pred$landuse_hii^2))/
     #(betas[17] * cov_stack_pred$elevation^2) +
   
terra::plot(preds)

#######################################################################
##
## 3. Bin predictions
##
##########################################################################

pred_vals <- terra::values(preds)

breaks <- quantile(pred_vals, probs = 0:10/10, na.rm = T)

bins <- cut(pred_vals, breaks, include.lowest=TRUE)

binned <- terra::classify(preds, rcl=as.vector(breaks))

terra::plot(binned)

levels(binned) <- 1:10



########################################################################
##
## 3. Internal validation (residents)
##
##########################################################################
locs <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_10-02-2023.csv")

residents <- locs %>% 
  filter(dispersal_status =="resident") %>% 
  select(-c(disp_date_nsd:dispersing))

eval_res <- residents %>% st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070)

#extract binned predicted probability of use values for each used location
loc_preds<- terra::extract(binned, eval_res)[,2] %>% as.numeric() %>% na.omit() 

#plot proportion of test points in each bin. 
ggplot() +
  geom_bar(aes(x=loc_preds, y = after_stat(prop))) 

n_bins <- length(unique(loc_preds))

#calculate pearson correlation between the bin number and the proportion of used locations in each bin
cor_cv_res <- cor.test(x = 1:n_bins, y = as.numeric(table(loc_preds))/(length(loc_preds)/n_bins), method = "spearman", exact = FALSE)$estimate

cor_cv_res 


########################################################################
##
## 3. External validation (dispersers)
##
##########################################################################

dispersers <- locs %>% 
  filter(dispersal_status =="disperser") %>% 
  filter(dispersing == TRUE)

eval_disp <- dispersers %>% st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070)

#extract binned predicted probability of use values for each used location
loc_preds_disp<- terra::extract(binned, eval_disp)[,2] %>% as.numeric() %>% na.omit() 

#plot proportion of test points in each bin. 
ggplot() +
  geom_bar(aes(x=loc_preds_disp, y = after_stat(prop))) 

n_bins <- length(unique(loc_preds_disp))

#calculate pearson correlation between the bin number and the proportion of used locations in each bin
cor_cv_res_disp <- cor.test(x = 1:n_bins, y = as.numeric(table(loc_preds_disp))/(length(loc_preds_disp)/n_bins), method = "spearman", exact = FALSE)$estimate

cor_cv_res_disp


