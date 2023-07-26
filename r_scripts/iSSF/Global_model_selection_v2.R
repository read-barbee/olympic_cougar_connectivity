#### All Subsets Model Selection v2 ####

# Author: Read Barbee

# Date:2023-07-19 

# Purpose:


################################ libraries #################################
library(tidyverse)
library(glmmTMB)
#library(MuMIn)
#library(INLA)
library(beepr)

#########################################################################
##
## 1. Import and format step data
##
##########################################################################
#no imputation
steps <- read_csv("data/Location_Data/Steps/2h_steps_unscaled_no_imp_7-12-2023.csv")


#set all negative elevations to 0 and filter dispersal tracks for post dispersal event
steps <- steps %>% 
  mutate(elevation = ifelse(elevation < 0, 0, elevation)) %>% 
  select(-c(aspect_deg, aspect_rad)) %>% 
  filter(is.na(dispersing) | dispersing==TRUE)

#collapse landcover categories
steps <- steps %>% mutate(land_cover_usfs_lumped = fct_collapse(land_cover_usfs, trees = c("trees", "tall_trees_shrubs", "gfh_tree_mix", "tree_shrub_mix",  "barren_tree_mix"),
                                                                shrubs = c("tall_shrubs", "shrubs", "gfh_shrub_mix", "barren_shrub_mix"),
                                                                gfh = c("gfh", "barren_gfh_mix"),
                                                                barren = c("barren_impervious", "snow_ice")),
                          land_use_usfs_lumped = fct_collapse(land_use_usfs, agriculture = c("agriculture", "rangeland_pasture")), .after = land_cover_usfs)

#scale and dummify covariates
steps_scaled <- steps %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), scale)) %>% 
  mutate(across(c(sl_, ta_, gpp:perc_nonveg, precip:tpi, roads_hii:power_hii), as.numeric)) %>% 
  select(-c(land_cover_usfs, land_use_usfs)) %>% 
  DataExplorer::dummify(select=c("land_use_usfs_lumped", "land_cover_usfs_lumped", "season", "hunting_season", "calving_season")) #%>% 
#select(case_, animal_id, step_id_, gpp:calving_season_yes) %>% 

#########################################################################
##
## 2. Fit global model and auto-dredge (not working)
##
##########################################################################

#steps_scaled_no_na <- steps_scaled %>% na.omit()

#Fit global glmmTMB model: takes about 1.7 hours
global <- glmmTMB(case_ ~ -1 +
                    #fixed effects
                    land_use_usfs_lumped_forest +
                    ndvi + I(ndvi^2) +
                    tree_cover_hansen +
                    gpp +
                    npp +
                    popdens_hii + I(popdens_hii^2) +
                    land_cover_usfs_lumped_water +
                    slope + I(slope^2) +
                    perc_tree_cover + I(perc_tree_cover^2) +
                    northing +
                    land_cover_usfs_lumped_gfh +
                    perc_nonveg + I(perc_nonveg^2) +
                    perc_nontree_veg + I(perc_nontree_veg^2) +
                    tri + I(tri^2) +
                    land_cover_usfs_lumped_shrubs +
                    landuse_hii + I(landuse_hii^2) +
                    easting + I(easting^2) +
                    roads_hii + I(roads_hii^2) +
                    dist_water + I(dist_water^2) +
                    rails_hii +
                    infra_hii +
                    #stratum-based intercept
                    (1|step_id_) +
                    #random slopes
                    (0 + land_use_usfs_lumped_forest | animal_id) +
                    (0 +  ndvi | animal_id) +
                    (0 +  tree_cover_hansen | animal_id) +
                    (0 +  gpp | animal_id) +
                    (0 +  npp | animal_id) +
                    (0 +  popdens_hii | animal_id) +
                    (0 +  land_cover_usfs_lumped_water | animal_id) +
                    (0 +   slope | animal_id) +
                    (0 +   perc_tree_cover | animal_id) +
                    (0 +  northing | animal_id) +
                    (0 +  land_cover_usfs_lumped_gfh | animal_id) +
                    (0 +  perc_nonveg | animal_id) +
                    (0 +  perc_nontree_veg | animal_id) +
                    (0 +  tri | animal_id) +
                    (0 +  land_cover_usfs_lumped_shrubs | animal_id) +
                    (0 +  landuse_hii | animal_id) +
                    (0 +  easting | animal_id) +
                    (0 +  roads_hii | animal_id) +
                    (0 +  dist_water | animal_id) +
                    (0 +  rails_hii | animal_id) +
                    (0 +  infra_hii | animal_id),
                  family=poisson,
                  data = steps_scaled,
                  doFit=FALSE); 
global$parameters$theta[1] <- log(1e3)
global$mapArg <- list(theta=factor(c(NA, 1:21)))

system.time(fit <- fitTMB(global)); beep("fanfare")

#dredging isn't computationally feasible
#system.time(all_fits <- dredge(fit))


#########################################################################
##
## 3. Manual backward stepwise selection
##
##########################################################################
#1. remove rails (p=0.73)

global2 <- glmmTMB(case_ ~ -1 +
                    #fixed effects
                    land_use_usfs_lumped_forest +
                    ndvi + I(ndvi^2) +
                    tree_cover_hansen +
                    #gpp +
                    #npp +
                    popdens_hii + I(popdens_hii^2) +
                    land_cover_usfs_lumped_water +
                    slope + #I(slope^2) +
                    perc_tree_cover + I(perc_tree_cover^2) +
                    northing +
                    land_cover_usfs_lumped_gfh +
                    perc_nonveg + I(perc_nonveg^2) +
                    perc_nontree_veg + I(perc_nontree_veg^2) +
                    #tri + I(tri^2) +
                    land_cover_usfs_lumped_shrubs +
                    #landuse_hii + I(landuse_hii^2) +
                    easting + #I(easting^2) +
                    roads_hii + I(roads_hii^2) +
                    #dist_water + I(dist_water^2) +
                    #rails_hii +
                    #infra_hii +
                    #stratum-based intercept
                    (1|step_id_) +
                    #random slopes
                    (0 + land_use_usfs_lumped_forest | animal_id) +
                    (0 +  ndvi | animal_id) +
                    (0 +  tree_cover_hansen | animal_id) +
                    #(0 +  gpp | animal_id) +
                    #(0 +  npp | animal_id) +
                    (0 +  popdens_hii | animal_id) +
                    (0 +  land_cover_usfs_lumped_water | animal_id) +
                    (0 +   slope | animal_id) +
                    (0 +   perc_tree_cover | animal_id) +
                    (0 +  northing | animal_id) +
                    (0 +  land_cover_usfs_lumped_gfh | animal_id) +
                    (0 +  perc_nonveg | animal_id) +
                    (0 +  perc_nontree_veg | animal_id) +
                    #(0 +  tri | animal_id) +
                    (0 +  land_cover_usfs_lumped_shrubs | animal_id) +
                    #(0 +  landuse_hii | animal_id) +
                    (0 +  easting | animal_id) +
                    (0 +  roads_hii | animal_id),
                    #(0 +  dist_water | animal_id),
                    #(0 +  rails_hii | animal_id) +
                    #(0 +  infra_hii | animal_id),
                  family=poisson,
                  data = steps_scaled,
                  doFit=FALSE); 
global2$parameters$theta[1] <- log(1e3)
global2$mapArg <- list(theta=factor(c(NA, 1:14)))

system.time(fit2 <- fitTMB(global2)); beep("fanfare")


#########################################################################
##
## 3. Auto dredge TMB
##
##########################################################################
vars <- c("land_use_usfs_lumped_forest + (0 + land_use_usfs_lumped_forest | animal_id)",
          "ndvi + I(ndvi^2) + (0 + ndvi | animal_id)",
          "tree_cover_hansen + (0 + tree_cover_hansen | animal_id)",
          "gpp + (0 + gpp | animal_id)",
          "npp + (0 + npp | animal_id)",
          "popdens_hii + I(popdens_hii^2) + (0 + popdens_hii | animal_id)",
          "land_cover_usfs_lumped_water + (0 + land_cover_usfs_lumped_water | animal_id)",
          "slope + I(slope^2) + (0 + slope | animal_id)",
          "perc_tree_cover + I(perc_tree_cover^2) + (0 + perc_tree_cover | animal_id)",
          "northing + (0 + northing | animal_id)",
          "land_cover_usfs_lumped_gfh + (0 + land_cover_usfs_lumped_gfh | animal_id)",
          "perc_nonveg + I(perc_nonveg^2) + (0 + perc_nonveg | animal_id)",
          "perc_nontree_veg + I(perc_nontree_veg^2) +  (0 + perc_nontree_veg | animal_id)",
          "tri + I(tri^2) + (0 + tri | animal_id)",
          "land_cover_usfs_lumped_shrubs + (0 + land_cover_usfs_lumped_shrubs | animal_id)",
          "landuse_hii + I(landuse_hii^2) +  (0 + landuse_hii| animal_id)",
          "easting + I(easting^2) + (0 + easting | animal_id)",
          "roads_hii + I(roads_hii^2) + (0 + roads_hii | animal_id)",
          "dist_water + I(dist_water^2) + (0 + dist_water | animal_id)",
          "rails_hii + (0 + rails_hii | animal_id)",
          "infra_hii + (0 + infra_hii | animal_id)")


all_comb <- do.call("c", lapply(seq_along(vars), function(i) combn(vars, i, FUN = list)))


#~ 10 min for 2,097,151 combinations of 21 covariates
forms <- list()

for (i in 1:length(all_comb)){
  var_i <- all_comb[[i]]
  forms[[i]] <- as.formula(paste("case_",  paste("-1", "(1|step_id_)", paste(var_i, collapse="+"), sep="+"), sep="~"))
}

mods <- list()
system.time(for (i in 1:length(forms)){
  form <- forms[[i]]
  mod <- glmmTMB(form, family=poisson, data=steps_scaled, doFit=FALSE)
  mod$parameters$theta[1] <- log(1e3)
  map_length <- 1:(length(mod$parameters$theta)-1)
  mod$mapArg <- list(theta=factor(c(NA, map_length)))
  
  mods[[i]] <- fitTMB(mod)
})


#########################################################################
##
## 3. Auto dredge (INLA)
##
##########################################################################

steps_scaled2 <- steps_scaled %>% 
  select(-c(disp_qual, disp_date_nsd, dispersing)) %>% 
  na.omit() %>% 
  mutate(y = as.numeric(case_),
         animal_id = as.numeric(factor(animal_id)),
         ndvi2 = ndvi^2,
         popdens_hii2 = popdens_hii^2,
         slope2 = slope^2,
         perc_tree_cover2 = perc_tree_cover^2,
         perc_nonveg2 = perc_nonveg^2,
         perc_nontree_veg2 = perc_nontree_veg^2,
         landuse_hii2 = landuse_hii^2,
         easting2 = easting^2,
         roads_hii2 = roads_hii^2,
         dist_water2 = dist_water^2) %>% 
  mutate( step_id = paste0(animal_id, step_id_, sep = "-"),
          id1 = animal_id,
          id2 = animal_id,
          id3 = animal_id,
          id4 = animal_id,
          id5 = animal_id,
          id6 = animal_id,
          id7 = animal_id,
          id8 = animal_id,
          id9 = animal_id,
          id10 = animal_id,
          id11 = animal_id,
          id12 = animal_id,
          id13 = animal_id,
          id14 = animal_id,
          id15 = animal_id,
          id16 = animal_id,
          id17 = animal_id,
          id18 = animal_id,
          id19 = animal_id,
          id20 = animal_id,
          id21 = animal_id,
          id22 = animal_id,
          id23 = animal_id,
          id24 = animal_id,
          id25 = animal_id,
          id26 = animal_id,
          id27 = animal_id,
          id28 = animal_id,
          id29 = animal_id)


vars <- c("land_use_usfs_lumped_forest + f(id1, land_use_usfs_lumped_forest, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "ndvi + f(id2, ndvi, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "ndvi2 + f(id3, ndvi2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "tree_cover_hansen + f(id4, tree_cover_hansen, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "gpp + f(id5, gpp, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "land_use_usfs_lumped_other + f(id6, land_use_usfs_lumped_other, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "grass + f(id6, grass, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "shrub + f(id7, shrub, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "wetland + f(id8, wetland, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "popbig + f(id9, popbig, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "riv_width + f(id10, riv_width, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "roads_int + f(id11, roads_int, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))")

all_comb <- do.call("c", lapply(seq_along(vars), function(i) combn(vars, i, FUN = list)))
response <- "y"
mods <- list()
dic_list <- list()
## Took 32273 seconds (~9hrs) on Legion for 11 covariates (2047 combinations)
system.time(
  for(i in 1:length(all_comb)){
    var_i <- all_comb[[i]]
    form <- as.formula(paste(response, paste("-1", "f(strat, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T)))", paste(var_i, collapse="+"), sep="+"), sep="~"))
    
    mods[[i]] <- inla(form, family ="Poisson", data=disp_inla,
                      control.fixed = list(
                        mean = mean.beta,
                        prec = list(default = prec.beta)), control.compute = list(dic = TRUE))
    
    dic_list[[i]] <- mods[[i]]$dic$dic
    
    print(i/length(all_comb))
  }
)

#########################################################################
##
## 3. Manual dredge (INLA)
##
##########################################################################

inla_global_form <-as.formula(y ~ -1 + 
                      f(step_id, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T))) +
                      land_use_usfs_lumped_forest + f(id1, land_use_usfs_lumped_forest, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      ndvi + f(id2, ndvi, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      ndvi2 + f(id3, ndvi2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      tree_cover_hansen + f(id4, tree_cover_hansen, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      gpp + f(id5, gpp, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      land_use_usfs_lumped_other + f(id6, land_use_usfs_lumped_other, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      popdens_hii + f(id7, popdens_hii, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      popdens_hii2 + f(id8, popdens_hii2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      slope + f(id9, slope, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      slope2 + f(id10, slope2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      perc_tree_cover + f(id11, perc_tree_cover, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      perc_tree_cover2 + f(id12, perc_tree_cover2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      northing + f(id13, northing, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      land_cover_usfs_lumped_gfh + f(id14, land_cover_usfs_lumped_gfh, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      perc_nonveg + f(id15, perc_nonveg, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      perc_nonveg2 + f(id16, perc_nonveg2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      perc_nontree_veg + f(id17, perc_nontree_veg, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      perc_nontree_veg2 + f(id18, perc_nontree_veg2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      land_cover_usfs_lumped_shrubs + f(id19, land_cover_usfs_lumped_shrubs, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      landuse_hii + f(id20, landuse_hii, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      landuse_hii2 + f(id21, landuse_hii2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      easting + f(id22, easting, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      easting2 + f(id23, easting2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      roads_hii + f(id24, roads_hii, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      roads_hii2 + f(id25, roads_hii2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      dist_water + f(id26, dist_water, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      dist_water2 + f(id27, dist_water2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      rails_hii + f(id28, rails_hii, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                      infra_hii + f(id29, infra_hii, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))))


inla_global_form2 <-as.formula(y ~ -1 + 
                                f(step_id, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T))) +
                                land_use_usfs_lumped_forest + f(id1, land_use_usfs_lumped_forest, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05)))) +
                                ndvi + f(id2, ndvi, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05)))) +
                                #ndvi2 + f(id3, ndvi2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                tree_cover_hansen + f(id4, tree_cover_hansen, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05)))) +
                                #gpp + f(id5, gpp, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                #land_use_usfs_lumped_other + f(id6, land_use_usfs_lumped_other, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                popdens_hii + f(id7, popdens_hii, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05)))) +
                                #popdens_hii2 + f(id8, popdens_hii2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                slope + f(id9, slope, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05)))) +
                                #slope2 + f(id10, slope2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                perc_tree_cover + f(id11, perc_tree_cover, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05)))) +
                                #perc_tree_cover2 + f(id12, perc_tree_cover2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                northing + f(id13, northing, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05)))) +
                                #land_cover_usfs_lumped_gfh + f(id14, land_cover_usfs_lumped_gfh, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                #perc_nonveg + f(id15, perc_nonveg, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                #perc_nonveg2 + f(id16, perc_nonveg2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                #perc_nontree_veg + f(id17, perc_nontree_veg, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                #perc_nontree_veg2 + f(id18, perc_nontree_veg2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                #land_cover_usfs_lumped_shrubs + f(id19, land_cover_usfs_lumped_shrubs, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                #landuse_hii + f(id20, landuse_hii, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                #landuse_hii2 + f(id21, landuse_hii2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                easting + f(id22, easting, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05)))) +
                                #easting2 + f(id23, easting2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                roads_hii + f(id24, roads_hii, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05)))) +
                                #roads_hii2 + f(id25, roads_hii2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                dist_water + f(id26, dist_water, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05)))) +
                                #dist_water2 + f(id27, dist_water2, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                #rails_hii + f(id28, rails_hii, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(1,0.05)))) +
                                infra_hii + f(id29, infra_hii, values=1:82,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05)))))

#mean.beta <- 1e-4
prec.beta <- 1e-4

#stopped manually at > 17 hours...
system.time(inla_global <- inla(inla_global_form2, 
                    family = "Poisson", 
                    data=steps_scaled2,
                    control.fixed = list(
                      mean = 0, #mean.beta
                      prec = list(default = prec.beta)), verbose = TRUE, control.compute = list(dic = TRUE)))







##### CV TEST ZONE
# cov_stack <- terra::rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2.tif")
# 
# names(cov_stack) <- c("tree_cover_hansen",
#                       "gpp",
#                       "infra_hii",
#                       "landuse_hii",
#                       "land_cover_usfs",
#                       "land_use_usfs",
#                       "npp",
#                       "popdens_hii",
#                       "power_hii",
#                       "precip",
#                       "rails_hii",
#                       "roads_hii",
#                       "elevation",
#                       "slope",
#                       "aspect",
#                       "tri",
#                       "tpi",
#                       "perc_tree_cover",
#                       "perc_nontree_veg",
#                       "perc_nonveg",
#                       "ndvi",
#                       "evi",
#                       "dist_water")
# 
# library(sf)
# 
# all_sf <- steps_scaled %>% st_as_sf(coords=c("x1_", "y1_"), crs = 5070)
# all_poly <- all_sf %>% sf::st_buffer(dist = 10000) %>%
#   sf::st_union() %>%
#   sf::st_convex_hull()
# cov_stack_cropped <- terra::crop(cov_stack, all_poly)






