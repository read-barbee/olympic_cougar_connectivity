library(tidyverse)
library(terra)
library(parallel)
library(data.table)

# #import covariate stack
# cov_stack <- rast("data/homerange_habitat_layers_JR_9-14-22/cov_stack1_jr_4-28-2023.tif")
# 
# #pull forest and elevation layers out of covariate stack (replace elvation with TRI when available)
# covs <- terra::subset(cov_stack, c(1,5))
# 
# #convert raster to dataframe with coordinates and vlaues
# cov_df <- as.data.frame(covs, xy=TRUE, cells=TRUE)

tri <- rast("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/homerange_habitat_layers_JR_9-14-22/assets/puma_tri_op.tif")


forest <- rast("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/homerange_habitat_layers_JR_9-14-22/assets/puma_forest_op.tif")

stack <- c(forest, tri)

stack_rp <- project(stack, y = "epsg: 5070", res = 100) 

#convert raster to dataframe with coordinates and vlaues
cov_df <- as.data.frame(stack_rp, xy=TRUE, cells=TRUE)

cov_df <- na.omit(cov_df)

cov_df <- data.table(cov_df)

#import locations
locs <- read_csv("data/Location Data/Source Files/locations_master/gps_locs_master_5-16-2023.csv", col_types = list(fix_type = col_character()))

#put locations in chronological order for each individual and remove duplicates
locs_arr <- locs %>% 
  group_by(animal_id) %>% 
  arrange(date_time_local, .by_group=TRUE) %>% 
  distinct(animal_id, date_time_local, .keep_all = TRUE)

# Number of cells in study area
n_cells <- nrow(cov_df)

# Number of fix attempts
fix_attempts <- nrow(locs_arr)

# Lag between observations.  Should be 1 when no locations are missing.
Lag = rep(1, times=fix_attempts)
for(i in 3:fix_attempts){
  if(is.na(locs_arr$latitude[i-1]) == T) Lag[i] = Lag[i-1] + 1
}

# Detected (Yes = 1; No = 0)
detected = !is.na(locs_arr$latitude)

rm(locs, stack, forest, tri)

# Matrices of habitat data
distance <- as.matrix(dist(cov_df[,puma_forest_op:puma_tri_op]))
forest_cover <- matrix(cov_df %>% select(puma_forest_op), nrow=n_cells, ncol=n_cells, byrow=T)
elevation <- matrix(cov_df %>% select(puma_tri_op), nrow=n_cells, ncol=n_cells, byrow=T)

test <- as.matrix(mclapply(cbind(cov_df$x, cov_df$y), dist, mc.cores = detectCores()))

distance = as.matrix(dist(cbind(cov_df$x, cov_df$y)))
prcnt.sage = matrix(cov_df$puma_forest_op, nrow=n_cells, ncol=n_cells, byrow=T)
elevation = matrix(cov_df$puma_tri_op, nrow=n_cells, ncol=n_cells, byrow=T)

# Link function for probability of detection
logistic <- function(x) {
  exp(x)/(1+exp(x)) 
}


# RSF Likelihood
RSF <- function( a ){
  # Matrix of movement probabilities (based on Discrete Choice Function)
  EXP = exp(a[1]*distance + a[2]*forest_cover + a[3]*elevation)
  sumEXP = matrix(apply(EXP, 1, sum), nrow=n_cells, ncol=n_cells, byrow=F)
  # Pr[detection] 
  P = logistic(a[4] + a[5]*forest_cover[1,])
  # Matrix of Pr[movement] x Pr[detection]
  D = (EXP / sumEXP) * matrix(P, nrow=n_cells, ncol=n_cells, byrow=T)
  # Matrix of Pr[movement] x (1 - Pr[detection])
  B = (EXP / sumEXP) * matrix(1-P, nrow=n_cells, ncol=n_cells, byrow=T)
  # Creating B^(m-2) when there are m-1 missing locations in a sequence
  if(max(Lag) > 1){
    B.s = vector("list", max(Lag-1))
    B.s[[1]] = diag(n_cells)
    B.s[[2]] = B
    if(max(Lag) > 2){
      for(j in 3:max(Lag-1)){
        output = B
        for(k in 3:j){
          output = output %*% B
        }
        B.s[[j]] = output
      }
    }
  }
  # Individual likelihood for each observation
  prob = rep(NA, times=fix_attempts)
  for(i in 2:fix_attempts){
    # If 2 consecutive scheduled GPS fixes were successful
    if(detected[i] & Lag[i] == 1){
      prob[i] = D[cov_df$cell[i-1], cov_df$cell[i]]
    }
    # If there is a gap
    if(detected[i] & Lag[i] > 1){
      prob[i] = B[cov_df$cell[i-Lag[i]],] %*%
        B.s[[Lag[i]-1]] %*%
        D[,cov_df$cell[i]]
    }
  }
  -sum(log(prob), na.rm=T)
}  


# Maximize likelihood to obtain estimates of coefficients
fit = nlminb(RSF, start=rep(0, 5))

# Standard errors obtained from the Fisher Information Matrix using the 
# following function.
F.2nd.deriv <- function(x, FUN){
  FUN <- match.fun( FUN )
  d <- length(x)   # number of dimensions
  hess <- matrix( 0, nrow=d, ncol=d )
  eps <- 10e-7
  h <- (eps^(0.25))*x
  for( i in 1:d ){
    ei <- rep(0,d)
    ei[i] <- 1
    # compute diagonal element
    hess[i,i] <- (-FUN(x+2*h*ei) + 16*FUN(x+h*ei) - 30*FUN(x) + 
                    16*FUN(x-h*ei) - FUN(x-2*h*ei)) / (12*h[i]*h[i])
    if( (i+1) <= d ){
      for( j in (i+1):d ){
        ej <- rep(0,d)
        ej[j] <- 1
        # compute off diagonal element
        hess[i,j] <- (FUN(x+h*ei+h*ej) - FUN(x+h*ei-h*ej) - 
                        FUN(x-h*ei+h*ej) + FUN(x-h*ei-h*ej)) / (4*h[i]*h[j])
        # Assume symetric
        hess[j,i] <- hess[i,j]
      }
    }
  }
  hess
}
# Function call
hessian = F.2nd.deriv(fit$par, RSF)
SEs = sqrt(diag(solve(hessian)))


# Print results to screen
cat(" RSF Model Coefficients (SE):", "\n", 
    "   Distance (km) from Previous Location:", round(fit$par[1], 4), 
    "(", round(SEs[1], 4), ")", "\n",
    "   Percent Wyoming Big Sage:", round(fit$par[2], 4), "(", round(SEs[2], 4), ")", "\n",
    "   Elevation (km):", round(fit$par[3], 4), "(", round(SEs[3], 4), ")", "\n", 
    "Probabiilty of Detection Model Coefficients (SE):", "\n", 
    "   Intercept:", round(fit$par[4], 4), "(", round(SEs[4], 4), ")", "\n",
    "   Percent Wyoming Big Sage:", round(fit$par[5], 4), "(", round(SEs[5], 4), ")", "\n")




################################ OG Code #################################

# Read in available_habitat.txt
habitat = read.table("/Users/tb201494/Library/CloudStorage/Box-Box/Reference Literature/Spatial data sampling and management/Correcting for habitat bias/Nielson 2009 Supplement/available_habitat.txt", header=T)


# Read in gps_locations.txt
locations = read.table("/Users/tb201494/Library/CloudStorage/Box-Box/Reference Literature/Spatial data sampling and management/Correcting for habitat bias/Nielson 2009 Supplement/gps_locations.txt", header=T)

# Number of cells in study area
n.cells = max(habitat$unit.id)


# Number of fix attempts
fix.attempts = length(locations$fix)


# Lag between observations.  Should be 1 when no locations are missing.
Lag = rep(1, times=fix.attempts)
for(i in 3:fix.attempts){
    if(is.na(locations$unit.id[i-1]) == T) Lag[i] = Lag[i-1] + 1
}


# Detected (Yes = 1; No = 0)
detected = !is.na(locations$unit.id)


# Matrices of habitat data
distance = as.matrix(dist(cbind(habitat$utmX, habitat$utmY)))
prcnt.sage = matrix(habitat$prcnt.sage, nrow=n.cells, ncol=n.cells, byrow=T)
elevation = matrix(habitat$elevation, nrow=n.cells, ncol=n.cells, byrow=T)*1000


# Link function for probability of detection
logistic <- function(x) {
    exp(x)/(1+exp(x)) 
    }


# RSF Likelihood
RSF <- function( a ){
    # Matrix of movement probabilities (based on Discrete Choice Function)
    EXP = exp(a[1]*distance + a[2]*prcnt.sage + a[3]*elevation)
    sumEXP = matrix(apply(EXP, 1, sum), nrow=n.cells, ncol=n.cells, byrow=F)
    # Pr[detection] 
    P = logistic(a[4] + a[5]*prcnt.sage[1,])
    # Matrix of Pr[movement] x Pr[detection]
    D = (EXP / sumEXP) * matrix(P, nrow=n.cells, ncol=n.cells, byrow=T)
    # Matrix of Pr[movement] x (1 - Pr[detection])
    B = (EXP / sumEXP) * matrix(1-P, nrow=n.cells, ncol=n.cells, byrow=T)
    # Creating B^(m-2) when there are m-1 missing locations in a sequence
    if(max(Lag) > 1){
        B.s = vector("list", max(Lag-1))
        B.s[[1]] = diag(n.cells)
        B.s[[2]] = B
        if(max(Lag) > 2){
            for(j in 3:max(Lag-1)){
                output = B
                for(k in 3:j){
                    output = output %*% B
                }
                B.s[[j]] = output
            }
        }
    }
    # Individual likelihood for each observation
    prob = rep(NA, times=fix.attempts)
    for(i in 2:fix.attempts){
        # If 2 consecutive scheduled GPS fixes were successful
        if(detected[i] & Lag[i] == 1){
            prob[i] = D[locations$unit.id[i-1], locations$unit.id[i]]
        }
        # If there is a gap
        if(detected[i] & Lag[i] > 1){
            prob[i] = B[locations$unit.id[i-Lag[i]],] %*%
                      B.s[[Lag[i]-1]] %*%
                      D[,locations$unit.id[i]]
        }
    }
    -sum(log(prob), na.rm=T)
}  


# Maximize likelihood to obtain estimates of coefficients
fit = nlminb(RSF, start=rep(0, 5))

# Standard errors obtained from the Fisher Information Matrix using the 
# following function.
F.2nd.deriv <- function(x, FUN){
    FUN <- match.fun( FUN )
    d <- length(x)   # number of dimensions
    hess <- matrix( 0, nrow=d, ncol=d )
    eps <- 10e-7
    h <- (eps^(0.25))*x
    for( i in 1:d ){
        ei <- rep(0,d)
        ei[i] <- 1
        # compute diagonal element
        hess[i,i] <- (-FUN(x+2*h*ei) + 16*FUN(x+h*ei) - 30*FUN(x) + 
                        16*FUN(x-h*ei) - FUN(x-2*h*ei)) / (12*h[i]*h[i])
        if( (i+1) <= d ){
            for( j in (i+1):d ){
                ej <- rep(0,d)
                ej[j] <- 1
                # compute off diagonal element
                hess[i,j] <- (FUN(x+h*ei+h*ej) - FUN(x+h*ei-h*ej) - 
                            FUN(x-h*ei+h*ej) + FUN(x-h*ei-h*ej)) / (4*h[i]*h[j])
                # Assume symetric
                hess[j,i] <- hess[i,j]
            }
        }
    }
    hess
}
# Function call
hessian = F.2nd.deriv(fit$par, RSF)
SEs = sqrt(diag(solve(hessian)))


# Print results to screen
cat(" RSF Model Coefficients (SE):", "\n", 
    "   Distance (km) from Previous Location:", round(fit$par[1], 4), 
    "(", round(SEs[1], 4), ")", "\n",
    "   Percent Wyoming Big Sage:", round(fit$par[2], 4), "(", round(SEs[2], 4), ")", "\n",
    "   Elevation (km):", round(fit$par[3], 4), "(", round(SEs[3], 4), ")", "\n", 
    "Probabiilty of Detection Model Coefficients (SE):", "\n", 
    "   Intercept:", round(fit$par[4], 4), "(", round(SEs[4], 4), ")", "\n",
    "   Percent Wyoming Big Sage:", round(fit$par[5], 4), "(", round(SEs[5], 4), ")", "\n")
