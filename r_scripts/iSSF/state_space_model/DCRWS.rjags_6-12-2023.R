
#require(chron)
require(rjags)
# dat4bugs(hseal.dat, 1) ... 1 = timestep of 1 day, this line must be changed to:
#		dat4bugs(gs617.dat, 2) to run for grey seal 617 with timestep=2 days, and similarly for grey seal 2986.

#x has one observation per regular interval
#y is all observations wtih NA rows for intervals with no observations

#Prep Data
#hseal.wb.dat <- dat4bugs(hseal.dat, 1)

#test$j <- rep(x = 1, times = length(idx))

y <- test$y
idx <- test$idx
x <- matrix(unlist(lapply(1:(length(idx)-1), function(i){y[idx[i],]})), length(idx)-1, 2, byrow = T)# if there are multiple observations in a regular interval, pull the last one from y
bmode <- rep(1, nrow(x) - 1)

#Define Initial Values
inits <- list(list(iSigma = matrix(c(1, 0, 0, 1), 2, 2), 
                   gamma = c(0.5, 0.5), 
                   theta = c(0, 0), 
                   alpha = c(0.5, 0.5), 
                   lambda = c(0.5, NA), 
                   bmode = bmode, 
                   x = x, 
                   psi = 1), 
              list(iSigma = matrix(c(1, 0, 0, 1), 2, 2), 
                   gamma = c(0.5, 0.5), 
                   theta = c(0, 0), 
                   alpha = c(0.5, 0.5), 
                   lambda = c(0.5, NA), 
                   bmode = bmode, 
                   x = x, 
                   psi = 1))

#make vector of parameters to track: add y
parameters <- c('Sigma', 'x', 'theta', 'gamma', 'alpha', 'bmode', 'psi')

#set path to model script
DCRWS.jags <- file.path("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/r_scripts/iSSF/state_space_model/DCRWS_rb.txt")

#run the model
testDCRWS.out <- jags.model(file = DCRWS.jags, 
                            data = test,
                            inits = inits, 
                            n.chains=2, 
                            n.adapt=5)
# clean up
rm(hseal.wb.dat, inits, parameters, y, idx, x, bmode)
hsealDCRWS.out<-rbind(hsealDCRWS.out[[1]],hsealDCRWS.out[[2]])
return(hsealDCRWS.out)


