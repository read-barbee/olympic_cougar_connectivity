require(chron)
require(rbugs)
require(dplyr)
# dat4bugs(hseal.dat, 1) ... 1 = timestep of 1 day, this line must be changed to:
#		dat4bugs(gs617.dat, 2) to run for grey seal 617 with timestep=2 days, and similarly for grey seal 2986.

hseal.dat <- read.delim("/Users/tb201494/Library/CloudStorage/Box-Box/Reference Literature/Spatial data sampling and management/Correcting for habitat bias/State Space Models/Jonsen Supplements/Supplement1/hseal.dat", header=TRUE, sep = ",")

hseal.wb.dat <- dat4bugs(hseal.dat, 1)
y <- hseal.wb.dat$y
idx <- hseal.wb.dat$idx
x <- matrix(unlist(lapply(1:(length(idx)-1), function(i){y[idx[i],]})), length(idx)-1, 2, byrow = T)
bmode <- rep(1, nrow(x) - 1)

inits <- list(list(iSigma = matrix(c(1, 0, 0, 1), 2, 2), gamma = c(0.5, 0.5), theta = c(0, 0), alpha = c(0.5, 0.5), lambda = c(0.5, NA), bmode = bmode, x = x, psi = 1), list(iSigma = matrix(c(1, 0, 0, 1), 2, 2), gamma = c(0.5, 0.5), theta = c(0, 0), alpha = c(0.5, 0.5), lambda = c(0.5, NA), bmode = bmode, x = x, psi = 1))

parameters <- c('Sigma', 'x', 'theta', 'gamma', 'alpha', 'bmode', 'psi')

DCRWS.bug <- file.path('~/data/SSM/ecology04Rev/DCRWS.txt')
hsealDCRWS.out <- rbugs(data = hseal.wb.dat, inits, parameters,
		DCRWS.bug, n.chains=2, n.iter=40000, n.burnin=20000, n.thin = 5, workingDir="/home/jonsen/.wine/winC/temp", bugsWorkingDir="c:/temp", bugs="c:/Program Files/WinBUGS14/WinBUGS14.exe", wine="/usr/local/bin/wine", useWine=TRUE, debug=F)
# clean up
rm(hseal.wb.dat, inits, parameters, y, idx, x, bmode)
hsealDCRWS.out<-rbind(hsealDCRWS.out[[1]],hsealDCRWS.out[[2]])
return(hsealDCRWS.out)


