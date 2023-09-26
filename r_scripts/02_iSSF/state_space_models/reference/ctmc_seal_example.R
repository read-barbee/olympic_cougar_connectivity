##
## Example of using a CTMC model for movement
##
## Steps:
##  1. Fit Quasi-Continuous Path Model to telemetry data (done using Buderman et al 2015)
##  2. Create covariate raster objects (the CTMC will be on the raster
##     grid cells)
##  3. Impute a quasi-continuous path
##  4. Turn quasi-continuous path into a CTMC discrete-space path
##  5. Turn discrete-space path into latent Poisson GLM format
##  6. Fit a Poisson GLM model to the data
##

#Import data

library(ctmcmove)

data(seal)
xyt=seal$locs[,3:1]
head(xyt)
plot(xyt[,1:2],type="b")
xy=xyt[,-3]
x=xyt[,1]
y=xyt[,2]
t=xyt[,3]

#########################################################################
##
## 1. Fit Quasi-Continuous Path Model
##
##########################################################################

library(fda)

## Define the knots of the spline expansion.
##
## Problems with fitting the functional movement model can often be fixed by
## varying the spacing of the knots.
knots = seq(min(t),max(t),by=1/4)
## create B-spline basis vectors used to approximate the path
b=create.bspline.basis(c(min(t),max(t)),breaks=knots,norder=3)
## define the sequence of times on which to sample the imputed path
tpred=seq(min(t),max(t),by=1/24/60)



## Fit latent Gaussian model using MCMC
out=mcmc.fmove(xy,t,b,tpred,QQ="CAR",n.mcmc=400,a=1,r=1,num.paths.save=30)
str(out)

## plot 3 imputed paths
plot(xy,type="b")
points(out$pathlist[[1]]$xy,col="red",type="l")
points(out$pathlist[[2]]$xy,col="blue",type="l")
points(out$pathlist[[3]]$xy,col="green",type="l")


#########################################################################
##
## 2. Creating rasters
##
##########################################################################

cov.df=seal$cov.df
str(cov.df)

#create sst raster
NN=sqrt(nrow(cov.df$X))
sst=matrix(seal$cov.df$X$sst,NN,byrow=TRUE)
sst=sst[NN:1,]
sst=raster(sst,xmn=min(seal$cov.df$X$x),xmx=max(seal$cov.df$X$x),
           ymn=min(seal$cov.df$X$y),ymx=max(seal$cov.df$X$y))


crs(sst)="+proj=longlat +datum=WGS84"
plot(sst)

#create chA raster
chA=matrix(seal$cov.df$X$chA,NN,byrow=TRUE)
chA=chA[NN:1,]
chA=raster(chA,xmn=min(seal$cov.df$X$x),xmx=max(seal$cov.df$X$x),
           ymn=min(seal$cov.df$X$y),ymx=max(seal$cov.df$X$y))
crs(chA)="+proj=longlat +datum=WGS84"

#create npp raster from pro covariate
pro=matrix(seal$cov.df$X$pro,NN,byrow=TRUE)
pro=pro[NN:1,]
npp=raster(pro,xmn=min(seal$cov.df$X$x),xmx=max(seal$cov.df$X$x),
           ymn=min(seal$cov.df$X$y),ymx=max(seal$cov.df$X$y))
crs(npp)="+proj=longlat +datum=WGS84"

#set all values of sst raster to 1
int=sst
values(int) <- 1

#create distance to rookery raster
d2r=int
rookery.cell=cellFromXY(int,xyt[1,1:2])
values(d2r)=NA
values(d2r)[rookery.cell]=0
d2r=distance(d2r)

#make raster stack
grad.stack=stack(sst,chA,npp,d2r)
names(grad.stack) <- c("sst","cha","npp","d2r")

plot(sst)
points(xyt[,1:2],type="b") #add points to sst raster

plot(grad.stack)


##########################################################################
##
## 3 Impute Quasi-Continuous Paths
##
##########################################################################

#number of paths to impute
P=20

plot(sst,col=grey.colors(100))
for(i in 1:P){
  points(out$pathlist[[i]]$xy,col=i,type="l",lwd=2)
}
points(xyt[,1:2],type="b",pch=20,cex=2,lwd=2)


##########################################################################
##
## 4. Turn continuous space path into a CTMC discrete space path
##
##########################################################################

path=out$pathlist[[1]]
ctmc=path2ctmc(path$xy,path$t,int,method="LinearInterp")
## alternate method, useful if you have impassible barriers, but slower
## ctmc=path2ctmc(path$xy,path$t,int,method="ShortestPath")

str(ctmc)


##########################################################################
##
## 5. Turn CTMC discrete path into latent Poisson GLM data
##
##########################################################################


loc.stack=stack(int,sst)
names(loc.stack) <- c("Intercept","sst.loc")

glm.list=list()
glm.list[[1]]=ctmc2glm(ctmc,loc.stack,grad.stack)

str(glm.list)

for(i in 2:P){
  cat(i," ")
  path=out$pathlist[[i]]
  ctmc=path2ctmc(path$xy,path$t,int,method="LinearInterp")
  glm.list[[i]]=ctmc2glm(ctmc,loc.stack,grad.stack)
}

## remove transitions that are nearly instantaneous
##  (These are essentially outliers in the following regression analyses)
for(i in 1:P){
  idx.0=which(glm.list[[i]]$tau<10^-5)
  if(length(idx.0)>0){
    glm.list[[i]]=glm.list[[i]][-idx.0,]
  }
  glm.list[[i]]$t=glm.list[[i]]$t-min(glm.list[[i]]$t)
}


##
## Stack the P imputations together
##

glm.data=glm.list[[1]]
for(i in 2:P){
  glm.data=rbind(glm.data,glm.list[[i]])
}

str(glm.data)

##########################################################################
##
## 6. Fit Poisson GLM
##    (here we are fitting all "M" paths simultaneously,
##     giving each one a weight of "1/M")
##
##########################################################################

fit.SWL=glm(z~cha+npp+sst+crw+d2r+sst.loc,
            weights=rep(1/P,nrow(glm.data)),family="poisson",offset=log(tau),data=glm.data)
summary(fit.SWL)

beta.hat.SWL=coef(fit.SWL)
beta.se.SWL=summary(fit.SWL)$coef[,2]



##########################################################################
##
## 6. Fit Poisson GLM
##    (here we are fitting using Multiple Imputation
##
##########################################################################

## Fit each path individually
glm.fits=list()
for(i in 1:P){
  glm.fits[[i]]=glm(z~cha+npp+sst+crw+d2r+sst.loc,
                    family="poisson",offset=log(tau),data=glm.list[[i]])
}

## get point estimates and sd estimates using Rubin's MI combining rules
beta.hat.mat=integer()
beta.se.mat=integer()
for(i in 1:P){
  beta.hat.mat=rbind(beta.hat.mat,coef(glm.fits[[i]]))
  beta.se.mat=rbind(beta.se.mat,summary(glm.fits[[i]])$coef[,2])
}

beta.hat.mat
beta.se.mat

## E(beta) = E_paths(E(beta|path))
beta.hat.MI=apply(beta.hat.mat,2,mean)
beta.hat.MI

## Var(beta) = E_paths(Var(beta|path))+Var_paths(E(beta|path))
beta.var.MI=apply(beta.se.mat^2,2,mean)+apply(beta.hat.mat,2,var)
beta.se.MI=sqrt(beta.var.MI)

cbind(beta.hat.MI,beta.se.MI)

##
## compare estimates from MI and Stacked Weighted Likelihood approach
##

## standardize regression coefficients by multiplying by the SE of the X matrix
sds=apply(model.matrix(fit.SWL),2,sd)
sds[1]=1

## plot MI and SWL regression coefficients
par(mfrow=c(1,2))
plot(beta.hat.MI*sds,beta.hat.SWL*sds,main="(a) Coefficient Estimates",
     xlab="Weighted Likelihood Coefficient",
     ylab="Multiple Imputation Coefficient",pch=20,cex=2)
abline(0,1,col="red")
plot(log(beta.se.MI),log(beta.se.SWL),
     main="(b) Estimated log(Standard Errors)",xlab="Weighted Likelihood log(SE)",
     ylab="Multiple Imputation log(SE)",pch=20,cex=2)
abline(0,1,col="red")


###########################################################################
##
## 6. (Alternate) We can use any software which fits Poisson glm data.
##    The following uses "gam" in package "mgcv" to fit a time-varying
##    effect of "d2r" using penalized regression splines.  The result
##    is similar to that found in:
##
##    Hanks, E.; Hooten, M.; Johnson, D. & Sterling, J. Velocity-Based
##    Movement Modeling for Individual and Population Level Inference
##    PLoS ONE, Public Library of Science, 2011, 6, e22795
##
###########################################################################

library(mgcv)

fit=gam(z~cha+npp+crw+sst.loc+s(t,by=-d2r),
        weights=rep(1/P,nrow(glm.data)),family="poisson",offset=log(tau),data=glm.data)
summary(fit)

plot(fit)
abline(h=0,col="red")

############################################################
##
## Overview Plot
##
############################################################


## pdf("sealfig.pdf",width=8.5,height=8.85)
par(mfrow=c(3,3))
##
plot(sst,col=(terrain.colors(30)),main="(a) Sea Surface Temperature")
points(xyt[1,1:2]-c(0,.05),type="p",pch=17,cex=2,col="red")
points(xyt[,1:2],type="b",pch=20,cex=.75,lwd=1)
##
plot(d2r/1000,col=(terrain.colors(30)),main="(b) Distance to Rookery")
points(xyt[1,1:2]-c(0,.05),type="p",pch=17,cex=2,col="red")
points(xyt[,1:2],type="b",pch=20,cex=.75,lwd=1)
##
image(sst,col=rev(terrain.colors(30)),main="(c) Imputed Functional Paths",xlab="",ylab="")
for(i in 1:5){
  ## points(out$pathlist[[i]]$xy,col=i+1,type="l",lwd=3)
  points(out$pathlist[[i]]$xy,col=i+1,type="l",lwd=2)
}
points(xyt[,1:2],type="p",pch=20,cex=.75,lwd=1)
##
ee=extent(c(188.5,190.5,58.4,59.1))
sst.crop=crop(sst,ee)
bg=sst.crop
values(bg)=NA
for(i in c(2)){
  values(bg)[cellFromXY(bg,out$pathlist[[i]]$xy)] <- 1
}
image(sst.crop,col=(terrain.colors(30)),xlim=c(188.85,190.2),
      ylim=c(58.5,59),main="(d) CTMC Path",xlab="",ylab="")
image(bg,col="blue",xlim=c(188.85,190.2),ylim=c(58.5,59),add=TRUE)
for(i in c(2)){
  points(out$pathlist[[i]]$xy,col=i,type="l",lwd=3)
}
points(xyt[,1:2],type="b",pch=20,cex=2,lwd=2)
##
image(sst.crop,col=(terrain.colors(30)),xlim=c(189.62,189.849),
      ylim=c(58.785,58.895),main="(e) CTMC Model Detail",xlab="",ylab="")
abline(v=189.698+res(sst)[1]*c(-1,0,1,2))
abline(h=58.823+res(sst)[2]*c(-1,0,1,2))
##
plot(fit,main="(f) Time-Varying Response to Rookery",shade=TRUE,
     shade.col="orange",lwd=3,rug=F,xlab="Day of Trip",
     ylab="Coefficient of Distance To Rookery")
abline(h=0,col="red")
##



###############################################
##
## Get UD (following Kenady et al 2017+)
##
###############################################

RR=get.rate.matrix(fit.SWL,loc.stack,grad.stack)
UD=get.UD(RR,method="lu")
ud.rast=sst
values(ud.rast) <- as.numeric(UD)
plot(ud.rast)


###############################################
##
## Get shortest path and current maps (following Brennan et al 2017+)
##
###############################################

library(gdistance)

## create a dummy transition layer from a raster.
## make sure the "directions" argument matches that used in path2ctmc
## also make sure to add the "symm=FALSE" argument
trans=transition(sst,mean,directions=4,symm=FALSE)
## now replace the transition object with the "rate" matrix
## so "conductance" values are "transition rates"
transitionMatrix(trans) <- RR
str(trans)

##
## now calculate least cost paths using "shortestPath" from gdistance
##

## pick start and end locations
plot(sst)
st=c(185,59.5)
en=c(190,57.3)

st.cell=cellFromXY(sst,st)
en.cell=cellFromXY(sst,en)

## shortest path
sp=shortestPath(trans,st,en,output="SpatialLines")
plot(sst,main="Shortest Path (SST in background)")
lines(sp,col="brown",lwd=7)



##
## Now calculate "current maps" that show space use of random walkers
## moving between two given locations.
##
## gdistance's "passage" function allows for asymmetric transition rates
##

passage.gdist=passage(trans,st,en,theta=.001,totalNet="net")
plot((passage.gdist))






