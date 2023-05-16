#SUPPLEMENTAL MATERIAL
#the script used to screen GPS location data in paper: Bjørneraas et al. 2010
#the script is made for the statistical software R (R Development Core Team 2009) and depends on the R-package Adehabitat (Calenge 2006)
#the animal locations must be defined by a projected coordinate system
#the function does not take into account the spherical surface of the earth, thus, the 
#script is not valid at a global scale, but at smaller spatial scales 
#we refer to the paper for further details on the screening procedure

#R Development Core Team. 2009. R: a language and environment for statistical
#computing. R Foundation for Statistical Computing, Vienna, Austria.
#Calenge, C. 2006. The package "adehabitat" for the R software: a tool for the analysis of
#space and habitat use by animals. Ecological Modelling 197:516–519.
###########################################################################

#the screening method depends on the 2 functions below: GPS.screening and 
#GPS.screening.wrp, which can be loaded by running the scripts in R

#the function GPS.screening applies the algorithm explained in the paper
#the screening occurs in 2 rounds:
#first, remove all points that are far off
#second, detect spikes: high outgoing and incoming speed and a sharp return

#the function GPS.screening.wrp prepares the data appropriately. It takes animal id's, 
#x-coordinates, y-coordinates and dates as a vector
#ensure that all these vectors are of the same length
#there is no restriction to the number of animals (i.e., different animal ids)

###format of those vectors has to be:
###	id: numeric, character or factor
###	x: numeric
###	y: numeric
###	date: POSIXct

#for the use of other criteria (see below and in the paper) you can change the values below:
#medcrit: detla, the cutoff value for the median, which is used in the first step from the first 	round
#meancrit: mu, the cutoff value for the mean, which is used in the second step from the 	first round
#spikesp: alpha, the criterion for the speed of the outgoing and incoming move, which is 	used in the second round
#spikecos: theta, the criterion for the cosine of the angle associated with the outgoing and 	incoming move, which is used in the second round

###########################################################################

### Define Parameters

medcrit = ""

meancrit = ""

spikesp = ""

spikecos = ""




###STEP 1: Run the function GPS.screening 

GPS.screening <- function(x, medcrit, meancrit, spikesp, spikecos){
  ###round one: #############
  ###remove points far off:
  RoundOne <- function(x, df, win, MM){
    if (x<=win){
      dfT <- df[c(1:2*win),]
    }
    if (x>(nrow(df)-win)){
      dfT <- df[c((nrow(df)-2*win):(nrow(df))),]
    }
    if (x>win & x<=(nrow(df)-win)){
      dfT <- df[c((x-win):(x+win)),]
    }
    if (MM=="mean"){
      temp <- sqrt((df$x[x]-mean(dfT$x, na.rm=T))^2+(df$y[x]-mean(dfT$y, na.rm=T))^2)
      #calculate the distance between the mean and the focal point
    }
    if (MM=="median"){
      temp <- sqrt((df$x[x]-median(dfT$x, na.rm=T))^2+(df$y[x]-median(dfT$y, na.rm=T))^2)
      #calculate the distance between the median and the focal point
    }
    return(temp)
  }
  
  ##get the points that are really far off, by using the median
  x$R1dmed <- unlist(lapply(c(1:nrow(x)), RoundOne, df=x, win=10, MM="median"))
  
  ##get the points that are rather far off, by using the mean
  x$R1dmean <- meancrit+1
  x$R1dmean[x$R1dmed<medcrit] <- unlist(lapply(c(1:nrow(x[x$R1dmed<medcrit,])),  RoundOne, df=x[x$R1dmed<medcrit,], win=10, MM="mean"))
  x$R1error <- ifelse(x$R1dmean>meancrit, TRUE, FALSE)
  
  ###round two:#############
  
  ##find spikes
  xT <- x[x$R1error==FALSE,]
  library(adehabitat)
  xT <- as.ltraj(xT[,c("x","y")], date=xT$date, id=xT$id)[[1]]
  
  x$R2error <- NA
  x$R2error[x$R1error==FALSE] <- ((xT$dist/xT$dt*3600)>spikesp & c(NA, (xT$dist/xT$dt*3600)[-nrow(xT)])>spikesp & cos(xT$rel.angle)<(spikecos))
  return(x)
}

###########################################################################

###STEP 2: Run the function GPS.screening.wrp 

GPS.screening.wrp <- function(id, x, y, da, medcrit, meancrit, spikesp, spikecos){
  mydata <- traj2df(as.traj(id = id, xy = data.frame(x=x,y=y), date = da))
  mydata <- split(mydata, mydata$id)
  mydata <- lapply(mydata, GPS.screening, medcrit=medcrit, meancrit=meancrit, spikesp=spikesp, spikecos=spikecos)
  mydata <- do.call("rbind", mydata)
  return(data.frame(id=mydata$id, x=mydata$x, y=mydata$y, date=mydata$date, R1dmed=mydata$R1dmed, R1dmean=mydata$R1dmean, R1error=mydata$R1error, R2error=mydata$R2error))
}

###########################################################################

###STEP 3: Screen the data

mydata <- GPS.screening.wrp(mydata$id, mydata$x, mydata$y, mydata$date, medcrit=100000, meancrit=10000, spikesp=1500, spikecos=(-0.97))

head(mydata)

### if an error message occur, the problem may be that the time zone should be set to UTC 	(i.e., tz = "UTC") for the column date

### the criteria are sat as:  = 100,000 m;  = 10,000 m;  = 1500 m/h; -0.97, 
###	which is the same as for the moose example in the paper, but can be changed. 

######################################################################################################################
###example: remove one "#" from the beginning of every sentence to see how it works

# ###create a trajectory
#set.seed(0)
#demo <- simm.crw(date=as.POSIXct(c(1:100)*3600, origin="1970-01-01"), h=10, r=0.2, id="bob")
#plot(demo)
#demo <- traj2df(ltraj2traj(demo))
#demo <- demo[, c("id", "x", "y", "date")]

# ###introduce some errors
# ##2 consecutive points real far off
#demo$x[5:6] <- demo$x[5:6]+500000
#demo$y[5:6] <- demo$y[5:6]+500000

# ##a spike
#demo$x[49] <- demo$x[50]
#demo$y[49] <- demo$y[50]
#demo$x[51] <- demo$x[50]
#demo$y[51] <- demo$y[50]
#demo$x[50] <- demo$x[50]+1500
#demo$y[50] <- demo$y[50]+1500

#plot(demo$x, demo$y)

# ##apply screening script:
#demo <- GPS.screening.wrp(demo$id, demo$x, demo$y, demo$date, medcrit=100000, 		meancrit=10000, spikesp=1500, spikecos=(-0.97))
#demo[c(1:9, 45:54),]

