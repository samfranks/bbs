###############################
##                           ##
##    GAMs AND TREND MAPS    ##
##                           ##
## Dario Massimino (2011/12) ##
##                           ##
###############################

# Specify the first and the last years
first.years<-c(1994, 1996)
last.years<-c(2007, 2009)

setwd("D:/Mapping trends 2011")
library(mgcv)
library(MASS)
library(maptools)

min.BBSsqs <- 400

species.list<-read.csv("data/Species' list.csv")
species.list<-subset(species.list, BBSsqs>=min.BBSsqs)
species.codes<-as.character(species.list$Code)

# Read land cover data
land<-read.csv("Data/Land_elev_UK.csv")
land<-subset(land, !is.na(county)&elev!=-9999, select=-X)

# Variable "island"
land$island<-"Mainland"
land$island[land$county=="Shetland"]<-"Shetland"
land$island[land$county=="Orkney"]<-"Orkney"
land$island[land$county=="Western Isles"]<-"WestIsl"
land$island[land$county=="Antrim"]<-"NorthIrl"
land$island[land$county=="Armagh"]<-"NorthIrl"
land$island[land$county=="Belfast"]<-"NorthIrl"
land$island[land$county=="Down"]<-"NorthIrl"
land$island[land$county=="Fermanagh"]<-"NorthIrl"
land$island[land$county=="Londonderry"]<-"NorthIrl"
land$island[land$county=="Tyrone"]<-"NorthIrl"
land$island[land$county=="Isle of Man"]<-"Man"
land$island[land$county=="Isles of Scilly"]<-"Scilly"
land$island[land$island=="Orkney"] <- "Orkney"
land$island[land$island=="Shetland"] <- "Shetland"

### GAMs

#for (species in species.codes) {
for (species in c("MG")) {

	datasetname <- paste("Data/Dataset_",species,".csv",sep="")
	dataset <- read.csv(datasetname)

	print(paste("Processing",species))
	print(Sys.time())

	# Normalise east, north, elev and year between 0 and 1
	dataset.01<-dataset
	dataset.01$east01<-(dataset.01$easting-min(dataset.01$easting))/(max(dataset.01$easting)-min(dataset.01$easting))
	dataset.01$north01<-(dataset.01$northing-min(dataset.01$northing))/(max(dataset.01$northing)-min(dataset.01$northing))
	dataset.01$elev01<-(dataset.01$elev-min(dataset.01$elev))/(max(dataset.01$elev)-min(dataset.01$elev))
	dataset.01$year01<-(dataset.01$year-min(dataset.01$year))/(max(dataset.01$year)-min(dataset.01$year))

	# Shetland and Orkney as same factor
	#dataset.01$island <- as.character(dataset.01$island)
	#dataset.01$island[dataset.01$island=="Orkney"] <- "NorthIsl"
	#dataset.01$island[dataset.01$island=="Shetland"] <- "NorthIsl"
	#dataset.01$island <- as.factor(dataset.01$island)

	# List islands where the species has not been detected in either periods
	island0 <- aggregate(dataset.01$count, by=list(dataset.01$island), sum)
	island0 <- as.character(island0$Group.1[island0$x==0])

	runif(1)	# Sometimes .Random.seed is not found so this line creates one

	model.t0<-gam(count~s(north01,east01,elev01, k=1)+MHB+BrlfWood+ConWood+ImpGrass+
		+SemNatGr+Arable+Human+island, family=quasipoisson, offset=log(fit.vec),
		weights=wt, subset=(year>=first.years[1]&year<=first.years[2]), data=dataset.01)

	model.t1<-gam(count~s(north01,east01,elev01, k=1)+MHB+BrlfWood+ConWood+ImpGrass+
		+SemNatGr+Arable+Human+island, family=quasipoisson, offset=log(fit.vec),
		weights=wt, subset=(year>=last.years[1]&year<=last.years[2]), data=dataset.01)

	r2 <- round(c(summary(model.t0)$r.sq,summary(model.t1)$r.sq),2)
	devexpl <- round(c(summary(model.t0)$dev.expl,summary(model.t1)$dev.expl),2)

	### Predicting values

	land.pred <- land
	land.pred$fit.vec <- 1

	# Normalise coordinates between 0 and 1
	land.pred$north01 <- (land$northing-min(dataset.01$northing))/(max(dataset.01$northing)-min(dataset.01$northing))
	land.pred$east01 <- (land$easting-min(dataset.01$easting))/(max(dataset.01$easting)-min(dataset.01$easting))
	land.pred$elev01 <- (land$elev-min(dataset.01$elev))/(max(dataset.01$elev)-min(dataset.01$elev))

	model<-list(model.t0,model.t1)
	for (t01 in c(1,2)) {
		pred <- predict.gam(model[[t01]], newdata=land.pred, type="response", se.fit=T)
		pred$fit <- as.numeric(signif(pred$fit,3))
		pred$se.fit <- as.numeric(signif(pred$se.fit,3))
		land.pred <- data.frame(land.pred, pred$fit, pred$se.fit)
	}

	names(land.pred)[names(land.pred)=="pred.fit"]<-"pred0"
	names(land.pred)[names(land.pred)=="pred.fit.1"]<-"pred1"
	names(land.pred)[names(land.pred)=="pred.se.fit"]<-"se0"
	names(land.pred)[names(land.pred)=="pred.se.fit.1"]<-"se1"
	
	# Convert into densities/km2 and calculate the difference
	land.pred$dens0 <- land.pred$pred0*2.5
	land.pred$se.dens0 <- land.pred$se0*2.5
	land.pred$dens1 <- land.pred$pred1*2.5
	land.pred$se.dens1 <- land.pred$se1*2.5
	land.pred$ddens<-land.pred$dens1-land.pred$dens0
	land.pred$relchange<-(land.pred$dens1-land.pred$dens0)/land.pred$dens0
	land.pred$maxdens<-pmax(land.pred$dens1,land.pred$dens0)
	land.pred$island0 <- as.numeric(!land.pred$island %in% island0)
	land.pred <- subset(land.pred, select=-c(Sea,Coastal,MHB,BrlfWood,ConWood,ImpGrass,SemNatGr,Arable,Human,Unclass,county,country,pred0,se0,pred1,se1))

	# Prepare data frame with 10-km resolution
	km10 <- subset(land.pred, select=c("easting","northing","InWater","dens0","dens1","se.dens0","se.dens1","ddens","relchange","maxdens","island","island0"))
	km10$easting10km <- round(km10$easting/10000)
	km10$northing10km <- round(km10$northing/10000)
	km10$count <- 1
	sums <- aggregate (km10$count, by=list(km10$easting10km, km10$northing10km), sum)
	means <- aggregate (subset(km10, select=-c(easting10km,northing10km,island)), by=list(km10$easting10km, km10$northing10km, km10$island), mean)	
	colnames(means)[1:3] <- c("easting10km","northing10km","island")
	land.pred.10km <- means[sums$x>=10 & means$InWater<50,]
	land.pred.10km$easting <- land.pred.10km$easting10km*10000
	land.pred.10km$northing <- land.pred.10km$northing10km*10000
	land.pred.10km <- subset(land.pred.10km, select=c("easting","northing","dens0","dens1","se.dens0","se.dens1","ddens","relchange","maxdens","island","island0"))

	# Append r2 and explained deviance to the csv file
	write.table(data.frame(species,r2[1],r2[2],devexpl[1],devexpl[2]), "Results/Model fit.txt", append=T, row.names=F, col.names=F)

	rm(km10, dataset, sums, means, pred, t01, r2, devexpl)

	filename <- paste("Data/",species," final model.Rdata",sep="")

	save.image(filename)
	#load(filename)
}

print("Done.")