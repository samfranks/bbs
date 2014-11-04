########################################################################
#            DISTANCE SAMPLING OF BBS LINE TRANSECT DATA TO            #
#                 OBTAIN DENSITY AND ABUNDANCE ESTIMATES               #
########################################################################

# The core program below was written by Ali Johnston 2009 with some original 
# code from Eric Rexstad and modified by Dario Massimino 2011
# The programs make use of the R library mrds
# The data formatting code was written by Ali Johnston 2009

# IMPORTANT! Please type the name of the species and the visit here:
species <- "SC"
Visit <- "M"	# It can be E (early), L (late), or M (maximum)
min.year <- 1994
max.year <- 2011	# Final year (2 digits)

setwd("D:/Whinchat")

# These source functions to read and manipulate the bbs data:
source("D:/R functions/get.data.per.yr.per.square function.txt")
source("D:/R functions/format.data.all.yr.square function.txt")
source("D:/R functions/get.data.per.yr function.txt")
source("D:/R functions/format.data.all.yr function.txt")
source("D:/R functions/get.zero.obs.weather function.txt")
library(mrds)


#######################
# MODEL DETECTABILITY

# READING IN BIRD DATA
# Need transect observations of all transects in which there WERE observations.

miny <- substr(min.year,3,4)
maxy <- substr(max.year,3,4)
filename <- ifelse(substr(species,2,2)==".",
	paste(substr(species,1,1),"_ ",miny,maxy," tran zF.txt", sep=""),
	paste(species," ",miny,maxy," tran zF.txt", sep=""))

if (length(list.files("Data",pattern=filename))==1) {
	bird.dat <- read.table(paste("Data/",filename,sep=""),sep=",")
} else {
	print(Sys.time())
	print(paste(filename,"has not been found"))
	print("Data are now being read from the BBS folder")
	print("This may take very long (several hours for the most abundant species!)")
	wd <- getwd()
	setwd("Z:/bbs")
	bird.dat.xx <- format.data.all.yr(spec.code=species, min.year=min.year,
		max.year=max.year, zero.obs.squares=F)
	# It often gives me this warnings, but the data frame seems ok:
	# In readLines(file.loc) : incomplete final line found on 'data05/ebird05'
	setwd(wd)
	write.table(bird.dat.xx,paste("Data/",filename,sep=""),sep=",")
	bird.dat <- bird.dat.xx
	print(Sys.time())
}

# Use the top one for transects and the bottom one for squares
bird.dat2 <- subset(bird.dat,select=c("site","observer","detected","year","distance","distbegin","distend","effort","eorl","habs1","tran.sec"))
visit <- ifelse(bird.dat2$eorl=="E",1,0)	# 1=Early;  0=Late
bird.dat2 <- cbind(bird.dat2,visit)

# Check for wrong effort.  Anything over 0.4 is wrong.  This is 2km*100m*2 = 2*0.1*2 km2
# Often the effort is wrong because the square is entered twice in the habitat file, so the real
# effort is half of the recorded effort.  Should check these, but for now is fudged:
sub1 <- subset(bird.dat2,effort>0.4)
sub2 <- subset(bird.dat2,effort<=0.4)
sub1$effort <- sub1$effort/2
bird.dat2 <- rbind(sub1,sub2)

# Values of the year as actual numbers (e.g. 1996) make convergence hard as they are 
# such big numbers.  So scale to be 1=min.year:
bird.dat2$year <- as.numeric(as.character(bird.dat2$year))
bird.dat2$year <- bird.dat2$year-min.year+1

# Create an 'object' column, with a unique index for each row.
bird.dat2$object <- seq(1,nrow(bird.dat2))

# The procedure in mrds assumes that the 'observer' field is used for double observer
# transects, so change the observer field to a vector of 1s, and put the observer code in
# a different field 'observer.code'.
bird.dat2$observer.code <- bird.dat2$observer
bird.dat2$observer <- rep(1,length(bird.dat2$observer))

# Remove any habitat types which don't have a big enough sample size of observations.
# Also remove any observations which don't have a habitat type recorded.
# Use table(bird.dat2$habs1) to check the sample sizes.

bird.dat2 <- bird.dat2[bird.dat2$habs1%in% c("A","B","C","D","E","F","G","H","I"),]
bird.dat2$habs1 <- bird.dat2$habs1[drop=T]

# Removing any with less than 20:
table.hab <- table(bird.dat2$habs1)
bird.dat2$habs1 <- as.character(bird.dat2$habs1)
bird.dat2$habs2 <- bird.dat2$habs1

if ((table.hab["H"]<=20|is.na(table.hab["H"])) | (table.hab["I"]<=20|is.na(table.hab["I"])) | (table.hab["G"]<=20|is.na(table.hab["G"]))) {
	for(i in 1:nrow(bird.dat2)){
		bird.dat2$habs2[i] <- ifelse(bird.dat2$habs1[i] %in% c("G","H","I"),"GHI",bird.dat2$habs1[i])
	}
}
	if ((table.hab["C"]<=20|is.na(table.hab["C"]<=20)) | (table.hab["D"]<=20|is.na(table.hab["D"]))) {
	for(i in 1:nrow(bird.dat2)){
		bird.dat2$habs2[i] <- ifelse(bird.dat2$habs1[i] %in% c("C","D"),"CD",bird.dat2$habs2[i])
	}
}
	if ((table.hab["A"]<=20|is.na(table.hab["A"]<=20)) | (table.hab["B"]<=20|is.na(table.hab["B"]<=20))) {
	for(i in 1:nrow(bird.dat2)){
		bird.dat2$habs2[i] <- ifelse(bird.dat2$habs1[i] %in% c("A","B"),"AB",bird.dat2$habs2[i])
	}
}

if (any(table(bird.dat2$habs2)<=20)) print ("There are still habitats <=20")

bird.dat2$habs1.old <- bird.dat2$habs1
bird.dat2$habs1 <- bird.dat2$habs2
bird.dat2$habs1 <- as.factor(bird.dat2$habs1)

bird.dat3 <- bird.dat2
bird.dat3$habs1 <- bird.dat3$habs1[drop=T]

bird.dat8 <- bird.dat3	# Jumped to bird.dat8 as we don't consider speed, etc.

#n<-function(x){ as.numeric(as.character(x))}
#bird.dat8$vis<-n(bird.dat8$vis)

bird.dat8 <- subset(bird.dat8, !is.na(habs1))

# Makes one line for each individual observation:

bird.dat9 <- bird.dat8
bird.dat8 <- bird.dat8[-c(1:nrow(bird.dat8)),]

bird.dat.null <- bird.dat8
print(Sys.time())
batch <- 0
if (nrow(bird.dat9)>10000) {
	for (batch in 1:(floor(nrow(bird.dat9)/10000))){
		print (paste("Processing batch number",batch,"of",floor(nrow(bird.dat9)/10000)+1))
		linecounter <- 1
		temp <- bird.dat.null
		for(i in 1:10000){	
			n <- bird.dat9$detected[(batch-1)*10000+i]
			temp[linecounter:(linecounter+n-1),] <- bird.dat9[(batch-1)*10000+i,]
			linecounter <- linecounter+n
		} # closes i
		bird.dat8 <- rbind(bird.dat8,temp)
	} # closes batch
} # closes if
batch <- batch+1
print (paste("Processing batch number",batch,"of",floor(nrow(bird.dat9)/10000)+1))
print (Sys.time())
linecounter <- 1
temp <- bird.dat.null
for(i in 1:(ifelse(nrow(bird.dat9)>10000, nrow(bird.dat9)%%((batch-1)*10000), nrow(bird.dat9)))){	
	n <- bird.dat9$detected[(batch-1)*10000+i]
	temp[linecounter:(linecounter+n-1),] <- bird.dat9[(batch-1)*10000+i,]
	linecounter <- linecounter+n
}
bird.dat8 <- rbind(bird.dat8,temp)
	
bird.dat8$detected <- rep(1,nrow(bird.dat8))
bird.dat8$size <- rep(1,nrow(bird.dat8))
bird.dat8$object <- seq(1,nrow(bird.dat8),by=1)
bird.dat8$habs1 <- as.character(bird.dat8$habs1)

# gets rid of 9s (code for missing values)
bird.dat8 <- subset(bird.dat8, visit!=9)


########################
# RUNNING THE MODELS
#
# As the purpose is running these models for a wide range of species, we need to keep
# computing time within reasonable limits.
# For this reason we decided to fit only one models with visit and habitat
# as predictors.

Sys.time()
ddf.visit.habs1 <- ddf(~mcds(key="hn",formula=~as.factor(habs1)+visit),
	data=bird.dat8, method="ds")
Sys.time()

final <- ddf.visit.habs1

save.image("Data/temporary rescue point.Rdata")
load("Data/temporary rescue point.Rdata")


###############################################################
# ESTIMATE DETECTABILITY AT THE SQUARE LEVEL ON ALL SQUARES
# Need the program "final" (a distance sampling model output) and the counts from
# each transect section of each square with the habitat for each transect section.
# This is used to sum the counts from each habitat, on each square.

# For the prediction step if the bird is a resident, only the first visit is used.
# And if it is a migrant, only the second visit is used.
# Therefore, if "visit" is a variable in the model "final", be careful to predict
# detecability in squares based only on the visit type which is going to be used
# for the abundance estimation.  e.g for residents predict detectability in the 
# first visit of the season.

filename <- ifelse(substr(species,2,2)==".",
	paste(substr(species,1,1),"_ ",miny,maxy," squr zT.txt", sep=""),
	paste(species," ",miny,maxy," squr zT.txt", sep=""))

if (length(list.files("Data",pattern=filename))==1) {
	bird.dat.xx.sz<-read.table(paste("Data/",filename,sep=""), sep=",", header=T)
} else {
	print(Sys.time())
	print(paste(filename,"has not been found"))
	print("Data are now being read from the BBS folder")
	print("This may take very long(several hours for the most abundant species!)")
	wd <- getwd()
	setwd("Z:/bbs")
	bird.dat.xx.sz <- format.data.all.yr.square(spec.code=species,
	min.year=min.year, max.year=max.year, zero.obs.squares=T)
	# It often gives me this warnings, but the data frame seems ok:
	# In readLines(file.loc) : incomplete final line found on 'data05/ebird05'
	setwd(wd)
	write.table(bird.dat.xx.sz,paste("Data/",filename,sep=""),sep=",")
	print(Sys.time())
}

# Run the following to get the square-visit-level information:
# Include zero squares for the prediction step
# Need to work out how many observation there were in each habitat per square:

ag.t <- aggregate(bird.dat$detect, by=list(bird.dat$site, bird.dat$year, bird.dat$eorl), mean)
sites <- as.character(ag.t$Group.1)
years <- ag.t$Group.2
eorls <- as.character(ag.t$Group.3)
habAdt<-habBdt<-habCdt<-habDdt<-habEdt<-habFdt<-habGdt<-habHdt<-habIdt<-vector(length=length(sites))
hab.det <- as.data.frame(cbind(habAdt,habBdt,habCdt,habDdt,habEdt,habFdt,habGdt,habHdt,habIdt,sites,eorls,years))
hab.det[,c(1:9)] <- matrix(nrow=length(sites),ncol=9,rep(0,9*length(sites)))
habs <- c("A","B","C","D","E","F","G","H","I")

for(i in 1:nrow(hab.det)){
	sub<-subset(bird.dat,site==sites[i] & eorl==eorls[i] & year==years[i])
	if(nrow(sub)>0){
		for(j in 1:length(habs)){
			g<-grep(habs[j],sub$habs1)
			if(length(g)>0){
				hab.det[i,j]<-sum(sub$detected[g])
			}	# close if
		}	# close j
	}	# close if
} 	# close i

hab.det$siteeorl <- paste(hab.det$sites,hab.det$eorl,hab.det$year,sep="")

# Match up the habitat counts with the square level information:

# Don't have the no of transects info for some squares:
bird.s2<-subset(bird.dat.xx.sz,!is.na(habA)) 

## Chose whether "E"arly or "L"ate visits:
## "E" if resident; "L" if migratory.
#if (Visit=="E"|Visit=="L")	bird.s3 <- subset(bird.s2,eorl==Visit)
#if (Visit=="M") {
#	bird.s2.1 <- aggregate(bird.s2$detected, by=list(bird.s2$site, bird.s2$year, bird.s2$dist.band), FUN=max)
#	names(bird.s2.1)<-c("site","year","dist.band","detected")
#	bird.s3<-bird.s2.1
#}

#bird.s3<-subset(bird.s3,dist.band==1)
bird.s3<-subset(bird.s2,dist.band==1)

# Merge the two datasets and delete the variable "detected" to avoid confusion
bird.s3$siteeorl <- paste(bird.s3$site, bird.s3$eorl, bird.s3$year, sep="")
bird.s4 <- merge(bird.s3,hab.det,by="siteeorl",all.x=T)
bird.s4 <- subset(bird.s4, select=-detected)

# Fill in zeros where no observations were made:
null.sites<-grep(TRUE,is.na(bird.s4$sites))
bird.s4[null.sites,c("habAdt","habBdt","habCdt","habDdt","habEdt","habFdt","habGdt","habHdt","habIdt")]<-0

# Amalgamate habitats which are from now on "the same":
# This will be different for each species!
# Also, be careful not to run more than once, otherwise the estimates are out.
#bird.s4$habH<-bird.s4$habH+bird.s4$habI
#bird.s4$habI <- 0
#bird.s4$habHdt<-bird.s4$habHdt+bird.s4$habIdt
#bird.s4$habIdt <- 0

bird.s5<-bird.s4

# Specify the pooled habitats:

#        A   B   C   D   E   F   G   H    I  
habs<-c("A","B","C","D","E","F","G","H","I")

if (any(names(table(bird.dat8$habs1))=="GHI")) {
	habs[c(7,8,9)]<-c("GHI","GHI","GHI")
	}
if (any(names(table(bird.dat8$habs1))=="CD")) {
	habs[c(3,4)]<-c("CD","CD")
	}
if (any(names(table(bird.dat8$habs1))=="AB")) {
	habs[c(1,2)]<-c("AB","AB")
	}

# Work out the total detections and total transects for each square:
st<-grep("habA",colnames(bird.s4))
st1<-st[1]-1
st2<-st[2]-1

tot.tran<-tot.det<-vector(length=nrow(bird.s5))
for(i in 1:nrow(bird.s5)){
  tot.tran[i]<-sum(bird.s5[i,c(st1+1:9)])
  tot.det[i]<-sum(bird.s5[i,c(st2+1:9)])
}

bird.s5<-cbind(bird.s5,tot.tran,tot.det)
bird.s5<-subset(bird.s5,tot.tran<11)

bird.s5<-subset(bird.s5,tot.tran>0)

# create visit and speed in the same way when we ran the model
bird.s5$visit<-bird.s5$speed<-NA
bird.s5$visit<-ifelse(bird.s5$eorl=="E", 1,0)

# Temporary rescue point 2
save.image("Data/temporary rescue point 2.Rdata")
#load("Data/temporary rescue point 2.Rdata")

# The following loop takes about 20 mins
habs1 <- as.vector(names(table(bird.dat8$habs1)))
fit.vec<-vector()
Sys.time()
for(i in 1:nrow(bird.s5)){
	sub <- bird.s5[i,]
	visit <- rep(sub$visit,length(habs1))
	nd <- as.data.frame(cbind(habs1,visit))
	nd$visit <- as.numeric(as.character(nd$visit))
	fitted <- predict.ds(final, newdata=nd, compute=T)$fitted[,1]
	fit <- as.data.frame(cbind(habs1,fitted))
	colnames(fit)[1]<-"habs"
	# Match up pooled habitats:
	ha <- as.data.frame(cbind(1:9,habs))
	ha2 <- merge(ha,fit,by="habs",all.x=T)
	fit.hab <- as.numeric(as.character(ha2$fitted))
	# Find which habitats are in the square and in what nos:
	g1 <- grep(TRUE,bird.s5[i,st1+1:9]>0)
	# Calculate the counts in each habitat and the weights:
	ct <- bird.s5[i,st2+g1]
	wt <- ct+bird.s5[i,st1+g1]
	# Calculate the square detectability:
	fit.vec[i] <- sum(fit.hab[g1]*wt)/sum(wt)
	if (i%%10000==0) print(paste(i,"of",nrow(bird.s5)))
}
bird.s6 <-cbind(bird.s5,fit.vec)

save.image("temporary rescue point whinchat.Rdata")

Sys.time()

# THIS MUST BE CHECKED
#if (Visit=="E"|Visit=="L") {
#	birds<-subset(bird.s6,eorl==Visit, select=c("site","year","eorl","tot.det","fit.vec","birds"))
#}	

if (Visit=="M") {
	#bird.s6$est.birds <- bird.s6$tot.det/bird.s6$fit.vec
	bird.s6$eorl <- as.character(bird.s6$eorl)
	bird.s6.E <- subset(bird.s6, eorl=="E")
	bird.s6.L <- subset(bird.s6, eorl=="L")
	bird.s7 <- merge(bird.s6.E,bird.s6.L,by=c("site","year"))
	bird.s7 <- subset(bird.s7, select=c("site","year","eorl.x","tot.det.x","fit.vec.x","eorl.y","tot.det.y","fit.vec.y"))
	# compare birds detected during early and late visit and set best to:
	# 0 if same number; 1 if more birds for early visit; 2 for late visit
	best <- ifelse(bird.s7$tot.det.x==bird.s7$tot.det.y, 0, 2)
	best <- ifelse(bird.s7$tot.det.x>bird.s7$tot.det.y, 1, best)
	birds <- bird.s7[,1:2]
	birds$eorl <- ifelse(best==0, "B", "L")
	birds$eorl <- ifelse(best==1, "E", birds$eorl)
	birds$tot.det <- ifelse(best==1, bird.s7$tot.det.x, bird.s7$tot.det.y)
	birds$fit.vec <- ifelse(best==0, (bird.s7$fit.vec.x+bird.s7$fit.vec.y)/2, bird.s7$fit.vec.y)
	birds$fit.vec <- ifelse(best==1, bird.s7$fit.vec.x, birds$fit.vec)
}

imagename <- paste("Data/",species," Step 1.Rdata",sep="")
save.image(imagename)
#load(imagename)

filename <- paste("Data/",species,"detect.csv",sep="")
write.csv(birds, filename, row.names=F)

print(paste("End of Step 1 for ",species,". Visit: ",Visit,sep=""))


