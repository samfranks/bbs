##################################################################
#        PREPARES THE DATASET TO BE USED FOR GLMs OR GAMs        #
##################################################################

# The program below was written by Ali Johnston 2009 and modified by
# Dario Massimino 2011

setwd("D:/Whinchat")

species <- "SC"

filename <- paste("Data/",species,"detect.csv",sep="")
birds <- read.csv(filename)

# This assumes Step 1 has already been run and a final model
# (ddf object) has been output, and is called "final".
# It also assumes that a final dataset called "birds" has been output.
# This should have a column called "fit.vec" which contains
# the (weighted) average detectability for the square.

memory.limit(size=4000)

print(paste("Processing",species))

spat1 <- subset(birds, select=c("site","eorl","year","tot.det","fit.vec"))
colnames(spat1)[colnames(spat1)=="tot.det"] <- "count"

# Read land cover and elevation
landa.elev<-read.csv("D:/Data/land_elev_uk.csv")
landa.elev<-subset(landa.elev, elev!=-9999&!is.na(county))
spat2<-merge(spat1,landa.elev, by="site")

# Calculate the weights by the BBS regions
rdasqs<-read.table("D:/Data/rdasqs.txt")
rdasqs<-rdasqs[,c(1,6)]
names(rdasqs)<-c("site","bbsreg")

spat3<-merge(spat2,rdasqs,by="site",all.x=T)

wt<-table(rdasqs$bbsreg)/table(spat3$bbsreg)
wts<-as.data.frame(cbind(names(wt),wt))
colnames(wts)<-c("bbsreg","wt")

spat4<-merge(spat3,wts,by="bbsreg",all.x=T)
spat4$wt<-as.numeric(as.character(spat4$wt))

spat5<-spat4

spat5$island<-"Mainland"
spat5$island[spat5$county=="Shetland"]<-"Shetland"
spat5$island[spat5$county=="Orkney"]<-"Orkney"
spat5$island[spat5$county=="Western Isles"]<-"WestIsl"
spat5$island[spat5$county=="Antrim"]<-"NorthIrl"
spat5$island[spat5$county=="Armagh"]<-"NorthIrl"
spat5$island[spat5$county=="Belfast"]<-"NorthIrl"
spat5$island[spat5$county=="Down"]<-"NorthIrl"
spat5$island[spat5$county=="Fermanagh"]<-"NorthIrl"
spat5$island[spat5$county=="Londonderry"]<-"NorthIrl"
spat5$island[spat5$county=="Tyrone"]<-"NorthIrl"
spat5$island[spat5$county=="Isle of Man"]<-"Man"
spat5$island[spat5$county=="Isles of Scilly"]<-"Scilly"
	
# Dataset ready for GAMs!

filename<-paste("Data/Dataset_",species,".csv", sep="")
write.csv(spat5, filename)

