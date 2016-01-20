#pedigree and kinship for Jays - for building kinship matrix for real (I hope)
#July 2015 

setwd("C:/Users/Erin/Dropbox/Jay survival analyses")

library(kinship2)

#upload file
demo.ped <- read.csv(file="Ped_withlastobs.csv")

#change column names
names(demo.ped)[names(demo.ped)=="Demo_Pedigree.Jays.USFWBand"] <- "ID"
names(demo.ped)[names(demo.ped)=="Jays_1.USFWBand"] <- "MomID"
names(demo.ped)[names(demo.ped)=="Jays_2.USFWBand"] <- "DadID"

#convert dates to date format
demo.ped$LastObsDate<- as.Date(demo.ped$LastObsDate, format = "%m/%d/%Y")

#add censorship status column 
demo.ped["censorship"] <- NA
demo.ped$censorship <- 1

#add 0's to those still alive 2015-03-11
demo.ped$censorship[which(demo.ped$LastObsDate=="2015-03-11")]<-0
demo.ped$censorship[which(demo.ped$LastObsDate=="2015-04-07")]<-0

#add sex indicator column
demo.ped["sexind"] <- 3
demo.ped$sexind[which(demo.ped$Sex=="M")] <- 1
demo.ped$sexind[which(demo.ped$Sex=="F")] <- 2

#add affected column (no meaning here, just needed for function)
demo.ped["aff"] <- 0

#have to do something with birds that have no parents listed (not in ID column
# but in dad or mom id column, aka founders)
demo.ped["missid"] <- "FALSE"
demo.ped$missid[which(demo.ped$JayID=="LL-S")] <- "TRUE"
demo.ped$missid[which(demo.ped$JayID=="SGY-")] <- "FALSE"
demo.ped$missid[which(demo.ped$JayID=="-ORS")] <- "TRUE"
demo.ped$missid[which(demo.ped$JayID=="SY-P")] <- "TRUE"
demo.ped$missid[which(demo.ped$JayID=="-SYW")] <- "TRUE"


jay.ped <- pedigree(id=demo.ped$ID, dadid=demo.ped$MBreeder, momid=demo.ped$FBreeder, sex=demo.ped$sexind, affected=demo.ped$aff, status=demo.ped$censorship, missid = demo.ped$missid)



