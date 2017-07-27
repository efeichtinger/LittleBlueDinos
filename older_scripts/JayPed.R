### Picking up June 2016 - although ideally I want this all on the same script file

#pedigree and kinship for Jays - for building kinship matrix for real (I hope)
#July 2015 

library(kinship2)

#upload file
demo.ped <- read.csv(file="jay_ped.csv", na.strings = c("","NA"))
demo.ped <- na.omit(demo.ped)

#change column names
colnames(demo.ped)[1] <- "ID"
colnames(demo.ped)[2] <- "IDband"
colnames(demo.ped)[6] <- "MomID"
colnames(demo.ped)[7] <- "Momband"
colnames(demo.ped)[8] <- "DadID"
colnames(demo.ped)[9] <- "Dadband"


#convert dates to date format
demo.ped$LastObsDate<- as.Date(demo.ped$LastObsDate, format = "%m/%d/%Y")

#add censorship status column 
demo.ped["censorship"] <- 1

#add 0's to those still alive 2016-4-12
demo.ped$censorship[which(demo.ped$LastObsDate=="2016-4-12")]<-0

#Subset the pedigree to get rid of birds with lastobsdate before 1981
#Not sure if this is actually what I want to do, but I did it and it 
#removed the NA's
#demo.ped <- subset(demo.ped, demo.ped$LastObsDate >="1981-05-10")

#add sex indicator column
demo.ped["sexind"] <- 0
demo.ped$sexind[which(demo.ped$Sex=="M")] <- 1
demo.ped$sexind[which(demo.ped$Sex=="F")] <- 2

#add affected column (no meaning here, just needed for function)
demo.ped["aff"] <- 0

demo.ped["missid"] <- "FALSE"

### Need to specify the founders 

jay.ped <- pedigree(id=demo.ped$IDband, dadid=demo.ped$Dadband, 
                    momid=demo.ped$Momband, sex=demo.ped$sexind, 
                    affected=demo.ped$aff, status=demo.ped$censorship, 
                    missid = demo.ped$missid)









