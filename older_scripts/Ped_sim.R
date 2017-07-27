# July 5 2016
# Pedigree simulation so I can learn the ins and outs
# Mostly to correctly code the founders
# The founders fields for dadid and momid have NA in these fields 

library(kinship2)
pdsim <- read.csv("Ped_sim.csv")
ped.df <- pdsim[-(12:185), ] 

ped.df["sexcode"] <- 1
ped.df["status"] <- 1
ped.df["affected"] <- 0

ped.df$sexcode[which(ped.df$Sex == "F")] <- 2

sim.ped <- pedigree(ped.df$JayID, ped.df$DadID, ped.df$MomID, ped.df$sexcode,
              ped.df$affected, ped.df$status)

k <- kinship(sim.ped)

#Easy enough to code for the founders but the real challenge is finding the founders


jped <- read.csv("Erin_Demo_Ped.csv")
jped[12] <- NULL

jped[jped == ""] <- NA
jped <- na.omit(jped)

## convert dates to date format
jped$FldgDate <- as.Date(jped$FldgDate, format = "%m/%d/%Y")
jped$LastObsDate <- as.Date(jped$LastObsDate, format = "%m/%d/%Y")


jped["sexcode"] <- 1
jped$sexcode[which(jped$Sex == "F")] <- 2


jped["status"] <- 1
jped$status[which(jped$LastObsDate == "2016-4-12")]<-0

jped["affected"] <- 0

#Error is that I need at least 2 founders 
#I need to find those founders, append rows for them and leave the 
#mom and dadid fields blank - then it knows these are the founders
ped.prac <- pedigree(jped$Jays.USFWBand, jped$Jays_1.USFWBand, jped$Jays_2.USFWBand,
        jped$sexcode, jped$affected, jped$status)




