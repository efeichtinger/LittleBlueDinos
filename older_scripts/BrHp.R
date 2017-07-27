## July 1 2016 (HOLY SHIT)
## Script for breeders and helpers  
## Source code from "JaySurv.R" 

library(survival)
library(car)
library(coxme)
library(kinship2)
library(SurvRegCensCov)
library(parfm)

############################################################
## Read in April Census Dates
aprilD <- read.csv("April_Census_Dates.csv")
## Convert to date object
aprilD$CensusDate <- as.Date(aprilD$CensusDate, format = "%m/%d/%Y")

## Read in data frame, all known age birds
## Note that this comes from query names "Erin_June_FY" 
bird.df <- read.csv("Erin_June_Birds.csv")
str(bird.df)

colnames(bird.df)[8] <- "Days"
bird.df["Yrs"] <- bird.df$Days/365.25

## convert dates to date format
bird.df$FldgDate <- as.Date(bird.df$FldgDate, format = "%m/%d/%Y")
bird.df$LastObsDate <- as.Date(bird.df$LastObsDate, format = "%m/%d/%Y")

## Remove birds that died in the first year so I'm left with helpers and breeders
bird.df <- subset(bird.df, bird.df$Yrs >= 1)
str(bird.df)

mean(bird.df$Yrs)
range(bird.df$Yrs)
sd(bird.df$Yrs)
median(bird.df$Yrs)

#remove duplicates 
#birds2 <- birds[!duplicated(birds),]


#1009 individuals which is consistent with the database
#after removing helpers (in database), left with 527 known age breeders
#need to figure out the birds that are in the 1009 dataframe but not
#in the 527 set, so 482 of the birds should be in both 
#also later need to add in unknown age breeders 

## The data called below includes only known age breeders, these same inds
#are in bird.df but the ones that died before becoming a breeder are not in 
#brdr.df


###########################################################
## Read in breeder data set with known-age breeders 
brdr.df <- read.csv("Erin_June_KA_Breeders.csv")
colnames(brdr.df)[8] <- "Days"
brdr.df["Yrs"] <- brdr.df$Days/365.25

brdr.df$FldgDate <- as.Date(brdr.df$FldgDate, format = "%m/%d/%Y")
brdr.df$LastObsDate <- as.Date(brdr.df$LastObsDate, format = "%m/%d/%Y")
brdr.df$BrdrDate <- as.Date(brdr.df$BrdrDate, format = "%m/%d/%Y")

###############################################
## censorship - 1 = dead, 0 = alive
brdr.df["Censor"] <- 1 
brdr.df$Censor[which(brdr.df$LastObsDate =="2016-4-12")] <- 0

str(brdr.df)

brdr.df[11] <- NULL
brdr.df[9] <- NULL
################################################################


################################################
#Matching breeders and helpers
#These birds did become breeders 
brdr.df["Censor"] <- 0

## Now how to find who is in the helper and breeder set (bird.df) and 
#who is only in the breeder set
#There are 482 birds that became helpers but did not make it to breeding 

#Birds that made it to 1 year old but did not become breeders
hlpr <- subset(bird.df, !(JayID %in% brdr.df$JayID))
str(hlpr)

#Need to add the birds in hlpr to the birds who did make it 
#Add to brdr.df but needs some restructuring 

#Add 1 to all rows because they all died before becoming breeders
hlpr["Censor"] <- 1

brhlp <- rbind(brdr.df, hlpr)

#But give 0 to birds still alive today! 
brhlp$Censor[which(brhlp$LastObsDate=="2016-4-12")] <- 0

brhlp <- subset(brhlp, brhlp$Yrs > 0)

#add censorship column, some helpers are still alive as of this past April
#hlpr["Censor"] <- 1
#hlpr$Censor[which(hlpr$LastObsDate=="2016-4-12")] <- 0

#hlpr <- subset(hlpr, hlpr$Days > 0)

#str(hlpr)

#Change dates back to numeric 
brhlp$FldgDate <- as.numeric(brhlp$FldgDate)
brhlp$LastObsDate <- as.numeric(brhlp$LastObsDate)

## Birds coded 0 are either still alive today or made it to the breeder stage
## Birds coded 1 died after turning 1 year old and 

##survival object for helpers
brhl.ob <- Surv(brhlp$Yrs, brhlp$Censor, type= c('right'))
km.fit <- survfit(brhl.ob ~ 1, conf.type = "log-log")
kmplot <- plot(km.fit, xlab="Time (years)", log = "y", 
            ylab = "Cumulative Survival (log)", main = "Helpers",
            ylim = c(0.45,1), xlim = c(1,5))
          

km.sex <- survfit(brhl.ob ~ brhlp$Sex, conf.type = "log-log")
sexplot <- plot(km.sex,col = c("navy","red"), xlab = "Time (years)", log = "y",
                ylab = "Cumulative Survival", main = "Helpers",
                ylim =c(0.40,1), xlim = c(1,5))
legend("topright", c("Females","Males"), col=c("navy","red"),
       lty = c(1,2),lwd=1)

cph.sex <- coxph(brhl.ob ~ Sex, data = brhlp)
summary(cph.sex)
cox.zph(cph.sex)
plot(cox.zph(cph.sex))


#mm <- coxme(brhl.ob ~ Sex + (1|NatalNest), data = brhlp)
#summary(mm)

####################################################
## July 21 2016 - breeders with both known and unknown age 

## Read back in for breeders only, need to add the unknown age ones too
unbr.df <- read.csv("Erin_UnknownAge_Breeders.csv")
colnames(unbr.df)[1] <- "ID"
colnames(unbr.df)[2] <- "Band"
colnames(unbr.df)[4] <- "BrdrDate"

unbr.df$LastObsDate <- as.Date(unbr.df$LastObsDate, format = "%m/%d/%Y")
unbr.df$BrdrDate <- as.Date(unbr.df$BrdrDate, format = "%m/%d/%Y")

unbr.df["Yrs"] <- (unbr.df$LastObsDate-unbr.df$BrdrDate)/365
unbr.df["Censor"] <- 1
unbr.df$Censor[which(unbr.df$LastObsDate =="2016-4-12")] <- 0

unbr.df <- subset(unbr.df, unbr.df$Yrs >0)


#####################################################
brdr.df <- read.csv("Erin_June_KA_Breeders.csv")
colnames(brdr.df)[8] <- "Days"
brdr.df["Yrs"] <- brdr.df$Days/365.25

brdr.df$FldgDate <- as.Date(brdr.df$FldgDate, format = "%m/%d/%Y")
brdr.df$LastObsDate <- as.Date(brdr.df$LastObsDate, format = "%m/%d/%Y")
brdr.df$BrdrDate <- as.Date(brdr.df$BrdrDate, format = "%m/%d/%Y")

#Add first year of breeding 
year <- as.POSIXlt(brdr.df$BrdrDate)$year+1900
brdr.df["FirstYear"] <- year
brdr.df$FirstYear <- as.factor(brdr.df$FirstYear)

## censorship - 1 = dead, 0 = alive
brdr.df["Censor"] <- 1 
brdr.df$Censor[which(brdr.df$LastObsDate =="2016-4-12")] <- 0

brdr.df <- subset(brdr.df, brdr.df$Yrs > 0 & brdr.df$Days > 0)

str(brdr.df)
#View(brdr.df)

#have to add in age at first breeding 
brdr.df <- brdr.df[,c(1,7,2,3,9,11,6)]

##Move after manipulation

#brdr.df$FldgDate <- as.numeric(brdr.df$FldgDate)
#brdr.df$LastObsDate <- as.numeric(brdr.df$LastObsDate)
#brdr.df$BrdrDate <- as.numeric(brdr.df$BrdrDate)

## Manipulate dfs to be in same format to bind
unbr.df <- unbr.df[,c(1,2,8,7,4,3,5,6,9,10)]





#######################################################
## Survival object for breeders (known age)
br.ob <- Surv(brdr.df$Yrs, brdr.df$Censor, type =c('right'))

## Kaplan Meier estimate
br.fit <- survfit(br.ob ~ 1, conf.type= "log-log")
br.plot <- plot(br.fit, xlab= "Time (years)", log = "y",
                ylab="Cumulative Survival (log)", main = "Breeders",
                xlim = c(0,15), ylim = c(0.01,1))

br.fit
str(br.fit)

br.sex <- survfit(br.ob ~ brdr.df$Sex, conf.type = "log-log")
bx.plot <- plot(br.sex,col = c("navy","red"), xlab = "Time (years)", log = "y",
                ylab = "Cumulative Survival", main = "Breeders",
                ylim =c(0.40,1), xlim = c(1,5))
legend("topright", c("Females","Males"), col=c("navy","red"),
       lty = c(1,2),lwd=1)

c1 <- coxph(br.ob ~ Sex, data = brdr.df)
summary(c1)
cox.zph(c1)
plot(cox.zph(c1))


