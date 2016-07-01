## July 1 2016 (HOLY SHIT)
## Script for breeders and helpers  
## Source code from "JaySurv.R" 

library(survival)
library(car)
library(coxme)
library(kinship2)
library(SurvRegCensCov)
library(parfm)

## Read in April Census Dates
aprilD <- read.csv("April_Census_Dates.csv")
## Convert to date object
aprilD$CensusDate <- as.Date(aprilD$CensusDate, format = "%m/%d/%Y")

## Read in data frame, all known age birds
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

#1009 individuals which is consistent with the database
#after removing helpers (in database), left with 527 known age breeders
#need to figure out the birds that are in the 1009 dataframe but not
#in the 527 set, so 482 of the birds should be in both 
#also later need to add in unknown age breeders 

## The data called below includes only known age breeders, these same inds
#are in bird.df but the ones that died before becoming a breeder are not in 
#brdr.df

## Read in breeder data set with known-age breeders 
brdr.df <- read.csv("Erin_June_KA_Breeders.csv")
colnames(brdr.df)[8] <- "Days"
brdr.df["Yrs"] <- brdr.df$Days/365.25

brdr.df$FldgDate <- as.Date(brdr.df$FldgDate, format = "%m/%d/%Y")
brdr.df$LastObsDate <- as.Date(brdr.df$LastObsDate, format = "%m/%d/%Y")
brdr.df$BrdrDate <- as.Date(brdr.df$BrdrDate, format = "%m/%d/%Y")

## censorship - 1 = dead, 0 = alive
brdr.df["Censor"] <- 1 
brdr.df$Censor[which(brdr.df$LastObsDate =="2016-4-12")] <- 0

str(brdr.df)

## Now how to find who is in the helper and breeder set (bird.df) and 
#who is only in the breeder set
#There are 482 birds that became helpers but did not make it to breeding 

#Birds that made it to 1 year old but did not become breeders
hlpr <- subset(bird.df, !(JayID %in% brdr.df$JayID))

#add censorship column, some helpers are still alive as of this past April
hlpr["Censor"] <- 1
hlpr$Censor[which(hlpr$LastObsDate=="2016-4-12")] <- 0

hlpr <- subset(hlpr, hlpr$Days > 0)

str(hlpr)

#Change dates back to numeric 
hlpr$FldgDate <- as.numeric(hlpr$FldgDate)
hlpr$LastObsDate <- as.numeric(hlpr$LastObsDate)

##survival object for hlpr
hlpr.ob <- Surv(hlpr$Yrs, hlpr$Censor, type= c('right'))
km.fit <- survfit(hlpr.ob ~ 1, conf.type = "log-log")
kmplot <- plot(km.fit, xlab="Time (years)", log = "y", 
            ylab = "Cumulative Survival (log)", main = "Helpers",
            ylim = c(0.01,1), xlim = c(1,5))
          

km.sex <- survfit(hlpr.ob ~ hlpr$Sex, conf.type = "log-log")
sexplot <- plot(km.sex,col = c("navy","red"), xlab = "Time (years)", log = "y",
                ylab = "Cumulative Survival", main = "Helpers",
                ylim =c(0.01,1), xlim = c(1,5))
legend("topright", c("Females","Males"), col=c("navy","red"),
       lty = c(1,2),lwd=1)

cph.sex <- coxph(hlpr.ob ~ Sex, data = hlpr)
summary(cph.sex)
cox.zph(cph.sex)
plot(cox.zph(cph.sex))


mm <- coxme(hlpr.ob ~ Sex + (1|NatalNest), data = hlpr)
summary(mm)


