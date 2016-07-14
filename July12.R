### New new script for modeling birds who made it to the yearling node 
### Prebreeders and breeders of known and unknown age
### Need slightly different approach then before 


library(survival)
library(car)
library(coxme)
library(kinship2)

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

## Read in breeder data set with known-age breeders 
brdr.df <- read.csv("Erin_June_KA_Breeders.csv")
colnames(brdr.df)[8] <- "Days"
brdr.df["Yrs"] <- brdr.df$Days/365.25

brdr.df$FldgDate <- as.Date(brdr.df$FldgDate, format = "%m/%d/%Y")
brdr.df$LastObsDate <- as.Date(brdr.df$LastObsDate, format = "%m/%d/%Y")
brdr.df$BrdrDate <- as.Date(brdr.df$BrdrDate, format = "%m/%d/%Y")

## Create new data frame for birds that made it 1 year but not to breeding
## at least not at ABS, could have died or dispersed from ABS 
hlpr.df <- subset(bird.df, !(JayID %in% brdr.df$JayID))

## remove natal nest field 
brdr.df[4] <- NULL

## get date became yearling (it's actually the April census date) but using 365 days
brdr.df["YrlgDate"] <- brdr.df$FldgDate + 365

## Remove fledge date column
brdr.df[4] <- NULL

## Change order of columns to get in desired form
brdr.df <- brdr.df[c(1,5,2,3,9,7,4,6,8)]
str(brdr.df)

## Add column for prebreeder duration 
brdr.df["PrbDur"] <- (brdr.df$BrdrDate - brdr.df$YrlgDate)/365
brdr.df[10] <- round(brdr.df$PrbDur, digits = 2)

brdr.df$PrbDur <- as.numeric(brdr.df$PrbDur)
hist(brdr.df$PrbDur, main = "Distribution of time in prebreeder state",
     xlab = "Time (years)")

mean(brdr.df$PrbDur)
sd(brdr.df$PrbDur)
range(brdr.df$PrbDur)

## Add 1 to get age 
brdr.df["Fromstart"] <- brdr.df$PrbDur + 1

hist(brdr.df$Fromstart, main = "Distribution of age at first breeding",
     xlab = "Age (years)")
#Total age from fledging 
mean(brdr.df$Fromstart)
sd(brdr.df$Fromstart)
range(brdr.df$Fromstart)

females <- subset(brdr.df, brdr.df$Sex == "F")
males <- subset(brdr.df, brdr.df$Sex == "M")

mean(females$Fromstart)
sd(females$Fromstart)
range(females$Fromstart)

hist(females$Fromstart, xlab = "Age (years)", main = "Age at first breeding females")

mean(females$PrbDur)
sd(females$PrbDur)
range(females$PrbDur)

hist(females$PrbDur, main = "Duration of prebreeder state",
     xlab = "Time (years)")


### males 
mean(males$Fromstart)
sd(males$Fromstart)
range(males$Fromstart)

hist(males$Fromstart, xlab = "Age (years)", main = "Age at first breeding males")

mean(males$PrbDur)
sd(males$PrbDur)
range(males$PrbDur)





####################################################
## Rearrange "hlpr.df" in the same way 
hlpr.df[4] <- NULL
hlpr.df["YrlgDate"] <- hlpr.df$FldgDate + 365


## Change order of columns to get in desired form
hlpr.df <- hlpr.df[c(1,6,2,3,4,9,7,5,10)]


str(hlpr.df)

hlpr.df$Days <- as.numeric(hlpr.df$Days)
hlpr.df["Timef1y"] <- hlpr.df$Days/365
hlpr.df["Timeff"] <- hlpr.df$Timef1y + 1

hlpr.df[9] <- NULL

hist(hlpr.df$Timeff)

## Could prove to be useful later - might want to add age at death, age at first breeding 
#year <- as.POSIXlt(bird.df$FldgDate)$year+1900
#yrlg.df["FYear"] <- year
#yrlg.df$FYear <- as.factor(yrlg.df$FYear)

## Read in file with unknown age breeders 
unbr.df <- read.csv("Erin_UnknownAge_Breeders.csv")

unbr.df <- unbr.df[!duplicated(unbr.df),]
unbr.df <- unbr.df[c(1,2,8,7,4,5,3,6)]
colnames(unbr.df)[1] <- "JayID"
colnames(unbr.df)[2] <- "USFWBand"
colnames(unbr.df)[3] <- "Sex"

unbr.df["YrlgDate"] <- NA

unbr.df <- unbr.df[c(1,2,3,4,9,5,6,7,8)]
unbr.df[8] <- NULL
unbr.df[9] <- NULL
colnames(unbr.df)[6] <- "BrdrDate"

jays <- rbind(hlpr.df, brdr.df)





