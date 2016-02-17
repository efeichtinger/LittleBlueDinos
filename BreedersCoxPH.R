## February 16, 2016
## Use SurvJan.R and refernce reports starting with FSJ_Surv.Rmd


#Breeders CoxPH - COX PH analyses for FSJ breeders 
#Summer 2015
#UPDATED September 2015

#use this on laptop
setwd("C:/Users/Erin/LittleBlueDinos")

#use this on school desktop
setwd("C:/Users/efeichtinger/LittleBlueDinos")

library(survival)
library(car)
library(kinship2)
df <- read.csv(file="Breeders_Fall.csv")
april.dates <- read.csv(file="April_Census_Dates.csv")
apdat <- as.list(april.dates)

#remove duplicate records (because I still have mulitple records for each individual)
brdrs <- df[!duplicated(df),]

#convert dates to date format
brdrs$FirstDateBred <- as.Date(brdrs$FirstDateBred, format = "%m/%d/%Y")
brdrs$LastObsDate <- as.Date(brdrs$LastObsDate, format = "%m/%d/%Y")

#get year only for YrClass
year <- as.POSIXlt(brdrs$LastObsDate)$year+1900
brdrs["YrDied"]<- year
brdrs$Yr <- as.numeric(brdrs$Yr)

#subtract dates to get number of days
date.diff<- brdrs$LastObsDate-brdrs$FirstDateBred

#add to data frame
brdrs["days"] <- date.diff
brdrs$days <-as.numeric(brdrs$days)

brdrs["yrs"] <- brdrs$days/365

brdrs$FirstDateBred <- as.numeric(brdrs$FirstDateBred)
brdrs$LastObsDate <- as.numeric(brdrs$LastObsDate)

#### Note: lines below NOT THE CORRECT WAY TO DO THIS!!!! For novice breeders

#add censorship, 0 = alive/right censored, 1 = dead after 1st year breeding, novice
#brdrs["censorshipnov"] <- NA
#brdrs$censorshipnov <- 1

#this is not quite right. I want the day of the April census the year 
#following becoming a breeder for the first time, so instead of 365, 
#something like April Census year of 1st breeding + 1? I don't know 
#the R snytax for this.
brdrs$censorshipnov[which(brdrs$days>=365)]<- 0

#censorship status for breeding span after first year, exper breeders
brdrs["censorship"] <- 1
#add 0's to those still alive 2015-06-17
brdrs$censorship[which(brdrs$LastObsDate=="2015-06-17")]<-0


#get rid of birds breeding before 1980
brdrs.new <- subset(brdrs, FirstYr > 1980, select=ID:censorship)
#very important piece of code for the model to work properly 
brdrs.new <- subset(brdrs.new, brdrs.new$days > 0 & brdrs.new$yrs > 0)

#I don't know if it is necessary to take out the NA's
brdrs.new <- na.omit(brdrs.new)

#Create survival object - have to find some way to account for birds with negative days
my.survyr <- Surv(brdrs.new$yrs, brdrs.new$censorship)
my.survdy <- Surv(brdrs.new$days, brdrs.new$censorship)
my.fityr <- survfit(my.survyr~1)
my.fitdy <- survfit(my.survdy~1)
plot(my.fitdy, xlim=c(0,4000), xlab="Days", ylab="Survival", main="Breeders")
plot(my.fityr, xlim=c(0,15), xlab="Years", ylab= "Survival", main="Breeders")
#summary(my.fitdy)

#Simple cox model for novice and experienced breeders 
jay.cox <- coxph(my.survyr~1, data= brdrs.new)
summary(jay.cox)
#year of first breeding
jay.yr <- coxph(my.survyr~brdrs.new$Yr, data=brdrs.new)
summary(jay.yr)
#age at first breeding (for some, minimum age)
jay.age <- coxph(my.survyr~brdrs.new$AgeFirstBreed)
summary(jay.age)

#
jay.exp <- survreg(my.survyr~1, dist="exponential")
summary(jay.exp)

jay.wb <- survreg(my.survyr~1, dist = "weibull")
summary(jay.wb)



#Survival object splitting up novice and experienced breeders
#Not sure if this is the correct way to do this 
nov.surv <- Surv(brdrs.new$yrs, brdrs.new$censorshipnov, type=c('right'))
nov.fit <- survfit(nov.surv~1)
plot(nov.fit, xlim=c(0,1), xlab="Years", ylab="Cumulative Survival", main="Novice and Exp Breeders")

#time variable not numeric it says 
#surv.rcens <- Surv(brdrs.new$LastObsDate~brdrs.new$FirstYr, brdrs.new$censorship)

#basic model
cox.fit <- coxph(my.survdy~1, data= brdrs.new)
null <- basehaz(cox.fit)
plot(null)
plot(null, log="xy", xlab = "Hazard (log)", ylab= "Time(log)")
summary(null)
summary(cox.fit)
cox.fit

yr.fit <- coxph(my.survdy~brdrs.new$Yr, data = brdrs.new)
yr.fit
summary(yr.fit)


sex.fit <- coxph(my.survdy~brdrs.new$Sex, data = brdrs.new)
sex.fit
summary(sex.fit)

#year of first breeding
jay.yr <- coxph(my.survyr~brdrs.new$Yr, data=brdrs.new)
summary(jay.yr)
#age at first breeding (for some, minimum age)
jay.age <- coxph(my.survyr~brdrs.new$AgeFirstBreed)
summary(jay.age)

#Anova in car package
