#Breeders CoxPH - COX PH analyses for FSJ breeders 
#Summer 2015
#UPDATED September 2015

#use this on laptop
setwd("C:/Users/Erin/Dropbox/Jay survival analyses")

#use this on school desktop
setwd("C:/Users/efeichtinger/Dropbox/Jay survival analyses")

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

#add censorship, 0 = alive/right censored, 1 = dead after 1st year breeding, novice
brdrs["censorshipnov"] <- NA
brdrs$censorshipnov <- 1

#change to 1 if bird died within 1 year of breeding, change to less than April of previous year?
brdrs$censorshipnov[which(brdrs$days>=365)]<- 0

#censorship status for breeding span after first year, exper breeders
brdrs["censorship"] <- 1
#add 0's to those still alive 2015-06-17
brdrs$censorship[which(brdrs$LastObsDate=="2015-06-17")]<-0


#get rid of birds breeding before 1980
brdrs.new <- subset(brdrs, FirstYr > 1980, select=ID:censorship)

#I don't know if it is necessary to take out the NA's
#brdrs.new <- na.omit(brdrs.new)

#Create survival object - have to find some way to account for birds with negative days
my.survyr <- Surv(brdrs.new$yrs, brdrs.new$censorship)
my.survdy <- Surv(brdrs.new$days, brdrs.new$censorship)
my.fityr <- survfit(my.survyr~1)
my.fitdy <- survfit(my.survdy~1)
plot(my.fitdy, xlim=c(0,4000), xlab="Days", ylab="Survival")
plot(my.fityr, xlim=c(0,15), xlab="Years", ylab= "Survival")
#summary(my.fitdy)

#time variable not numeric it says 
surv.rcens <- Surv(brdrs.new$LastObsDate~brdrs.new$FirstYr, brdrs.new$censorship)

#basic model
cox.fit <- coxph(my.survdy~1, data= brdrs.new)
null <- basehaz(cox.fit)
plot(null)
summary(null)
summary(cox.fit)
cox.fit

#error
yr.fit <- coxph(my.survdy~brdrs.new$Yr, data = brdrs.new)
yr.fit
summary(yr.fit)


sex.fit <- coxph(my.survdy~brdrs.new$Sex, data = brdrs.new)
sex.fit
summary(sex.fit)

#Anova in car package
