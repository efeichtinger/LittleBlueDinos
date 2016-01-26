### JANUARY 25 2016 
### Fitting survival models to the data
### Although the Cox PH models may very well be the best, have to fit basic
#survival models first 

#Much of this code is taken from "BreedersCoxPH.R" - KEEP THAT FILE AROUND
#It's the breeder only file 

## Start with all birds at first where they are NOT split up by stage 
## Start with Kaplan Meier and by reading Gordon's chapter 

#use this on laptop
setwd("C:/Users/Erin/LittleBlueDinos")
#use this on school desktop
setwd("C:/Users/efeichtinger/LittleBlueDinos")

library(survival)
library(car)
library(kinship2)

#read in CSV of all birds 
birds <- read.csv("Erin_Surv_All.csv")
str(birds)

#remove duplicates 
birds2 <- birds[!duplicated(birds),]
str(birds2)
#2322 observations of 6 variables - known age birds only, so birds that were
#hatched on the study tract (fledge date knowm), from Jan 1981 to Dec 2015 

colnames(birds2)[7] <- "LastObsDate"


#convert dates to date format
birds2$FldgDate <- as.Date(birds2$FldgDate, format = "%m/%d/%Y")
birds2$LastObsDate <- as.Date(birds2$LastObsDate, format = "%m/%d/%Y")

#subtract dates to get number of days
date.diff<- birds2$LastObsDate-birds2$FldgDate

#add to data frame - survival period in days
birds2["days"] <- date.diff
birds2$days <-as.numeric(birds2$days)
#and survival period in years 
birds2["yrs"] <- birds2$days/365

#very important piece of code for the model to work properly 
birds2 <- subset(birds2, birds2$days > 0)

#add column for censorship status 
birds2["censorship"] <- 1
#If last observed date = 10/14/2015, 0 for still alive today
birds2$censorship[which(birds2$LastObsDate=="2015-10-14")]<-0


year <- as.POSIXlt(birds2$FldgDate)$year+1900
birds2["FYear"] <- year

#change back to numeric for survival object 
birds2$FldgDate <- as.numeric(birds2$FldgDate)
birds2$LastObsDate <- as.numeric(birds2$LastObsDate)
birds2$days <- as.numeric(birds2$days)
birds2$yrs <- as.numeric(birds2$yrs)


#Create survival object based off Gordon's Cactus Finch example
<<<<<<< HEAD
survobj <- Surv(birds2$days, birds2$censorship, type =c('right'))
jay.lifetab <- survfit(survobj~1)
jay.fit <- plot(jay.lifetab, xlab = "Time (days)", 
      ylab = "Cumulative Survival", main = "FL Scrub Jay survival")
jay.fitlog <- plot(jay.lifetab, log= "xy", xlab = "Time (years)", 
                  main = "FL Scrub Jay survival")
jay.cox <- coxph(survobj~1, data= birds2)
#jay.coxsex <- coxph(survobj ~ birds2$sex, data = birds2)

#Exponential model
#First fit with just intercept
## error of invalid survival times for this distribution - mitigated by removing 0's
jay.int <- survreg(Surv(birds2$days, birds2$censorship)~1, dist="exponential")
summary(jay.int)

jay.weib <- survreg(survobj~1, dist = "weibull")
summary(jay.weib)
survobj <- Surv(birds2$yrs, birds2$censorship)
jay.lifetab <- survfit(survobj~1)
jay.fit <- plot(jay.lifetab, xlab = "Time (years)", 
      ylab = "Cumulative Survival", main = "FL Scrub Jay survival")

jay.cox <- coxph(survobj~1, data= birds2)
jay.coxsex <- coxph(survobj ~ birds2$sex, data = birds2)

#Exponential model
#First fit with just intercept 
jay.int <- survreg(Surv(birds2$days, birds2$censorship)~1, dist="exponential")
## error of invalid survival times for this distribution

jay.weib <- survreg(survobj~1, dist = "weibull")

