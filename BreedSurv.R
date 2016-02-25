##February 21 2016
##Script for survival models - breeders for now
##USE THIS ERIN ANY TIME AFTER FEBRUARY 21
###Refer to SurvJan.R and BreedersCoxPH.R for code help
## These data have the start date as the minimum date that a bird is
# known to have laid or sired an egg, this is when a bird is considered to be
# a breeder

#use this on laptop
setwd("C:/Users/Erin/LittleBlueDinos")
#use this on school desktop
setwd("C:/Users/efeichtinger/Dropbox/Jay_data_nogithub")

library(survival)
library(car)
library(kinship2)
library(SurvRegCensCov)

##Read in CSV file of male and female breeders with mulitple rows for each bird
bird.df <- read.csv("Breeders_Fall.csv")
str(bird.df)

#remove duplicates  - for years where there was more than one nest in a year
jay.df<- bird.df[!duplicated(bird.df),]
str(jay.df)

#colnames(jay.df)[1] <- "ID"
#colnames(jay.df)[2] <- "Band"
colnames(jay.df)[5] <- "MinDate"

#convert dates to date format
jay.df$MinDate <- as.Date(jay.df$MinDate, format = "%m/%d/%Y")
jay.df$LastObsDate <- as.Date(jay.df$LastObsDate, format = "%m/%d/%Y")

#subtract dates to get number of days
date.diff<- jay.df$LastObsDate-jay.df$MinDate

#and survival period in years, account for leap year 
jay.df["Yrs"] <- as.numeric(date.diff/365.25)

jay.df$FirstYr <- as.factor(jay.df$FirstYr)

#very important piece of code for the model to work properly, remove any 
#weird entries like birds that have negative years of experience or a negative
#survival interval 
jay.df <- subset(jay.df, jay.df$Days > 0)


#add column for censorship status, in survival package - 0=alive, 1=dead
jay.df["censorship"] <- 1
#If last observed date = 10/14/2015, 0 for still alive today
jay.df$censorship[which(jay.df$LastObsDate=="2015-10-14")]<-0

#Check for correct structure
str(jay.df)

#How many males and females?
sum(jay.df$Sex == "M")
sum(jay.df$Sex == "F")

#change back to numeric for survival object 
jay.df$MinDate <- as.numeric(jay.df$MinDate)
jay.df$LastObsDate <- as.numeric(jay.df$LastObsDate)

#Create survival object - IS THIS CORRECT?? 
jay.ob <- Surv(jay.df$Yrs, jay.df$censorship, type =c('right'))
jay.lifetab <- survfit(jay.ob~1, conf.type = "log-log")
jay.fit <- plot(jay.lifetab, xlab = "Time (years)", 
          ylab = "Cumulative Survival", main = "FL Scrub Breeder survival",
          pin = c(5,5))
#Log scale
jay.fitlog <- plot(jay.lifetab, xlab = "Time (years)",
                   log = "y", ylim = c(0.001,2), ylab = "Cumulative Survival", 
                   main = "FL Scrub Breeder survival Log Scale")
#Grouping by sex - following example "Cox Regression in R" J. Fox
km.sex <- survfit(jay.ob ~ jay.df$Sex, conf.type = "log-log")
km.fit <- plot(km.sex, col=c("dodgerblue2","red"),lty = c(1,2), lwd = 2, xlab = "Time (years)", 
               ylab = "Survival", main = "Survival by Sex", pin = c(20,20))
legend("topright", c("Females","Males"), col = c("dodgerblue2","red"),
       lty = c(1,2), lwd = 2)
sex.log <- plot(km.sex, col = c("navy","red"), log = "y", ylim = c(0.001,2),
          lty  = c(1,2), xlab = "Time (years)",ylab = "Cumulative Survival", 
                main = "Survival by Sex Log Scale")
legend("topright", c("Females","Males"), col = c("navy","red"), 
       lty = c(1,2), lwd =1)


#Summary statistics
mean(jay.df$Yrs)
sd(jay.df$Yrs)
median(jay.df$Yrs)
range(jay.df$Yrs)

males <- subset(jay.df, jay.df$Sex == "M")
females <- subset(jay.df, jay.df$Sex == "F")

#Table of means & SD of males and females for age at first breeding
ma <- cbind(mean(na.omit(females$AgeFirstBreed)), mean(na.omit(males$AgeFirstBreed)))
sa <- cbind(sd(na.omit(females$AgeFirstBreed)), sd(na.omit(males$AgeFirstBreed)))
age.table <- rbind(ma,sa)
colnames(age.table) <- c("Females", "Males")
rownames(age.table) <- c("Mean", "SD")
age.table

#table of means & SD for males and females of years survived once becoming a breeder
means <- cbind(mean(females$Yrs), mean(males$Yrs))
sd <- cbind(sd(females$Yrs), sd(males$Yrs))
colnames(means) <- c("Females", "Males")
m.sd <- rbind(means, sd)
rownames(m.sd) <- c("Mean", "SD")
m.sd

#What does the range of ages at first breeding look like?
jay <- na.omit(jay.df$AgeFirstBreed)

range(jay)
mean(jay)
sd(jay)
median(jay)

#simple plots to look at age and years of breeding experience
#Mental gymnastics trying to figure out which one should be y
#Can't do these plots on data frame with one record per individual 
#plot(jay.df$YrsExp, jay.df$CurrentAge, xlab = "Years Experience",
     #ylab = "Current Age")
#plot(jay.df$CurrentAge, jay.df$YrsExp, xlab = "Current Age",
     #ylab = "Years Experience")

#Look at distribution of age of first time breeders 
#Bar plot? Bubble plot? 
#ggplot2




#First Cox Models

#Age at first breeding 
cox1 <- coxph(jay.ob ~ AgeFirstBreed, data = jay.df)
summary(cox1)
#Check for violation of proportional hazard 
res.cox1 <- cox.zph(cox1)
res.cox1
plot(res.cox1)
#They are proportional? "A non-zero slope is evidence against proportionality"

#Sex as predictor
cox2 <- coxph(jay.ob ~ Sex, data = jay.df)
summary(cox2)
res.cox2 <- cox.zph(cox2)
res.cox2
plot(res.cox2)
#First year of breeding as predictor
cox3 <- coxph(jay.ob ~ FirstYr, data = jay.df)
summary(cox3)
res.cox3 <- cox.zph(cox3)
res.cox3
plot(res.cox3)

#Age and sex
cox4 <- coxph(jay.ob ~ AgeFirstBreed + Sex, data = jay.df)
summary(cox4)

#Hazard ratio - B1 + B2
hr.4 <- cox4$coefficients[1] + cox4$coefficients[2]
hr.4

#Include years experience, age, and sex, no interactions
cox5 <- coxph(jay.ob ~ AgeFirstBreed + Sex + FirstYr, data = jay.df)
summary(cox5)
res.cox5 <- cox.zph(cox5, transform = "km")
res.cox5
plot(res.cox45)

cox6 <- coxph(jay.ob ~ Sex + FirstYr, data = jay.df)
summary(cox6)

cox7 <- coxph(jay.ob ~ AgeFirstBreed + FirstYr, data = jay.df)
summary(cox7)

#AFT model with Weibull distribution and years of experience
AFT.weibull <- survreg(jay.ob ~ AgeFirstBreed, data = jay.df, dist = "weibull")
summary(AFT.weibull)

#AFT model with Weibull distribution and sex
AFT.weibull2 <- survreg(jay.ob ~ Sex, data = jay.df, dist = "weibull")
summary(AFT.weibull2)

#AFT model with Weibull distribution and sex
AFT.weibull3 <- survreg(jay.ob ~ FirstYr, data = jay.df, dist = "weibull")
summary(AFT.weibull3)

#All 3 covariates 
AFT.weibull4 <- survreg(jay.ob ~ AgeFirstBreed + Sex + FirstYr, 
      data = jay.df, dist = "weibull")

#Compare AFT model (weibull) with Cox PH model that has 3 covariates 
summary(cox4)
summary(AFT.weibull4)

#Analysis of deviance
anova(cox4)


###################################################################
#Just for fun, let's look at survival of all known age birds in the population

#Known age birds

#read in CSV of all known-age birds 
birds <- read.csv("Erin_Surv_All.csv")
#str(birds)

#remove duplicates 
birds2 <- birds[!duplicated(birds),]
#str(birds2)

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
birds2["yrs"] <- birds2$days/365.25

#very important piece of code for the model to work properly
#remove any zero or negative values in days and years 
birds2 <- subset(birds2, birds2$days > 0 & birds2$yrs > 0)

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
birds2$FYear <- as.factor(birds2$FYear)

birds2 <- birds2[-which(birds2$Sex == ""),]


#Create survival object based off Gordon's Cactus Finch example
survobj <- Surv(birds2$yrs, birds2$censorship, type =c('right'))

all.lifetab <- survfit(survobj~1)
all.fit <- plot(jay.lifetab, xlab = "Time (years)", 
ylab = "Cumulative Survival", main = "All known-age birds",
pin = c(5,5))

all.log <- plot(jay.lifetab, log = "y", ylim=c(0.001,2),
     xlab =  "Time (years)", ylab = "Cumulative Survival", 
      main = "Survival of Known Age Birds - Log Scale")
        

all.sex <- survfit(survobj~birds2$Sex, conf.type = "log-log")
sex.kmall <- plot(all.sex, col = c("blue", "red"),
    xlab = "Time (years)", ylab = "Survival", main = "Survival by Sex")
legend("topright", c("Females","Males"), col = c("blue","red"),
       lty = c(1,2), lwd = 2)

all.logsex <- plot(all.sex, col = c("blue", "red"), log = "y", 
            ylim = c(0.001,2), xlab = "Time (years)",
    ylab = "Cumulative Survival", main = "Survival by Sex - Log Scale")
legend("topright", c("Females","Males"), col = c("blue","red"),
       lty = c(1,2), lwd = 2)

#Doesn't work right - warning message: X matrix deemed to be singular
cox <- coxph(survobj ~ birds2$Sex, data = birds2)
cox1 <- coxph(survobj ~ birds2$FYear, data = birds2)

#AFT model for sex
AFT.sex <- survreg(survobj ~ Sex, data = birds2, dist = "weibull")
summary(AFT.sex)
