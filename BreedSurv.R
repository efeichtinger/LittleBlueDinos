##February 21 2016
##Script for survival models - breeders for now
##USE THIS ERIN ANY TIME AFTER FEBRUARY 21
###Refer to SurvJan.R and BreedersCoxPH.R for code help
## These data have the start date as the minimum date that a bird is
# known to have laid or sired an egg, this is when a bird is considered to be
# a breeder

#use this on laptop
setwd("C:/Users/Erin/Dropbox/Jay_data_nogithub")
#use this on school desktop
setwd("C:/Users/efeichtinger/Dropbox/Jay_data_nogithub")

library(survival)
library(car)
library(kinship2)

##Read in CSV file of male and female breeders with mulitple rows for each bird
bird.df <- read.csv("Erin_Breeders_All_Years.csv")
str(bird.df)

#remove duplicates  - for years where there was more than one nest in a year
jay.df<- bird.df[!duplicated(bird.df),]
str(jay.df)

colnames(jay.df)[1] <- "ID"
colnames(jay.df)[2] <- "Band"

#convert dates to date format
jay.df$MinDate <- as.Date(jay.df$MinDate, format = "%m/%d/%Y")
jay.df$LastObsDate <- as.Date(jay.df$LastObsDate, format = "%m/%d/%Y")

#subtract dates to get number of days
date.diff<- jay.df$LastObsDate-jay.df$MinDate

#and survival period in years, account for leap year 
jay.df["Yrs"] <- date.diff/365.25

#very important piece of code for the model to work properly, remove any 
#weird entries like birds that have negative years of experience or a negative
#survival interval 
jay.df <- subset(jay.df, jay.df$YrsExp >= 0 & jay.df$Yrs > 0)
str(jay.df)

#add column for censorship status, in survival package - 0=alive, 1=dead
jay.df["censorship"] <- 1
#If last observed date = 10/14/2015, 0 for still alive today
jay.df$censorship[which(jay.df$LastObsDate=="2015-10-14")]<-0


#change back to numeric for survival object 
jay.df$MinDate <- as.numeric(jay.df$MinDate)
jay.df$LastObsDate <- as.numeric(jay.df$LastObsDate)
jay.df$Yrs <- as.numeric(jay.df$Yrs)

#Create survival object - IS THIS CORRECT?? 
jay.ob <- Surv(jay.df$Yrs, jay.df$censorship, type =c('right'))
jay.lifetab <- survfit(jay.ob~1)
jay.fit <- plot(jay.lifetab, xlab = "Time (years)", 
          ylab = "Cumulative Survival", main = "FL Scrub Breeder survival")

#Summary statistics
mean(jay.df$Yrs)
sd(jay.df$Yrs)
median(jay.df$Yrs)
range(jay.df$Yrs)

#simple plots to look at age and years of breeding experience
#Mental gymnastics trying to figure out which one should be y
plot(jay.df$YrsExp, jay.df$CurrentAge, xlab = "Years Experience",
     ylab = "Current Age")
plot(jay.df$CurrentAge, jay.df$YrsExp, xlab = "Current Age",
     ylab = "Years Experience")

#First Cox Models

#Years of breeder experience as predictor
cox1 <- coxph(jay.ob ~ YrsExp, data = jay.df)
#Sex as predictor
cox2 <- coxph(jay.ob ~ Sex, data = jay.df)
#Age as predictor
cox3 <- coxph(jay.ob ~ CurrentAge, data = jay.df)

#Include years experience, age, and sex, no interactions
cox4 <- coxph(jay.ob ~ YrsExp + Sex + CurrentAge, data = jay.df)

#AFT model with exponential distribution and years of experience 
AFT.exp1 <- survreg(jay.ob ~ YrsExp, data = jay.df, dist = "exponential")
#AFT model with Weibull distribution and years of experience
AFT.weibull <- survreg(jay.ob ~ YrsExp, data = jay.df, dist = "weibull")

#AFT model with exponential distribution and sex
AFT.exp2 <- survreg(jay.ob ~ Sex, data = jay.df, dist = "exponential")
#AFT model with Weibull distribution and sex
AFT.weibull2 <- survreg(jay.ob ~ Sex, data = jay.df, dist = "weibull")

#AFT model with exponential distribution and sex
AFT.exp3 <- survreg(jay.ob ~ CurrentAge, data = jay.df, dist = "exponential")
#AFT model with Weibull distribution and sex
AFT.weibull3 <- survreg(jay.ob ~ CurrentAge, data = jay.df, dist = "weibull")

#All 3 covariates 
AFT.weibull4 <- survreg(jay.ob ~ YrsExp + Sex + CurrentAge, data = jay.df,
                        dist = "weibull")

#Compare AFT model (weibull) with Cox PH model that has 3 covariates 
summary(cox4)
summary(AFT.weibull4)

