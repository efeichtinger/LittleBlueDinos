## August 31 2016
## Separate scripts for each age/stage based method 
## All known age birds

library(survival)
library(car)
library(coxme)
library(kinship2)
library(SurvRegCensCov)
library(parfm)


bird.df <- read.csv("Erin_June_All_Birds.csv")
bird.df<- bird.df[!duplicated(bird.df),]
bird.df$ClutchNum <- as.numeric(bird.df$ClutchNum)


str(bird.df)
## 2370 individuals 

## Change column name to Days for days lived and add a year column
colnames(bird.df)[8] <- "Days"
bird.df["Yrs"] <- bird.df$Days/365.25

## convert dates to date format
bird.df$FldgDate <- as.Date(bird.df$FldgDate, format = "%m/%d/%Y")
bird.df$LastObsDate <- as.Date(bird.df$LastObsDate, format = "%m/%d/%Y")


## Add column to indicate whether or not bird was alive after one year 
## 0 = alive, 1 = dead (same convention as censorship indicator)
bird.df["Censor"] <- 1
## If the date is the April 2016 census date, consider still alive as of then
bird.df$Censor[which(bird.df$LastObsDate == "2016-4-12")] <- 0

## Change back to numeric for survival object 
bird.df$FldgDate <- as.numeric(bird.df$FldgDate)
bird.df$LastObsDate <- as.numeric(bird.df$LastObsDate)

bird.df <- subset(bird.df, bird.df$Days > 0 & bird.df$Yrs > 0)

## Create survival object - only birds that are right censored are still alive
jay.ob <- Surv(bird.df$Yrs, bird.df$Censor, type =c('right'))

## Kaplan-Meier estimate
survfit(jay.ob~1)

my.fit<- survfit(jay.ob~1)
#summary(my.fit)$surv
## Get information on the KM estimate - returns a list
str(my.fit)

## Plot KM
p1 <- plot(my.fit, main="Kaplan-Meier estimate with 95% CI",log = "y",
           xlab="Years", ylab="survival (log)",xlim = c(1,15),
           ylim = c(0.001,1))

p1 <- plot(my.fit, main="Kaplan-Meier estimate with 95% CI",
           xlab="Years", ylab="survival (log)",xlim = c(1,15))

## KM for sex 
my.fit2 <- survfit(jay.ob ~ bird.df$Sex)
my.fit2
str(my.fit2)

## Plot KM estimate by sex
p2 <- plot(my.fit2, main = "Kaplan-Meier estimate with 95% CI",log = "y",
           xlab="Years", ylab="surv function (log)",ylim = c(0.001,2))

# Cox model assuming time constant covariates 
mod1 <- coxph(jay.ob ~ Sex, data = bird.df)
mod1
summary(mod1)

# estimate the distribution of survival times 
p3 <- plot(survfit(mod1), xlab = "years", ylab = "Cumulative Survival",
           main="All known age birds")

#Tests for the PH assumption - tests based on scaled Schoenfeld residuals
cox.zph(mod1)
plot(cox.zph(mod1))

myob <- survfit(mod1)
plot(myob)
