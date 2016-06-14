## Summer 2016
## Moving on to the "real" analyses from the preliminary ones this spring
## Code adapted from BreedSurv.R

setwd("C:/Users/efeichtinger/LittleBlueDinos/Data")

library(survival)
library(car)
library(coxme)
library(kinship2)
library(SurvRegCensCov)

## Read in April Census Dates
aprilD <- read.csv("April_Census_Dates.csv")

## Convert to date object
aprilD$CensusDate <- as.Date(aprilD$CensusDate, format = "%m/%d/%Y")

## First data table to deal with, known age birds to figure out who died in the first year
## Read in CSV file 
bird.df <- read.csv("Erin_June_FY.csv")
str(bird.df)
## 2371 individuals 

## Change column name to Days for days lived and add a year column
colnames(bird.df)[8] <- "Days"
bird.df["Yrs"] <- bird.df$Days/365.25

## Subset those birds who more than 365 days to find those who became helpers
prebrdr <- subset(bird.df, Days >= 365)
str(prebrdr)
## 1009 individuals 

## How to the April Census Dates that are the year following fledge year?? 
## Can R even do this? Well, regardless, first step is to convert to dates

## convert dates to date format
bird.df$FldgDate <- as.Date(bird.df$FldgDate, format = "%m/%d/%Y")
bird.df$LastObsDate <- as.Date(bird.df$LastObsDate, format = "%m/%d/%Y")

## Add column to indicate whether or not bird was alive after one year 
## 0 = alive, 1 = dead (same convention as censorship indicator)
bird.df["DeadAlive"] <- 0
## If years lived is less than 1, dead
bird.df$DeadAlive[which(bird.df$Yrs < 1)] <- 1

## Change back to numeric for survival object 
bird.df$FldgDate <- as.numeric(bird.df$FldgDate)
bird.df$LastObsDate <- as.numeric(bird.df$LastObsDate)

#very important piece of code for the model to work properly, remove any 
#weird entries like birds that have negative years of experience or a negative
#survival interval 
bird.df <- subset(bird.df, bird.df$Days > 0)

## I don't think this is truly the way to do this, is there a way to restrict
## this to one year?? Because the birds are not truly censored if they are 
## alive after one year, so I don't know how to do this 
bird.ob <- Surv(bird.df$Days, bird.df$DeadAlive, type =c('right'))

jay.lifetab <- survfit(bird.ob~1, conf.type = "log-log")
jay.fit <- plot(jay.lifetab, xlab = "Time (years)", 
                ylab = "Cumulative Survival")

