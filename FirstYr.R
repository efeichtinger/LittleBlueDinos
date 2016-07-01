## July 1 2016 (HOLY SHIT)
## Script for first year birds 
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

## This is the full data set with all known age birds 
bird.df <- read.csv("Erin_June_All_Birds.csv")
bird.df<- bird.df[!duplicated(bird.df),]
bird.df$ClutchNum <- as.numeric(bird.df$ClutchNum)

str(bird.df)

## Change column name to Days for days lived and add a year column
colnames(bird.df)[8] <- "Days"
bird.df["Yrs"] <- bird.df$Days/365.25

## convert dates to date format
bird.df$FldgDate <- as.Date(bird.df$FldgDate, format = "%m/%d/%Y")
bird.df$LastObsDate <- as.Date(bird.df$LastObsDate, format = "%m/%d/%Y")

## add column for censorship status, in survival package - 0=alive, 1=dead
bird.df["Censor"] <- 1
#Birds that are right censored (indicated by 0) are still alive after 1 year
bird.df$Censor[which(bird.df$Yrs > 1)]<-0

year <- as.POSIXlt(bird.df$FldgDate)$year+1900
bird.df["FYear"] <- year
bird.df$FYear <- as.factor(bird.df$FYear)

#subset to get rid of years less than 0
yrlg.df <- subset(bird.df, bird.df$Yrs > 0 & bird.df$Days > 0)

#change back to numeric for survival object 
yrlg.df$FldgDate <- as.numeric(yrlg.df$FldgDate)
yrlg.df$LastObsDate <- as.numeric(yrlg.df$LastObsDate)

#survival object 
yrlg.ob <- Surv(yrlg.df$Yrs, yrlg.df$Censor, type = c('right'))

my.fit <- survfit(yrlg.ob~1, conf.type = "log-log")
my.fit
str(my.fit)

## Plot survival curve from KM estimate (my.fit)
plot.fit <- plot(my.fit, xlab = "Time (years)",
        log = "y", ylim = c(0.4, 1),xlim=c(0,1), ylab = "Cumulative Survival", 
        main = "Yearlings")

## KM for sex 
my.fit2 <- survfit(yrlg.ob ~ yrlg.df$Sex, conf.type = "log-log")
my.fit2
str(my.fit2)

plot.sex <- plot(my.fit2, col = c("navy","red"),
                 xlab="years", log="y",  ylim = c(0.4, 1),xlim=c(0,1),
                 ylab = "Cumulative survival", main = "Yearlings by sex")
legend("topright", c("Females","Males"), col=c("navy","red"),
       lty = c(1,2),lwd=1)


## Cox PH Models

# Model 1 with sex as only predictor
cx.yrl <- coxph(yrlg.ob ~ Sex, data = yrlg.df)
summary(cx.yrl)

## Check PH assumption
cox.zph(cx.yrl)
plot(cox.zph(cx.yrl))

# Model 2 with clutch size as single predictor 
cx.yrl2 <- coxph(yrlg.ob ~ ClutchNum, data = yrlg.df)
summary(cx.yrl2)

## Check PH assumption
cox.zph(cx.yrl2)
plot(cox.zph(cx.yrl2))

# Model 3 with clutch size and sex, no interaction
cx.yrl3 <- coxph(yrlg.ob ~ Sex + ClutchNum, data = yrlg.df)
summary(cx.yrl3)

## Check PH assumption
cox.zph(cx.yrl3)
plot(cox.zph(cx.yrl3))

# Model 4 interaction between clutch size and sex
cx.yrl4 <- coxph(yrlg.ob ~ ClutchNum + Sex + ClutchNum*Sex, data = yrlg.df)
cx.yrl4

## Check PH assumption
cox.zph(cx.yrl4)
plot(cox.zph(cx.yrl4))

# Model 5 cohort year as factor
cx.yrl5 <- coxph(yrlg.ob ~ FYear, data = yrlg.df)

# Model 6 cohort year and clutch size 
cx.yrl6 <- coxph(yrlg.ob ~ FYear + ClutchNum, data = yrlg.df)

# Model 7 cohort year, clutch size and sex
cx.yrl7 <- coxph(yrlg.ob ~ FYear + ClutchNum + Sex, data = yrlg.df)


## Get AIC for models 1-7
extractAIC(cx.yrl)
extractAIC(cx.yrl2)
extractAIC(cx.yrl3)
extractAIC(cx.yrl4)
extractAIC(cx.yrl5)
extractAIC(cx.yrl6)
extractAIC(cx.yrl7)

#Mixed effects cox model with year as randome effect
mm1 <- coxme(yrlg.ob ~ ClutchNum + (1|FYear), data = yrlg.df)
mm2 <- coxme(yrlg.ob ~ Sex + (1|FYear), data = yrlg.df)
mm3 <- coxme(yrlg.ob ~ Sex + (ClutchNum|FYear), data = yrlg.df)
#mm4 does not work 
mm4 <- coxme(yrlg.ob ~ ClutchNum + (Sex|FYear), data = yrlg.df)

#Frailty models using functions in the survival package 
f1 <- coxph(yrlg.ob ~ Sex + ClutchNum + 
              frailty(FYear, dist='gamma'), data = yrlg.df)
f2 <- coxph(yrlg.ob ~ ClutchNum + frailty(FYear, dist='gamma'), data = yrlg.df)
summary(f2)
extractAIC(f1)
extractAIC(f2)




