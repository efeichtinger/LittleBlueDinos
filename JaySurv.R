## Summer 2016
## Moving on to the "real" analyses from the preliminary ones this spring
## Code adapted from BreedSurv.R


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

## First data table to deal with, known age birds to figure out who died in the first year
## Read in CSV file
## June 28 It says all birds but really it can be used for first year
## June 28 added clutch size 
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

#### Code above is necessary to start for any stage

########################################################
## Yearlings 
## add column for censorship status, in survival package - 0=alive, 1=dead
bird.df["Censor"] <- 1
#Birds that are right censored (indicated by 0) are still alive after 1 year
#post fledgling 
bird.df$Censor[which(bird.df$Yrs > 1)]<-0
bird.df<- subset(bird.df, bird.df$Yrs > 0)
#subset to get rid of years less than 0
yrlg.df <- subset(bird.df, bird.df$Yrs > 0)

year <- as.POSIXlt(bird.df$FldgDate)$year+1900
yrlg.df["FYear"] <- year
yrlg.df$FYear <- as.factor(yrlg.df$FYear)


#change back to numeric for survival object 
yrlg.df$FldgDate <- as.numeric(yrlg.df$FldgDate)
yrlg.df$LastObsDate <- as.numeric(yrlg.df$LastObsDate)

#survival object 
yrlg.ob <- Surv(yrlg.df$Yrs, yrlg.df$Censor, type = c('right'))

my.fit <- survfit(yrlg.ob~1, conf.type = "log-log")
my.fit
str(my.fit)

plot.fit <- plot(my.fit, xlab = "Time (years)",
      log = "y", ylim = c(0.4, 1),xlim=c(0,1), ylab = "Cumulative Survival", 
                   main = "Yearlings")

## KM for sex 
my.fit2 <- survfit(yrlg.ob ~ yrlg.df$Sex, conf.type = "log-log")
my.fit2
str(my.fit2)

plot.sex <- plot(my.fit2, col = c("navy","red"),
                 xlab="years", log="y",  ylim = c(0.4, 1),xlim=c(0,1),
                 ylab = "survival", main = "Yearlings by sex")
legend("topright", c("Females","Males"), col=c("navy","red"),
       lty = c(1,2),lwd=1)

### Appears to be no difference between sexes for this first year, which 
# isn't surprising, but right before 1 year old there is a difference

## Cox PH model with sex as a predictor 
cx.yrl <- coxph(yrlg.ob ~ Sex, data = yrlg.df)
summary(cx.yrl)
extractAIC(cx.yrl)
#Check PH assumption
cox.zph(cx.yrl)
plot(cox.zph(cx.yrl))

#Add clutch size as a predictor 
cx.yrl2 <- coxph(yrlg.ob ~ ClutchNum, data = yrlg.df)
summary(cx.yrl2)
extractAIC(cx.yrl2)
#Check PH assumption
cox.zph(cx.yrl2)
plot(cox.zph(cx.yrl2))

#Clutch size and sex
cx.yrl3 <- coxph(yrlg.ob ~ Sex + ClutchNum, data = yrlg.df)
summary(cx.yrl3)
extractAIC(cx.yrl3)
#Check PH assumption
cox.zph(cx.yrl3)
plot(cox.zph(cx.yrl3))

#Interaction between clutch size and sex
cx.yrl4 <- coxph(yrlg.ob ~ ClutchNum + Sex + ClutchNum*Sex, data = yrlg.df)
cx.yrl4
extractAIC(cx.yrl4)

#Cohort year, approach where it is a factor - trying frailty models next
cx.yrl4 <- coxph(yrlg.ob ~ FYear, data = yrlg.df)
extractAIC(cx.yrl4)

#Cohort year and clutch size
cx.yrl5 <- coxph(yrlg.ob ~ FYear + ClutchNum, data = yrlg.df)
extractAIC(cx.yrl5)

#cohort year, clutch size and sex
cx.yrl6 <- coxph(yrlg.ob ~ FYear + ClutchNum + Sex, data = yrlg.df)
extractAIC(cx.yrl6)
plot(cox.zph(cx.yrl6))

#parfm package 
#frailty model 
fr.mod <- parfm(yrlg.ob ~ ClutchNum + Sex, cluster = "FYear", 
                data = yrlg.df, dist = "weibull", frailty = "gamma")
#error - no rows to aggregate

#Mixed effects cox model with year as randome effect
mm1 <- coxme(yrlg.ob ~ ClutchNum + (1|FYear), data = yrlg.df)
mm2 <- coxme(yrlg.ob ~ Sex + (1|FYear), data = yrlg.df)



######################################################################
## All known age birds 
## June 21 2016 
## Using the known-age bird data subset to understand how to get all the 
## important stuff from hazard models including the survival function
## Following "Survival Analysis in R" David M Diez

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
           xlab="Years", ylab="survival (log)",ylim = c(0.001,2))

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
p3 <- plot(survfit(mod1), xlab = "years")

#Tests for the PH assumption - tests based on scaled Schoenfeld residuals
cox.zph(mod1)
plot(cox.zph(mod1))

myob <- survfit(mod1)
plot(myob)


###########################################################################
## Subset those birds who more than 365 days to find those who became helpers
## June 21 - note that this includes birds that did go on to be breeders
## See notes on the challenges of getting the right data

prebrdr <- subset(bird.df, Days >= 365)
str(prebrdr)
## 1009 individuals 

## How to the April Census Dates that are the year following fledge year?? 
## Can R even do this? Well, regardless, first step is to convert to dates

## convert dates to date format
prebrdr$FldgDate <- as.Date(prebrdr$FldgDate, format = "%m/%d/%Y")
prebrdr$LastObsDate <- as.Date(prebrdr$LastObsDate, format = "%m/%d/%Y")

## Add column to indicate whether or not bird was alive after one year 
## 0 = alive, 1 = dead (same convention as censorship indicator)
prebrdr["Censor"] <- 1
## If years lived is less than 1, dead
prebrdr$Censor[which(prebrdr$LastObsDate =="2016-4-12")] <- 0

## Change back to numeric for survival object 
prebrdr$FldgDate <- as.numeric(prebrdr$FldgDate)
prebrdr$LastObsDate <- as.numeric(prebrdr$LastObsDate)

pre.ob <- Surv(prebrdr$Yrs, prebrdr$Censor, type =c('right'))

#KM curve 
pre.fit <- survfit(pre.ob~1, conf.type = "log-log")
jay.fit <- plot(pre.fit,log="y", xlim = c(1,14), ylim=c(0.1,1), xlab = "Time (years)", 
                ylab = "Cumulative Survival")
#Km curve by sex
sex.fit <- survfit(pre.ob~ Sex, data = prebrdr)
sex.plot <- plot(sex.fit, log="y", xlim = c(1,14), ylim=c(0.01,1), xlab = "Time (years)", 
                 ylab = "Cumulative Survival")

model1 <- coxph(pre.ob ~ Sex, data = prebrdr)
summary(model1)

cox.zph(model1)
plot(cox.zph(model1))
