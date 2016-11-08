## 9 13 2016 Major modifications

## 7 1 2016 - First version  
## Script for first year birds 
## Source code from "JaySurv.R" 

library(survival)
library(ggplot2)
library(car)
library(coxme)
library(kinship2)
library(SurvRegCensCov)
library(parfm)
library(plyr)
library(corrplot)

##########################################################################

# Block 1 - Input data, manipulate, make survival object

## Read in April Census Dates
aprilD <- read.csv("April_Census_Dates.csv")

## Convert to date object
aprilD$CensusDate <- as.Date(aprilD$CensusDate, format = "%m/%d/%Y")

## This is the full data set with all known age birds 
#bird.df <- read.csv("Erin_June_All_Birds.csv")
#bird.df <- read.csv("Erin_June_All_BirdsMass.csv")
bird.df <- read.csv("Erin_Oct_FY.csv")
bird.df<- bird.df[!duplicated(bird.df),]


## Finding number of NA's
#sum(is.na(bird.df$HatchNum))

##Remove NA's for hatch number 
#bird.df <- bird.df[!is.na(bird.df$HatchNum),]

str(bird.df)


################################################################################

#Mean and sd number hatchlings
mean(bird.df$HatchNum)
sd(bird.df$HatchNum)

#Mean and sd number fledglings 
mean(bird.df$FldgNum)
sd(bird.df$FldgNum)

hmean <- mean(bird.df$HatchNum)
fmean <- mean(bird.df$FldgNum)
hsd <- sd(bird.df$HatchNum)
fsd <- sd(bird.df$FldgNum)

hfmeans <- data.frame(hmean, fmean)
sds <- c(hsd, fsd)
mnsd <- rbind(hfmeans, hsd)
colnames(mnsd) <- c("Hatchlings","Fledglings")
rownames(mnsd) <- c("mean","sd")

mnsd
############################################################################

## convert dates to date format
bird.df$FldgDate <- as.Date(bird.df$FldgDate, format = "%m/%d/%Y")
bird.df$LastObsDate <- as.Date(bird.df$LastObsDate, format = "%m/%d/%Y")

bird.df$MeasDate <- as.Date(bird.df$MeasDate, format = "%m/%d/%Y")
bird.df$HatchDate <- as.Date(bird.df$HatchDate, format = "%m/%d/%Y")
bird.df["MeasDiff"] <- bird.df$MeasDate - bird.df$HatchDate

bird.df["Days"] <- bird.df$LastObsDate - bird.df$FldgDate
bird.df["Yrs"] <- bird.df$Days/365.25
bird.df$Days <- as.numeric(bird.df$Days)
bird.df$Yrs <- as.numeric(bird.df$Yrs)
bird.df$Days <- round(bird.df$Days, digits = 2)
bird.df$Yrs <- round(bird.df$Yrs, digits = 2)

bird.df <- subset(bird.df, bird.df$MeasDiff == 11)

## add column for censorship status, in survival package - 0=alive, 1=dead
bird.df["Censor"] <- 1
#Birds that are right censored (indicated by 0) are still alive after 1 year
bird.df$Censor[which(bird.df$Yrs >= 1)]<-0

#year <- as.POSIXlt(bird.df$FldgDate)$year+1900
#bird.df["Cohort"] <- year
#bird.df$Cohort <- as.factor(bird.df$Cohort)

#subset to get rid of years less than 0
yrlg.df <- subset(bird.df, bird.df$Yrs > 0 & bird.df$Days > 0)

##problem with the 2015 cohort since my condition was to be still alive at
##one year, however, the 2015 cohort 
yrlg.df$Censor[which(yrlg.df$LastObsDate == "2016-04-12")] <- 0

##Just to check out 2015 cohort for correctness 
#cohort15 <- subset(yrlg.df, yrlg.df$Cohort == "2015")

#change back to numeric for survival object 
yrlg.df$FldgDate <- as.numeric(yrlg.df$FldgDate)
yrlg.df$LastObsDate <- as.numeric(yrlg.df$LastObsDate)

yrlg.df$NestYear <- as.factor(yrlg.df$NestYear)

#survival object 
yrlg.ob <- Surv(yrlg.df$Yrs, yrlg.df$Censor, type = c('right'))

my.fit <- survfit(yrlg.ob~1, conf.type = "log-log")
my.fit
str(my.fit)

my.fit$surv

#stepfunction
kms <- survfit(yrlg.ob ~ 1)
survest <- stepfun(kms$time, c(1, kms$surv))

## Plot survival curve from KM estimate (my.fit)
plot.fit <- plot(my.fit, xlab = "Time (years)",
        log = "y", ylim = c(0.4, 1),xlim=c(0,1), ylab = "Cumulative Survival", 
        main = "Fledge to 1Yr")

## KM for sex 
my.fit2 <- survfit(yrlg.ob ~ yrlg.df$Sex, conf.type = "log-log")
my.fit2
str(my.fit2)

plot.sex <- plot(my.fit2, col = c("navy","red"),
                 xlab="years", log="y",  ylim = c(0.4, 1),xlim=c(0,1),
                 ylab = "Cumulative survival", main = "Yearlings by sex")
legend("topright", c("Females","Males"), col=c("navy","red"),
       lty = c(1,2),lwd=1)

#Year to year estimates 
fit.yr <- survfit(yrlg.ob ~ yrlg.df$NestYear)
fit.yr
p.yr <- plot(fit.yr, xlab = "Time", log= "y", xlim=c(0,1), ylim=c(0.1,1))

#Life table construction 
deaths <- c(2,5,1,7,2,2,0,9,90,4,59,3,12,7,10,6,4,2,55,
            38,33,59,28,67,60,74,30,98,39,91,122,15,47,72,49)
inds <- fit.yr$n
inds <- cbind(inds)
deaths <- cbind(deaths)
life.table <- data.frame(inds, deaths)
colnames(life.table) <- c("Inds", "Deaths")
life.table$Deaths <- as.numeric(life.table$Deaths)
life.table$Inds <- as.numeric(life.table$Inds)

life.table["p"] <- 1-(life.table$Deaths/life.table$Inds)
life.table$p <- round(life.table$p, 2)
life.table["Cohort"] <- c(1981,1982,1983,1984,1985,1986,1987,1988,1989,
                  1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,
                  2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,
                  2010,2011,2012,2013,2014,2015)
life.table$Cohort <- as.factor(life.table$Cohort)

##Fairly nice figures - change color later 
p1 <- ggplot(life.table, aes(x=Cohort, y=p)) +
 geom_bar(stat="identity", fill="black") +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x = "Cohort", y = "Proportion of Survivors")

p2 <- ggplot(life.table, aes(x=Cohort, y=Inds))+ 
  geom_bar(stat="identity", fill="gray47") +
  theme(axis.text.x=element_text(angle=300,hjust=0, size=12)) +
  theme(axis.ticks.x=element_blank()) +
  labs(x = "Cohort", y= "Number of Fledglings Produced")

grid.arrange(p1, p2)

m.cumuls <- mean(life.table$p)
var.cumuls <- var(life.table$p)
sd.cumuls <- sd(life.table$p)


m.cumuls
var.cumuls
sd.cumuls

surv.stats <- cbind(m.cumuls,var.cumuls, sd.cumuls)
colnames(surv.stats) <- c("Mean","Variance","Standard Deviation")

surv.stats[,1:3] <- round(surv.stats[,1:3], 3)
surv.stats

range(life.table$p)

#########################################################################

# Block 2 - fit Cox PH and AFT models 

## Cox PH Models
cox.intercept <- coxph(yrlg.ob ~ 1, data = yrlg.df)
cox.intercept

#Model 1 with sex as only predictor
cox.sex <- coxph(yrlg.ob ~ Sex, data = yrlg.df)
summary(cox.sex)
#Check PH assumption
cox.zph(cox.sex)
plot(cox.zph(cox.sex))


#Model 2 with number of hatchlings in nest 
cox.hatch <- coxph(yrlg.ob ~ HatchNum, data = yrlg.df)
summary(cox.hatch)
#Check PH assumption
cox.zph(cox.hatch)
plot(cox.zph(cox.hatch))


#Model 3 with number of fledglings from nest 
cox.flg <- coxph(yrlg.ob ~ FldgNum, data= yrlg.df)
summary(cox.flg)
#Check PH assumption
cox.zph(cox.flg)
plot(cox.zph(cox.flg))

#Model 4 with day 11/nestling mass as predictor
#nestling mass (day 11) 
cox.mass <- coxph(yrlg.ob ~ Weight, data = yrlg.df)
summary(cox.mass)
#Check PH assumption
cox.zph(cox.mass)
plot(cox.zph(cox.mass))

#Model 5 with cohort year as predictor - not the best way to do this
cox.cohort <- coxph(yrlg.ob ~ NestYear, data = yrlg.df)
cox.cohort$loglik

#Model Loglik converged before variable 6. beta may be infinite
#So the beta estimates are not reliable but it does show first year surv
#is influenced by year, almost each year coefficient has p < 0.05 
#But I think frailty models will be better here for this 

##Model without cohort because it is not the best way to do it 
cox.mult<- coxph(yrlg.ob ~ Sex + HatchNum + FldgNum + Weight, 
                  data = yrlg.df)
summary(cox.mult)
cox.mult$loglik
anova(cox.mult)
anova(cox.mult, type=("II"))
anova(cox.mult, type=("III"))

cox.full <- coxph(yrlg.ob ~ Sex + HatchNum + FldgNum + Weight + Cohort, 
                   data = yrlg.df)
summary(cox.full)
cox.full$loglik

cox.mc <- coxph(yrlg.ob ~ Weight + NestYear, data = yrlg.df)
cox.mc$loglik
anova(cox.mc)

anova(cox.full)
anova(cox.cohort)

extractAIC(cox.sex)
extractAIC(cox.mass)
extractAIC(cox.hatch)
extractAIC(cox.flg)
extractAIC(cox.cohort)
extractAIC(cox.mult)
extractAIC(cox.full)

#best of these simple models is the cohort one according to AIC 

anova(cox.sex) #adds pretty much nothing to model (sex not good predictor)
anova(cox.mass) #reduction in loglik but not by much 
anova(cox.hatch)
anova(cox.flg)
anova(cox.cohort) #the biggest reduction in loglik

#################################################################
#Old as of 9 13 2016, using incorrect data and have adjusted 

# Model 3 with clutch size and sex, no interaction
cx.yrl3 <- coxph(yrlg.ob ~ Sex + HatchNum, data = yrlg.df)
summary(cx.yrl3)

## Check PH assumption
cox.zph(cx.yrl3)
plot(cox.zph(cx.yrl3))

# Model 4 interaction between clutch size and sex
cx.yrl4 <- coxph(yrlg.ob ~ ClutchNum + Sex + ClutchNum*Sex, data = yrlg.df)
summary(cx.yrl4)

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

########################################################################

# Block 3 - frailty models, mixed effects Cox models 
mm.int <- coxme(yrlg.ob ~ 1 + (1|Cohort), data = yrlg.df)
mm.int

#Mixed effects cox model with year as random effect
#First model - nest size fixed effect 
mm1 <- coxme(yrlg.ob ~ FldgNum + (1|Cohort), data = yrlg.df)
mm1

#Natal nest as random effect, does this make sense with clutch size as fixed?
mm2 <- coxme(yrlg.ob ~ FldgNum + (1|NatalNest), data = yrlg.df)
mm2

#Day 11 mass as fixed effect
mm3 <- coxme(yrlg.ob ~ Weight + (1|Cohort), data = yrlg.df)
mm3

mm4 <- coxme(yrlg.ob ~ Weight + (1|NatalNest), data = yrlg.df)
mm4


#Frailty models using functions in the survival package 
f1 <- coxph(yrlg.ob ~ Sex + HatchNum + 
              frailty(Cohort, dist='gamma'), data = yrlg.df)
f2 <- coxph(yrlg.ob ~ HatchNum + frailty(Cohort, dist='gamma'), data = yrlg.df)
summary(f1)
summary(f2)
extractAIC(f1)
extractAIC(f2)

#AFT models
AFT.weibull <- survreg(yrlg.ob ~ HatchNum + frailty(Cohort, dist='gamma'), 
                        data = yrlg.df, dist = "weibull")
AFT.weibull
summary(AFT.weibull)

AFT.int <- survreg(yrlg.ob ~ 1 + frailty(Cohort, dist='gamma'),
              data = yrlg.df, dist = "weibull")
summary(AFT.int)

#Mass and cohort as frailty
AFT.mass <- survreg(yrlg.ob ~ Weight + frailty(Cohort, dist='gamma'),
                    data = yrlg.df, dist ="weibull")
summary(AFT.mass)


#Mass as fixed, natal nest as the frailty term
AFT.mnest <- survreg(yrlg.ob ~ Weight + frailty(NatalNest, dist='gamma'),
                     data=yrlg.df, dist="weibull")
summary(AFT.mnest)

########################################################################
# Block 4 - Input of Vegetation, fire, and terr size data
# Manipulation for use in models

dom.veg <- read.csv("dom_veg.csv")
tsf <- read.csv("tsf_terr.csv")
terr <- read.csv("terr_size.csv")

## include terrs with zero scrub and zero patches in 1-9 fire window?
## extract too
#

#Subset data to only keep scrub veg types

#create object to store charc strings corresponding to scrub types
keep <- c("RSh", "SFi", "SFl", "SFx", "SSo", "SSr")
#creat new data frame with only the scrub types in "keep"
vegdf <- dom.veg[dom.veg$Dom.Veg %in% keep, ]

#Is this correct?? - seems so 
#Creating another new data frame from "vegdf" with summarise argument
#Summing the values for cell counts of each scrub type by terryr
scrub.terr <- ddply(vegdf, .(Terryr), summarise, 
                    Count.Dom.Veg=sum(Count.Dom.Veg))
colnames(scrub.terr) <- c("TerrYr", "ScrubCount")

#TSF data
#Create object for the numbers, same logic as with veg data
keep2 <- c(1,2,3,4,5,6,7,8,9)
firedf <- tsf[tsf$TSF_years %in% keep2, ]
tsf.terr <- ddply(firedf, .(TERRYR), summarise, CellCount=sum(CellCount))
colnames(tsf.terr) <- c("TerrYr", "FireCount")


#Territory size, cell counts by terryr
#terr - already has one count per terryr
str(terr)
terr.size <- terr

## Need to add a terryear variable, can do this with string manipulation
## With this, I can add group size, territory size (cells), cells of oak scrub
## cells of tsf in 1-9 
## The final piece is the pedigree 