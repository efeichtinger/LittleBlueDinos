### 1981 to 2015 

## Started 10-18-2016

library(survival)
library(ggplot2)
library(car)
library(coxme)
library(kinship2)
library(plyr)
library(corrplot)

#1981-2015
aprilD <- read.csv("April_Census_Dates.csv")
aprilD$CensusDate <- as.Date(aprilD$CensusDate, format = "%m/%d/%Y")

jays <- read.csv("Erin_HY_1981.csv")

jays$NestYear <- as.factor(jays$NestYear)

#Convert dates to date format 
jays$FldgDate <- as.Date(jays$FldgDate, format = "%m/%d/%Y")
jays$LastObsDate <- as.Date(jays$LastObsDate, format = "%m/%d/%Y")
jays$HatchDate <- as.Date(jays$HatchDate, format = "%m/%d/%Y")
jays$MeasDate <- as.Date(jays$MeasDate, format = "%m/%d/%Y")

jays["Days"] <- jays$LastObsDate - jays$FldgDate
jays["Years"] <- round(jays$Days/365.25, digits = 2)

str(jays)

#Add censorship column, 1 if dead before 1 yr, 0 if yr > 1 or 4/12/12016
jays["Censor"] <- 1
jays$Censor[which(jays$Years >= 1)]<-0
yrlg.df <- subset(jays, jays$Years > 0 & jays$Days > 0)
yrlg.df$Censor[which(yrlg.df$LastObsDate == "2016-04-12")] <- 0

#Jay ID not unique but the band number is 
#yrlg.df[,2] <- NULL

str(yrlg.df)

yrlg.df["Day11"] <- yrlg.df$MeasDate - yrlg.df$HatchDate
yrlg.df <- subset(yrlg.df, yrlg.df$Day11 <= 13)

yrlg.df[,8] <- NULL
yrlg.df[,15] <- NULL
colnames(yrlg.df)[7] <- "Mass"

group.size <- read.csv("groupsize1981.csv")

#String function and regular expressions to get a "TerrYr" variable to match
#Remember this is just for 1999 to 2015, the territory data goes to 1981
new.col <- gsub(".$","",yrlg.df$NatalNest)
colTY <- as.vector(new.col)
colTY <- cbind(colTY)
yrlg.df["TerrYr"] <- colTY

yrlg.df <- merge(yrlg.df, group.size, by="TerrYr")

## Add in territory info
#Scrub data
dom.veg <- read.csv("dom_veg.csv")
dom.veg <- subset(dom.veg, InABS == TRUE)
#Subset data to only keep scrub veg types

#create object to store charc strings corresponding to scrub types
keep <- c("RSh", "SFi", "SFl", "SFx", "SSo", "SSr")
#creat new data frame with only the scrub types in "keep"
vegdf <- dom.veg[dom.veg$Dom.Veg %in% keep, ]

scrub.terr <- ddply(vegdf, .(Terryr), summarise, 
                    Count.Dom.Veg=sum(Count.Dom.Veg))
colnames(scrub.terr) <- c("TerrYr", "Scrub")

no.scr <- subset(dom.veg, !(Terryr %in% scrub.terr$TerrYr))
no.scr["scrb.count"] <- 0

#Keep only terryr and scrb.count
vars <- c("Terryr","scrb.count")
no.scr <- no.scr[vars]

#remove duplicate rows 
no.scr <- no.scr[!duplicated(no.scr),]
colnames(no.scr)[1] <- "TerrYr"
colnames(scrub.terr)[2] <- "scrb.count"

#This includes terryears from 1981 to 2015, have to add year
#String operations? 
scr.ct <- rbind(scrub.terr, no.scr)

#Time since fire data
tsf <- read.csv("tsf_terr.csv")
tsf <- subset(tsf, InABS == TRUE)

#TSF data
#Create object for the numbers, same logic as with veg data
keep2 <- c(1,2,3,4,5,6,7,8,9)
firedf <- tsf[tsf$TSF_years %in% keep2, ]
tsf.terr <- ddply(firedf, .(TERRYR), summarise, CellCount=sum(CellCount))
colnames(tsf.terr) <- c("TerrYr", "FireCount")

no.tsf1 <- subset(tsf, !(TERRYR %in% tsf.terr$TerrYr))
no.tsf1["tsf.count"] <- 0

no.tsf1 <- no.tsf1[,c(1,8)]
colnames(no.tsf1)[1] <- "TerrYr"

colnames(tsf.terr)[2] <- "tsf.count"

#All TerrYrs including counts of 0
tsf.ct <- rbind(tsf.terr,no.tsf1)

##territory size
terr <- read.csv("terr_size.csv")
terr <- subset(terr, InABS == TRUE)

#Keep only terryr and scrb.count
vars1 <- c("TERRYR","Count")
terr <- terr[vars1]
colnames(terr) <- c("TerrYr", "TerrSize")

#territory size 
veg.size <- merge(scr.ct,terr)
#info on territory quality, cell count of scrub, size, tsf in 1-9 window
terr.info <- merge(veg.size, tsf.ct)

#remove duplicate rows 
terr.info <- terr.info[!duplicated(terr.info),]

yrlg.df <- merge(yrlg.df, terr.info, by="TerrYr")

###Add code to see who drops out when the territory info is added to check
#for bias by year (i.e. more drop out as prop of total in some years)
#Who is not in the new df?


## Rearrange columns 
colnames(yrlg.df)[16] <- "GroupSize"
colnames(yrlg.df)[18] <- "OakScrub"
colnames(yrlg.df)[20] <- "TSF"

yrlg.df["stdscr"] <- scale(yrlg.df$OakScrub, center = FALSE, scale = TRUE)
yrlg.df["stdtsf"] <- scale(yrlg.df$TSF, center = FALSE, scale = TRUE)
yrlg.df["stdsize"] <- scale(yrlg.df$TerrSize, center= FALSE, scale = TRUE)


# Data frame of covariates for correlation matirx
covars.stad <- yrlg.df[,c(8,11,12,16,17,21,22,23)]
covars.no <- yrlg.df[,c(8,11,12,16,17,18,19,20)]


corrs <- cor(covars.stad, use="complete.obs")
corrplot(corrs, method="pie", type="lower")

corrs2 <- cor(covars.no, use="complete.obs")
corrplot(corrs2, method="pie", type="lower")


#Change to numeric for survival object 
yrlg.df$FldgDate <- as.numeric(yrlg.df$FldgDate)
yrlg.df$LastObsDate <- as.numeric(yrlg.df$LastObsDate)

yrlg.df$Years <- as.numeric(yrlg.df$Years)
yrlg.df$Days <- as.numeric(yrlg.df$Days)

yrlg.ob <- Surv(yrlg.df$Years, yrlg.df$Censor, type = c('right'))
my.fit <- survfit(yrlg.ob~1, conf.type = "log-log")

plot.fit <- plot(my.fit, xlab = "Time (years)", conf.int=TRUE,
        log = "y", ylim = c(0.4, 1),xlim=c(0,1), ylab = "Cumulative Survival", 
        main = "Fledge to 1Yr - 1981 to 2015")

fit.yr <- survfit(yrlg.ob ~ yrlg.df$NestYear)
plot.yr <- plot(fit.yr,xlab = "Time (years)",
                log = "y", ylim = c(0.1, 1),xlim=c(0,1), ylab = "Cumulative Survival", 
                main = "Curves for each year 1981 - 2015")
fit.yr

### Make a life table 

## Models 

#Cohort year
cox1 <- coxph(yrlg.ob ~ NestYear, data= yrlg.df)
cox1
anova(cox1)

#Day 11 Mass
cox.wt <- coxph(yrlg.ob ~ Mass, data = yrlg.df)

#Cohort and Mass
cox1b <- coxph(yrlg.ob ~ Mass + NestYear, data=yrlg.df)
anova(cox1b)

#Mixed effects with year as random effect
cox2 <- coxme(yrlg.ob ~ Mass + (1|NestYear), data= yrlg.df)
cox2
anova(cox2)

#Mixed effects with nestID as random effect 
cox3 <- coxme(yrlg.ob ~ Mass + (1|NatalNest), data =yrlg.df)
cox3
anova(cox3)

#Interestingly, there is a lot of variation from year to year, but the effect
#of cohort year has a relatively small variance 
#So I suppose there could be a lot of variation in p during the first year but 
#the main source of variation may not be the stochastic (and other) effects
#represented by using cohort year identity 
#High variance for natal nest ID as random term 


#Checking other social factors
coxgrp <- coxph(yrlg.ob ~ GroupSize, data=yrlg.df)
coxgrp
anova(coxgrp)

coxgm <- coxph(yrlg.ob ~ Mass + GroupSize, data=yrlg.df)
coxgm
anova(coxgm)

coxterr <- coxph(yrlg.ob ~ OakScrub + TerrSize + GroupSize + GroupSize:TerrSize, 
                 data = yrlg.df)

coxterr
anova(coxterr)
