#### Breeders
#### July 21 2016

library(survival)
library(ggplot2)
library(car)
library(coxme)
library(kinship2)
library(SurvRegCensCov)
library(parfm)
library(plyr)

#################################################################

# Block 1 - Input data, manipulate, make survival object, fix basic Cox Model

## Read in file of breeders - known and unknown age
## One record per individual 

brd <- read.csv("Breeders.csv")

#remove duplicate rows 
brd <- brd[!duplicated(brd),]

colnames(brd)[1] <- "ID"
colnames(brd)[2] <- "Band"
colnames(brd)[3] <- "FY"
colnames(brd)[4] <- "Fbreed"


brd <- brd[,c(1,2,8,7,6,3,4,5)]

#minimum age at first breeding 
colnames(brd)[5] <- "MinAgeFBr"

#Set all blanks to NA
brd[brd==""] <- NA
#remove na
brd <- na.omit(brd)

brd$Fbreed <- as.Date(brd$Fbreed, format ="%m/%d/%Y")
brd$LastObsDate <- as.Date(brd$LastObsDate, format = "%m/%d/%Y")

#Years spent as a breeder, 365.25 to account for leap years 
brd["Yrs"] <- (brd$LastObsDate - brd$Fbreed)/365.25
brd$Yrs <- round(brd$Yrs, digits = 1)

#Add censorship column
brd["Censor"] <- 1
brd$Censor[which(brd$LastObsDate =="2016-4-12")] <- 0

brd <- subset(brd, brd$Yrs > 0)
brd <- subset(brd, brd$FY >= "1981")
brd$FY <- as.factor(brd$FY)



males <- subset(brd, Sex == "M")
females <- subset(brd, Sex == "F")
mean(males$Yrs)
sd(males$Yrs)
mean(females$Yrs)
sd(females$Yrs)

cat("No. of Known Age Birds = ", nrow(brd[brd$AgeKnown == "Y", ]))
cat("No. of Unknown Age Birds = ", nrow(brd[brd$AgeKnown == "N", ]))

cat("No. of Females = ", nrow(brd[brd$Sex == "F", ]))
cat("No. of Males = ", nrow(brd[brd$Sex == "M", ]))

cat("No. of unknown age females = ", nrow(females[females$AgeKnown == "N", ]))
cat("No. of unknown age males = ", nrow(males[males$AgeKnown == "N", ]))

cat("No. of known age females = ", nrow(females[females$AgeKnown == "Y", ]))
cat("No. of known age males = ", nrow(males[males$AgeKnown == "Y", ]))

#Change dates back to numeric for surv object 
brd$Fbreed <- as.numeric(brd$Fbreed)
brd$LastObsDate <- as.numeric(brd$LastObsDate)
brd$Yrs <- as.numeric(brd$Yrs)


#summary stats for breeder lifespan, i.e. time spent as a breeder
#not total lifespan 
mean(brd$Yrs)
var(brd$Yrs)
sd(brd$Yrs)

#mean age at first breeding 
mean(brd$MinAgeFBr)
sd(brd$MinAgeFBr)
var(brd$MinAgeFBr)

males <- subset(brd, Sex == "M")
females <- subset(brd, Sex == "F")
mean(males$Yrs)
sd(males$Yrs)
mean(females$Yrs)
sd(females$Yrs)

mean(males$MinAgeFBr)
sd(males$MinAgeFBr)

mean(females$MinAgeFBr)
sd(females$MinAgeFBr)

#mean age at first breeding 
mean(males$MinAgeFBr)
sd(males$MinAgeFBr)
mean(females$MinAgeFBr)
sd(females$MinAgeFBr)

ggplot(brd, aes(MinAgeFBr, fill = Sex)) +
  geom_histogram(binwidth = 0.2) +
  xlab("Minimum Age at First Breeding") +
  ylab("Frequency of Individuals") +
  labs(title = "1981 to 2016")

ggplot(brd, aes(MinAgeFBr, fill = Sex)) +
  geom_bar(position="dodge") +
  xlab("Minimum Age at First Breeding") +
  ylab("Frequency") +
  labs(title = "Distribution of Ages at First Breeding")


#qplot(Yrs, data=brd, geom="histogram")
#Histogram of breeding life span, color by sex
ggplot(brd, aes(Yrs, fill = Sex)) +
  geom_histogram(binwidth = 0.08) +
  xlab("Number of Years Bred") +
  ylab("Frequency of Individuals") + 
  labs(title = "1981 to 2016")

ggplot(brd, aes(Yrs, fill = Sex)) +
  geom_bar(position="dodge")

ggplot(males, aes(Yrs)) + 
  geom_histogram(binwidth = 0.1) +
  xlab("Number of Years Bred") +
  ylab("Frequency of Individuals") +
  labs(title = "Males") 

ggplot(females, aes(Yrs)) + 
  geom_histogram(binwidth = 0.1) +
  xlab("Number of Years Bred") +
  ylab("Frequency of Individuals") +
  labs(title = "Females")

brd.ob <- Surv(brd$Yrs, brd$Censor,type= c('right'))
brd.fit <- survfit(brd.ob ~ 1, conf.type = "log-log")
brd.fit

year.fit <- survfit(brd.ob ~ brd$FY, conf.type="log-log")
year.fit
plot(year.fit, log="y")

# adjust x axis it's breeding span 
kmplot <- plot(brd.fit, xlab="Breeding time span (years)", log = "y", 
               ylab = "Cumulative Survival (log)", main = "Breeders",
               ylim = c(0.01,1), xlim = c(0,15))

brd.sex <- survfit(brd.ob ~ brd$Sex, conf.type="log-log")
sxplot <- plot(brd.sex, conf.int=TRUE, col = c("darkblue","darkorange3"), 
               xlab = "Breeding time span (years)", log="y", 
               ylab = "Cumulative Survival (log)", main = "Breeders",
               ylim = c(0.01,1), xlim = c(0,15))
legend("topright", c("Females","Males"), col=c("darkblue","darkorange3"),
       lty = c(1,1),lwd=1)

yrplot <- plot(year.fit, xlab ="Breeding time span (years)", log = "y",
               ylab = "Cumulative Survival (log)", main = "Breeders", 
               ylim = c(0.02,1), xlim = c(0,15))

cox.null <- coxph(brd.ob ~ 1, data = brd)
summary(cox.null)

cx1 <- coxph(brd.ob ~ Sex, data = brd)
summary(cx1)
cox.zph(cx1)
plot(cox.zph(cx1))
anova(cx1)

cx2 <- coxph(brd.ob ~ MinAgeFBr, data = brd)
summary(cx2)
cox.zph(cx2)
plot(cox.zph(cx2))
anova(cx2)

cx3 <- coxph(brd.ob ~ FY, data = brd)
anova(cx3)
cx4 <- coxph(brd.ob ~ Sex + MinAgeFBr, data=brd)
anova(cx4)

cx5 <- coxph(brd.ob ~ Sex + MinAgeFBr + FY, data=brd)
anova(cx5)
cx6 <- coxph(brd.ob ~ Sex + MinAgeFBr + FY + Sex*MinAgeFBr, data = brd)
anova(cx6)


extractAIC(cx1)
extractAIC(cx2)
extractAIC(cx3)
extractAIC(cx4)
extractAIC(cx5)
extractAIC(cx6)

#Deviance 
anova(cx1)
anova(cx2)
anova(cx3)
anova(cx4)
anova(cx5)
anova(cx6)

anova(cx2,cx4)



#Analysis of deviance table comparing first three Cox PH models 
dev.compare <- anova(cx1, cx2, cx3, cx4, cx5, cx6, test="Chisq")
dev.compare

dev2 <- anova(cx3,cx5,cx6)
dev2

#Calculate residuals for Coxph fit 
resd <- residuals(cx1, type="deviance", collapse=brd$ID)

#Frailty models where cohort year is a random effect 
frail1 <- coxme(brd.ob ~ MinAgeFBr + (1|FY), data = brd)
summary(frail1)
frail2 <- coxme(brd.ob ~ Sex + (1|FY), data= brd)
summary(frail2)
frailnull <- coxme(brd.ob ~ (1|FY), data= brd)

## AFT frailty models 

AFT1 <- survreg(brd.ob ~ MinAgeFBr + frailty(FY, dist='gamma'), 
                       data = brd, dist = "weibull")
summary(AFT1)

AFT2 <- survreg(brd.ob ~ Sex + MinAgeFBr + frailty(FY, dist='gamma'),
                data=brd, dist = "weibull")
summary(AFT2)


####################################################################
# Block 2 - Input of Vegetation, fire, terr size, group size
# Time-varying covarites 
# Manipulation for use in models

dom.veg <- read.csv("dom_veg.csv")
dom.veg <- subset(dom.veg, InABS == TRUE)


##############################################################################

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
colnames(scrub.terr) <- c("TerrYr", "Dom.Veg")

##Use scrub.terr and dom.veg to find terryrs that dropped out? 
#I used this to find birds who made it 1 yr but not to breeding 
#Comparing two data frames 
#hlpr <- subset(bird.df, !(JayID %in% brdr.df$JayID))
#str(hlpr)

#I think this gives me the terryrs that have no oak scrub 
no.scr <- subset(dom.veg, !(Terryr %in% scrub.terr$TerrYr))
no.scr["scrb.count"] <- 0

#Keep only terryr and scrb.count
vars <- c("Terryr","scrb.count")
no.scr <- no.scr[vars]

#remove duplicate rows 
no.scr <- no.scr[!duplicated(no.scr),]
colnames(no.scr)[1] <- "TerrYr"
colnames(scrub.terr)[2] <- "scrb.count"

## add no.scr to scrub.terr

#Object with territories in ABS with cell count of oak scrub by terryr
#Six types of oak scrub, also included 19 terryrs with 0 cell count
scr.ct <- rbind(scrub.terr, no.scr)



########################################################################
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
#Do I want zeros or cell count
#I think GF said not to use proportional data
#In that case it would be zeros if no part of terryr has patches in the 1-9
#tsf window

head(tsf.terr)
head(no.tsf1)

names(tsf.terr)
names(no.tsf1)

#Keep only terryr and scrb.count
no.tsf1 <- no.tsf1[,c(1,8)]
colnames(no.tsf1)[1] <- "TerrYr"

colnames(tsf.terr)[2] <- "tsf.count"


tsf.ct <- rbind(tsf.terr,no.tsf1)

##########################################################################

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


#####################################################################
# Block 3 - Input of data for each year bred 
## Read in file for breeders all years, multiple records per individual
## Time varying Cox Models 

# Input of data with helpers at each terr year with breeder IDs
hlp.byterr <- read.csv("Erin_Helpers_ByTerryr.csv")
# Input of data with breeders all years 
brd.allyrs <- read.csv("breeders_allyears.csv")

colnames(hlp.byterr)[1] <- "NestYear"
colnames(brd.allyrs)[1] <- "JayID"

#Change dates to dates and calculate time lived 

brd.allyrs$BrdrDate <- as.Date(brd.allyrs$BrdrDate, format = "%m/%d/%Y")
brd.allyrs$LastObsDate <- as.Date(brd.allyrs$LastObsDate,format = "%m/%d/%Y")

brd.allyrs["Days"] <- brd.allyrs$LastObsDate - brd.allyrs$BrdrDate
brd.allyrs$Days <- as.numeric(brd.allyrs$Days)
brd.allyrs["Years"] <- brd.allyrs$Days/365.25

#Add censorship column
brd.allyrs["Censor"] <- 1
brd.allyrs$Censor[which(brd.allyrs$LastObsDate =="2016-4-12")] <- 0

all.brd <- merge(brd.allyrs, hlp.byterr)
#Duplicate records for multiple nests
#This seems to have worked = I checked a few jays to make sure matches
#are correct 

#remove duplicate rows mulitple nests per year
all.brd <- all.brd[!duplicated(all.brd),]

#losing rows but I don't understand why, I think bc no breeding at terr
#Data frame with most information needed 
all <- merge(all.brd, terr.info)
### seems to have worked.......

#reorganize data frame 
#two columns for sex and some other unnecessary cols 
#change order to better facilitate analyses

#brd <- brd[,c(1,2,8,7,6,3,4,5)]    #to rearrange cols
names(all)
all.info <- all[,c(2,4,8,9,5,6,7,18,11,12,14,15,16,1,17,21,22,23,19,20)]
names(all.info)
colnames(all.info)[16] <- "OakScrub"
colnames(all.info)[18] <- "TSF"

#Change dates back to numeric for survival object

all.info$BrdrDate <- as.numeric(all.info$BrdrDate)
all.info$LastObsDate <- as.numeric(all.info$LastObsDate)

#####################################################################

# Block 4 - Time varying Models 


#Covariates change value each year 