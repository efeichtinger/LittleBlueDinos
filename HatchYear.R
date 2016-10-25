# Hatch Year/First Year Birds 

# Updated 10/14/2016
# Source code from "FirstYr.R" 

# Subsample 1999 - 2015 to look for effects of sex 
# 1981- 2015 full data set but most first year birds in the 80's were not sexed 

library(survival)
library(ggplot2)
library(car)
library(coxme)
library(kinship2)
library(plyr)
library(data.table)
library(stringr)

################################################################################
### Block 1 - data input, 1999-2015 and initial analyses from KM to Coxme models

#April Census Dates, but I don't know how to match them to the birds 1 year later
aprilD <- read.csv("April_Census_Dates.csv")
aprilD$CensusDate <- as.Date(aprilD$CensusDate, format = "%m/%d/%Y")

#Problems yesterday so I just pulled out the birds from 1999 on
#bird.df <- read.csv("Erin_1999_FY.csv")
bird.df <- read.csv("Erin_Oct_FY_mass_sex.csv")

group.size <- read.csv("groupsize.csv")

bird.df$NestYear <- as.factor(bird.df$NestYear)

#Convert dates to date format 
bird.df$FldgDate <- as.Date(bird.df$FldgDate, format = "%m/%d/%Y")
bird.df$LastObsDate <- as.Date(bird.df$LastObsDate, format = "%m/%d/%Y")
bird.df$HatchDate <- as.Date(bird.df$HatchDate, format = "%m/%d/%Y")
bird.df$MeasDate <- as.Date(bird.df$MeasDate, format = "%m/%d/%Y")

bird.df["Days"] <- bird.df$LastObsDate - bird.df$FldgDate
bird.df["Years"] <- round(bird.df$Days/365.25, digits = 2)

str(bird.df)

#Add censorship column, 1 if dead before 1 yr, 0 if yr > 1 or 4/12/12016
bird.df["Censor"] <- 1
bird.df$Censor[which(bird.df$Years >= 1)]<-0
yrlg.df <- subset(bird.df, bird.df$Years > 0 & bird.df$Days > 0)
yrlg.df$Censor[which(yrlg.df$LastObsDate == "2016-04-12")] <- 0

#Jay ID not unique but the band number is 
#yrlg.df[,2] <- NULL

str(yrlg.df)

yrlg.df["Day11"] <- yrlg.df$MeasDate - yrlg.df$HatchDate
yrlg.df <- subset(yrlg.df, yrlg.df$Day11 <= 13)

#Data frame includes the birds from 1981 on but it needs to be filtered since
#sex is included 
yrlg.df <- subset(yrlg.df, yrlg.df$FldgDate > "1999-01-01")

##year to year sex ratio 
#this is some nifty code basic but useful
cat("No. of Females = ", nrow(yrlg.df[yrlg.df$Sex == "F", ]))
cat("No. of Males = ", nrow(yrlg.df[yrlg.df$Sex == "M", ]))

DT <- data.table(yrlg.df)
nums <-DT[, .N, by = list(NestYear, Sex)]

##Probably some fancy way to do this, but I'll just use brute force
fls <- subset(nums, Sex == "F")
mls <- subset(nums, Sex == "M")
all <- cbind(fls,mls)
all[,4] <- NULL
names(all)[3] <- "Females"
names(all)[5] <- "Males"
setkey(all, NestYear)
ratio <- all$Males/all$Females
col1 <- cbind(ratio)
Year <- all$NestYear
#Males to females, male/female
sex.ratios <- data.frame(Year,col1)
colnames(sex.ratios) <- c("Year","SexRatio")
#sex.ratios

#visualize sex ratios of fledglings 
#Add number of fledglings produced, sampling error?
qplot(Year, data=sex.ratios, geom="bar", weight=SexRatio,
  ylab = "Fledgling Ratio of Males to Females", xlab = "Year") +
  theme(axis.text.x=element_text(angle=300,hjust=0, size=12))




#Change to numeric for survival object 
yrlg.df$FldgDate <- as.numeric(yrlg.df$FldgDate)
yrlg.df$LastObsDate <- as.numeric(yrlg.df$LastObsDate)

yrlg.df$Years <- as.numeric(yrlg.df$Years)
yrlg.df$Days <- as.numeric(yrlg.df$Days)

#Create survival object 
yrlg.ob <- Surv(yrlg.df$Years, yrlg.df$Censor, type = c('right'))
my.fit <- survfit(yrlg.ob~1, conf.type = "log-log")

plot.fit <- plot(my.fit, xlab = "Time (years)", conf.int=TRUE,
                 log = "y", ylim = c(0.4, 1),xlim=c(0,1), ylab = "Cumulative Survival", 
                 main = "Fledge to 1Yr - 1999 to 2015")

##Year 
fit.yr <- survfit(yrlg.ob ~ yrlg.df$NestYear)
plot.yr <- plot(fit.yr,xlab = "Time (years)",
                log = "y", ylim = c(0.2, 1),xlim=c(0,1), ylab = "Cumulative Survival", 
                main = "Curves for each year 1999 - 2015")
fit.yr

#Sex -  actually looks alright, shows that the CI's for males and females 
#overlapp
sex.fit <- survfit(yrlg.ob ~ yrlg.df$Sex, conf.type = "log-log")
plot.sex <- plot(sex.fit, conf.int = FALSE,col = c("navy","red"), xlab = "Time (years)", log = "y",
                   ylim = c(0.4,1), xlim=c(0,1), ylab = "Cumulative Survival", 
                   main = "Fledge to 1yr, 1999 - 2015")
legend("topright", c("Females","Males"), col=c("navy","red"),
       lty = c(1,1),lwd=1)

plot.sexci <- plot(sex.fit, conf.int = TRUE,col = c("navy","red"), xlab = "Time (years)", log = "y",
                 ylim = c(0.4,1), xlim=c(0,1), ylab = "Cumulative Survival", 
                 main = "Fledge to 1yr, 1999 - 2015")
legend("topright", c("Females","Males"), col=c("navy","red"),
       lty = c(1,1),lwd=1)


#without color 
plot.sex2 <- plot(sex.fit, conf.int =TRUE, xlab= "Time (yrs)", 
              ylim = c(0.3,1), xlim = c(0,1), 
              ylab = "Cumulative Surivival", 
              main = "Fledge to 1yr, 1999 to 2015")

start <- as.numeric(fit.yr$n)
deaths <- c(47,36,30,54,23,64,56,71,28,96,39,72,116,11,41,65,48)
life.table <- cbind(start, deaths)
life.table <- as.data.frame(life.table)
life.table["p"] <- 1 - (round(life.table$deaths/life.table$start, digits = 3))
life.table["Year"] <- c(1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009, 
                      2010,2011,2012,2013,2014,2015)
str(life.table)

range(life.table$p)
mean(life.table$p)
sd(life.table$p)

#Proportion of survivors at 1 yr for each cohort (year)
qplot(Year, data=life.table, geom="bar", weight=p,
      ylab = "Proportion of survivors at 1 yr", xlab = "Year") +
  theme(axis.text.x=element_text(angle=300,hjust=0, size=12))


#Now some models 
cox.sex <- coxph(yrlg.ob ~ Sex, data = yrlg.df)
cox.sex
anova(cox.sex)

cox.mass <- coxph(yrlg.ob ~ Weight, data = yrlg.df)
cox.mass
anova(cox.mass)

mass.sex <- coxph(yrlg.ob ~ Sex + Weight, data = yrlg.df)
mass.sex
anova(mass.sex)

ms.sexint <- coxph(yrlg.ob ~ Sex + Weight + Sex:Weight, data = yrlg.df)
ms.sexint
anova(ms.sexint)
#lol this adds nothing to the model 
#I think it's safe to say sex does not substanitally change survival between sexes
#At least not with the sample size we have 

extractAIC(cox.sex)
extractAIC(cox.mass)
extractAIC(mass.sex)
extractAIC(ms.sexint)
#based on this sex really doesn't add anything 

flg.num <- coxph(yrlg.ob ~ FldgNum, data=yrlg.df)
flg.num

hatch.num <- coxph(yrlg.ob ~ HatchNum, data=yrlg.df)
hatch.num

fld.hatch <- coxph(yrlg.ob ~ HatchNum + FldgNum, data=yrlg.df)
fld.hatch

fhmass <- coxph(yrlg.ob ~ HatchNum + FldgNum + Weight, data=yrlg.df)
fhmass

#Mass and flg num interaction
mass.fl <- coxph(yrlg.ob ~ FldgNum + Weight + FldgNum:Weight, data = yrlg.df)
mass.fl
anova(mass.fl)

extractAIC(mass.fl)
#AIC very close to that of the model with day 11 mass alone 

##something weird is happening because this worked yesterday and now 
##there seems to be a problem with perfect classification 
cox.year <- coxph(yrlg.ob ~ NestYear, data = yrlg.df)
cox.year
anova(cox.year)
extractAIC(cox.year)

cox3 <- coxph(yrlg.ob ~ NestYear + Sex, data = yrlg.df)
anova(cox3)
extractAIC(cox3)

#### Mixed effects models
#NestYear/Cohort Year 
mm.year <- coxme(yrlg.ob ~ Weight + (1|NestYear), data = yrlg.df)
mm.year

mm.nest <- coxme(yrlg.ob ~ Weight + (1|NatalNest), data = yrlg.df)
mm.nest
anova(mm.year)
anova(mm.nest)


#########################################################################
#Block 2 - territory quality data 

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

### Write to csv file 
#write.csv(terr.info, file ="TerrInfo.csv", sep =" ", col.names=TRUE)


#String function and regular expressions to get a "TerrYr" variable to match
#Remember this is just for 1999 to 2015, the territory data goes to 1981
new.col <- gsub(".$","",yrlg.df$NatalNest)
colTY <- as.vector(new.col)
colTY <- cbind(colTY)

yrlg.df["TerrYr"] <- colTY
yrlg.df <- merge(yrlg.df, terr.info, by="TerrYr")

yrlg.df <- merge(yrlg.df, group.size, by="TerrYr")
colnames(yrlg.df)[22] <- "GroupSize"
colnames(yrlg.df)[19] <- "oakscrub"
colnames(yrlg.df)[21] <- "tsf"
### Make a new survival object? 
new.ob <- Surv(yrlg.df$Years, yrlg.df$Censor, type = c('right'))

names(yrlg.df)

# Add columns for centered variables for territory metrics
#Studentized (I think) - standard deviation 
yrlg.df["stdscr"] <- scale(yrlg.df$oakscrub, center = FALSE, scale = TRUE)
#Centered - mean
yrlg.df["centscr"] <- scale(yrlg.df$oakscrub, center = TRUE, scale =FALSE)

yrlg.df["stdtsf"] <- scale(yrlg.df$tsf, center = FALSE, scale = TRUE)
yrlg.df["centtsf"] <- scale(yrlg.df$ts, center = TRUE, scale = FALSE)

yrlg.df["stdsize"] <- scale(yrlg.df$TerrSize, center= FALSE, scale = TRUE)
yrlg.df["centsize"] <- scale(yrlg.df$TerrSize, center = TRUE, scale = FALSE)



## Need to think about the modeling procedure using terr metrics 
#Model with studentized data for oak scrub cell count and tsf 1-9 year cell count
#Units of standard deviation from the mean (I think)

#Oak scrub
#Untransformed 
oak <- coxph(new.ob ~ oakscrub, data = yrlg.df)
#Studentized
oak.cx <- coxph(new.ob ~ stdscr, data = yrlg.df)
oak.cx
anova(oak.cx)
#Centered
oak.cx2 <- coxph(new.ob ~ centscr, data = yrlg.df)
oak.cx2

#TSF
#Untransformed
tsf.cx1 <- coxph(new.ob ~ tsf, data = yrlg.df)
tsf.cx1
#Studentized
tsf.cx2 <- coxph(new.ob ~ stdtsf, data = yrlg.df)
tsf.cx2
anova(tsf.cx2)
#Centered
tsf.cx3 <- coxph(new.ob ~ centtsf, data = yrlg.df)
tsf.cx3

#terr size only why noy
terrsize.cx <- coxph(new.ob ~ stdsize, data = yrlg.df)
terrsize.cx
anova(terrsize.cx)

#Oak and TSF
oak.tsf <- coxph(new.ob ~ stdscr + stdtsf, data= yrlg.df)
oak.tsf
anova(oak.tsf)

oak.tsf.int <- coxph(new.ob ~ stdscr + stdtsf + stdscr:stdtsf, data = yrlg.df)
oak.tsf.int
anova(oak.tsf.int)

tsf.size <- coxph(new.ob ~ stdtsf + stdsize, data = yrlg.df)
tsf.size
anova(tsf.size)

oak.tsf.size <- coxph(new.ob ~ stdscr + stdtsf + stdsize, data = yrlg.df)
oak.tsf.size
anova(oak.tsf.size)

#What about untransformed data for oak and terr size, interaction?
oak.size <- coxph(new.ob ~ oakscrub + TerrSize + 
              scrb.count:TerrSize, data = yrlg.df)
oak.size
anova(oak.size)

oak.size2 <- coxph(new.ob ~ stdscr + stdsize, data = yrlg.df)
oak.size2
oak.size3 <- coxph(new.ob ~ stdscr + stdsize + stdscr:stdsize, data = yrlg.df)
oak.size3
anova(oak.size3)

oak.size4 <- coxph(new.ob ~ stdscr + stdsize +stdtsf +stdscr:stdsize
                   + stdtsf:stdsize, data= yrlg.df)
oak.size4
anova(oak.size4)


##### ADD IN GROUP SIZE!!!!! Then look for interactions with that and territory
### After that, start combining predictors that have turned out to be significant

# Natal Nest ID, day 11 mass, oak scrub, pedigree, group size, year

##Group size 

grp.sz <- coxph(new.ob ~ GroupSize, data=yrlg.df)
grp.sz
anova(grp.sz)

grp.oak <- coxph(new.ob ~ GroupSize + stdscr, data = yrlg.df)
grp.oak
anova(grp.oak)

grp.oakx <- coxph(new.ob ~ GroupSize + stdscr + GroupSize:stdscr, data=yrlg.df)
grp.oakx
anova(grp.oakx)

#Oak scrub, terr size, group size 

otg <- coxph(new.ob ~ GroupSize + stdscr + stdsize, data = yrlg.df)
otg
anova(otg)

#I think this might be the best one 
otgx <- coxph(new.ob ~ GroupSize + stdscr + stdsize + stdscr:stdsize, data=yrlg.df)
otgx
anova(otgx)

otgfx <- coxph(new.ob ~ GroupSize + stdscr + stdsize + stdtsf +
                 stdscr:stdsize + stdtsf:stdsize, data= yrlg.df)
otgfx
anova(otgfx)
