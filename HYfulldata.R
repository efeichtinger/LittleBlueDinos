## Jan 18 2017

# Hatch year birds with full data set 1981-2015 
# HatchYear.R is code and analyses for the subsample of 1999-2015 when birds
# were sexed using blood samples 

#my_env = new.env(hash = TRUE, parent = .GlobalEnv)
#saveRDS(my_env, file = "HYApril19.RData")

library(survival)
library(ggplot2)
library(car)
library(coxme)
library(kinship2)
library(plyr)
library(data.table)
library(stringr)
library(simPH)
library(AICcmodavg)

##########################################################################
#Block 1 - upload jay data and tidying up 

#April Census Dates, but I don't know how to match them to the birds 1 year later
aprilD <- read.csv("April_Census_Dates.csv")
aprilD$CensusDate <- as.Date(aprilD$CensusDate, format = "%m/%d/%Y")

#Data file with HY bird data, 1981-2016
hybrd <- read.csv("Erin_Jan_FY.csv")
#Data file with helpers by TerrYr
hlpty <- read.csv("Helpers_TerrYr.csv")

group.size <- read.csv("groupsize.csv")

hybrd$NestYear <- as.factor(hybrd$NestYear)

#Convert dates to date format 
hybrd$FldgDate <- as.Date(hybrd$FldgDate, format = "%m/%d/%Y")
hybrd$LastObsDate <- as.Date(hybrd$LastObsDate, format = "%m/%d/%Y")
hybrd$HatchDate <- as.Date(hybrd$HatchDate, format = "%m/%d/%Y")
hybrd$MeasDate <- as.Date(hybrd$MeasDate, format = "%m/%d/%Y")

hybrd["Days"] <- hybrd$LastObsDate - hybrd$FldgDate
hybrd["Years"] <- round(hybrd$Days/365.25, digits = 2)

names(hybrd)

#Create data frame to see repeat band numbers 
reprows <- ddply(hybrd, .(hybrd$JayID), nrow)

colnames(hybrd)[7] <- "Mass"

#Add censorship column, 1 if dead before 1 yr, 0 if yr > 1 or 4/12/12016
hybrd["Censor"] <- 1

#this is the line where I could change to April census + 1
hybrd$Censor[which(hybrd$Days >= 365.25)]<-0

yrlg.df <- subset(hybrd, hybrd$Days > 0)
yrlg.df$Censor[which(yrlg.df$LastObsDate == "2016-04-12")] <- 0

#Jay ID not unique but the band number is 
#yrlg.df[,2] <- NULL

yrlg.df["Expr1003"] <- NULL
yrlg.df["Day11"] <- yrlg.df$MeasDate - yrlg.df$HatchDate


yrlg.df <- subset(yrlg.df, yrlg.df$Day11 > 10 & yrlg.df$Day11 <= 20)

str(yrlg.df)

yrlg.df["SocialClass"] <- NULL


#rearrange columns 
yrlg.df <- yrlg.df[,c(1,2,3,11,12,8,4,5,6,7,9,10,13,14,15,16)]

#create TerrYr
yrlg.df["TerrYr"] <- paste(yrlg.df$Terr, str_sub(yrlg.df$NestYear, start= -2),
                       sep ="")
#Change NestYear to Year to match other data files 
colnames(yrlg.df)[5] <- "Year"
str(yrlg.df)

yrlg.df1 <- merge(yrlg.df, hlpty)

str(yrlg.df1)

length(unique(yrlg.df1$USFWBand))
#3145 "levels" of band numbers, does this sound like a reasonable sample size?
#JayID not unique at nestling stage so there are 2261 levels of JayID
#Unique color band at day 75 assuming a bird makes it there

yrlg.df1 <- yrlg.df1[!duplicated(yrlg.df1$USFWBand),]
length(unique(yrlg.df1$USFWBand))

#filter out 2016 birds 
yrlg.df1 <- subset(yrlg.df1, yrlg.df1$FldgDate <= "2015-12-31")
length(unique(yrlg.df1$USFWBand))

yrlg.df1 <- na.omit(yrlg.df1)

#### April 8 2017 
# Just now discovered a slight problem, birds were dropped due to the measurement filtering
#N = 2492 now

###################################################################
####See if you can get census date in
col.year <- cbind(as.character(yrlg.df1$FldgDate))

## Year Month Date
colnames(col.year) <- "FldgDate"
col.year <- gsub("-","", col.year)
col.april <- cbind(as.character(aprilD$CensusDate))
col.april <- gsub("-","",col.april)
colnames(col.april) <- "AprilDate"

#Fledge years
new.col <- substr(col.year, 1, 4)
#April census years
new.april <- substr(col.april, 1, 4)

#Match FldgDate form new.col (year) to first 4 characters in col.april
d$match <- map$c2[match(d$c1,map$c1)]

new <- new.col[match(new.april, new.col)]
new1 <- new.col[match(new.col, new.april)]
new1 <- cbind(new.col[match(new.col, new.april)])

test <- cbind(new.col, new1)

##################################################################
yrlg.df1$FldgDate <- as.numeric(yrlg.df1$FldgDate)
yrlg.df1$LastObsDate <- as.numeric(yrlg.df1$LastObsDate)

yrlg.df1$Years <- as.numeric(yrlg.df1$Years)
yrlg.df1$Days <- as.numeric(yrlg.df1$Days)


str(yrlg.df1)

# 3 3 2017 leave these variables as integers for now, can change later
# I *think* helper number is a count so it should be an integer and not
# a category
# need this data type for the simPH program to work 
#yrlg.df1$Help <- as.factor(yrlg.df1$Help)
#yrlg.df1$HelpNum <- as.factor(yrlg.df1$HelpNum)

length(unique(yrlg.df1$USFWBand))

## 6 12 2017 add density 
#2492 now, updated April 8 2017

# 6 15 2017 new df with density, don't overwirte yrlg.df1 bc 
#it messes things up downstream 
den.df <- merge(yrlg.df1, jay.den, by = ("Year"), all.x=TRUE)



#basic survival object and KM plots
hy.ob <- Surv(yrlg.df1$Days, yrlg.df1$Censor, type = c('right'))
hy.ob1 <- Surv(den.df$Days, den.df$Censor, type =c('right'))
hy.fit <- survfit(hy.ob ~ 1, conf.type = "log-log")
hy.fit1 <- survfit(hy.ob1 ~ 1, conf.type = "log-log")
summary(hy.fit)
hy.fit
str(summary(hy.fit))
plot(hy.fit$time, hy.fit$surv, xlim = c(0,365), ylim = c(0,1), type = "l")
lines(hy.fit$time, hy.fit$upper, lty = 2)
lines(hy.fit$time, hy.fit$lower, lty = 2)



yr.fit <- survfit(hy.ob ~ yrlg.df1$Year, conf.type="log-log")
yr.fit1 <- survfit(hy.ob1 ~ den.df$Year, conf.type="log-log")
yr.fit
summary(yr.fit)
str(summary(yr.fit))
#plot(yr.fit, log="y", xlim = c(0,1))

new11 <- coxph(hy.ob ~ Mass + Help, data = yrlg.df1)
new1 <- coxph(hy.ob ~ Mass + Help + Year, data = yrlg.df1)
new2 <- coxme(hy.ob ~ Mass + Help + (1|Year), data =yrlg.df1)
new3 <- coxme(hy.ob ~ Mass + Help + (1|USFWBand) + (1|Year),
              varlist=coxmeMlist(kins), data = yrlg.df1)

#6 4 2017
#Changed df to den.df, same as yrlg.df1 but with density 
#messes things up downstream if I overwrite yrlg.df1
yr.full <- coxph(hy.ob1 ~ Year, data = den.df)
ym <- coxph(hy.ob1 ~ Mass + Year, data = den.df)
summary(ym)
mbb <- coxph(hy.ob1 ~ Mass + HatchNum, data = den.df)
mbf <- coxph(hy.ob1 ~ Mass + FldgNum, data = den.df)
mb <- coxph(hy.ob1 ~ Mass + HatchNum + FldgNum, data= den.df)
extractAIC(mb)
ymh <- coxph(hy.ob1 ~ Mass + Help + Year, data = den.df)
ymhi <- coxph(hy.ob1 ~ Mass +  Year*Help, data = den.df)
ymhii <- coxph(hy.ob1 ~ Mass + Help + Year + Year*Help, data = den.df)
summary(ymhii)

# 6 12 2017 density 
dens <- coxph(hy.ob1 ~ Mass, data = den.df)
dens1 <- coxph(hy.ob1 ~ Mass + Den, data = den.df)
dens2 <- coxph(hy.ob1 ~ Mass + Den + Help , data = den.df)
dens3 <- coxph(hy.ob1 ~ Mass + Den + Help + Year, data= den.df)
extractAIC(dens)
extractAIC(dens1)
extractAIC(dens2)
extractAIC(dens3)
summary(dens1)

myb <- coxph(hy.ob1 ~ Mass +Help + HatchNum + Year, data = den.df)
mybd <- coxph(hy.ob1 ~ Mass + Den + Help + HatchNum + Year, data = den.df)
myf <- coxph(hy.ob1 ~ Mass  + Help + FldgNum + Year, data = den.df)

Cand.modelsE <- list()
Cand.modelsE[[1]] <- dens1
Cand.modelsE[[2]] <- dens2
Cand.modelsE[[3]] <- dens3
Cand.modelsE[[4]] <- ym
Cand.modelsE[[5]] <- ymh
Cand.modelsE[[6]] <- ymhii
Cand.modelsE[[7]] <- myb
Cand.modelsE[[8]] <- mybd
ModnamesE <- paste("Model", 1:length(Cand.modelsE), sep="")
tableE <- aictab(Cand.modelsE, ModnamesE, sort = TRUE, second.ord=FALSE)
tableE

extractAIC(myb)
extractAIC(myf)
extractAIC(ymh)
extractAIC(ym)
Cand.modelsB <- list()
Cand.modelsB[[1]] <- mass.only
Cand.modelsB[[2]] <- ym
Cand.modelsB[[3]] <- ymh
Cand.modelsB[[4]] <- myb
Cand.modelsB[[5]] <- myf
Cand.modelsB[[6]] <- desn2
ModnamesB <- paste("Model", 1:length(Cand.modelsB), sep="")
tableB <- aictab(Cand.modelsB, ModnamesB, sort = TRUE, second.ord=FALSE)
tableB


# 6 7 2017 Just to check estimates w/o year
mass.only <- coxph(hy.ob ~ Mass, data = yrlg.df1)
summary(mass.only)
help.only <- coxph(hy.ob ~ Help, data = yrlg.df1)
summary(help.only)
helpn <- coxph(hy.ob ~ HelpNum, data = yrlg.df1)
summary(helpn)
mh <- coxph(hy.ob ~ Mass + Help, data = yrlg.df1)
summary(mh)
mhn <- coxph(hy.ob ~ Mass + HelpNum, data = yrlg.df1)
summary(mhn)
extractAIC(mh)
extractAIC(mhn)
#AIC same for either form of helper 

extractAIC(mass.only)
extractAIC(mh)
extractAIC(mbb)
extractAIC(mb)
extractAIC(mbf)

## June 21 - figures of coefficients
#use den.df
#hy.ob1

## I didn't scale mass and helper number

fig1 <- coxph(hy.ob1 ~ Mass + HelpNum, data = den.df)

fig33.df <- data.frame(predict(fig1, newdata=transform(den.df,Mass=mean(Mass)),
              type="risk", se.fit=TRUE))
fig1.df <- data.frame(predict(fig1, newdata=den.df),
                              type="risk", se.fit=TRUE)
fig1.df$HelpNum <- den.df$HelpNum
#fig1.df$Mass <- den.df$Mass
pf1 <- with(fig1.df, data.frame(HelpNum, 
                  hazard = fit, lwr = fit- 1.96*se.fit,
                  upr = fit+ 1.96*se.fit))
fig1A <- ggplot(pf1, aes(HelpNum, hazard)) +
  geom_point() +
  xlab("Number of Helpers") +
  ylab("Predicted RelatiVe Risk") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 
  #geom_ribbon(data = pf1, aes(ymin=lwr, ymax=upr), alpha = 0.2) 
fig1A

tst <- data.frame(predict(fig1, type="risk", se.fit=TRUE))


trts1 <- expand.grid(HelpNum = den.df$HelpNum, 
                      Mass = mean(den.df$Mass))
dat.temp <- trts1
dat.temp$Risk <- predict(fig1, type='risk')



#density
fig2 <- coxph(hy.ob1 ~ Den, data = den.df)
summary(fig2)
fig2.df <- data.frame(predict(fig2, type="risk", se.fit=TRUE))
fig2.df$Den <- den.df$Den
pf2 <- with(fig2.df, data.frame(Den, hazard = fit, lwr = fit - 1.96*se.fit,
                    upr = fit + 1.96*se.fit))

fig2A <- ggplot(pf2, aes(Den, hazard)) + 
  geom_line() +
  xlab("Population Density") +
  ylab("Predicted Relative Risk") + 
  theme_bw() +
  geom_hline(yintercept = 1.0, linetype = 2) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  geom_ribbon(data = pf2, aes(ymin=lwr, ymax=upr), alpha = 0.2) 
  
fig2A


Cand.modelsA <- list()
Cand.modelsA[[1]] <- mass.only
Cand.modelsA[[2]] <- mh
Cand.modelsA[[3]] <- mbb
Cand.modelsA[[4]] <- mb
Cand.modelsA[[5]] <- mbf
ModnamesA <- paste("Model", 1:length(Cand.modelsA), sep="")
tableA <- aictab(Cand.modelsA, ModnamesA, sort = TRUE, second.ord=FALSE)
tableA

Cand.modelsC <- list()
Cand.modelsC[[1]] <- mass.only
Cand.modelsC[[2]] <- mh
Cand.modelsC[[3]] <- mbb
Cand.modelsC[[4]] <- mb
Cand.modelsC[[5]] <- mbf
Cand.modelsC[[6]] <- ym
Cand.modelsC[[7]] <- ymh
Cand.modelsC[[8]] <- myb
Cand.modelsC[[9]] <- myf
Cand.modelsC[[10]] <- ymhii
ModnamesC <- paste("Model", 1:length(Cand.modelsC), sep="")
tableC <- aictab(Cand.modelsC, ModnamesC, sort = TRUE, second.ord=FALSE)
tableC


#Cand.modsA <- list("Mass" = mass.only, "Mass + Helper" = mh, 
     # "Mass + HN" = mbb, "Mass + HN + FN" = mb, "Mass + FN")
#tableA <- aictab(Cand.modsA, sort = TRUE, second.ord=FALSE)

summary(yr.full)
summary(ym)
summary(ymh)
summary(ymhi)
#summary(ymhii)

extractAIC(yr.full)
extractAIC(ym)
extractAIC(ymh)
extractAIC(ymhi)
extractAIC(ymhii)

library(AICcmodavg)
Cand.models2 <- list()
Cand.models2[[1]] <- yr.full
Cand.models2[[2]] <- ym
Cand.models2[[3]] <- ymh
Cand.models2[[4]] <- ymhi
Modnames2 <- paste("Model", 1:length(Cand.models2), sep="")
#with second.ord=FALSE, will use AIC (not AICc)
hymods2 <- aictab(Cand.models2, Modnames2, sort = TRUE, second.ord=FALSE)
hymods2


### Extract number of fledglings and number of events/deaths 
no.flg <- yr.fit$n
no.events <- c(29,24,35,54,20,26,30,23,39,9,51,13,35,
                   55,34,41,70,22,45,33,32,42,28,64,60,
                   73,30,91,38,71,106,15,49,72,50)
start <- as.Date("1981-01-01")
end <- as.Date("2015-12-31")
years <- seq(start, end, "years")
years1 <- substring(years, 1,4)
kind.oflt <- data.frame(years1,no.flg, no.events)
colnames(kind.oflt) <- c("Year","NumFldg", "NumDeaths")
kind.oflt$d <- round(kind.oflt$NumDeaths/kind.oflt$NumFldg, digits = 2)
kind.oflt$p <- 1 - kind.oflt$d
kind.oflt$propend <- kind.oflt$p*100


#### need # fledgs per breeding pair
num.pairs <- read.csv("PairsByYearCount.csv")
prs1981 <- subset(num.pairs, num.pairs$NestYear >= 1981)
prs1981 <- subset(prs1981, prs1981$NestYear <= 2015)
tmp.pairs <- cbind(kind.oflt, prs1981)
tmp.pairs$MPP <- tmp.pairs$NumFldg/tmp.pairs$CountOfNestYear



###Number of fledglings produced, looks decent 
#Suppl. info
inputs <- ggplot(tmp.pairs, aes(Year, MPP)) +
  geom_bar(stat='identity', color = 'gray', fill='black', width = 0.9)
inputs  +
  scale_x_discrete(breaks = c("1981","1985","1990","1995","2000","2005",
                              "2010","2015")) +
  ylab("Mean Number of Fledglings") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  scale_y_continuous(expand=c(0,0), limits = c(0,2)) +
  theme(axis.text.x = element_text(angle=30)) +
  #geom_hline(yintercept = (median(kind.oflt$NumFldg))) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

##Proportion alive each year 
p.year <- ggplot(kind.oflt, aes(Year, p)) +
  geom_bar(stat='identity', color = 'gray', fill='black', width = 0.9) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  scale_y_continuous(expand=c(0,0), limits = c(0,1)) + 
  scale_x_discrete(breaks =c("1981","1985","1990","1995","2000","2005",
                   "2010","2015")) +
  theme(axis.text.x = element_text(angle=30)) +
  ylab("Proportion of Cohort") +
  geom_hline(yintercept = 0.4) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

p.year


library(GGally)
#KM curve for all years 
all.yrs <- ggsurv(yr.fit, plot.cens=FALSE, 
                  back.white=TRUE) +
xlim(0,365) +
  ylim(0,1) +
  theme(legend.position="none") +
  ylab("Cumulative Survival") +
  xlab("Time (Days)") + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 

all.yrs

all.yrs.bw <- ggsurv(yr.fit, plot.cens=FALSE, 
                     back.white=TRUE, surv.col="black") +
  xlim(0,365) +
  ylim(0,1) +
  theme(legend.position="none") +
  ylab("Cumulative Survival") +
  xlab("Time (Days)") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 
all.yrs.bw

days.plot2 <- ggsurv(hy.fit, plot.cens = FALSE, lty.est = 1, lty.ci =2, 
                     back.white=TRUE, xlab = "Time (Days)",
                     ylab = "Cumulative Survival") +
  ylim(0,1)  +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  scale_x_continuous(limits=c(0,365),
                     breaks=c(0,50,100,150,200,250,300,350)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
days.plot2


m.surv <- mean(kind.oflt$p)
sd.surv <- sd(kind.oflt$p)
m.surv
sd.surv
range(kind.oflt$p)
median(kind.oflt$p)

plot(yr.fit$time, yr.fit$sur)

mean(yr.fit$n)
sd(yr.fit$n)
range(yr.fit$n)
median(yr.fit$n)
plot(yr.fit$n)

grid.arrange(days.plot2, p.year, ncol=1)
multiplot(days.plot2, p.year, cols=2)
mulitplot(days.plot2, p.year, cols=1)

gA <- ggplotGrob(days.plot2)
gB <- ggplotGrob(p.year)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
grid.arrange(gA, gB, ncol=2)
#####################################################################
#4 8 2017 
#code below is test

#4 4 2017
##Plot options for publication using base plot
par(mfrow = c(1,2), mar = c(4,4,0.5,0.5), oma = c(0,0,0,0))
kmcurve <- plot(hy.fit, xlab="Time (Days)", 
               ylab = "Proportion Surviving", 
               xlim = c(0,1), ylim = c(0,1))

yr.curves <- plot(yr.fit, xlim = c(0,365), ylim = c(0,1),
     xlab = "Time (Days)")

### I'd like to use ggplot, but if I can't get it figured out 
#just have to use base plot (above)
library(GGally)
surv.plot1 <- ggsurv(hy.fit, plot.cens = FALSE, xlab = "Time (Days)",
                     ylab = "Proportion Surviving")  
  xlim(0,1) + 
  surv.plot1 + theme_bw() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 
surv.plot1
  
try <- ggsurv(hy.fit, plot.cens = FALSE, lty.est = 1, lty.ci = 2,
              back.white=TRUE, xlab = "Time (Years)", 
              ylab = "Cumulative Survival")  
try + ylim(0,1) + xlim(0,1)

days.ob <- Surv(yrlg.df1$Days, yrlg.df1$Censor, type = c('right'))
days.fit <- survfit(days.ob ~ 1, conf.type = "log-log")

days.plot <- ggsurv(days.fit, plot.cens = FALSE, lty.est = 1, lty.ci =2, 
                    back.white=TRUE, xlab = "Time (Days)",
                    ylab = "Cumulative Survival")
days.plot + ylim(0,1) + xlim(0,365)


#### Go from here, try to work out the axis limits issue 
days.plot2 <- ggsurv(days.fit, plot.cens = FALSE, lty.est = 1, lty.ci =2, 
                     back.white=TRUE, xlab = "Time (Days)",
                     ylab = "Cumulative Survival")
days.plot2  +ylim(0,1)  +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  scale_x_continuous(limits=c(0,365),
                     breaks=c(0,50,100,150,200,250,300,350))
  
#days yr.fit
days.yr <- survfit(days.ob ~ yrlg.df1$Year)
days.plot3 <- ggsurv(days.yr, plot.cens = FALSE,
                     back.white=TRUE, xlab= "Time (Days)",
                     ylab = "Cumulative Survival")
days.plot3 + ylim(0,1) + scale_x_continuous(limits=c(0,365),
              breaks=c(0,50,100,150,200,250,300,350)) +
  theme(legend.position="none")


## All years
ggsurv(yr.fit, plot.cens =FALSE, xlab = "Time (Years)", 
       ylab="Proportion Surviving") +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  theme(legend.position="none")


#store summary of survfit object 
stuff1 <- summary(hy.fit)
str(stuff1)
cols.st <- lapply(c(2,6,8:11), function(x) stuff1[x])
store.it <- do.call(data.frame, cols.st)
#select columns by index from str(res)
st.plot <- ggplot(store.it, aes(time, surv)) +
  geom_line(linetype=1, size=1, color = "black") +
  geom_line(aes(time, lower), linetype=2, color= "black", size = 1) +
  geom_line(aes(time, upper), linetype=2, color = "black", size=1) +
  theme(axis.text = element_text(size=10)) +
  theme(axis.text.x = element_text(size=10, face="bold"))+
  theme(axis.text.y =element_text(size=10, face="bold")) +
  theme(axis.title.x = element_text(face = "bold")) +
  theme(axis.title.y = element_text(face = "bold")) +
  xlab("Time (Years)") +
  ylab("Proportional Survival") +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  xlim(0,1) +
  ylim(0.38,1)
st.plot



stuff2 <- summary(yr.fit)
str(stuff2)
cols.st2 <- lapply(c(2,6,9,10,11,8), function(x) stuff2[x])
store.it2 <- do.call(data.frame, cols.st2)
eightone <- subset(store.it2, strata = "yrlg.df1$Year=1981")
eighttwo <- subset(store.it2, strata = "yrlg.df1$Year=1982")
eightthree <- subset(store.it2, strata = "yrlg.df1$Year=1983")

ggplot(store.it2, aes(time, surv)) +
  geom_line(eightone, aes(time, surv)) +
  geom_line(eighttwo, aes(time, surv)) +
  geom_line(eightthree, aes(time, surv))

library(survminer)
ggsurvplot(yr.fit, data = yrlg.df1, risk.table = FALSE)

#######################################################################

#Correlation plots
library(corrplot)
covars <- yrlg.df1[,c(10,12,13,18,19)]
covars$Help <- as.numeric(covars$Help)
covars$HelpNum <- as.numeric(covars$HelpNum)
corrs <- cor(covars, use="complete.obs")
corrs
corrplot(corrs, method="pie", type= "lower")

####################################################################
#EDA plots

ggplot(yrlg.df1, aes(Mass)) +
  geom_histogram(binwidth = 0.3, colour = "darkblue", fill = "darkblue") +
  xlab("Day 11 Mass (g)") +
  ylab("Frequency")

mass.plot <- ggplot(yrlg.df1, aes(Mass)) +
  geom_histogram(binwidth = 0.3, colour = "black", fill = "black") +
  xlab("Day 11 Mass (g)") +
  ylab("Frequency") + 
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 
  #geom_vline(xintercept = mean(yrlg.df1$Mass), color = "red")
mass.plot
  

cat("Mean nestling mass= ", mean(yrlg.df1$Mass))
cat("Standard deviation of mass", sd(yrlg.df1$Mass))
cat("Median of mass", median(yrlg.df1$Mass))

survs <- subset(yrlg.df1, yrlg.df1$Censor == 0)
notsurvs <- subset(yrlg.df1, yrlg.df1$Censor == 1)

mean(survs$Mass)
sd(survs$Mass)
mean(notsurvs$Mass)
sd(notsurvs$Mass)

#Fledglings 
names(HYsurv)
prop.plot <- ggplot(HYsurv, aes(Year, N)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45)) +
  xlab("Year") +
  ylab("Number of Fledglings produced") +
  scale_x_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank())
  
prop.plot

fldgs.plot <- ggplot(HYsurv, aes(Year, N)) +
  geom_bar(width = 0.75, stat = "identity", color = "black", fill = "red") +
  theme(axis.text.x = element_text(angle = 45)) + 
  xlab("Year") +
  ylab("Count of Individuals") +
  geom_bar(aes(Year, D), stat = "identity", width = 0.75,
           color= "black", fill = "blue")
fldgs.plot
#can I overlay the number of survivors?? Basically adding bars on top
#showing the # of survivors (not proportion)


# 3 6 2017 yanked from the R cookbook to combine plots 
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


multiplot(surv.plot1, fldgs.plot, cols=2)
########################################################################
#Pedigree and kinship matrix here
names(yrlg.df1)
str(yrlg.df1)

#upload pedigree
jayped <- read.csv("Demo_Pedigree_2017.csv")
names(jayped)
jayped$SexIND[which(jayped$Sex=="F")] <- 2
jayped$SexIND[which(jayped$Sex=="M")] <- 1
jayped$SexIND[which(jayped$Sex=="")] <- 3
jayped[,11:12] <- NULL


HYped <- merge(yrlg.df1, jayped, intersect = c("USFWBand", "NatalNest"))
names(HYped)
str(HYped)

wtf <- subset(yrlg.df1, !(yrlg.df1$USFWBand %in% HYped$USFWBand))

HYped <- HYped[,c(2,1,6,21,23,24,25,26,16)]


names(HYped)
names(pars)

#Use this if you have added the territory data
#HYped <- HYped[,c(1,2,6,30,32,33,34,35,16)]
str(HYped)
HYped$Year <- as.character(HYped$Year)
HYped$Year <- as.numeric(HYped$Year)


#So it appears that everyone in the ID list must appear in parent list
#Birds who have parents who are not in the list get NA for parents 

length(unique(HYped$FUSFWBand))  #398 moms, there is one with a missing band#
length(unique(HYped$FBreeder))      
length(unique(HYped$MUSFWBand))   #383 ids, 2 missing band #
length(unique(HYped$MBreeder))


parents <- data.frame(HYped$FBreeder, HYped$FUSFWBand, HYped$MBreeder, HYped$MUSFWBand)
colnames(parents) <- c("FID", "Fband", "MID", "Mband")
library(plyr)
mom.counts <- count(parents, "FID")
dad.counts <- count(parents, "MID")

short.ped <- data.frame(HYped$JayID, HYped$FBreeder, HYped$FUSFWBand, 
                        HYped$MBreeder, HYped$MUSFWBand)
colnames(short.ped) <-c("JayID", "Mom", "MomUS", "Dad", "DadUS")
moms <- data.frame(short.ped$Mom, short.ped$MomUS)
colnames(moms) <- c("Mom", "MomUS")
dads <- data.frame(short.ped$Dad, short.ped$DadUS)
colnames(dads) <- c("Dad", "DadUS")
kids <- data.frame(short.ped$JayID)
colnames(kids) <- "kids"

#moms not in the focal bird list (fledglings)
notmoms <- subset(short.ped, !(short.ped$Mom %in% short.ped$JayID), 
                  select=JayID:MomUS)
length(unique(notmoms$Mom)) # 271 mamas not in the fledge list
notmoms[,1] <- NULL
notmoms <- data.frame((notmoms[!duplicated(notmoms$Mom),]))
#Moms not in focal ID list (not found in HY bird individuals)
#Color band combo
colnames(notmoms) <- c("Mom", "MomUS")

notdads <- subset(short.ped, !(short.ped$Dad %in% short.ped$JayID))
notdads[,2] <- NULL
notdads[,1] <- NULL
notdads[,1] <- NULL
notdads <- data.frame((notdads[!duplicated(notdads$Dad),]))
#Dads not in focal ID list (not found in HY bird individuals)
colnames(notdads) <- c("Dad", "DadUS")

#Any birds with parent(s) in this list get 0 for parent ID
notdads
notmoms

length(unique(notdads$Dad)) # 195 dads not in the fledge list 
count.notf <- count(notmoms, "Mom")
count.notm <- count(notdads, "Dad")

yesmoms <- subset(short.ped, (short.ped$Mom %in% short.ped$JayID))
yesmoms <- yesmoms[!duplicated(yesmoms$Mom),]
yesmoms[1] <- NULL
yesmoms[3] <- NULL
yesmoms[3] <- NULL
yesdads <- subset(short.ped, (short.ped$Dad %in% short.ped$JayID))
yesdads <- yesdads[!duplicated(yesdads$Dad),]
yesdads[1] <- NULL
yesdads$Mom <- NULL
yesdads$MomUS <- NULL
#Founders need to be any bird that has a mom or dad in the data frame above

#Mom and Dad ID's need to be changed to NA for any bird that has parents
#on list from "count.notf" and "count.notm" - those birds are not in fledge list
#There are some birds that were fledglings then became parents later 

#2 14 2017 
#Add in birds that are not on focal id list but have kids on it
names(HYped)
notmoms["SexIND"] <- 2
notdads["SexIND"] <- 1
colnames(notmoms) <- c("JayID", "USFWBand", "SexIND")
colnames(notdads) <- c("JayID", "USFWBand", "SexIND")

#founders
pars <- rbind(notmoms, notdads)
pars["Year"] <- ""
pars["FBreeder"] <- NA
pars["FUSFWBand"] <- NA
pars["MBreeder"] <- NA
pars["MUSFWBand"] <- NA
pars["Censor"] <- ""

names(HYped)
names(pars)

pars <- pars[,c(1,2,4,3,5,6,7,8,9)]

newPed <- rbind(HYped, pars)
newPed$FBreeder <- as.character(newPed$FBreeder)
newPed$FUSFWBand <- as.character(newPed$FUSFWBand)
newPed$MBreeder <- as.character(newPed$MBreeder)
newPed$MUSFWBand <- as.character(newPed$MUSFWBand)

newPed[is.na(newPed)] <- 0

newPed$FUSFWBand[which(newPed$FBreeder == "6893")] <- "1"
newPed$FUSFWBand[which(newPed$FBreeder == "-LQ")] <- "2"
newPed$MUSFWBand[which(newPed$MBreeder == "RLC-U")] <- "3"
newPed$MUSFWBand[which(newPed$MBreeder == "-SW_")] <- "4"

##optional but need as.character
##removes dash in band #
#newPed$USFWBand <- gsub("-", "", newPed$USFWBand)
#newPed$USFWBand <- as.numeric(newPed$USFWBand)
#newPed$FUSFWBand <- gsub("-", "", newPed$FUSFWBand)
#newPed$MUSFWBand <- gsub("-", "", newPed$MUSFWBand)
#newPed$FUSFWBand <- as.numeric(newPed$FUSFWBand)
#newPed$MUSFWBand <- as.numeric(newPed$MUSFWBand)

newPed$USFWBand <- as.character(newPed$USFWBand)
newPed$FUSFWBand <- as.character(newPed$FUSFWBand)
newPed$MUSFWBand <- as.character(newPed$MUSFWBand)

newPed$USFWBand[which(newPed$JayID =="6893")] <- "1"
newPed$USFWBand[which(newPed$Jay=="-LQ")] <- "2"
newPed$USFWBand[which(newPed$JayID=="RLC-U")] <- "3"
newPed$USFWBand[which(newPed$JayID =="-SW_")] <- "4"

str(newPed)

### 2 15 2017 have to add in a year for those birds I think
founders1 <- subset(newPed, newPed$FUSFWBand==0 & newPed$MUSFWBand==0)
# Run BreedersJuly script through line 98 for the breeder data set to get year
str(brd)
str(founders1)
#founders1 - parent info for founders 
founders1$Year <- NULL
founders1$Censor <- NULL
found.year <- brd[,1:6]
colnames(found.year)[1] <- "JayID"
found.year[,4:5] <- NULL
colnames(found.year)[2] <- "USFWBand"
found.year["SexIND"] <- 0
found.year$SexIND[which(found.year$Sex == "F")] <- 2
found.year$SexIND[which(found.year$Sex == "M")] <- 1

#JayID and year of breeding, I actually don't think this step is necessary
#I added it thinking that this was the problem with the kinship matrix 
found.year$USFWBand <- NULL
found.year$SexIND <- NULL
found.year$Sex <- NULL

#founders with year of first breeding 
founders <- merge(founders1, found.year, by = "JayID")
founders["yearind"] <- founders$Year - 1
founders$Year <- NULL
colnames(founders)[8] <- "Year"

#Who is in founders1 and not in founders? 
notin <- subset(founders1, !(founders1$JayID %in% founders$JayID))
notin <- notin[,c(2,1,3,4,5,6,7)]
#need the year
notin["Year"] <- 1978
founders["MinAgeFBr"] <- NULL
founders["yearind"] <- NULL

birds <- rbind(founders, notin)
birds <- birds[,c(2,1,3,4,5,6,7,8)]
birds["Censor"] <- 0

### Have to add pedigree info to the original data set 
### yrlg.df1 + birds
birds2 <- subset(newPed, !newPed$FUSFWBand == 0)
birds2$Year <- NULL
birds2$Censor <- NULL
names(birds2)


#re-arrange dataframe with covariates to merge
yrlg.df1 <- yrlg.df1[,c(2,3,4,1,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)]
yrlg.df1$USFWBand <- as.character(yrlg.df1$USFWBand)

#called please as in "PLEASE WORK!" 
#dataframe with parent info and all covariates
please <- merge(yrlg.df1, birds2, by = c("USFWBand","JayID"))

birds["NatalNest"] <- NA
birds["TerrYr"] <- NA
birds["Terr"] <- NA
birds["HatchDate"] <- NA
birds["FldgDate"] <- NA
birds["LastObsDate"] <- NA
birds["Mass"] <- NA
birds["MeasDate"] <- NA
birds["HatchNum"] <-  NA
birds["FldgNum"] <- NA
birds["Days"] <- NA
birds["Years"] <- NA
birds["Day11"] <- NA
birds["Help"] <-  NA
birds["HelpNum"] <- NA
birds["Censor"] <- NA

birds <- birds[,c(1,2,10,11,12,8,13,14,15,17,18,19,20,21,22,9,23,
                  24,16,3,4,5,6,7)]

new.df <- rbind(please, birds)
str(new.df)

## 2 20 2017 need to add family id
ids <- makefamid(new.df$USFWBand, new.df$MUSFWBand, new.df$FUSFWBand)
famid <- cbind(ids)
colnames(famid) <- "famid"
new.df <- cbind(new.df, famid)
str(new.df)


## 3 5 2017 Don't run this code 
#new.df$USFWBand <- gsub("-", "", new.df$USFWBand)
#new.df$USFWBand <- as.numeric(new.df$USFWBand)
#new.df$USFWBand <- as.integer(new.df$USFWBand)
#new.df$FUSFWBand <- gsub("-", "", new.df$FUSFWBand)
#new.df$MUSFWBand <- gsub("-", "", new.df$MUSFWBand)
#new.df$FUSFWBand <- as.numeric(new.df$FUSFWBand)
#new.df$MUSFWBand <- as.numeric(new.df$MUSFWBand)
#new.df$FUSFWBand <- as.integer(new.df$FUSFWBand)
#new.df$MUSFWBand <- as.integer(new.df$MUSFWBand)

#
#new.df$USFWBand <- as.numeric(new.df$USFWBand)
#new.df$FUSFWBand <- as.numeric(new.df$FUSFWBand)
#new.df$MUSFWBand <- as.numeric(new.df$MUSFWBand)
#new.df$gid <- paste(new.df$famid, new.df$JayID, sep="/")
#still duplicated 

#jay.ped <- with(new.df, pedigree(id = JayID, dadid=MBreeder,
#momid = FBreeder, sex = SexIND, famid = famid, missid = 0))

jay.ped <- with(new.df, pedigree(id = USFWBand, dadid=MUSFWBand, 
                                 momid = FUSFWBand, sex = SexIND, famid =famid, missid = 0))

#family 1 is the largest because it's everyone who is connected
#the other families are disconnected subfamilies 
ped1 <- jay.ped[1]
plot(ped1)
ped2 <- jay.ped[2]
plot(ped2)
ped3 <- jay.ped[3]
plot(ped3)
ped4 <- jay.ped[4]
plot(ped4)
 
ped5 <- jay.ped[5]
plot(ped5)
#this is good
ped6 <- jay.ped[6]
plot(ped6)
ped8 <- jay.ped[8]
plot(ped8)
ped9 <- jay.ped[9]
plot(ped9)
ped10 <- jay.ped[10]
plot(ped10)
ped15 <- jay.ped[15]
plot(ped15)
ped20 <- jay.ped[20]
plot(ped20)
ped25 <- jay.ped[25]
plot(ped25)
ped30 <- jay.ped[30]
plot(ped30)
ped45 <- jay.ped[45]
plot(ped45)

#estimate kinship matrix
kins <- kinship(jay.ped)

allFY <- Surv(new.df$Days, new.df$Censor, type = 'right')
allFY1 <- coxph(allFY ~ Mass + Help, data = new.df)
allFY2 <- coxph(allFY ~ Mass + Help + Year, data = new.df)
allFY3 <- coxme(allFY ~ Mass + Help + (1|Year), data = new.df)
allFY4 <- coxme(allFY ~ Mass + Help + (1|USFWBand) + (1|Year), 
                varlist=coxmeMlist(kins), data=new.df)
summary(allFY4)



all.FYdata <- coxme(allFY ~ Mass + Help + HatchNum + (1|USFWBand) +
                      (1|Year), varlist = coxmeMlist(kins), data =new.df)

##################################################################
## 3 27 2017 Add in parent age 
##need parent ages to add to yrlg.df1
#"ages" is from new.df, which focuses on HY birds with parental info
#got rid of NA's bc those are the founders, new.df has the ped info
ages <- na.omit(new.df)
#Year is nest year (I'm pretty sure)
ages <- ages[,c(1,2,3,4,6,21,22,23,24)]
ages.m <- ages[,c(1,2,4,5,6,7)]
ages.d <- ages[,c(1,2,4,5,8,9)]
colnames(ages.m)[4] <- "NestYear"
colnames(ages.d)[4] <- "NestYear"

#These are the breeders, above the focal birds are HY
parent.age <- all.brd[,c(1,2,6,10,11,12,18)]
names(parent.age)
#parent.age <- breeders.real[ ,c(1,5,14,15,18,19)]
#parent.age1 <- new.realdf[ ,c(2,7,16,17,20,1)]
#colnames(parent.age)[2] <- "NestYear"
#this is the band number of parents 
#colnames(parent.age)[1] <- "USFWBand"
mom.age <- subset(parent.age, sex == "F")
names(mom.age)
colnames(mom.age)[1] <- "FBreeder"
colnames(mom.age)[2] <- "FUSFWBand"

dad.age <- subset(parent.age, sex == "M")
names(dad.age)
colnames(dad.age)[1] <- "MBreeder"
colnames(dad.age)[2] <- "MUSFWBand"

# 6 2 2017
#bizzare problem with merge
#ages.m and ages.d are the focal birds/fledglings 
#mom.age and dad.age are the parents to be matched 
names(mom.age)
names(dad.age)
names(ages.m)
names(ages.d)

###Problem of NA's with yrs experience comes here - it's because some terr
#years are missing
new.d <- merge(ages.d, dad.age, by = c("MUSFWBand","NestYear"), all.x=TRUE)
new.m <- merge(ages.m, mom.age, by = c("FUSFWBand","NestYear"), all.x=TRUE)


new.m <- new.m[,c(1:9)]
names(new.m)
colnames(new.m)[6] <- "momid"
colnames(new.m)[1] <- "momband"
colnames(new.m)[8] <- "momexp"
colnames(new.m)[9] <- "momage"
#6 2 2017 something is messed up, for some reason merge code adding more cols
new.m[7] <- NULL
colnames(new.m)[5] <- "TerrYr"

#new.d <- new.d[,c(3,4,5,6,1,7,9,10,11)]
new.d <- new.d[,c(1:9)]
names(new.d)
colnames(new.d)[1] <- "dadband"
colnames(new.d)[6] <- "dadid"
colnames(new.d)[8] <- "dadexp"
colnames(new.d)[9] <- "dadage"
#6 2 2017 something weird happening 
new.d[7] <- NULL
colnames(new.d)[5] <- "TerrYr"

#new.d <- new.d[,c(1,2,5,7,8,9,10)]
#new.m <- new.m[,c(1,2,5,7,8,9,10)]

new1 <- merge(yrlg.df1, new.d, by = c("USFWBand","JayID"), all.x=TRUE)
new2 <- merge(new1, new.m, by = c("USFWBand","JayID"), all.x = TRUE)

yrlg.df1 <- new2

#6 2 2017 - duplicate cols for some reason now
#yrlg.df1 <- yrlg.df1[,c(1,2,3,4,5,6,7,8,9,10,
                        #11,12,13,14,15,16,17,18,19,22,23,25,
                        #29,30)]


#yrlg.df1[,20] <- NULL
#yrlg.df1[,22] <- NULL
#colnames(yrlg.df1)[4] <- "TerrYr"


#how to plot...?
#mexp0 <- subset(agexp.count, momexp == 0)
#mexp0 <- na.omit(mexp0)
#qplot(momage, freq, data = mexp0, geom="bar")
#ggplot(mexp0, aes(x=momage, y=freq)) +
#geom_bar(stat = "identity")


#counts <- na.omit(counts)
#gplot(counts, aes(x=momage, y=freq), fill=momexp) +
#geom_bar(colour="blue", stat = "identity")

#a <- ggplot(counts, aes(x=momage, y=freq))
#a + geom_point(aes(colour=factor(momexp), y=freq))

#b <- ggplot(counts, aes(x=momexp, y=freq))
#b + geom_point(aes(colour=factor(momage), y=freq))

#b+ geom_bar(aes(fill = factor(momage), y=freq), stat="identity")
  

#A few duplicates because I matched by TerrYr, there were some breakups
#of breeding pairs
#reps1 <- ddply(new1, .(new1$USFWBand), nrow)
#reps2 <- ddply(new2, .(new2$USFWBand), nrow)

sum(is.na(yrlg.df1$dadexp))
sum(is.na(yrlg.df1$momexp))

# 6 2 2017 - run models without terr data first so I can exclude NAs
library(tidyr)
#remove records with NA's for parent age and experience 
june2 <- yrlg.df1 %>% drop_na(24,25,30,31)
str(june2)

june2 %>% count(momage)
june2 %>% count(momexp)
mom.counts <- june2 %>% count(momexp, momage)

a <- ggplot(mom.counts, aes(x=momage, y=n))
a + geom_bar(aes(fill=factor(momexp), y=n), stat="identity")

## add in density 
june2 <- merge(june2, jay.den, by = ("Year"), all.x=TRUE)

#survival object 
june2.mob <- Surv(june2$Days, june2$Censor, type=c('right'))

## 6 2 2017 Models here for everything (including mxe)
mass2 <- coxph(june2.mob ~ Mass, data = june2)
hatch.num2 <- coxph(june2.mob ~ HatchNum, data=june2)
flg.num2 <- coxph(june2.mob ~ FldgNum, data=june2)
cox.year2 <- coxph(june2.mob ~ Year, data = june2)
summary(mass2)
summary(hatch.num2)
summary(flg.num2)
summary(cox.year2)

help1 <- coxph(june2.mob ~ Help, data = june2)
help2 <- coxph(june2.mob ~ HelpNum, data = june2)
summary(help1)
summary(help2)
extractAIC(help1)
extractAIC(help2)

my <- coxph(june2.mob ~ Mass + Year, data = june2)
summary(my)

myh <- coxph(june2.mob ~ Mass + Help + Year, data = june2)
summary(myh)
myhn <- coxph(june2.mob ~ Mass + HelpNum + Year, data = june2)
summary(myhn)
extractAIC(myh)
extractAIC(myhn)
myhd <- coxph(june2.mob ~ Mass + Help + Den + Year, data = june2)
myhdn <- coxph(june2.mob ~ Mass + HelpNum + Den + Year, data=june2)
extractAIC(myh)
extractAIC(myhd)
mhyi <- coxph(june2.mob ~ Mass + Help + Year + Year*Help, data = june2)
summary(mhyi)
extractAIC(mhyi)
test11 <- coxph(june2.mob ~ Mass + Year + Help + momexp + FldgNum, 
                data = june2)
test12 <- coxph(june2.mob ~ Mass + Year + Help + momexp + HatchNum, 
                data = june2)
test13 <- coxph(june2.mob ~ Mass + dadexp, data=june2)
test14 <- coxph(june2.mob ~ Mass + momexp, data = june2)


myha <- coxph(june2.mob ~ Mass + Year + Help + momexp, data = june2)
myhad <- coxph(june2.mob ~ Mass + Year + Help + momexp + Den, data=june2)
#age or exp, doesn't matter one isn't better or worse 
myhage <- coxph(june2.mob ~ Mass + Help + Year + momage, data =june2)
myhaii <- coxph(june2.mob ~ Mass + Year + Help + momexp + Year*Help, 
                data = june2)

extractAIC(hatch.num2)
extractAIC(flg.num2)
extractAIC(cox.year2)
extractAIC(help1)
extractAIC(help2)
extractAIC(mass2)

library(AICcmodavg)
Cand.models1 <- list()
Cand.models1[[1]] <- test11
Cand.models1[[2]] <- test12
Cand.models1[[3]] <- cox.year2
Cand.models1[[4]] <- help1
Cand.models1[[5]] <- help2
Cand.models1[[6]] <- mass2
Cand.models1[[7]] <- my
Cand.models1[[8]] <- myh
Cand.models1[[9]] <- mhyi
Cand.models1[[10]] <- myha
Cand.models1[[11]] <- myhad
Cand.models1[[12]] <- myhd
Modnames1 <- paste("Model", 1:length(Cand.models1), sep="")
hymods1 <- aictab(Cand.models1, Modnames1, sort = TRUE, second.ord = FALSE)
hymods1

Cand.models3 <- list()
Cand.models3[[1]] <- my
Cand.models3[[2]] <- myh
Cand.models3[[3]] <- mhyi
Cand.models3[[4]] <- myha
Cand.models3[[5]] <- myhaii
Cand.models3[[6]] <- test11
Cand.models3[[7]] <- test12
Modnames3 <- paste("Model", 1:length(Cand.models3), sep="")
hymods3 <- aictab(Cand.models3, Modnames3, sort = TRUE, second.ord=FALSE)
hymods3

extractAIC(my)
extractAIC(myh)
extractAIC(mhyi)
extractAIC(myha)
extractAIC(myhaii)

#Do same models but with year as a random term 
mxx1 <- coxme(june2.mob ~ Mass + (1|Year), data=june2)
mxx2 <- coxme(june2.mob ~ Mass + Help + (1|Year), data=june2)
mxx3 <- coxme(june2.mob ~ Mass + Help + momexp + (1|Year), data=june2)
mxx4 <- coxme(june2.mob ~ Mass + momexp + (1|Year), data=june2)
summary(mxx1)
summary(mxx2)
summary(mxx3)
summary(mxx4)

#mxx5 <- coxme(june2.mob ~ Mass + momexp + (1|Year) + (Help|Year) , 
              #data=june2)
mxden <- coxme(june2.mob ~ Mass + Den + (1|USFWBand) + (1|Year), 
              varlist = coxmeMlist(kins), data = june2)
#summary(mxden)
mxden1 <- coxme(june2.mob ~ Mass + Help + Den + (1|USFWBand) + (1|Year),
              varlist=coxmeMlist(kins), data = june2)
mxden2 <- coxme(june2.mob ~ Mass + Help + momexp + Den + (1|USFWBand) + 
                (1|Year), varlist=coxmeMlist(kins), data = june2)
#summary(mxden2)
mxden3 <- coxme(june2.mob ~ Mass + Help + momexp + (1|USFWBand) + (1|Year),
                varlist = coxmeMlist(kins), data=june2)
#summary(mxden3)
mxdenAB <- coxme(june2.mob ~ Mass + Help + (1|USFWBand) + (1|Year),
            varlist = coxmeMlist(kins), data = june2)

####################################################################
mxdennA <- coxme(june2.mob ~ Mass + HelpNum + Den + (1|USFWBand) + (1|Year),
                 varlist=coxmeMlist(kins), data=june2)
summary(mxdennA)
mxdennB <- coxme(june2.mob ~ Mass + HelpNum + momexp + (1|USFWBand) + (1|Year),
                 varlist=coxmeMlist(kins), data=june2)
summary(mxdennB)
mxdennC <- coxme(june2.mob ~ Mass + HelpNum + momexp + Den + (1|USFWBand) + 
                   (1|Year), varlist=coxmeMlist(kins), data=june2)
summary(mxdennC)

mxdenn <- coxme(june2.mob ~ Mass + HelpNum + (1|USFWBand) + (1|Year), 
                varlist = coxmeMlist(kins), data=june2)
summary(mxdenn)
Cand.modelsD <- list()
Cand.modelsD[[1]] <- mxdenn
Cand.modelsD[[2]] <- mxdennA
Cand.modelsD[[3]] <- mxdennB
Cand.modelsD[[4]] <- mxdennC
ModnamesD <- paste("Model", 1:length(Cand.modelsD), sep="")
tableD <- aictab(Cand.modelsD, ModnamesD, sort = TRUE, second.ord=FALSE)
tableD
#####################################################################
mxdennD <- coxme(june2.mob ~ Mass + Help + Den + momexp + HatchNum
      + (1|USFWBand) + (1|Year), varlist = coxmeMlist(kins), data = june2)
mxdennE <- coxme(june2.mob ~ Mass + Help + momexp + HatchNum 
              + (1|USFWBand) + (1|Year), varlist = coxmeMlist(kins),
              data = june2)
mxdennF <- coxme(june2.mob ~ Mass + Help + HatchNum + (1|USFWBand) +
                (1|Year), varlist = coxmeMlist(kins), data=june2)
summary(mxdennF)

Cand.modelsAAA <- list()
Cand.modelsAAA[[1]] <- mxdennD
Cand.modelsAAA[[2]] <- mxdennE
Cand.modelsAAA[[3]] <- mxdennF
Cand.modelsAAA[[4]] <- mxdenAB
Cand.modelsAAA[[5]] <- mxden1
Cand.modelsAAA[[6]] <- mxden2
Cand.modelsAAA[[7]] <- mxden3
ModnamesAAA <- paste("Model", 1:length(Cand.modelsAAA), sep="")
tableAAA <- aictab(Cand.modelsAAA, ModnamesAAA, sort = TRUE, second.ord=FALSE)
tableAAA


Cand.modelsF <- list()
Cand.modelsF[[1]] <- mxden
Cand.modelsF[[2]] <- mxden1
Cand.modelsF[[3]] <- mxden2
Cand.modelsF[[4]] <- mxden3

ModnamesF <- paste("Model", 1:length(Cand.modelsF), sep="")
tableF <- aictab(Cand.modelsF, ModnamesF, sort = TRUE, second.ord=FALSE)
tableF


mxc1 <- coxme(june2.mob ~ Mass + (1|USFWBand) + (1|Year),
          varlist = coxmeMlist(kins), data = june2)
mxc2 <- coxme(june2.mob ~ Mass + Help + (1|USFWBand) + (1|Year),
              varlist = coxmeMlist(kins), data = june2)
mxc3 <- coxme(june2.mob ~ Mass + Help + momexp + (1|USFWBand) + (1|Year),
              varlist = coxmeMlist(kins), data = june2)
mxc4 <- coxme(june2.mob ~ Mass + momexp + (1|USFWBand) + (1|Year), 
              varlist = coxmeMlist(kins), data = june2)
mxc5 <- coxme(june2.mob ~ Mass + momage + (1|USFWBand) + (1|Year), 
              varlist = coxmeMlist(kins), data = june2)
mxc6 <- coxme(june2.mob ~ Mass + Help + momage + (1|USFWBand) + (1|Year),
              varlist = coxmeMlist(kins), data = june2)
mxc7 <- coxme(june2.mob ~ Mass + Help + momexp + HatchNum +
                (1|USFWBand) + (1|Year), varlist=coxmeMlist(kins), 
              data = june2)
mxc8 <- coxme(june2.mob ~ Mass + Help + momexp + FldgNum + (1|USFWBand) 
              + (1|Year), varlist=coxmeMlist(kins), data= june2)
#mxc9 <- coxme(june2.mob ~ Mass + Help + momexp + Year*Help + (1|USFWBand),
              #varlist=coxmeMlist(kins), data=june2)
########################################################
summary(mxc1)
summary(mxc2)
summary(mxc3)
summary(mxc4)
summary(mxc5)
summary(mxc6)
########################################################

Cand.models4 <- list()
Cand.models4[[1]] <- mxc1
Cand.models4[[2]] <- mxc2
Cand.models4[[3]] <- mxc3
Cand.models4[[4]] <- mxc4
#Cand.models4[[5]] <- mxc5
#Cand.models4[[6]] <- mxc6
Cand.models4[[5]] <- mxc7
Cand.models4[[6]] <- mxc8
#Cand.models4[[7]] <- my
#Cand.models4[[10]] <- myh
#Cand.models4[[11]] <- mhy
#Cand.models4[[12]] <- myha
#Cand.models4[[13]] <- myhaii
#Cand.models4[[14]] <- test11
#Cand.models4[[15]] <- test12
Modnames4 <- paste("Model", 1:length(Cand.models4), sep="")
table4 <- aictab(Cand.models4, Modnames4, sort = TRUE, second.ord=FALSE)
table4


#doesn't work
#mxc7 <- coxme(june2.mob ~ Mass + Help + momexp + (1|USFWBand), varlist = coxmeMlist(kins), data = june2)


##########################################################################
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

# 3 5 2017
yrlg.df1 <- merge(yrlg.df1, terr.info, by = "TerrYr", all.x = TRUE)


#yrlg.df1 <- merge(yrlg.df1, terr.info, by="TerrYr")
names(yrlg.df1)
str(yrlg.df1)

#yrlg.df1[is.na(yrlg.df1)] <- 4
length(unique(yrlg.df1$USFWBand))

# Add columns for centered variables for territory metrics
#Studentized (I think) - standard deviation 
yrlg.df1["stdscr"] <- scale(yrlg.df1$scrb.count, center = FALSE, scale = TRUE)
#Centered - mean
yrlg.df1["centscr"] <- scale(yrlg.df1$scrb.count, center = TRUE, scale =FALSE)

yrlg.df1["stdtsf"] <- scale(yrlg.df1$tsf.count, center = FALSE, scale = TRUE)
yrlg.df1["centtsf"] <- scale(yrlg.df1$tsf.count, center = TRUE, scale = FALSE)

yrlg.df1["stdsize"] <- scale(yrlg.df1$TerrSize, center= FALSE, scale = TRUE)
yrlg.df1["centsize"] <- scale(yrlg.df1$TerrSize, center = TRUE, scale = FALSE)


#Change to numeric for survival object 
yrlg.df1$FldgDate <- as.numeric(yrlg.df1$FldgDate)
yrlg.df1$LastObsDate <- as.numeric(yrlg.df1$LastObsDate)

yrlg.df1$Years <- as.numeric(yrlg.df1$Years)
yrlg.df1$Days <- as.numeric(yrlg.df1$Days)
yrlg.df1$stdsize <- as.numeric(yrlg.df1$stdsize)

names(yrlg.df1)
yrlg.df1[,c(21,22,26,28,29,33)] <- NULL
names(yrlg.df1)

library(corrplot)
covars1 <- yrlg.df1[,c(20,21,22)]
#covars1["oakfrac"] <- covars1$scrb.count/covars1$TerrSize
#covars1["tsffrac"] <- covars1$tsf.count/covars1$TerrSize
colnames(covars1) <- c("Oak Scrub", "Territory Size", "2-9 year TSF")
corrs1 <- cor(covars1, use="complete.obs")

# new colors
my.cols <- colorRampPalette(c("black", "gray50", "gray87"),bias=0.5)

corrplot(corrs1, method="pie", type= "lower")
corrs1
#######################################################################
# PCA analysis for terr size, oak scrub, and tsf

#princomp
terrvars <- data.frame(yrlg.df1$TerrYr, yrlg.df1$TerrSize, yrlg.df1$scrb.count,
                       yrlg.df1$tsf.count)
colnames(terrvars) <- c("TerrYr","TerrSize", "OakScrub", "TSF")
terrvars <- na.omit(terrvars)
terrvars["reloak"] <- terrvars$OakScrub/terrvars$TerrSize
terrvars["reltsf"] <- terrvars$TSF/terrvars$TerrSize
terr.info1 <- terrvars[,2:4]
#terr.info2 <- terrvars[,5:6]
terrs <- terrvars[, 1]
#I don't know which one to use...?
#I think the absolute area, not the relative area
vars.pca1 <- prcomp(terr.info1, center = TRUE, scale. = TRUE)
#vars.pca2 <- prcomp(terr.info2, center = TRUE, scale. = TRUE)
#prints the rotation, or loadings 
print(vars.pca1)
#plot(vars.pca)
summary(vars.pca1)
#print(vars.pca2)
#summary(vars.pca2)


#plots variances of the PCs
plot(vars.pca1, type = "l")
biplot(vars.pca1)

plot(vars.pca2, type = "l")

terrs.reg <- terrvars[,2:4]
vrreg.pca <- prcomp(terrs.reg, center = TRUE, scale. =TRUE)
print(vrreg.pca)
summary(vrreg.pca)

library(devtools)
install_github("ggbiplot", "vqv")
library(ggbiplot)
#doesn't work
#g <- ggbiplot(vars.pca, obs.scale = 1, var.scale = 1, groups = terrvars$TerrYr,
#ellipse = TRUE, circle = TRUE)

theta <- seq(0,2*pi, length.out=100)
circle <- data.frame(x = cos(theta), y = sin(theta))
p <- ggplot(circle, aes(x,y)) + geom_path()

loadings <- data.frame(vars.pca1$rotation, 
                       .names = row.names(vars.pca1$rotation))
p + geom_text(data = loadings, mapping = aes(x = PC1, y = PC2, 
                                             label = .names, colour = .names)) + 
  coord_fixed(ratio=1) +
  labs(x = "PCI", y = "PC2") +
  theme(panel.grid = element_blank()) +
  theme_bw()

## Cox models with PCI and II as predictors for territory quality 
#PCA object
names(vars.pca1)
# x is what we want
require(caret)

#get scores from PC1 and 2, put in data frame 
pc.scores <- data.frame(vars.pca1$x)
pc.scores <- pc.scores[,1:2]
str(pc.scores)

#subset data frame to include birds with terr data
library(tidyr)
flg.sub <- yrlg.df1
#flg.sub1 <- flg.sub %>% drop_na(24,25,26,27,28,29,30,31,32)
#June 2
flg.sub1 <- flg.sub %>% drop_na(22,23,26,27,28,29,30,31,32,33,34,35,36)

#should be 1529, more now 4 15 2017
flg.sub1 <- cbind(flg.sub1, pc.scores)
plot(flg.sub1$PC1, flg.sub1$stdsize)
plot(flg.sub1$PC2, flg.sub1$stdscr)
plot(flg.sub1$PC2, flg.sub1$stdtsf)

##untransformed data 
plot(flg.sub1$PC1, flg.sub1$TerrSize)
plot(flg.sub1$PC1, flg.sub1$scrb.count)
plot(flg.sub1$PC1, flg.sub1$tsf.count)

plot(flg.sub1$PC2, flg.sub1$scrb.count)
plot(flg.sub1$PC2, flg.sub1$tsf.count)


str(flg.sub1)
names(flg.sub1)
flg.sub2 <- flg.sub1 %>% drop_na(22,23)
flg.sub2 <- flg.sub2[!duplicated(flg.sub2$USFWBand),]

#yearling territory data object 
yrlg.tob <- Surv(flg.sub2$Years, flg.sub2$Censor, type = c('right'))

#####################################################################
#4 3 2017
#Imputation of terr data 
library(mice)
library(VIM)
library(mi)

#should have all info plus mom and dad age/exp
str(yrlg.df1)
names(yrlg.df1)

new.yrlgdf <- yrlg.df1[,1:26]
str(new.yrlgdf)
names(new.yrlgdf)
colnames(new.yrlgdf)[24] <- "OakScrub"
colnames(new.yrlgdf)[26] <- "TSF"
new.yrlgdf$TerrSize <- as.numeric(new.yrlgdf$TerrSize)

#Impute terr data and parent ages
mice(new.yrlgdf)

########################################################################
## Models - survival object 

# 3 21 2017 I;m rearranging this script so it makes more sense
# I adjusted the data frame to include the territory info for birds 
# that have it AND to include birds that don't (NAs)
# makes it much easier because then I just have to build 1 data frame 

#Survival object for all models 
#yrlg.mob <- Surv(yrlg.df1$Years, yrlg.df1$Censor, type = c('right'))
#days
#Use days because there are some that lived 1 day (which R says is 0 years)
yrlg.dob <- Surv(yrlg.df1$Days, yrlg.df1$Censor, type = c('right'))

## June 5 2017 
#Do models here 
#june5 <- yrlg.df1 %>% drop_na(27:34)
june5 <- yrlg.df1 %>% drop_na(32:40)
june5.ob <- Surv(june5$Days, june5$Censor, type = c('right'))

#June 19 - some repeated band #s
rep.juv <- ddply(june5, .(june5$USFWBand), nrow)
june5a <- june5[!duplicated(june5$USFWBand),]
june5a.ob <- Surv(june5a$Days, june5a$Censor, type = c('right'))


## corrplot 




terr1 <- coxph(june5a.ob ~ Mass + stdscr + Year, data = june5a)
summary(terr1)
terr2 <- coxph(june5a.ob ~ Mass + stdsize + Year, data = june5a)
summary(terr2)
terr3 <- coxph(june5a.ob ~ Mass + stdtsf + Year, data = june5a)
summary(terr3)
terr4 <- coxph(june5a.ob ~ Mass + Year, data = june5a)
extractAIC(terr4)
extractAIC(terr1)
extractAIC(terr2)
extractAIC(terr3)

terrh <- coxph(june5a.ob ~ Mass + Help + Year, data = june5a)
terrh1 <- coxph(june5a.ob ~ Mass + Help + Year + stdsize, data = june5a)
terrh1i <- coxph(june5a.ob ~ Mass + Help + Year + stdsize + Year*Help,
                 data = june5a)
extractAIC(terrh)
extractAIC(terrh1)

Cand.models4 <- list()
Cand.models4[[1]] <- terr1
Cand.models4[[2]] <- terr2
Cand.models4[[3]] <- terr3
Cand.models4[[4]] <- terr4
Cand.models4[[5]] <- terrh
Cand.models4[[6]] <- terrh1
Cand.models4[[7]] <- terrh1i
Modnames4 <- paste("Model", 1:length(Cand.models4), sep="")
hymods4 <- aictab(Cand.models4, Modnames4, sort = TRUE, second.ord = FALSE)
hymods4

## June 21 u

fig3 <- coxph(june5a.ob ~ stdsize, data = june5a)
fig3.df <- data.frame(predict(fig3, type="risk", se.fit=TRUE))
fig3.df$stdsize <- june5a$stdsize
pf3 <- with(fig3.df, data.frame(stdsize, hazard = fit,
            lwr = fit - 1.96*se.fit, upr = fit + 1.96*se.fit))


fig3A <- ggplot(pf3, aes(stdsize, hazard)) + 
  geom_line() +
  xlab("Territory Size") +
  ylab("Predicted Relative Risk") + 
  theme_bw() +
  geom_hline(yintercept = 1.0, linetype = 2) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  geom_ribbon(data = pf3, aes(ymin=lwr, ymax=upr), alpha = 0.2) +
  scale_x_continuous(breaks=c(0,1,2,3))
fig3A

fig4 <- coxph(june5a.ob ~ TerrSize, data = june5a)
fig4.df <- data.frame(predict(fig4, type="risk", se.fit=TRUE))
fig4.df$TerrSize <- june5a$TerrSize
pf4 <- with(fig4.df, data.frame(TerrSize, hazard = fit,
            lwr = fit - 1.96*se.fit, upr = fit + 1.96*se.fit))

fig4A <- ggplot(pf4, aes(TerrSize, hazard)) + 
  geom_line() +
  xlab("Territory Size") +
  ylab("Predicted Relative Risk") + 
  theme_bw() +
  geom_hline(yintercept = 1.0, linetype = 2) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  geom_ribbon(data = pf4, aes(ymin=lwr, ymax=upr), alpha = 0.2) 
fig4A


#Mixed models
mxb1 <- coxme(june5a.ob ~ Mass + Help + stdsize + (1|Year), data=june5a)
mxb2 <- coxme(june5a.ob ~ Mass + Help + (1|Year), data = june5a)
mxb3 <- coxme(june5a.ob ~ Mass + (1|Year), data=june5a)

summary(mxb1)
summary(mxb2)
summary(mxb3)

den.dfAA <- merge(june5a, jay.den, by = ("Year"), all.x=TRUE)
june5b.ob <- Surv(den.dfAA$Days, den.dfAA$Censor, type='right')

mxb4 <- coxme(june5b.ob ~ Mass + (1|USFWBand) + (1|Year), 
      varlist = coxmeMlist(kins), data = den.dfAA)
mxb5 <- coxme(june5b.ob ~ Mass + Help + stdsize + (1|USFWBand) + (1|Year),
              varlist = coxmeMlist(kins), data = den.dfAA)
mxb6 <- coxme(june5b.ob ~ Mass + Help + (1|USFWBand) + (1|Year),
              varlist = coxmeMlist(kins), data = den.dfAA)
mxb7 <- coxme(june5b.ob ~ Mass + Help + Den + (1|USFWBand) + (1|Year),
              varlist = coxmeMlist(kins), data = den.dfAA)
mxb8 <- coxme(june5b.ob ~ Mass + Help + HatchNum + (1|USFWBand) + (1|Year),
              varlist = coxmeMlist(kins), data = den.dfAA)
mxb9 <- coxme(june5b.ob ~ Mass + Help + Den + stdsize + (1|USFWBand) +
                (1|Year), varlist = coxmeMlist(kins), data = den.dfAA)
summary(mxb9)

mxb10 <- coxme(june5b.ob ~ Mass + HelpNum + stdsize + HelpNum*stdsize +
                 (1|USFWBand) + (1|Year), varlist = coxmeMlist(kins), 
               data = den.dfAA)
mxb11 <- coxme(june5b.ob~ Mass + HelpNum + stdsize + (1|USFWBand) + 
              (1|Year), varlist= coxmeMlist(kins), data = den.dfAA)

# 6 9 2017 
Cand.models5 <- list()
Cand.models5[[1]] <- mxb4
Cand.models5[[2]] <- mxb5
Cand.models5[[3]] <- mxb6
Cand.models5[[4]] <- mxb7
Cand.models5[[5]] <- mxb8
Cand.models5[[6]] <- mxb9
Modnames5 <- paste("Model", 1:length(Cand.models5), sep="")
hymods5 <- aictab(Cand.models5, Modnames5, sort = TRUE, second.ord = FALSE)
hymods5







#######################################################################
## Cox PH Models 
#Brood size 
hatch.num <- coxph(yrlg.mob ~ HatchNum, data=yrlg.df1)
hatch.num
anova(hatch.num)

#Number of fledglings from nest
flg.num <- coxph(yrlg.mob ~ FldgNum, data=yrlg.df1)
flg.num
anova(flg.num)

#4 27 2017 fledge date? probably not the right way to do it 
flg.date <- coxph(yrlg.mob ~ FldgDate, data = yrlg.df1)
summary(flg.date)

#Year as a factor
cox.year <- coxph(yrlg.mob ~ Year, data = yrlg.df1)
cox.year
anova(cox.year)

#5 31 2017
yrb.h <- coxph(yrlg.mob ~ Mass + HelpNum + Year + HelpNum*Year, data=yrlg.df1)
#logistic regression 
lg.yh <- glm(Censor ~ Mass + HelpNum + Year:HelpNum, data = yrlg.df1, 
      family = binomial(link=logit))
lg.yh1 <- glm(Censor ~ Mass + HelpNum + Year, data = yrlg.df1, 
             family=binomial(link=logit))
summary(lg.yh1)

#Compare the first 3, brood size and fledge # add nothing 
dev.compare1 <- anova(hatch.num, flg.num, cox.year, test="Chisq")
dev.compare1

#Parent age and experience
#just a simple plot
plot(yrlg.df1$momage, yrlg.df1$Years)
plot(yrlg.df1$Years, yrlg.df1$momage)
plot(yrlg.df1$momexp, yrlg.df1$Years)

mom.age <- coxph(yrlg.mob ~ momage, data = yrlg.df1)
summary(mom.age)
cox.zph(mom.age, global=TRUE)
plot(cox.zph(mom.age))
anova(mom.age)

dadage <- coxph(yrlg.mob ~ dadage, data = yrlg.df1)
plot(cox.zph(dadage))
summary(dadage)
anova(mom.age)
anova(dadage)

##polynomials
momage.py <- coxph(yrlg.mob ~ momage + I(momage^2), data = yrlg.df1)
summary(momage.py)
anova(mom.age, momage.py)

momexp.mod <- coxph(yrlg.mob ~ momexp, data = yrlg.df1)
cox.zph(momexp.mod)
summary(momexp.mod)
plot(cox.zph(momexp.mod))
anova(momexp.mod)


dadexp.mod <- coxph(yrlg.mob ~ dadexp, data = yrlg.df1)
summary(dadexp.mod)
anova(dadexp.mod)

sim.moage <- coxsimLinear(mom.age, b = "momage", Xj = 1:15)
simGG(sim.moage)
sim.mexp <- coxsimLinear(momexp.mod, b = "momexp", Xj = 0:15)
simGG(sim.mexp)

#Plot of juvenile survival as a function of mom age at each experience level
#We want the predicted juvenile survival probability as a function 
#of mom age and experience 
#Given that my mom has 2 years of experience, it is better for my mom
#to be 4 than say 7? So at any given experience level, want to know
#the effect of age
temp.mom <- yrlg.df1
temp.momA <- temp.mom
temp.momA <- temp.momA %>% drop_na(23)
ob.A <- Surv(temp.momA$Years, temp.momA$Censor, type = c('right'))
zillion <- coxph(ob.A ~ Mass + momage, data = temp.momA)

sickofthis <- data.frame(predict(zillion, type="risk", se.fit=TRUE))
sickofthis$momage <- temp.momA$momage
predframe.age <- with(sickofthis, data.frame(momage, 
        hazard = fit, lwr = fit - 1.96*se.fit,
        upr = fit + 1.96*se.fit))

age.plot <- ggplot(predframe.age, aes(momage, hazard)) +
  geom_point() + 
  xlab("Age of Mother (Years)") +
  ylab("Predicted Relative Risk") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  scale_y_continuous(breaks = c(0.8,1.0,1.2,1.4,1.6,1.8,2.0)) +
  scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15))
age.plot


temp.mom <- temp.mom %>% drop_na(22,23)
temp.ob <- Surv(temp.mom$Years, temp.mom$Censor, type = c('right'))
temp.cox <- coxph(temp.ob ~ momage, data = temp.mom)
summary(temp.cox)
anova(temp.cox)
temp.cox1 <- coxph(temp.ob ~ momexp, data = temp.mom)
summary(temp.cox1)
anova(temp.cox1)
temp.cox2 <- coxph(temp.ob ~ momage + momexp, data = temp.mom)
summary(temp.cox2)
anova(temp.cox2)

temp.cox3 <- coxph(temp.ob ~ momage + momexp + momage*momexp, 
                   data = temp.mom)
summary(temp.cox3)
anova(temp.cox2, temp.cox3)

temp.mom$pred.age <- predict(temp.cox, temp.mom, type="risk")
temp.mom$pred.exp <- predict(temp.cox1, temp.mom, type="risk")
temp.mom$pred.both <- predict(temp.cox2, temp.mom, type="expected")
temp.mom$expsur <- exp(-temp.mom$pred.both)                    

plot(temp.mom$momage, temp.mom$momexp)

p1 <- ggplot(temp.mom, aes(x=momage, y=momexp))
p1 + geom_point()

str(summary(temp.cox2))

trts <- expand.grid(momexp=temp.mom$momexp, momage=temp.mom$momage)
trts$Risk <- predict(temp.cox2, trts, type="risk")
trts$RiskAge <- predict(temp.cox, trts, type="risk")
trts$RiskExp <- predict(temp.cox1, trts, type="risk")
trts$survae <- predict(temp.cox1, trts, type="expected")
trts$Inter <- predict(temp.cox3, trts, type = "risk")
trts <- trts[!duplicated(trts),]
trts <- subset(trts, momexp < momage)
trts <- trts[order(trts[,2,1]),]

ggplot(trts, aes(momage, RiskAge)) +
  geom_point()
ggplot(trts, aes(momexp, RiskExp)) +
  geom_point()

p.momint <- ggplot(trts, aes(momexp, Inter)) + geom_point()
p.momint + facet_grid(.~momage) + 
xlab("Mom Experience") +
ylab("Relative Risk")

## 4 15 2017 adding to dissertation chapter/manuscript 
## 3 29 2017 Getting closer with the plot below 
p.mom1 <- ggplot(trts, aes(momage, Risk)) + geom_point() 
p.mom1 + facet_grid(.~momexp) +
  xlab("Mom Age") +
  ylab("Predicted Relative Risk") +
  theme_bw()+
theme(panel.grid.major.y = element_blank(), 
      panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  scale_y_continuous(breaks = c(0.8, 1.0, 1.2, 1.4))


ggplot(trts, aes(x=momage,y=momexp)) +
  geom_point(aes(size=Risk))

ggplot(trts, aes(momage, momexp)) + 
  geom_point(aes(colour=Risk))

qplot(momage, Risk, data = trts, colour = momexp)
qplot(momage, momexp, data=trts, colour= Risk)
##only bottom half of plot makes sense
##email Gordon about this 

#Now use each one (exp and age separately by dropping NA's sep)
temp.mom2 <- yrlg.df1
temp.mom2 <- temp.mom2[!duplicated(temp.mom2$USFWBand),]

temp.mom2 <- temp.mom2 %>% drop_na(23)
length(unique(temp.mom2$USFWBand))
tm2.ob <- Surv(temp.mom2$Days, temp.mom2$Censor, type = c('right'))

tm2.exp <- coxph(tm2.ob ~ Mass + momexp, data = temp.mom2)
summary(tm2.exp)
anova(tm2.exp)
tm2.exp1 <- coxph(tm2.ob ~ Mass + momexp + HelpNum, data = temp.mom2)
tm2.exp2 <- coxph(tm2.ob ~ Mass + momexp +Help, data = temp.mom2)
summary(tm2.exp1)
summary(tm2.exp2)
anova(tm2.exp1)
anova(tm2.exp2)
tm2.exp3 <- coxph(tm2.ob ~ Mass + Help, data = temp.mom2)
tm2.exp4 <- coxph(tm2.ob ~ Mass + HelpNum, data = temp.mom2)
anova(tm2.exp3)
anova(tm2.exp4)

fit.exp <- data.frame(predict(tm2.exp, type="risk",se.fit=TRUE))
fit.exp$momexp <- temp.mom2$momexp
predframe.exp <- with(fit.exp, data.frame(momexp, 
                      hazard = fit, lwr = fit- 1.96*se.fit,
                      upr = fit+ 1.96*se.fit))

exp.plot <- ggplot(predframe.exp, aes(momexp, hazard)) +
  geom_point() +
  xlab("Experience of Mother (Years)") +
  ylab("Predicted Relative Hazard") +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14)) +
  scale_y_continuous(breaks = c(0.8,1.0,1.2,1.4,1.6))
  
exp.plot

###Try AFT model to get survival probabilites 
#temp.ob - surv object
momAFT <- survreg(temp.ob ~ momage, data = temp.mom, dist="weibull")
summary(momAFT)
names(momAFT)

#pct <- 1:98/100
#ptime <- predict(momAFT, type='quantile', p=c(.1,.5,.9))
#with(temp.mom, plot(momage, Years, xlab = "momage", ylav = "Years",))
#matlines(1:65, ptime)

#Mass
mass.full <- coxph(yrlg.mob ~ Mass, data = yrlg.df1)
mass.full <- coxph(yrlg.dob ~ Mass, data = yrlg.df1)
summary(mass.full)
anova(mass.full)
temp <- subset(yrlg.df1, yrlg.df1$Year == 2005)
temp1 <- subset(temp, temp$Mass == 45.5)
temp2 <- subset(temp, temp$Mass == 33.7)
temp3 <- subset(temp, temp$Mass == 56.2)
sfit <- survfit(mass.full, temp)
sfit1 <- survfit(mass.full, temp1)
sfit2 <- survfit(mass.full, temp2)
sfit3 <- survfit(mass.full, temp3)
plot(sfit1, xlim = c(0,365))
lines(sfit2)
lines(sfit3)

#compare mass within year
tempA <- subset(yrlg.df1, yrlg.df1$Year == 1997)
tempA1 <- subset(tempA, tempA$Mass == 45)
tempA2 <- subset(tempA, tempA$Mass == 59.8)
tempA3 <- subset(tempA, tempA$Mass == min(tempA$Mass))
sfitA1 <- survfit(mass.full, tempA1)
sfitA2 <- survfit(mass.full, tempA2)
sfitA3 <- survfit(mass.full, tempA3)
plot(sfitA2, xlim = c(0,365))
lines(sfitA1, lty = 2)
lines(sfitA2, lty = 3)


# 4 27 2017
#from internet example
tab <- data.frame(table(yrlg.df1[yrlg.df1$Censor == 1, "Days"]))
y <- as.numeric(levels(tab[, 1]))[tab[, 1]]
d <- tab[, 2]
H0 <- basehaz(mass.full, centered = FALSE)
H0 <- H0[H0[, 2] %in% y, ]
betaHat <- mass.full$coef
#h0 <- rep(NA, length(y))
#for(i in 1:length(y)) {
 # h0[i] <- d[i] / sum(exp(yrlg.df1[yrlg.df1$time >= y[1], "Mass"] * betaHat))
  #print(h0)
  #
#H0$surv <- (-integrate(exp(H0$hazard)))


#cbind(H0, cumsum(h0))

#str(summary(survfit(mass.full), time = 0.5))
#str(summary(survfit(mass.full), time = 1))
#str(summary(survfit(mass.full), time = 5))


#plot of scaled Schoenfield residuals to check PH assumption 
plot(cox.zph(mass.full))
bsl.haz <- basehaz(mass.full)
plot(bsl.haz$time, bsl.haz$hazard, type = 'l', xlim = c(0,1))
lines(bsl.haz$time, 0.979*bsl.haz$hazard)
#simPH method - I don't understand what it is plotting here, this 
#does not match with the predict funciton...? 
sim.mass <- coxsimLinear(mass.full, b = "Mass", qi = "Relative Hazard", 
                         Xj = c(15, 65), spin = TRUE)
simGG(sim.mass,type="ribbons", alpha = 0.35)

sim.mass2 <- coxsimLinear(mass.full, b = "Mass", 
                      qi = "Hazard Ratio", Xj = c(15,65))
simGG(sim.mass2)

mass.poly <- coxph(yrlg.mob ~ Mass + I(Mass^2), data = yrlg.df1)
summary(mass.poly)
preds.mass <- predict(mass.poly, type="risk", se.fit=TRUE)
preds.mass <- as.data.frame(preds.mass)
preds.mass <- data.frame(preds.mass, yrlg.df1$Mass)
colnames(preds.mass) <- c("fit", "se", "Mass")

coefMass <- coef(mass.poly)
mean.mass <- mean(yrlg.df1$Mass)
rMean.mass <- exp(coefMass["Mass"])
r1234.mass <- exp(coefMass["Mass"]*(as.numeric(yrlg.df1[1:4, "Mass"])))
#relative risk of first 4 individuals divided by mean relative risk
r1234.mass/rMean.mass
relRisk.mass <- predict(mass.poly, yrlg.df1, type="risk")
risks2 <- data.frame(relRisk.mass[1:50], yrlg.df1$Mass[1:50])
colnames(risks2) <- c("RelRisk", "Mass")
plot(risks2$Mass, risks2$RelRisk)

fit.mass <- data.frame(predict(mass.poly, type="risk", se.fit=TRUE))
fit.mass$Mass <- yrlg.df1$Mass
colnames(fit.mass)[3] <- "Mass"
predframe.mass <- with(fit.mass, data.frame(Mass, hazard=fit, 
                                  lwr = fit - 1.96*se.fit, upr = fit + 1.96*se.fit))

##Mass poly - adapted from code written for mass with no polynomial
pB <- ggplot(predframe.mass, aes(Mass, hazard)) +
  geom_line()
pB + ylim(0,2) +
  xlab("Nestling Mass (g)") +
  ylab("Predicted Relative Risk") + 
  theme_bw() +
  geom_hline(yintercept = 1.0, linetype = 2) +
  geom_ribbon(data = predframe.mass, aes(ymin=lwr, ymax=upr), alpha = 0.2) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 
#### Looking at this graph I don't think there is much support for 
# a non-linear relationship
# it's approximately linear 

sim.poly <- coxsimPoly(mass.poly, b = "Mass", qi ="Hazard Rate", pow=2, 
                       Xj = c(19,66), ci = 0.95)
simGG(sim.poly)

####### 4 10 2017 Use this 

#predicted values 
preds <- predict(mass.full, type="risk", se.fit=TRUE)
preds <- as.data.frame(preds)
preds <- data.frame(preds, yrlg.df1$Mass)
colnames(preds) <- c("fit","se","Mass")

#From John Fox's appednix in applied regression book 
#Get coefficients 
coefCPH <- coef(mass.full)
mean.mass <- mean(yrlg.df1$Mass)
rMean <- exp(coefCPH["Mass"])
r1234 <- exp(coefCPH["Mass"]*(as.numeric(yrlg.df1[1:4, "Mass"])))
#relative risk of first 4 individuals divided by mean relative risk
r1234/rMean
relRisk <- predict(mass.full, yrlg.df1, type="risk")
risks1 <- data.frame(relRisk[1:50], yrlg.df1$Mass[1:50])
colnames(risks1) <- c("RelRisk", "Mass")
plot(risks1$Mass, risks1$RelRisk)



###use this an change to ggplot2
plot(preds$yrlg.df1.Mass, preds$fit)
lines(preds$yrlg.df1.Mass, preds$se.fit)

  
## Example online for nlme 
fit <- data.frame(predict(mass.full, type="risk", se.fit=TRUE))
fit$Mass <- yrlg.df1$Mass
colnames(fit)[3] <- "Mass"
predframe <- with(fit, data.frame(Mass, hazard=fit, 
        lwr = fit - 1.96*se.fit, upr = fit + 1.96*se.fit))

##4 11 2017 as long as this is correct, then this is nice 
pA <- ggplot(predframe, aes(Mass, hazard)) +
  geom_line()
pA + ylim(0,2) +
  xlab("Nestling Mass (g)") +
  ylab("Predicted Relative Risk") + 
  theme_bw() +
  geom_hline(yintercept = 1.0, linetype = 2) +
  geom_ribbon(data = predframe, aes(ymin=lwr, ymax=upr), alpha = 0.2) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 

#survival curve estimated from cox model
plot(survfit(mass.full, se.fit=TRUE, conf.int=.95), xlim=c(0,1))

## Make same plot but add helper number?
## Use fill or something like that for helpers 


#### Not correct but this is the way, integration of the hazard 
# function 
h <- function(m) exp((basehaz(mass.full)) + 0.979*m)
m <- yrlg.df1$Mass
t <- hy.fit$time
exp(-integrate(h, 0, t)$value)


test.plot <- termplot(mass.full, data = yrlg.df1,
                  se=TRUE, terms = "residuals")

#Helper presence/absence 1,0
cx.hlpr <- coxph(yrlg.mob ~ Help, data=yrlg.df1)
plot(cox.zph(cx.hlpr))
summary(cx.hlpr)
anova(cx.hlpr)
#simPH method
sim.help <- coxsimLinear(cx.hlpr, b = "Help", 
                        Xj = 0:1)
#3 6 2017 Figure for manuscript
helpplot <- simGG(sim.help, xlab = "Helper", method = "lm", type = "points", alpha = 0.8)
helpplot
plota <- helpplot + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank()) +
        scale_x_continuous(name="Helper Presence", breaks = c(0,1))
plota      
#predict function
predvals <- predict(cx.hlpr, type="risk")
#plot(yrlg.df1$Help, predvals)

#Helper number as a count 
cx.group <- coxph(yrlg.mob ~ HelpNum, data=yrlg.df1)
summary(cx.group)
anova(cx.group)
#hlphaz <- basehaz(cx.group)
#plot(hlphaz$time, hlphaz$hazard, type = "l", xlim = c(0,1), ylim = c(0,2))
#lines(hlphaz$time, exp(0.95059)*hlphaz$hazard, col = "green")
#simPH method
#Helper number is an integer 
sim.help2 <- coxsimLinear(cx.group, b = "HelpNum", qi = "Relative Hazard",
                      Xj = c(0,6))
helpplot2 <- simGG(sim.help2, xlab = "Helper Number",
                   ylab = "Relative Hazard", alpha = 0.4)
plotb <- helpplot2 + ylab("Simulated Relative Hazard") +
  theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
                  scale_x_continuous(breaks =c(0,1,2,3,4,5,6),
                  scale_y_continuous(expand=c(0,0))+
                  labels=c("0","1","2","3","4","5","6")) +
                  theme(axis.title.y = element_blank()) 
                 
                  
plotb

#non-linear
help.poly <- coxph(yrlg.mob ~ HelpNum + I(HelpNum^2), data = yrlg.df1)
summary(help.poly)
anova(cx.group, help.poly)
#no support for this term 


#adjust scaling
multiplot(plota, plotb, cols = 2)
#This figure is close to what I want, needs some aes. adj
sim.help3 <- coxsimLinear(cx.group, b ="HelpNum",
                  qi = "Hazard Ratio", Xj = c(0,6))
simGG(sim.help3, alpha = 0.4, xlab = "Helper Number")
#predict
predvals1 <- predict(cx.group, type="risk")
plot(yrlg.df1$HelpNum, predvals1,
     xlab = "Natal Territory Helper Number", ylab="Predicted Risk Scores")


#Mass and Helper (0,1)
hlp.mass <- coxph(hy.ob ~ Mass + Help, data = yrlg.df1)
summary(hlp.mass)
anova(hlp.mass)
predvals2 <- predict(hlp.mass, type="risk")
plot(yrlg.df1$Help, predvals2, xlab = "Helper Presence", 
     ylab = "Hazard Ratio", main = "Predicated Values")

#mass, helper 0:1, and predicted risk scores 
#treatments <- expand.grid(Helper = yrlg.df1$Help, risks = predvals2) #too large 

tmp <- data.frame(yrlg.df1$Help, yrlg.df1$HelpNum, yrlg.df1$Mass, predvals2)
colnames(tmp)[1] <- "Help"
colnames(tmp)[2] <- "HelpNum"
colnames(tmp)[3] <- "Mass"

tmp$Help <- as.factor(tmp$Help)
tmp$HelpNum <- as.factor(tmp$HelpNum)

#qplot(Help, predvals2, data = tmp)
#qplot(HelpNum, predvals2, data = tmp)
#qplot(HelpNum, predvals2, data = tmp, colour = yrlg.df1$Mass)
#qplot(yrlg.df1$Mass, predvals2, data = tmp, colour=HelpNum)
#qplot(HelpNum, yrlg.df1$Mass, data = tmp, colour=predvals2)
#Plot of predicted values given a mass of risk/hazard ratio as a function of
#helpers
treatments <- expand.grid(Help = levels(tmp$Help), 
            Mass = tmp$Mass)
temp <- treatments
temp$Risk <- predict(hlp.mass, temp, interval = c("confidence"), type="risk")
qplot(Mass, Risk, data = temp, colour = Help)

#get standard errors of predicted values 
prdse <- predict(hlp.mass, temp, interval = c("confidence"), 
                     type = "Relative Risk", se.fit=TRUE)

temp$SEfit <- prdse$se.fit

#2 15 2017 add error bars or some measure of uncertainty

# 3 3 2017 Code for publishable figure 
#4 15 2017

ggplot(data = temp, aes(x = Mass, y = Risk, group=Help), colour=Help) +
  geom_line(size = 1) +
  xlab("Nestling Mass (g)") + 
  ylab("Predicted Relative Risk")
  scale_color_manual(values = c("darkgrey", "black")) +
  #theme_bw() + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  scale_y_continuous(breaks = c(0.6,0.8,1.0,1.2,1.4,1.6,1.8)) +
    geom_hline(yintercept=1)
  
### 4 17 2017
ggplot(data = temp, aes(Mass, Risk, group=Help, colour=Help)) +
  geom_line(size=1) +
  scale_color_manual(values= c("darkgrey", "black")) +
  ylab("Predicted Relative Risk") +
  xlab("Nestling Mass (g)") +
  theme_bw() +
  scale_y_continuous(breaks = c(0.6,0.8,1,1.2,1.4,1.6,1.8)) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  geom_hline(yintercept=1, lty=2)
  

ggplot(data = temp, aes(Mass, Risk, group=Help, colour=Help)) +
  geom_point(size = 1) +
  scale_color_manual(values= c("blue", "darkorange"))


#Mass and Helper number
hlpn.mass <- coxph(hy.ob ~ Mass + HelpNum, data = yrlg.df1)
summary(hlpn.mass)
anova(hlpn.mass)
predvals3 <- predict(hlpn.mass, type = "risk")
plot(yrlg.df1$HelpNum, predvals3)

anova(mass.full, hlpn.mass)
anova(mass.full, hlp.mass)

# 4 15 2017

fit.helpers <- data.frame(predict(hlpn.mass, type = "risk", se.fit=TRUE))
fit.helpers$HelpNum <- yrlg.df1$HelpNum
fit.helpers$Mass <- yrlg.df1$Mass
predframeH <- with(fit.helpers, data.frame(HelpNum, Mass, 
                  hazard = fit, 
                  lwr = fit - 1.96*se.fit, 
                  upr = fit + 1.96*se.fit))
hlp.plot <- ggplot(predframeH, aes(Mass, hazard), color=HelpNum) +
  geom_point() + 
  xlab("Mass(g)") + 
  ylab("Predicted Relative Hazard") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  scale_y_continuous(breaks = c(0.8,1.0,1.2,1.4,1.6))
hlp.plot

#Plot for showing the predicted values of hazard ratio from model of
#risk as a function of mass and helper number 
treats <- expand.grid(HelpNum = levels(tmp$HelpNum), 
                      Mass = tmp$Mass)
treats$Risk <- predict(hlpn.mass, type='risk')
#how to get SE??

# 4 17 2017 FUCK FUCK FUCK FUCK FUCK FUCK 
ggplot(data = treats, aes(Mass, Risk, group=HelpNum, colour=HelpNum)) +
  geom_point() + 
  ylab("Predicted Relative Risk") +
  xlab("Nestling Mass (g)") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  geom_hline(yintercept=1, lty=2)


#temp2 <- treats
treats2 <- data.frame(predict(hlpn.mass, type = 'risk', se.fit=TRUE))
treats.se <- data.frame(predict(hlpn.mass, type = 'risk', se.fit = TRUE))
#temp2$Risk <- predict(hlpn.mass, temp2, type = "risk")
#qplot(Mass, Risk, data = temp2, colour=HelpNum)

str(treats)


#interaction between helper number and mass
hlpmni <- coxph(yrlg.mob ~ Mass*HelpNum, data = yrlg.df1)
summary(hlpmni)
anova(hlpmni)
predvals4 <- predict(hlpmni, type="risk")
plot(yrlg.df1$Mass, predvals4)
plot(yrlg.df1$HelpNum, predvals4)
anova(mass.full, hlpn.mass, hlpmni)

#interaction between helpers and mass, does this even make biological sense?
hlpmi <- coxph(yrlg.mob ~ Mass*Help, data=yrlg.df1)
summary(hlpmi)
anova(hlpmi)
anova(mass.full, hlpmi)

#mom age 
mha <- coxph(yrlg.mob ~ Mass + momage + HelpNum, data= yrlg.df1)
summary(mha)
anova(mha)

#Help num or age? They are correlated
new1 <- coxph(yrlg.mob ~ Mass + momage, data = yrlg.df1)
anova(new1)
newb <- coxph(yrlg.mob ~ Mass + momexp, data = yrlg.df1)
anova(newb)
new2 <- coxph(yrlg.mob ~ Mass + HelpNum, data = yrlg.df1)
summary(new2)

anova(mass.full, hlp.mass, hlpn.mass, hlpmi, hlpmni)

extractAIC(mass.full)
extractAIC(hlp.mass)
extractAIC(hlpn.mass)
extractAIC(hlpmi)
extractAIC(hlpmni)

##AIC table using AICcmodavg package
Cand.models <- list()
Cand.models[[1]] <- mass.full
Cand.models[[2]] <- hlp.mass
Cand.models[[3]] <- hlpn.mass
Cand.models[[4]] <- hlpmi
Cand.models[[5]] <- hlpmni

Modnames <- paste("Model", 1:length(Cand.models), sep="")

#AIC
hymods <- aictab(Cand.models, Modnames, sort = TRUE)
#AICc                       
hymods2 <- aictab(Cand.models, Modnames, sort = TRUE,
                  second.ord = TRUE)
hymods
hymods2

fuck <- coxph(yrlg.mob ~ Mass + momage + HelpNum, data = yrlg.df1)
anova(new1, fuck)

########################################################################
#Territory variables 
#use yrlg.tob

#PC1 and 2
cx.pc <- coxph(yrlg.tob ~ PC1 + PC2, data = flg.sub2)
summary(cx.pc)
cx.pc1 <- coxph(yrlg.tob ~ PC1, data = flg.sub2)
summary(cx.pc1)

anova(cx.pc1, cx.pc)

mpc1 <- coxph(yrlg.tob ~ Mass + PC1, data = flg.sub2)
summary(mpc1)

fit.pc1 <- data.frame(predict(cx.pc1, type="risk", se.fit=TRUE))
fit.pc1$PC1 <- flg.sub2$PC1
colnames(fit.pc1)[3] <- "PC1"
predframe.PC1 <- with(fit.pc1, data.frame(PC1, hazard=fit, 
      lwr = fit - 1.96*se.fit, upr = fit + 1.96*se.fit))

##4 17 2017 as long as this is correct, then this is nice 
pA <- ggplot(predframe.PC1, aes(PC1, hazard)) +
  geom_line()
pA + ylim(0,2) +
  xlab("PC1") +
  ylab("Predicted Relative Risk") + 
  theme_bw() +
  geom_hline(yintercept= 1.0, linetype = 2) +
  geom_ribbon(data = predframe, aes(ymin=lwr, ymax=upr), alpha = 0.2) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 

mmpc <- coxph(yrlg.tob ~ Mass + PC1 + momexp, data = flg.sub2)
summary(mmpc)
anova(mmpc)
bb <-coxsimLinear(mmpc, "momexp", qi = "Relative Hazard", Xj = c(1,15))
simGG(bb)
pc1 <- coxsimLinear(mmpc, "PC1", qi = "Relative Hazard", 
                    Xj = c(-2.8, 7.6))
simGG(pc1)

mpch <- coxph(yrlg.tob ~ Mass + PC1 + HelpNum, data = flg.sub2)
summary(mpch)
anova(mpch)

mmpch <- coxph(yrlg.tob ~ Mass + PC1 + momexp + HelpNum, data = flg.sub2)
summary(mmpch)

mmphint <- coxph(yrlg.tob ~ Mass + PC1 + HelpNum + PC1*HelpNum, data = flg.sub2)
summary(mmphint)

anova(mpc1,mmphint, mmpc, mmpch)



##########################################################################
# Regression with each terr variable, need to use pcs as predictors 
#scrb.count unttransformed
cx.scb <- coxph(yrlg.mob ~ scrb.count, data = yrlg.df1)
summary(cx.scb)
anova(cx.scb)
sim.oak1 <- coxsimLinear(cx.scb, b = "scrb.count", qi ="Relative Hazard",
                         Xj = c(0,11547))
simGG(sim.oak1, xlab = "Amount of Oak Scrub", alpha = 0.3)

#Standardized
cx.oak <- coxph(yrlg.mob ~ stdscr, data = yrlg.df1)
summary(cx.oak)
anova(cx.oak)
plot(cox.zph(cx.oak))
sim.oak <- coxsimLinear(cx.oak, b = "stdscr", qi = "Relative Hazard",
                        Xj = c(0, 3.4))
simGG(sim.oak, xlab = "Standardized Amount of Oak Scrub", alpha = 0.3)

#Amount of cells in 2-9 year fire patches, standardized 
cx.tsf <- coxph(yrlg.mob ~ stdtsf, data = yrlg.df1)
summary(cx.tsf)
anova(cx.tsf)
plot(cox.zph(cx.tsf))

#Terrsize
cx.size <- coxph(yrlg.mob ~ stdsize, data = yrlg.df1)
summary(cx.size)
anova(cx.size)
plot(cox.zph(cx.size))

anova(cx.oak, cx.tsf, cx.size)

#All three
cx.terr <- coxph(yrlg.mob ~ stdsize + stdscr + stdtsf, data=yrlg.df1)
summary(cx.terr)
anova(cx.terr)


## 3 21 2017 this won't matter because I'm using the PCs as regression
#coefficients 
#Terr size and helpers - binary 
sz.hlp <- coxph(yrlg.mob ~ stdsize + Help + stdsize*Help, data = yrlg.df1)
summary(sz.hlp)
anova(sz.hlp)

#Helpers as count, not binary
sz.hlp2 <- coxph(yrlg.mob ~ stdsize + HelpNum + stdsize*HelpNum, data = yrlg.df1)
summary(sz.hlp2)
anova(sz.hlp2)

scr.hlp <- coxph(yrlg.mob ~ stdscr + Help + stdscr*Help, data = yrlg.df1)
summary(scr.hlp)
anova(scr.hlp)

anova(sz.hlp2, scr.hlp)
anova(cx.size, sz.hlp2)


#######################################################################
## Mixed effects cox models

yrlg.dfA <- subset(yrlg.df1, yrlg.df1$USFWBand %in% new.df$USFWBand)
mob.A <- Surv(yrlg.dfA$Days, yrlg.dfA$Censor, type = c('right'))

yrlg.dfA$HelpNum <- as.numeric(yrlg.dfA$HelpNum)
## 5 31 2017 
yr.hlp <- coxme(mob.A ~ Mass + HelpNum + (HelpNum|Year),data=yrlg.dfA)

#ID as random effect, no fixed effect - KEEP 3/21/2017
kin1 <- coxme(mob.A ~ (1|USFWBand), data = yrlg.dfA, 
              varlist = coxmeMlist(kins))
summary(kin1)
kin2 <- coxme(mob.A ~ (1|USFWBand) + (1|Year), data = yrlg.dfA, 
              varlist = coxmeMlist(kins))
summary(kin2)
anova(kin1, kin2)
#Add mass as fixed effect
#kin.mass <- coxme(mob.A ~ Mass + (1|USFWBand), data = yrlg.dfA,
                 #varlist = coxmeMlist(kins))
flg.sub1 <- flg.sub %>% drop_na(24,25,26,27,28,29,30,31,32)
yrlg.dfA <- yrlg.dfA %>% drop_na(22)
mob.sub <- Surv(yrlg.dfA$Years, yrlg.dfA$Censor, type =c('right'))
kin.massyr <- coxme(mob.sub ~ Mass + (1|USFWBand) + (1|Year), data = yrlg.dfA, 
                    varlist = coxmeMlist(kins))
summary(kin.massyr)
anova(kin2, kin.massyr)
anova(kin1, kin2, kin.massyr)
#anova(kin.mass)
## Add mom age
kin.mage <- coxme(mob.A ~ Mass + momage + (1|USFWBand) + (1|Year), 
                  data = yrlg.dfA, varlist = coxmeMlist(kins))
summary(kin.mage)
anova(kin.massyr, kin.mage)

kin.mexp <- coxme(mob.sub ~ Mass + momexp + (1|USFWBand) + (1|Year),
                  data = yrlg.dfA, varlist = coxmeMlist(kins))
summary(kin.mexp)

#### Decent fitting model 4 7 2017
kin.mah <- coxme(mob.A ~ Mass + momage + HelpNum + (1|USFWBand) +
                   (1|Year), data = yrlg.dfA, varlist = coxmeMlist(kins))
summary(kin.mah)
anova(kin.mage, kin.mah)
#Helpers add ever so slightly better fit, just below the general rule
#for interpreting these types of statistics 

#### Decent fitting model 4 7 2017
kin.e2 <- coxme(mob.sub ~ Mass + momexp + HelpNum + (1|USFWBand) +
                  (1|Year), data = yrlg.dfA, varlist = coxmeMlist(kins))
summary(kin.e2)
anova(kin.e2)

anova(kin.massyr, kin.mexp, kin.e2)
#anova(kin.e2)
anova(kin.mexp, kin.e2)

anova(kin.mexp, kin.e2)

#not run 4 7 2017
#not convinced interaction adds anything 4 7 2017
kin.int <- coxme(mob.A ~ Mass + momage + HelpNum + momage:HelpNum +
          (1|USFWBand)  + (1|Year), data = yrlg.dfA, 
                 varlist = coxmeMlist(kins))
summary(kin.int)
#no support for interaction 

###############################################################################
###### Want one or the other in these models, colinearity here 4 7 2017
## There are more or less equivalent, use experience I think 
kin.maep <- coxme(yrlg.mob ~ Mass + momage + momexp + momage*momexp+
          (1|USFWBand) + (1|Year), data = yrlg.df1, varlist = coxmeMlist(kins))
summary(kin.maep)


kin.age.mom <- coxme(yrlg.mob ~ Mass + momage + momexp + HelpNum +
                (1|USFWBand) + (1|Year), data = yrlg.df1,
                varlist = coxmeMlist(kins))

kin.mass.exp <- coxme(yrlg.mob ~ Mass + momage + momexp + HelpNum +
                  Mass:momexp + (1|USFWBand) + (1|Year), data = yrlg.df1, 
                  varlist = coxmeMlist(kins))
anova(kin.mah, kin.maep, kin.age.mom, kin.mass.exp)

#coxme(yrlg.mob ~ Mass + momexp + HelpNum)

###########################################################################
#kinship, year, mass and principal component 1 (terr size)
#remember that this is a smaller data set (N=1529) and cannot be compared 
#to models with the larger data set (N = 2323)

#### April 8 2017
kin.t.mass <- coxme(yrlg.tob ~ Mass + (1|USFWBand) + (1|Year),
                    data = flg.sub2, varlist = coxmeMlist(kins))
summary(kin.t.mass)

kin.terr <- coxme(yrlg.tob ~ Mass + PC1 + (1|USFWBand) + (1|Year), 
                  data = flg.sub2, varlist = coxmeMlist(kins))
summary(kin.terr)

kin.terrh <- coxme(yrlg.tob ~ Mass + PC1 + HelpNum + (1|USFWBand) + (1|Year),
                   data = flg.sub2, varlist = coxmeMlist(kins))
summary(kin.terrh)
anova(kin.t.mass, kin.terr, kin.terrh)

kin.mom.t <- coxme(yrlg.tob ~ Mass + PC1 + momage + (1|USFWBand) + (1|Year),
                   data = flg.sub2, varlist= coxmeMlist(kins))
summary(kin.mom.t)
anova(kin.t.mass,kin.terr, kin.mom.t)

kin.mom.exp <- coxme(yrlg.tob ~ Mass + PC1 + momexp +
                       (1|USFWBand) + (1|Year), 
                     data =flg.sub2, varlist=coxmeMlist(kins))
summary(kin.mom.exp)
anova(kin.t.mass, kin.terr, kin.mom.exp)

kin.mom.th <- coxme(yrlg.tob ~ Mass + PC1 + momage + HelpNum + 
                      (1|USFWBand) + (1|Year), data= flg.sub2, 
                    varlist= coxmeMlist(kins))
summary(kin.mom.th)

kin.mom.eh <- coxme(yrlg.tob ~ Mass + PC1 + momexp + HelpNum + 
                      (1|USFWBand) + (1|Year), data = flg.sub2, 
                    varlist = coxmeMlist(kins))
summary(kin.mom.eh)

anova(kin.terr, kin.terrh, kin.mom.exp, kin.mom.th, kin.mom.eh)
anova(kin.terr, kin.mom.exp)
anova(kin.terr, kin.mom.exp, kin.mom.eh, int.pc)

library(corrplot)
covars2 <- flg.sub2[,c(19,22,23,33)]
colnames(covars2) <- c("HelpNum","MomExp","MomAge","PC1")
corrs2 <- cor(covars2, use="complete.obs")
corrplot(corrs2, method="pie", type= "lower")

#Interaction of helper number and PC1
int.pc <- coxme(yrlg.tob ~ Mass + PC1 + momexp + HelpNum + 
      HelpNum*PC1 + (1|USFWBand) + (1|Year), data = flg.sub2,
      varlist = coxmeMlist(kins))

anova(kin.terr,kin.terrh,kin.mom.t)
anova(kin.terr, kin.mom.t)
anova(kin.terr, kin.terrh)
anova(kin.terr, kin.mom.th)

anova(kin.terr, kin.mom.exp, kin.mom.eh)
####################################################################
# 4 8 2017 old 
summary(kin.mom.th)
anova(kin.terrh, kin.mom.t, kin.mom.th)

kin.h <- coxme(yrlg.tob ~ Mass + HelpNum + (1|USFWBand) + (1|Year),
               data = flg.sub2, varlist = coxmeMlist(kins))
summary(kin.terr)
summary(kin.h)
summary(kin.maep)
#asses fit of adding mass 
anova(kin1, kin.mass, test = "chisq")
#model with mass better



april7 = new.env(hash = TRUE, parent = .GlobalEnv)
saveRDS(april7, file = "ErinApril7.RData")

#Add helper (0,1)
kin.help <- coxme(yrlg.mob ~ Mass + Help + (1|USFWBand), data = yrlg.df1,
                  varlist = coxmeMlist(kins))
summary(kin.help)

#compare mods 1,2,3 - random effect only, add mass, add mass + helper
anova(kin1, kin.mass, kin.help, test = "chisq")

#Helper number with mass as fixed effect
kin.hlnum <- coxme(yrlg.mob ~ Mass + HelpNum + (1|USFWBand), data= yrlg.df1, 
                   varlist = coxmeMlist(kins))
summary(kin.hlnum)

anova(kin1, kin.mass, kin.hlnum, test = "chisq")

kin.year <- coxme(yrlg.mob ~ Mass + HelpNum + (1|USFWBand) + (1|Year), 
                  data = yrlg.df1, varlist = coxmeMlist(kins))

#Mixed effects or frailty Cox models are better fits than Cox PH models
#In the Cox PH models, the model with the interaction between helpers and
#and mass had the best support of the estimated models 
kin.help1 <- coxme(yrlg.mob ~ Mass + HelpNum + Mass:HelpNum + (1|USFWBand),
                data = yrlg.df1, varlist = coxmeMlist(kins))
summary(kin.help1)

kin.my <- coxme(yrlg.mob ~ Mass + (1|USFWBand) + (1|Year), 
                data = yrlg.df1, varlist = coxmeMlist(kins))
summary(kin.my)

#same as above but with another random effect of year
#Does it make sense to "nest" any of the random effects? 
kin.help2 <- coxme(yrlg.mob ~ Mass + HelpNum + Mass:HelpNum + (1|USFWBand)
                   + (1|Year), data = yrlg.df1, varlist = coxmeMlist(kins))
summary(kin.help2)

anova(kin1, kin.mass, kin.my, kin.hlnum, kin.help2, test = "chisq")
anova(kin.my, kin.year, kin.help2)

#kinship, year, mass and principal component 1 (terr size)
#remember that this is a smaller data set (N=1529) and cannot be compared 
#to models with the larger data set (N = 2323)
kin.terr <- coxme(yrlg.tob ~ Mass + PC1 + (1|USFWBand) + (1|Year), 
                  data = flg.sub2, varlist = coxmeMlist(kins))
kin.h <- coxme(yrlg.tob ~ Mass + HelpNum + (1|USFWBand) + (1|Year),
               data = flg.sub2, varlist = coxmeMlist(kins))
kin.terrh <- coxme(yrlg.tob ~ Mass + PC1 + HelpNum + (1|USFWBand) + (1|Year),
              data = flg.sub2, varlist = coxmeMlist(kins))
kin.mom.t <- coxme(yrlg.tob ~ Mass + PC1 + momage + (1|USFWBand) + (1|Year),
                    data = flg.sub2, varlist= coxmeMlist(kins))
kin.mom.th <- coxme(yrlg.tob ~ Mass + PC1 + momage + HelpNum + 
                  (1|USFWBand) + (1|Year), data= flg.sub2, 
                  varlist= coxmeMlist(kins))

summary(kin.terr)
summary(kin.h)
summary(kin.terrh)
summary(kin.mom.t)
summary(kin.mom.th)

anova(kin.mom.t, kin.mom.th)
anova(kin.terr, kin.mom.t)

anova(kin.terr, kin.h, kin.terrh, kin.mom.t, kin.mom.th, test = "chisq")

anova(kin.terr, kin.h, kin.terrh, test = "chisq")
anova(kin.terr, kin.terrh, test = "chisq")

kin.terr2 <- coxme(yrlg.tob ~ Mass + stdsize + (1|USFWBand) + (1|Year), 
                  data = flg.sub1, varlist = coxmeMlist(kins))
summary(kin.terr2)

kin.terr3 <- coxme(yrlg.tob ~ Mass + PC1 + HelpNum + PC1*HelpNum +
  (1|USFWBand) + (1|Year), data = flg.sub1, varlist = coxmeMlist(kins))
summary(kin.terr3)

anova(kin.terr, kin.terr2, kin.terr3)
anova(kin.terr, kin.terrh, kin.terr3)

mod.compare <- coxph(yrlg.tob ~ Mass + PC1 + HelpNum, data = flg.sub1)
anova(mod.compare, kin.terrh)


#### Using family id instead of lineage (lineage is the better method)
test <- coxme(suv.ob ~ Mass + Help + (1|Year) + (1|famid), data =yrlg.df1)
summary(test)
test2 <- coxme(suv.ob ~ Mass + HelpNum + (1|Year) + (1|famid), data= yrlg.df1)
summary(test2)

famid.fit <- coxme(suv.ob ~ Mass + (1|famid), data = yrlg.df1)
summary(famid.fit)

me2 <- coxme(suv.ob ~ Mass + (1|Year), data = new.df)
summary(me2)

anova(kin1, kin.fit, kin.help, test, kin.hlnum, kin.help2, kin.help1)
anova(kin.help, kin.hlnum)

anova(kin.help1, kin.help2)

anova(test, kin.help2)
#Year only
me1<- coxme(suv.ob ~ (1|Year), data = yrlg.df1)
summary(me1)


#compare model with just mass to model with mass and year
anova(mass.full, me2)
anova(me1, me2)

#mass and helpers as fixed effects
me3 <- coxme(suv.ob ~ Mass + Help + (1|Year), data = yrlg.df1)
summary(me3)

anova(me1, me2, me3)

me4 <- coxme(suv.ob ~ Mass + Help + Mass*Help + (1|Year), 
            data = yrlg.df1)
summary(me4)

anova(me1, me2, me3, me4)

me5 <- coxme(suv.ob ~ Mass + HelpNum + (1|Year), data = yrlg.df1)
summary(me5)

me6 <- coxme(suv.ob ~ Mass + HelpNum + Mass:HelpNum + (1|Year),
             data = yrlg.df1)
summary(me6)
summary(kin.help2)

anova(me6, kin.help2)

anova(me5,me6)
anova(me3,me4,me4,me6)

#Remember that this random effect is nested within year
#So the var estimate will be higher than year alone
me7 <- coxme(suv.ob ~ Mass + HelpNum + Mass:HelpNum + (1|NatalNest),
             data= yrlg.df1)
anova(me7, kin.help2)
summary(me7)
anova(me6, me7)
#as expected me7 is a slightly better fit by reduction of deviance 
#but not by much and I suspect pedigree will be ever more 

anova(me3, me5, kin.fit, kin.help, kin.help2)
summary(kin.help2)

#Model 5 - mixed effect model with year as random/frailty term
mm.year <- coxme(yrlg.ob ~ Mass + (1|Year), data = yrlg.df1)
mm.year

anova(cox.mass, mm.year)
#me model better fit by LRT chi sq and reduction in deviance 

#Model 6 - mixed effect model with natal nest as random term
#Realize that year is included in this variable 
mm.nest <- coxme(yrlg.ob ~ Mass + (1|NatalNest), data = yrlg.df1)
mm.nest

anova(mm.year, mm.nest)
#Nest is a slightly better fit, not suprising since natal nest is nested 
#within year 

## Model with oak scrub, mass, helpers, interactions, random term
mxe1 <- coxme(yrlg.ob ~ stdscr + (1|Year), data= yrlg.df1)
mxe2 <- coxme(yrlg.ob ~ stdsize + (1|Year), data = yrlg.df1)
#compare this to ph model
anova(cx.oak, mxe1)
#much better fit with random term

mxe3 <- coxme(yrlg.ob ~ Mass + stdscr + (1|Year), data = yrlg.df1)
summary(mxe3)
anova(mxe1, mxe3)
#model with mass and terr info better fit 

oak2 <- coxph(yrlg.ob ~ Mass + stdscr + Mass:stdscr, data = yrlg.df1)
summary(oak2)
anova(mxe3, oak2)
#way worse of a model than simpler one with mass, terr info, and year, no
#interaction is better model 

#compare models using oak scrub and size with year term
anova(mxe1, mxe2)
#similar fit 

#full model with helper and mass plus interaction, compare to model with terr
#size only and year as random effect
anova(sz.hlp2, mxe2)
#Simpler model better!

oak.mx <- coxme(yrlg.ob ~ stdsize + Mass + HelpNum + (1|Year), 
                data = yrlg.df1)
oak.mx2 <- coxme(yrlg.ob ~ stdsize + Mass + (1|Year), data =yrlg.df1)
summary(oak.mx)
summary(oak.mx2)
anova(mxe2, oak.mx)
anova(oak.mx2, oak.mx)

# 3 5 2017
mods <- coxme(yrlg.ob ~ Mass + stdscr + (1|USFWBand), 
            data = yrlg.df1 ,varlist = coxmeMlist(kins))
summary(mods)
mods2 <- coxme(yrlg.ob ~ Mass + HelpNum + stdscr + (1|USFWBand), 
               data = yrlg.df1, varlist = coxmeMlist(kins))
summary(mods2)
mods3 <- coxme(yrlg.ob ~ Mass + HelpNum + stdscr + HelpNum:stdscr +
                 (1|USFWBand), data = yrlg.df1, varlist = coxmeMlist(kins))
summary(mods3)

anova(mods, mods2, mods3)

#problem
mxe <- coxme(yrlg.ob ~ Mass + stdscr + HelpNum + Mass:HelpNum, 
            + stdscr:HelpNum, Mass:stdscr + (1|Year), 
            data = yrlg.df1)


#Stuff below is test code and the building blocks to this script 
#Not needed in many body but I didn't want to just delete it 

################################################################


##############################################################################
# 3 4 2017 I moved the code below to the beginning before any Cox models 
#Pedigree using kinship2
#Estimating covars between individuals

#example
data("sample.ped")

pedAll <- pedigree(id=sample.ped$id, dadid=sample.ped$father, momid=sample.ped$mother, sex=sample.ped$sex, famid=sample.ped$ped)
print(pedAll)
#subsetting by family, family id = 1
ped1basic <- pedAll['1']
#family id = 2
ped2basic <- pedAll['2']
#returns description of ped for family 2, number of subjects
print(ped2basic)
#plot pedigree
plot(ped2basic)
print(ped1basic)
plot(ped1basic)

#calculate kinship matrix
kin2 <- kinship(ped2basic)
kin2

names(yrlg.df1)
str(yrlg.df1)

#upload pedigree
jayped <- read.csv("Demo_Pedigree_2017.csv")
names(jayped)
jayped$SexIND[which(jayped$Sex=="F")] <- 2
jayped$SexIND[which(jayped$Sex=="M")] <- 1
jayped$SexIND[which(jayped$Sex=="")] <- 3
jayped[,11:12] <- NULL


HYped <- merge(yrlg.df1, jayped, intersect = c("USFWBand", "NatalNest"))
names(HYped)
str(HYped)

HYped <- HYped[,c(2,1,6,21,23,24,25,26,16)]


names(HYped)
names(pars)

#Use this if you have added the territory data
#HYped <- HYped[,c(1,2,6,30,32,33,34,35,16)]
str(HYped)
HYped$Year <- as.character(HYped$Year)
HYped$Year <- as.numeric(HYped$Year)


#So it appears that everyone in the ID list must appear in parent list
#Birds who have parents who are not in the list get NA for parents 

length(unique(HYped$FUSFWBand))  #398 moms, there is one with a missing band#
length(unique(HYped$FBreeder))      
length(unique(HYped$MUSFWBand))   #383 ids, 2 missing band #
length(unique(HYped$MBreeder))


parents <- data.frame(HYped$FBreeder, HYped$FUSFWBand, HYped$MBreeder, HYped$MUSFWBand)
colnames(parents) <- c("FID", "Fband", "MID", "Mband")
library(plyr)
mom.counts <- count(parents, "FID")
dad.counts <- count(parents, "MID")

short.ped <- data.frame(HYped$JayID, HYped$FBreeder, HYped$FUSFWBand, 
                        HYped$MBreeder, HYped$MUSFWBand)
colnames(short.ped) <-c("JayID", "Mom", "MomUS", "Dad", "DadUS")
moms <- data.frame(short.ped$Mom, short.ped$MomUS)
colnames(moms) <- c("Mom", "MomUS")
dads <- data.frame(short.ped$Dad, short.ped$DadUS)
colnames(dads) <- c("Dad", "DadUS")
kids <- data.frame(short.ped$JayID)
colnames(kids) <- "kids"

#moms not in the focal bird list (fledglings)
notmoms <- subset(short.ped, !(short.ped$Mom %in% short.ped$JayID), 
               select=JayID:MomUS)
length(unique(notmoms$Mom)) # 271 mamas not in the fledge list
notmoms[,1] <- NULL
notmoms <- data.frame((notmoms[!duplicated(notmoms$Mom),]))
#Moms not in focal ID list (not found in HY bird individuals)
#Color band combo
colnames(notmoms) <- c("Mom", "MomUS")

notdads <- subset(short.ped, !(short.ped$Dad %in% short.ped$JayID))
notdads[,2] <- NULL
notdads[,1] <- NULL
notdads[,1] <- NULL
notdads <- data.frame((notdads[!duplicated(notdads$Dad),]))
#Dads not in focal ID list (not found in HY bird individuals)
colnames(notdads) <- c("Dad", "DadUS")

#Any birds with parent(s) in this list get 0 for parent ID
notdads
notmoms

length(unique(notdads$Dad)) # 195 dads not in the fledge list 
count.notf <- count(notmoms, "Mom")
count.notm <- count(notdads, "Dad")

yesmoms <- subset(short.ped, (short.ped$Mom %in% short.ped$JayID))
yesmoms <- yesmoms[!duplicated(yesmoms$Mom),]
yesmoms[1] <- NULL
yesmoms[3] <- NULL
yesmoms[3] <- NULL
yesdads <- subset(short.ped, (short.ped$Dad %in% short.ped$JayID))
yesdads <- yesdads[!duplicated(yesdads$Dad),]
yesdads[1] <- NULL
yesdads$Mom <- NULL
yesdads$MomUS <- NULL
#Founders need to be any bird that has a mom or dad in the data frame above

#Mom and Dad ID's need to be changed to NA for any bird that has parents
#on list from "count.notf" and "count.notm" - those birds are not in fledge list
#There are some birds that were fledglings then became parents later 

#2 14 2017 
#Add in birds that are not on focal id list but have kids on it
names(HYped)
notmoms["SexIND"] <- 2
notdads["SexIND"] <- 1
colnames(notmoms) <- c("JayID", "USFWBand", "SexIND")
colnames(notdads) <- c("JayID", "USFWBand", "SexIND")

#founders
pars <- rbind(notmoms, notdads)
pars["Year"] <- ""
pars["FBreeder"] <- NA
pars["FUSFWBand"] <- NA
pars["MBreeder"] <- NA
pars["MUSFWBand"] <- NA
pars["Censor"] <- ""

pars <- pars[,c(2,1,4,3,5,6,7,8,9)]

newPed <- rbind(HYped, pars)
newPed$FBreeder <- as.character(newPed$FBreeder)
newPed$FUSFWBand <- as.character(newPed$FUSFWBand)
newPed$MBreeder <- as.character(newPed$MBreeder)
newPed$MUSFWBand <- as.character(newPed$MUSFWBand)

newPed[is.na(newPed)] <- 0

newPed$FUSFWBand[which(newPed$FBreeder == "6893")] <- "1"
newPed$FUSFWBand[which(newPed$FBreeder == "-LQ")] <- "2"
newPed$MUSFWBand[which(newPed$MBreeder == "RLC-U")] <- "3"
newPed$MUSFWBand[which(newPed$MBreeder == "-SW_")] <- "4"

##optional but need as.character
##removes dash in band #
#newPed$USFWBand <- gsub("-", "", newPed$USFWBand)
#newPed$USFWBand <- as.numeric(newPed$USFWBand)
#newPed$FUSFWBand <- gsub("-", "", newPed$FUSFWBand)
#newPed$MUSFWBand <- gsub("-", "", newPed$MUSFWBand)
#newPed$FUSFWBand <- as.numeric(newPed$FUSFWBand)
#newPed$MUSFWBand <- as.numeric(newPed$MUSFWBand)



newPed$USFWBand <- as.character(newPed$USFWBand)
newPed$FUSFWBand <- as.character(newPed$FUSFWBand)
newPed$MUSFWBand <- as.character(newPed$MUSFWBand)

newPed$USFWBand[which(newPed$JayID =="6893")] <- "1"
newPed$USFWBand[which(newPed$Jay=="-LQ")] <- "2"
newPed$USFWBand[which(newPed$JayID=="RLC-U")] <- "3"
newPed$USFWBand[which(newPed$JayID =="-SW_")] <- "4"

### 2 15 2017 have to add in a year for those birds I think
founders1 <- subset(newPed, newPed$FUSFWBand==0 & newPed$MUSFWBand==0)
# Run BreedersJuly script through line 98 for the breeder data set to get year
str(brd)
str(founders1)
#founders1 - parent info for founders 
founders1$Year <- NULL
founders1$Censor <- NULL
found.year <- brd[,1:6]
colnames(found.year)[1] <- "JayID"
found.year[,4:5] <- NULL
colnames(found.year)[2] <- "USFWBand"
found.year["SexIND"] <- 0
found.year$SexIND[which(found.year$Sex == "F")] <- 2
found.year$SexIND[which(found.year$Sex == "M")] <- 1

#JayID and year of breeding, I actually don't think this step is necessary
#I added it thinking that this was the problem with the kinship matrix 
found.year$USFWBand <- NULL
found.year$SexIND <- NULL
found.year$Sex <- NULL

#founders with year of first breeding 
founders <- merge(founders1, found.year, by = "JayID")
founders["yearind"] <- founders$Year - 1
founders$Year <- NULL
colnames(founders)[8] <- "Year"

#Who is in founders1 and not in founders? 
notin <- subset(founders1, !(founders1$JayID %in% founders$JayID))
notin <- notin[,c(2,1,3,4,5,6,7)]
#need the year
notin["Year"] <- 1978


birds <- rbind(founders, notin)
birds <- birds[,c(2,1,3,4,5,6,7,8)]
birds["Censor"] <- 0

### Have to add pedigree info to the original data set 
### yrlg.df1 + birds
birds2 <- subset(newPed, !newPed$FUSFWBand == 0)
birds2$Year <- NULL
birds2$Censor <- NULL
names(birds2)


#re-arrange dataframe with covariates to merge
yrlg.df1 <- yrlg.df1[,c(2,3,4,1,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)]
yrlg.df1$USFWBand <- as.character(yrlg.df1$USFWBand)

#called please as in "PLEASE WORK!" 
#dataframe with parent info and all covariates
please <- merge(yrlg.df1, birds2, by = c("USFWBand","JayID"))

birds["NatalNest"] <- NA
birds["TerrYr"] <- NA
birds["Terr"] <- NA
birds["HatchDate"] <- NA
birds["FldgDate"] <- NA
birds["LastObsDate"] <- NA
birds["Mass"] <- NA
birds["MeasDate"] <- NA
birds["HatchNum"] <-  NA
birds["FldgNum"] <- NA
birds["Days"] <- NA
birds["Years"] <- NA
birds["Day11"] <- NA
birds["Help"] <-  NA
birds["HelpNum"] <- NA
birds["Censor"] <- NA

birds <- birds[,c(1,2,10,11,12,8,13,14,15,17,18,19,20,21,22,9,23,
                  24,16,3,4,5,6,7)]
                
new.df <- rbind(please, birds)
str(new.df)

## 2 20 2017 need to add family id
ids <- makefamid(new.df$USFWBand, new.df$MUSFWBand, new.df$FUSFWBand)
famid <- cbind(ids)
colnames(famid) <- "famid"
new.df <- cbind(new.df, famid)
str(new.df)

new.df$USFWBand <- gsub("-", "", new.df$USFWBand)
new.df$USFWBand <- as.numeric(new.df$USFWBand)
new.df$USFWBand <- as.integer(new.df$USFWBand)
new.df$FUSFWBand <- gsub("-", "", new.df$FUSFWBand)
new.df$MUSFWBand <- gsub("-", "", new.df$MUSFWBand)
new.df$FUSFWBand <- as.numeric(new.df$FUSFWBand)
new.df$MUSFWBand <- as.numeric(new.df$MUSFWBand)
new.df$FUSFWBand <- as.integer(new.df$FUSFWBand)
new.df$MUSFWBand <- as.integer(new.df$MUSFWBand)

#
#new.df$USFWBand <- as.numeric(new.df$USFWBand)
#new.df$FUSFWBand <- as.numeric(new.df$FUSFWBand)
#new.df$MUSFWBand <- as.numeric(new.df$MUSFWBand)
str(new.df)
#new.df$gid <- paste(new.df$famid, new.df$JayID, sep="/")
#still duplicated 

#jay.ped <- with(new.df, pedigree(id = JayID, dadid=MBreeder,
            #momid = FBreeder, sex = SexIND, famid = famid, missid = 0))

jay.ped <- with(new.df, pedigree(id = USFWBand, dadid=MUSFWBand, 
              momid = FUSFWBand, sex = SexIND, famid =famid, missid = 0))

ped1 <- jay.ped[1]
#plot(ped1)
ped2 <- jay.ped[2]
plot(ped2)
ped5 <- jay.ped[5]
plot(ped5)
ped8 <- jay.ped[8]
plot(ped8)
ped45 <- jay.ped[45]
plot(ped45)

kins <- kinship(jay.ped)

#Doesn't run, or I guess it does but just sits there until it crashes
#t.fit <- coxme(Surv(new.df$Years, new.df$Censor) ~ Mass + (1|USFWBand), 
               #data = new.df, varlist = coxmeMlist(kins))
suv.ob <- Surv(new.df$Years, new.df$Censor, type=c('right'))

famid.fit <- coxme(suv.ob ~ Mass + (1|famid), data = new.df)
summary(famid.fit)
#exp(0.27) = 1.31
#According to T. Therneau coxme pdf, 


fit1 <- coxme(suv.ob ~ (1|new.df$USFWBand), 
              data = new.df, varlist = coxmeMlist(kins))

system.time(fit1 <-coxme(suv.ob ~ (1|new.df$USFWBand), 
                         data = new.df, varlist = coxmeMlist(kins)))
fit3 <- coxme(suv.ob ~ Mass + (1|USFWBand), 
              data = new.df, varlist = coxmeMlist(kins))


data(minnbreast)
mped <- with(minnbreast, pedigree(id, fatherid, motherid, sex, 
              affected=cancer, famid=famid, status=proband))
minnfemale <- minnbreast[minnbreast$sex == 'F' & !is.na(minnbreast$sex),]
str(minnfemale)
kmat <- kinship(mped)
fit2 <- coxme(Surv(endage, cancer) ~ I(parity>0) + (1|id),
              data=minnfemale, varlist = coxmeMlist(2*kmat, rescale=F), 
              subset=(proband==0))
system.time(fit2<-coxme(Surv(endage, cancer) ~ I(parity>0) + (1|id),
                        data=minnfemale, varlist = coxmeMlist(2*kmat, rescale=F), 
                        subset=(proband==0)))

                
####################################################################
newPed1 <- rbind(birds2, birds)
str(newPed1)
newPed1$Year <- as.factor(newPed1$Year)

pedtest1 <- with(yrlg.df1, pedigree(id = newPed1$USFWBand, sex = newPed1$SexIND, 
              dadid = newPed1$MUSFWBand, momid = newPed1$FUSFWBand, 
              missid = 0))

birdkin <- kinship(pedtest1, dadid = newPed$MUSFWBand, 
            momid = newPed$FUSFWBand)
test.fit <- coxme(yrlg.obfll ~ Mass + (1|JayID), data= yrlg.df1,
                  varlist=coxmeMlist(birdkin))

### Breast cancer study example 
data("minnbreast")
mped <- with(minnbreast, pedigree(id, fatherid, motherid, sex,
          affected=cancer, famid=famid, status=proband))

class(mped)
minnfemale <- minnbreast[minnbreast$sex == 'F' & !is.na(minnbreast$sex),]
fit0 <- coxph(Surv(endage, cancer) ~ I(parity > 0), minnfemale,
              subset=(proband==0))
summary(fit0)

kmat <- kinship(mped)


#######################################################################
#Somehow have to code parents not in ID list to 0 or NA
#birds with both parents not on list 
test <- subset(HYped, (!HYped$FBreeder %in% HYped$JayID) & (!HYped$MBreeder %in% HYped$JayID))

test$FBreeder <- NA
test$FUSFWBand <- NA
test$MBreeder <- NA
test$MUSFWBand <- NA
#birds with both parents on the list
test4 <- subset(HYped, (HYped$FBreeder %in% HYped$JayID) & (HYped$MBreeder %in% HYped$JayID))

#bind df's with birds with both parents either on list or not
checkcheck <- rbind(test, test4)
#subset 

#birds with one parent in list, one not
birdies <- subset(HYped, (!HYped$JayID %in% checkcheck$JayID))
birdies$FUSFWBand <- 0
birdies$MUSFWBand <- 0

#Now need birds that are in HYped but not in test or test4
cb <- rbind(checkcheck, birdies)
#whoelse <- subset(HYped, (HYped$JayID %in% cb$JayID))
whoelse1 <- subset(HYped, !(HYped$USFWBand %in% cb$USFWBand))
whoelse1$FUSFWBand <- 0
whoelse1$MUSFWBand <- 0

#cb$FBreeder[which(cb$FBreeder %in% notmoms$Mom)] <- NA
#cb$MBreeder[which(cb$MBreeder %in% notdads$Dad)] <- NA
#cb$FUSFWBand[which(cb$FUSFWBand %in% notmoms$MomUS)] <- NA
#cb$MUSFWBand[which(cb$MUSFWBand %in% notdads$DadUS)] <- NA

#whoelse1$FBreeder[which(whoelse1$FBreeder %in% notmoms$Mom)] <- NA
#whoelse1$MBreeder[which(whoelse1$MBreeder %in% notdads$Dad)] <- NA
#whoelse1$FUSFWBand[which(whoelse1$FUSFWBand %in% notmoms$MomUS)] <- NA
#whoelse1$MUSFWBand[which(whoelse1$MUSFWBand %in% notdads$DadUS)] <- NA

#This has all birds, need to put blanks for birds with 1 parent missing 
real.ped <- rbind(cb, whoelse1)
real.ped[is.na(real.ped)] <- 0

real.ped$USFWBand <- gsub("-", "", real.ped$USFWBand)
real.ped$USFWBand <- as.numeric(real.ped$USFWBand)

real.ped$FUSFWBand <- gsub("-", "", real.ped$FUSFWBand)
real.ped$MUSFWBand <- gsub("-", "", real.ped$MUSFWBand)

real.ped$FUSFWBand <- as.numeric(real.ped$FUSFWBand)
real.ped$MUSFWBand <- as.numeric(real.ped$MUSFWBand)
##################################################################################
##########################################################################
pedtest <- pedigree(id = real.ped$USFWBand, sex = real.ped$SexIND, dadid = real.ped$MUSFWBand,
          momid = real.ped$FUSFWBand, status = real.ped$Censor, missid = 0)
pedtest
plot(pedtest, cex = 0.1)

#I guess it does work
birdkin <- kinship(pedtest, dadid = real.ped$MUSFWBand, momid = real.ped$FUSFWBand)
test.fit <- coxme(yrlg.obfll ~ Mass + (1|Year), data= yrlg.df1,
                  varlist=coxmeMlist(birdkin))

yrlg.df1$Year <- as.character(yrlg.df1$Year)
yrlg.df1$Year <- as.numeric(yrlg.df1$Year)
new <- subset(yrlg.df1, yrlg.df1$Year >= 2005)
new.ob <- Surv(new$Years, new$Censor, type =c('right'))

test.fit2 <- coxme(new.ob ~ Mass + (1|NatalNest), data = new,
                  varlist=coxmeMlist(birdkin))

#None of this matters because I just kept all birds when I added the territory
#info so some birds just have blanks 
#birds2 is the parental info for fledglings 
names(birds2)

#Works now that band # all characters, I didn't change to integer
ter.df <- merge(yrlg.df1, birds2, by = "USFWBand")
#ter.df[,29:46] <- NULL
names(ter.df)
str(ter.df)
#colnames(ter.df)<- c("USFWBand", "JayID", "NatalNest", "TerrYr", "Terr", 
#"Year","HatchDate", "FldgDate","LastObsDate", "Mass", "MeasDate",
#"HatchNum", "FldgNum", "Days", "Years", "Censor", 
#"Day11", "Help", "HelpNum", "JayID", "NatalNest", "TerrYr", 
#"Terr", "Year", "HatchDate", "FldgDate", "LastObsDate", "Mass")

names(ter.df)
ter.df["JayID.y"] <- NULL
names(ter.df)
colnames(ter.df)[3] <- "JayID"

names(ter.df)
names(birds)

# 3 5 2017 I left them as characters 
##change band numbers to digits only 

#ter.df$USFWBand <- gsub("-", "", ter.df$USFWBand)
#ter.df$USFWBand <- as.integer(ter.df$USFWBand)
#ter.df$FUSFWBand <- gsub("-", "", ter.df$FUSFWBand)
#ter.df$MUSFWBand <- gsub("-", "", ter.df$MUSFWBand)
#ter.df$FUSFWBand <- as.integer(ter.df$FUSFWBand)
#ter.df$MUSFWBand <- as.integer(ter.df$MUSFWBand)


birds["scrb.count"] <- ""
birds["TerrSize"] <- ""
birds["tsf.count"] <- ""
birds["stdscr"] <- ""
birds["centscr"] <- ""
birds["stdtsf"] <- ""
birds["centtsf"] <- ""
birds["stdsize"] <- ""
birds["centsize"] <- ""

names(birds)
names(ter.df)

ter.df <- ter.df[,c(1,3,4,2,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
                    22,23,24,25,26,27,28,29,30,31,32,33)]

birds <- birds[,c(1,2,3,4,5,6,7,8,9,19,10,11,12,13,14,16,15,17,18,
                  25,26,27,28,29,30,31,32,33,20,21,22,23,24)]
names(ter.df)
names(birds)

test.df <- rbind(ter.df, birds)



#names(ter.tmp)
#colnames(ter.tmp)[2] <- "JayID"
#ter.tmp["JayID.y"] <- NULL
#birds <- birds[,c(1,2,3,4,5,6,7,8,9,11,12,13,14,15,17,16,18,19,20,21,22,
#23,24)]
#names(birds)
#names(ter.tmp)
#ter.tmp["MeasDate"] <- NULL

#ter.ped <- rbind(ter.tmp, birds)

#make family ids
ter.ids <- makefamid(ter.ped$USFWBand, ter.ped$MUSFWBand, ter.ped$FUSFWBand)
ter.fam <- cbind(ter.ids)
colnames(ter.fam) <- "famid"
ter.df <- cbind(ter.ped, ter.fam)

#have to add founders 

#make pedigree
ter.ped <- with(ter.df, pedigree(id = USFWBand, dadid=MUSFWBand, 
                                 momid = FUSFWBand, sex = SexIND, famid =famid, missid = 0))

#########################################################################