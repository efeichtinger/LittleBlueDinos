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
library(stringr)
library(simPH)
library(lubridate)
library(corrplot)
library(tidyverse)

#################################################################

# Block 1 - Input data, manipulate, make survival object, fix basic Cox Model

## Read in file of breeders - known and unknown age
## One record per individual

#Original file, had to revise from database to add territory name to get terryear
#brd1 <- read.csv("Breeders.csv")
brd <- read.csv("Breeders2.csv")

#fecundity data
jays.fec <- read.csv("QRY_TerrCh_FIYB_Aggr.csv")

#change NAs to zero 
jays.fec[is.na(jays.fec)] <- 0

##5 2 2017 - density, first file is just count of birds, need area 
jay.pop <- read.table("QRY_TC2_GetHelperCount.txt", header =TRUE)
jay.pop$Terr <- str_sub(jay.pop$TerrYr, start = 1, end = 4)

#sum #jays per year
jay.num <- ddply(jay.pop, .(Year), summarize, 
                 CountofJayID = sum(CountOfJayID))
terr.num <- ddply(jay.pop, .(Year), summarize, 
                  NumTerr = length(TerrYr))
jay.den <- cbind(jay.num, terr.num)
jay.den[,3] <- NULL
jay.den$Den <- jay.den$CountofJayID/jay.den$NumTerr
  
#group size
group.size <- read.csv("GroupsizeMarch.csv")
str(group.size)

#remove duplicate rows 
brd <- brd[!duplicated(brd),]

#brd1 <- brd1[!duplicated(brd1),]

brd$F_AccessID <- NULL
brd$Expr1 <- NULL

colnames(brd)[1] <- "JayID"
colnames(brd)[2] <- "USFWBand"
colnames(brd)[3] <- "Year"
colnames(brd)[4] <- "Fbreed"

#colnames(brd1)[1] <- "ID"
#colnames(brd1)[2] <- "Band"
#colnames(brd1)[3] <- "Year"
#colnames(brd1)[4] <- "Fbreed"


#brd <- brd[,c(1,2,8,7,6,3,4,5)] #corresponds to original df
#minimum age at first breeding 
colnames(brd)[6] <- "MinAgeFBr"
#colnames(brd1)[6] <- "MinAgeFBr"

#Set all blanks to NA
brd[brd==""] <- NA
#brd1[brd1==""] <- NA
#remove na
#brd <- na.omit(brd)
#brd1 <- na.omit(brd1)

brd$Fbreed <- as.Date(brd$Fbreed, format ="%m/%d/%Y")
brd$LastObsDate <- as.Date(brd$LastObsDate, format = "%m/%d/%Y")

#brd1$Fbreed <- as.Date(brd1$Fbreed, format ="%m/%d/%Y")
#brd1$LastObsDate <- as.Date(brd1$LastObsDate, format = "%m/%d/%Y")

#Years spent as a breeder, 365.25 to account for leap years 
brd["Yrs"] <- (brd$LastObsDate - brd$Fbreed)/365.25
brd$Yrs <- round(brd$Yrs, digits = 3)

#brd1["Yrs"] <- (brd1$LastObsDate - brd1$Fbreed)/365.25
#brd1$Yrs <- round(brd1$Yrs, digits = 1)

#Add censorship column
brd["Censor"] <- 1
brd$Censor[which(brd$LastObsDate =="2016-4-12")] <- 0

#brd1["Censor"] <- 1
#brd1$Censor[which(brd1$LastObsDate =="2016-4-12")] <- 0

brd <- subset(brd, brd$Yrs > 0)
#breeders who started for first time before 1981
pre1981 <- subset(brd, brd$Year <= 1980)
brd <- subset(brd, brd$Year >= 1981)
brd$FY <- as.factor(brd$Year)


brd$JayID <- as.character(brd$JayID)
brd$USFWBand <- as.character(brd$USFWBand)

brd$USFWBand[which(brd$JayID == "RLC-U")] <- "3"
brd$USFWBand[which(brd$JayID == "-SW_")] <- "4"
brd$USFWBand[which(brd$JayID == "6962")] <- "5"

brd <- brd[!duplicated(brd$JayID),]

#brd1 <- subset(brd1, brd1$Yrs > 0)
#brd1 <- subset(brd1, brd1$Year >= "1981")
#brd1$FY <- as.factor(brd1$Year)

#combine territory name with last two digits of year for "terryr" 
brd["TerrYr"] <- paste(brd$Terr, str_sub(brd$Year, start= -2),
                          sep ="")

brd["SexInd"] <- 0
brd$SexInd[which(brd$Sex == "M")] <- 1

#read in breeders by territory and year
#brdrty <- read.csv("QRY_Brdrs_TerrYr.csv")
#str(brdrty)


########################################################################

#plots and summary stats for male and female breeders 

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

## Mean years survived as breeder per cohort 
cohorts <- ddply(brd, .(FY), summarise, 
      mean.yrs = round(mean(Yrs), digits = 2), 
                 sd.yrs = round(sd(Yrs), digits = 2),
      var.yrs = round(var(Yrs), digits = 2))
cohorts$FY <- as.character(cohorts$FY)
cohorts$FY <- as.numeric(cohorts$FY)
cohorts$N <- c(8,21,19,23,16,19,22,29,23,22,31,19,35,25,31,42,33,19,32,13,24,30,21,
               41,29,32,21,57,16,51,37,13,68,16,37,15)

mean(cohorts$N)
sd(cohorts$N)
### plot breeding lifespans for some years with line for the mean 

closed.cohorts <- subset(cohorts, cohorts$FY < 2004)
#mean of a mean
mean(closed.cohorts$mean.yrs)
mean(brd$Yrs)
sd(brd$Yrs)

#write.table(cohorts,"cohorts.txt",sep="\t",row.names=FALSE)
#write.table(closed.cohorts, "closedcohorts.txt", sep="\t", row.names=FALSE)

#standard error of the mean
se.yrs <- (sd(brd$Yrs))/(sqrt(990))
se.yrs
se.yrs*1.96

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

sef <- (sd(females$Yrs))/(sqrt(478))
sef
sem <- (sd(males$Yrs))/(sqrt(480))
sem

mean(males$MinAgeFBr)
sd(males$MinAgeFBr)

mean(females$MinAgeFBr)
sd(females$MinAgeFBr)

#mean age at first breeding 
mean(males$MinAgeFBr)
sd(males$MinAgeFBr)
mean(females$MinAgeFBr)
sd(females$MinAgeFBr)

library(RColorBrewer)

# 4 6 2017 this does look good, use this in paper 
#remove gridlines 
#white background
# 3 1 2017 Use this one 
ggplot(brd, aes(MinAgeFBr, fill = Sex)) +
  geom_bar(position="dodge") +
  scale_fill_grey(start=0.6, end = 0.1) +
  xlab("Minimum Age at First Breeding") +
  ylab("Number of Individuals") +
  #labs(title = "Distribution of Ages at First Breeding (Years)") +
  #theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,250)) +
  scale_x_continuous(expand = c(0,0), breaks = c(1,2,3,4,5,6,7,8,9,10)) +
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank()) +
  theme_bw() +
  theme(panel.grid.major.x =element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank())
  #theme(panel.grid = element_blank(), panel.border = element_blank()) +
  #coord_cartesian(xlim = c(1,12), ylim = c(0,250))
  
  
#ggplot(brd, aes(MinAgeFBr, fill = Sex)) +
  #geom_histogram(binwidth = 0.2) +
  #xlab("Minimum Age at First Breeding") +
  #ylab("Frequency of Individuals") +
  #labs(title = "1981 to 2016")



#qplot(Yrs, data=brd, geom="histogram")
#Histogram of breeding life span, color by sex
#ggplot(brd, aes(Yrs, fill = Sex)) +
  #geom_histogram(binwidth = 0.08) +
  #xlab("Number of Years Bred") +
  #ylab("Frequency of Individuals") + 
  #labs(title = "1981 to 2016")

#ggplot(brd, aes(Yrs, fill = Sex)) +
  #geom_bar(position="dodge")

#ggplot(males, aes(Yrs)) + 
  #geom_histogram(binwidth = 0.1) +
  #xlab("Number of Years Bred") +
  #ylab("Frequency of Individuals") +
  #labs(title = "Males") 

#ggplot(females, aes(Yrs)) + 
  #geom_histogram(binwidth = 0.1) +
  #xlab("Number of Years Bred") +
  #ylab("Frequency of Individuals") +
  #labs(title = "Females")

##########################################################################
######## 4 15 2017 Don't use anymore, must move this later due to 
#my adding of birds from the all bird list
#Survival and hazard estimates
brd$Yrs <- as.numeric(brd$Yrs)

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

##ggplot2 version
library(GGally)
#KM curve for all years 
brds.all <- ggsurv(brd.fit, plot.cens=FALSE, 
                  back.white=TRUE) +
  theme(legend.position="none") +
  ylab("Cumulative Survival") +
  xlab("Time (Years)") + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank())
brds.all

brd.sex <- survfit(brd.ob ~ brd$Sex, conf.type="log-log")
sxplot <- plot(brd.sex, conf.int=TRUE, col = c("black","gray50"), 
               xlab = "Breeding time span (years)", 
               ylab = "Cumulative Survival",
               ylim = c(0.01,1), xlim = c(0,15), lwd = 2)
legend("topright", c("Females","Males"), col=c("black","gray50"),
       lty = c(1,1),lwd=2)

brds.sex <- ggsurv(brd.sex, plot.cens = FALSE, 
 back.white=TRUE,lty.est=1, surv.col = c("black","darkgray")) +
  xlab("Time (Years)") +
  ylab("Cumulative Survival") +
theme(panel.grid.major.y = element_blank(), 
      panel.grid.major.x=element_blank()) +
theme(panel.grid.minor.x = element_blank()) +
theme(panel.grid.minor.y = element_blank())
brds.sex 

str(brd.sex)
new.info <- data.frame(brd.sex$time, brd.sex$n.risk, brd.sex$n.event, 
        brd.sex$surv, brd.sex$upper, brd.sex$lower)
  

###Plot for Sam using ggsurv
#######################################################################
ggsurv1 <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                   cens.col = 'red', lty.est = 1, lty.ci = 2,
                   cens.shape = 3, back.white = F, xlab = 'Time',
                   ylab = 'Survival', main = ''){
  
  library(ggplot2)
  strata <- ifelse(is.null(s$strata) ==T, 1, length(s$strata))
  stopifnot(length(surv.col) == 1 | length(surv.col) == strata)
  stopifnot(length(lty.est) == 1 | length(lty.est) == strata)
  
  ggsurv.s <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = ''){
    
    dat <- data.frame(time = c(0, s$time),
                      surv = c(1, s$surv),
                      up = c(1, s$upper),
                      low = c(1, s$lower),
                      cens = c(0, s$n.censor))
    dat.cens <- subset(dat, cens != 0)
    
    col <- ifelse(surv.col == 'gg.def', 'black', surv.col)
    
    pl <- ggplot(dat, aes(x = time, y = surv)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(col = col, lty = lty.est)
    
    pl <- if(CI == T | CI == 'def') {
      pl + geom_step(aes(y = up), color = col, lty = lty.ci) +
        geom_step(aes(y = low), color = col, lty = lty.ci)
    } else (pl)
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  
  ggsurv.m <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = '') {
    n <- s$strata
    
    groups <- factor(unlist(strsplit(names
                                     (s$strata), '='))[seq(2, 2*strata, by = 2)])
    gr.name <-  unlist(strsplit(names(s$strata), '='))[1]
    gr.df <- vector('list', strata)
    ind <- vector('list', strata)
    n.ind <- c(0,n); n.ind <- cumsum(n.ind)
    for(i in 1:strata) ind[[i]] <- (n.ind[i]+1):n.ind[i+1]
    
    for(i in 1:strata){
      gr.df[[i]] <- data.frame(
        time = c(0, s$time[ ind[[i]] ]),
        surv = c(1, s$surv[ ind[[i]] ]),
        up = c(1, s$upper[ ind[[i]] ]),
        low = c(1, s$lower[ ind[[i]] ]),
        cens = c(0, s$n.censor[ ind[[i]] ]),
        group = rep(groups[i], n[i] + 1))
    }
    
    dat <- do.call(rbind, gr.df)
    dat.cens <- subset(dat, cens != 0)
    
    pl <- ggplot(dat, aes(x = time, y = surv, group = group)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(aes(col = group, lty = group))
    
    col <- if(length(surv.col == 1)){
      scale_colour_manual(name = gr.name, values = rep(surv.col, strata))
    } else{
      scale_colour_manual(name = gr.name, values = surv.col)
    }
    
    pl <- if(surv.col[1] != 'gg.def'){
      pl + col
    } else {pl + scale_colour_discrete(name = gr.name)}
    
    line <- if(length(lty.est) == 1){
      scale_linetype_manual(name = gr.name, values = rep(lty.est, strata))
    } else {scale_linetype_manual(name = gr.name, values = lty.est)}
    
    pl <- pl + line
    
    pl <- if(CI == T) {
      if(length(surv.col) > 1 && length(lty.est) > 1){
        stop('Either surv.col or lty.est should be of length 1 in order
             to plot 95% CI with multiple strata')
      }else if((length(surv.col) > 1 | surv.col == 'gg.def')[1]){
        pl + geom_step(aes(y = up, color = group), lty = lty.ci) +
          geom_step(aes(y = low, color = group), lty = lty.ci)
      } else{pl +  geom_step(aes(y = up, lty = group), col = surv.col) +
          geom_step(aes(y = low,lty = group), col = surv.col)}
    } else {pl}
    
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  pl <- if(strata == 1) {ggsurv.s(s, CI , plot.cens, surv.col ,
                                  cens.col, lty.est, lty.ci,
                                  cens.shape, back.white, xlab,
                                  ylab, main)
  } else {ggsurv.m(s, CI, plot.cens, surv.col ,
                   cens.col, lty.est, lty.ci,
                   cens.shape, back.white, xlab,
                   ylab, main)}
  pl
}

ggsurv1(brd.sex)
######################################################################
library(GGally)
#store summary of survfit object 
res <- summary(brd.sex)
str(res)
#select columns by index from str(res)
cols <- lapply(c(2:4,6:10), function(x) res[x])
#make data.frame
tbl <- do.call(data.frame, cols)
fems <- subset(tbl, strata == "brd$Sex=F")
mals <- subset(tbl, strata == "brd$Sex=M")
fems["strata"] <- "F"
mals["strata"] <- "M"

test.plot <- ggplot(fems, aes(time, surv)) +
  geom_line(colour="red", size = 1) 
test.plot

p <- ggsurv(brd.sex, CI=TRUE)
p + geom_ribbon(ymin = brd.sex$lower, ymax = brd.sex$upper)


survcurve <- data.frame(brd.sex$time, brd.sex$n.event,
                        brd.sex$surv, brd.sex$lower, brd.sex$upper)
colnames(survcurve) <- c("Time","Events","Surv","Lower","Upper")
ggplot(survcurve, aes(x=Time, y=Surv)) +
  geom_line()

###################################################################
yrplot <- plot(year.fit, xlab ="Breeding time span (years)", log = "y",
               ylab = "Cumulative Survival (log)", main = "Breeders", 
               ylim = c(0.02,1), xlim = c(0,15))

all.yrs.brds<- ggsurv(year.fit, plot.cens=FALSE, 
                     back.white=TRUE, surv.col="black") +
  
  theme(legend.position="none") +
  ylab("Cumulative Survival") +
  xlab("Time (Years)") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 
all.yrs.brds

cox.null <- coxph(brd.ob ~ 1, data = brd)
summary(cox.null)

cx1 <- coxph(brd.ob ~ Sex, data = brd)
summary(cx1)
cox.zph(cx1)
plot(cox.zph(cx1))
anova(cx1)

brd$SexInd <- as.integer(brd$SexInd)

#simPH functions 
sex2 <- coxph(brd.ob ~ SexInd, data= brd)
sim2 <- coxsimLinear(sex2, b = "SexInd", Xj = 0:1)
## 3 14 2017 add some stylistic modifications and adjust it so it just shows 
# 1 sex 
simGG(sim2, xlab = "Sex", method = "lm", type = "points",
      theme(axis.text.x = element_blank(), axis.ticks=element_blank())) 


sim2b <- coxsimLinear(sex2, b = "SexInd", Xj = c(0,1),
                      qi = "Hazard Rate")
simGG(sim2b)
cumhaz2 <- basehaz(sex2)
plot(cumhaz2$time, cumhaz2$hazrad, type = "l", 
     xlab = "Time", ylab = "Cumulative Hazard")
lines(cumhaz2$time, exp(0.9418)*cumhaz2$hazard, col="blue")

#get cumulative hazards
cumhaz <- basehaz(cx1)
## is this correct? Females with a higher cumulative hazard
plot(cumhaz$time, cumhaz$hazard, type = "l", lwd = 3, xlab = "Time", ylab = "Cumulative Hazard")
lines(cumhaz$time, exp(0.94177)*cumhaz$hazard, col="blue", lwd = 3)

#plot survival function
bslf <- survfit(cx1)
plot(bslf, xlab = "Breeding span (years)", ylab = "Cumulative Survival")

cx2 <- coxph(brd.ob ~ MinAgeFBr, data = brd)
summary(cx2)
cox.zph(cx2)
plot(cox.zph(cx2))
anova(cx2)
bslf2 <- survfit(cx2)
plot(bslf2, xlab = "Breeding span (years)", ylab = "Cumulative Survival")

sim3 <- coxsimLinear(cx2, b = "MinAgeFBr", Xj = c(1,12))
simGG(sim3)

cx3 <- coxph(brd.ob ~ FY, data = brd)
anova(cx3)
bslf3 <- survfit(cx3)
plot(bslf3, xlab = "Breeding span (years)", ylab = "Cumulative Survival")

sim4 <- coxsimLinear(cx3, b = "FY", Xj = c(1981,2016))

cx4 <- coxph(brd.ob ~ Sex + MinAgeFBr, data=brd)
anova(cx4)
bslf4 <- survfit(cx4)
plot(bslf4, xlab = "Breeding span (years)", ylab = "Cumulative Survival")

cxb <- coxph(brd.ob ~ MinAgeFBr + FY, data =brd)
anova(cxb)
summary(cxb)

cx5 <- coxph(brd.ob ~ Sex + MinAgeFBr + FY, data=brd)
anova(cx5)
cx6 <- coxph(brd.ob ~ Sex + MinAgeFBr + FY + Sex*MinAgeFBr, data = brd)
anova(cx6)
summary(cx6)

cx.noyear <- coxph(brd.ob ~ MinAgeFBr + Sex*MinAgeFBr, data = brd)
summary(cx.noyear)
anova(cx.noyear)

## AIC table
library(AICcmodavg)
Cand.models1 <- list()
Cand.models1[[1]] <- cx1
Cand.models1[[2]] <- cx2
Cand.models1[[3]] <- cx3
Cand.models1[[4]] <- cx4
Cand.models1[[5]] <- cx5
Cand.models1[[6]] <- cx6
Cand.models1[[7]] <- cxb
Modnames1 <- paste("Model", 1:length(Cand.models1), sep="")
hymodsA <- aictab(Cand.models1, Modnames1, sort = TRUE)
hymodsA
hymodsB <- aictab(Cand.models1, Modnames1, sort = FALSE)
hymodsB
#I think there is a "penalty" for large numbers of parameters
#The AIC scores are the exact opposite from the deviance reduction
#I suspect it's because of the # of parameters for year

anova(cx3, cxb, cx5, cx6, test = "chisq")

extractAIC(cx1)
extractAIC(cx2)
extractAIC(cx3)
extractAIC(cx4)
extractAIC(cx5)
extractAIC(cx6)
extractAIC(cxb)

#Deviance 
anova(cx1)
anova(cx2)
anova(cx3)
anova(cx4)
anova(cx5)
anova(cx6)

anova(cx2,cx4)



#Analysis of deviance table comparing first three Cox PH models 
dev.compare <- anova(cx1, cx2, cx3, cx4, cx5, cxb, cx6, test="Chisq")
dev.compare

dev2 <- anova(cx3,cxb, cx5,cx6)
dev2

#Calculate residuals for Coxph fit 
resd <- residuals(cx1, type="deviance", collapse=brd$ID)
#Ask GF about model diagnostics 
plot(resd)

#Frailty models where cohort year is a random effect 
frail1 <- coxme(brd.ob ~ MinAgeFBr + (1|FY), data = brd)
summary(frail1)
frail2 <- coxme(brd.ob ~ Sex + (1|FY), data= brd)
summary(frail2)
frailnull <- coxme(brd.ob ~ (1|FY), data= brd)

anova(frail2, frailnull)

frail3 <- coxme(brd.ob ~ MinAgeFBr + Sex*MinAgeFBr + (1|FY), data = 
                  brd)
summary(frail3)
anova(cx2, cx.noyear, frail3)
anova(cx2, frail3)
anova(cx2, frail1)
anova(cx2, frail1, frail3)

## AFT frailty models 

AFT1 <- survreg(brd.ob ~ MinAgeFBr + frailty(FY, dist='gamma'), 
                       data = brd, dist = "weibull")
summary(AFT1)
#shape parameter a
shape.a <- 1/AFT1$scale
shape.a
#Good it's over one, meaning failure rate increases with time
#As it should be for jays

#Try Gaussian distribution 
fr.gaus <- survreg(brd.ob ~ MinAgeFBr + frailty(FY, dist='gaussian'), 
                   data=brd, dist = "weibull")
summary(fr.gaus)


AFT2 <- survreg(brd.ob ~ Sex + MinAgeFBr + frailty(FY, dist='gamma'),
                data=brd, dist = "weibull")
summary(AFT2)

fit.age <- data.frame(predict(cx2, type = "risk", se.fit=TRUE))
age <- brd$MinAgeFBr
age <- na.omit(age)
fit.age$age <- age
colnames(fit.age)[3] <- "Age"


predframe.age <- with(fit.age, data.frame(Age, hazard=fit, 
          lwr = fit - 1.96*se.fit, upr = fit + 1.96*se.fit))

##Mass poly - adapted from code written for mass with no polynomial
pAge <- ggplot(predframe.age, aes(Age, hazard)) +
  geom_line()
pAge + ylim(0,3) +
  xlab("Age at First Breeding") +
  ylab("Predicted Relative Hazard") + 
  theme_bw() +
  geom_hline(yintercept = 1.0, linetype = 2) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  geom_ribbon(data = predframe.age, aes(ymin=lwr, ymax=upr), alpha = 0.2) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 

####################################################################
# Block 2 - Input of Vegetation, fire, terr size, group size
# Time-varying covarites 
# Manipulation for use in models

dom.veg <- read.csv("dom_veg.csv")
dom.veg <- subset(dom.veg, InABS == TRUE)


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

#Time since fire data

tsf <- read.csv("tsf_terr.csv")
tsf <- subset(tsf, InABS == TRUE)

#TSF data
#Create object for the numbers, same logic as with veg data
keep2 <- c(2,3,4,5,6,7,8,9)
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

#Weighted mean of fire patches

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
#Want to add the terr information to the one record data for first year
brd.terr <- merge(brd, terr.info, by = "TerrYr")
#brd.terr <- brd.terr[,c(2,3,4,5,6,7,8,9,10,11,12,13,1,14,15,16,17)]

#standardized for interpretation of coefficients, model fit the same
brd.terr["stdscrb"] <- scale(brd.terr$scrb.count, 
                        center = FALSE, scale = TRUE)
brd.terr["stdtsf"] <- scale(brd.terr$tsf.count, 
                        center = FALSE, scale = TRUE)
brd.terr["stdsize"] <- scale(brd.terr$TerrSize, 
                        center = FALSE, scale = TRUE)

colnames(brd.terr)[15] <- "OakScrub"
colnames(brd.terr)[17] <- "TSF"


str(brd.terr)

brd.terr$Fbreed <- as.numeric(brd.terr$Fbreed)
brd.terr$LastObsDate <- as.numeric(brd.terr$LastObsDate)
brd.terr$Yrs <- as.numeric(brd.terr$Yrs)

#corrplot
bcovs <- brd.terr[,c(15,16,17)]
bcorr <- cor(bcovs, use="complete.obs")
corrplot(bcorr, method = "pie", type = "lower")

######################################################################
#PCA 
bd.tvars <- data.frame(brd.terr$OakScrub, brd.terr$TerrSize, brd.terr$TSF)
colnames(bd.tvars) <- c( "OakScrub","TerrSize", "TSF")
bd.pca <- prcomp(bd.tvars, center = TRUE, scale. = TRUE)
print(bd.pca)
summary(bd.pca)
plot(bd.pca)

pc.scores1 <- data.frame(bd.pca$x)
pc.scores1 <- pc.scores1[,1:2]

brd.terr <- cbind(brd.terr, pc.scores1)

#get pc scores

#####################################################################
#Survival object
brdtr.ob <- Surv(brd.terr$Yrs, brd.terr$Censor, type=c('right'))
fit1 <- survfit(brdtr.ob ~ 1, conf.type="log-log")
fit1
plot(fit1, log="y")
year.fit <- survfit(brdtr.ob ~ brd.terr$FY, conf.type="log-log")
year.fit
plot(year.fit, log="y")
plot(year.fit)

### PC's s
cox.pc <- coxph(brdtr.ob ~ PC1 + PC2, data = brd.terr)
summary(cox.pc)
anova(cox.pc)

#Cox model with oak scrub, is PH assumption violated?
cox.oak <- coxph(brdtr.ob ~ stdscrb, data = brd.terr)
sim1 <- coxsimLinear(cox.oak, b = "stdscrb", Xj = c(0,3.6))
simGG(sim1)

plot(cox.zph(cox.oak))

summary(cox.oak)
plot(cox.zph(cox.oak))
anova(cox.oak)

#2-9 year tsf
cox.tsf <- coxph(brdtr.ob ~ stdtsf, data = brd.terr)
summary(cox.tsf)
sim2 <- coxsimLinear(cox.tsf, b = "stdtsf", Xj = c(0, 3.8))
simGG(sim2)
anova(cox.tsf)

plot(cox.zph(cox.tsf))

cox.size <- coxph(brdtr.ob ~ stdsize, data = brd.terr)
summary(cox.size)
plot(cox.zph(cox.size))

cox.int <- coxph(brdtr.ob ~ 1, data= brd.terr)
summary(cox.int)

anova(cox.int, cox.oak, cox.tsf, cox.size)

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

length(unique(brd.allyrs$band))
length(unique(brd.allyrs$JayID))

#Change dates to dates and calculate time lived 

brd.allyrs$BrdrDate <- as.Date(brd.allyrs$BrdrDate, format = "%m/%d/%Y")
brd.allyrs$LastObsDate <- as.Date(brd.allyrs$LastObsDate,format = "%m/%d/%Y")

brd.allyrs["Days"] <- brd.allyrs$LastObsDate - brd.allyrs$BrdrDate
brd.allyrs$Days <- as.numeric(brd.allyrs$Days)
brd.allyrs["Years"] <- brd.allyrs$Days/365.25

#Add censorship column
brd.allyrs["Censor"] <- 1
brd.allyrs$Censor[which(brd.allyrs$LastObsDate =="2016-4-12")] <- 0

brd.allyrs <- subset(brd.allyrs, brd.allyrs$Days > 0 & brd.allyrs$Years > 0)
length(unique(brd.allyrs$JayID))
length(unique(brd.allyrs$band))

all.brd <- merge(brd.allyrs, hlp.byterr)
#Duplicate records for multiple nests
#This seems to have worked = I checked a few jays to make sure matches
#are correct 

all.brd["Status"] <- 0 

length(unique(all.brd$band))
length(unique(all.brd$JayID))

#### Jan 30 2017 At this point, the data frame "all.brd" has the largest N
# N = 993
all.brd$band <- as.character(all.brd$band)
all.brd$JayID <- as.character(all.brd$JayID)


all.brd$band[which(all.brd$JayID == "RLC-U")] <- "3"
all.brd$band[which(all.brd$JayID == "-SW_")] <- "4"
all.brd$band[which(all.brd$JayID == "6962")] <- "5"

#remove duplicate rows mulitple nests per year
all.brd <- all.brd[!duplicated(all.brd),]

length(unique(all.brd$band))
length(unique(all.brd$JayID))

#rows removed for multiple nest attempts in a year if one failed but on
#succeed, I suppose I could keep them in and add nest attempt if 
#we think that sort of analyses is worthwhile, and doesn't take a zillion
#hours 

### Models with these data on helpers? The sample size is greater

#Still have to rearrange data

str(all.brd)
## N = 993, Females N = 497 Males N = 495
## Jan 30 2017 
## Work on transforming these data to be compatible with time dep models
names(all.brd)

#rearrange columns 
all.brd <- all.brd[ ,c(1,3,4,5,6,2,7,8,9,10,11,12,13,14,15,21,16,17,18,19,20)]
names(all.brd)
colnames(all.brd)[2] <- "USFWBand"

#all.brd <- subset(all.brd, all.brd$Days > 0 & all.brd$Years > 0)
length(unique(all.brd$USFWBand))
length(unique(all.brd$JayID))
# N = 968 Jay ID (a few breeders dont have a band number)
# F = 482
# M = 483


##################################################################
# 4 15 2017 FUCK MY FUCKING LIFE I'M A PIECE OF SHIT
#Who is in all birds but not the original list? 
whoru <- subset(all.brd, !(all.brd$USFWBand %in% brd$USFWBand))
whoru2 <- whoru[!duplicated(whoru$JayID),]
length(unique(whoru2$JayID))
names(whoru2)
names(brd)

colnames(brd)[4] <- "BrdrDate"
colnames(brd)[3] <- "FirstYr"

whoru2["Help"] <- NULL
whoru2["Status"] <- NULL
whoru2["sex"] <- NULL
whoru2["YrsExp"] <- NULL
whoru2["HelpNum"] <- NULL
whoru2["NestYear"] <- NULL
whoru2["Days"] <- NULL
whoru2["CurrentAge"] <- NULL
whoru2["MinimumAge"] <- NULL

whoru2 <- whoru2[, c(1,2,5,3,4,7,6,10,11,8,9,12)]
names(whoru2)
brd["FY"] <- NULL
whoru2["SexInd"] <- 0
colnames(whoru2)[6] <- "MinAgeFBr"
colnames(whoru2)[10] <- "Yrs"


#this is for the original data set with one record per individual 
brd.April15 <- rbind(brd, whoru2)
brd.April15$FirstYr <- as.factor(brd.April15$FirstYr)
brd.April15$SexInd[which(brd.April15$Sex == "F")] <- 0
brd.April15$SexInd[which(brd.April15$Sex == "M")] <-1
brd.April15$Sex[which(brd.April15$SexInd == "0")] <- "F"
brd.April15$Sex[which(brd.April15$SexInd == "1")] <- "M"

brd.April15$SexInd <- as.factor(brd.April15$SexInd)

length(unique(brd.April15$JayID))
length(unique(brd.April15$USFWBand))

##### 4 15 2017 Do basic models HERE 
brd.April15$Yrs <- as.numeric(brd.April15$Yrs)

brd.ob <- Surv(brd.April15$Yrs, brd.April15$Censor,type= c('right'))
brd.fit <- survfit(brd.ob ~ 1, conf.type = "log-log")
brd.fit

year.fit <- survfit(brd.ob ~ brd.April15$FirstYr, conf.type="log-log")
year.fit
plot(year.fit, log="y")

## Mean years survived as breeder per cohort 
cohorts1 <- ddply(brd.April15, .(FirstYr), summarise, 
                 mean.yrs = round(mean(Yrs), digits = 2), 
                 sd.yrs = round(sd(Yrs), digits = 2),
                 var.yrs = round(var(Yrs), digits = 2))
cohorts1$FY <- as.character(cohorts1$FY)
cohorts1$FY <- as.numeric(cohorts1$FY)
cohorts1$N <- c(year.fit$n)
  
#write.table(cohorts1, "cohorts1.txt", sep="\t", row.names = FALSE)  
  
# adjust x axis it's breeding span 
kmplot <- plot(brd.fit, xlab="Breeding time span (years)", log = "y", 
               ylab = "Cumulative Survival (log)", main = "Breeders",
               ylim = c(0.01,1), xlim = c(0,15))

##ggplot2 version
library(GGally)
#KM curve for all years 
brds.all <- ggsurv(brd.fit, plot.cens=FALSE, 
                   back.white=TRUE) +
  theme(legend.position="none") +
  ylab("Cumulative Survival") +
  xlab("Time (Years)") + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank())
brds.all

brd.sex <- survfit(brd.ob ~ brd.April15$Sex, conf.type="log-log")
#brd.sex <- survfit(brd.ob ~ brd.April15$SexInd, conf.type="log-log")
sxplot <- plot(brd.sex, conf.int=TRUE, col = c("firebrick3","dodgerblue1"), 
               xlab = "Breeding time span (years)", 
               ylab = "Cumulative Survival",
               ylim = c(0.01,1), xlim = c(0,15), lwd = 2)
legend("topright", c("Females","Males"), col=c("firebrick3","dodgerblue1"),
       lty = c(1,1),lwd=2)

legend.title <- "Sex"

brds.sex <- ggsurv(brd.sex, plot.cens = FALSE,
  back.white=TRUE,lty.est=c(1,2), surv.col = c("black","black")) +
  xlab("Time (Years)") +
  ylab("Cumulative Survival") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  labs(colour = "Sex")
brds.sex

brds.sex + ggplot2::guides(linetype=FALSE) + 
  ggplot2::scale_colour_discrete(name = 'Sex', breaks =c(1,2),
                                 labels = c('Male','Female'))

str(brd.sex)
new.info <- data.frame(brd.sex$time, brd.sex$n.risk, brd.sex$n.event, 
                       brd.sex$surv, brd.sex$upper, brd.sex$lower)

yrplot <- plot(year.fit, xlab ="Breeding time span (years)", log = "y",
               ylab = "Cumulative Survival (log)", main = "Breeders", 
               ylim = c(0.02,1), xlim = c(0,15))

all.yrs.brds<- ggsurv(year.fit, plot.cens=FALSE, 
                      back.white=TRUE, surv.col="black") +
  
  theme(legend.position="none") +
  ylab("Cumulative Survival") +
  xlab("Time (Years)") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 
all.yrs.brds

cox.null <- coxph(brd.ob ~ 1, data = brd.April15)
summary(cox.null)

cx1 <- coxph(brd.ob ~ SexInd, data = brd.April15)
summary(cx1)
cox.zph(cx1)
plot(cox.zph(cx1))
anova(cx1)

#brd$SexInd <- as.integer(brd.April15$SexInd)


#simPH functions 
sex2 <- coxph(brd.ob ~ SexInd, data= brd.April15)
sim2 <- coxsimLinear(sex2, b = "SexInd", Xj = 0:1)
## 3 14 2017 add some stylistic modifications and adjust it so it just shows 
# 1 sex 
simGG(sim2, xlab = "Sex", method = "lm", type = "points",
      theme(axis.text.x = element_blank(), axis.ticks=element_blank())) 


sim2b <- coxsimLinear(sex2, b = "SexInd", Xj = c(0,1),
                      qi = "Hazard Rate")
simGG(sim2b)
cumhaz2 <- basehaz(sex2)
plot(cumhaz2$time, cumhaz2$hazrad, type = "l", 
     xlab = "Time", ylab = "Cumulative Hazard")
lines(cumhaz2$time, exp(0.9418)*cumhaz2$hazard, col="blue")

#get cumulative hazards
cumhaz <- basehaz(cx1)
## is this correct? Females with a higher cumulative hazard
plot(cumhaz$time, cumhaz$hazard, type = "l", lwd = 3, xlab = "Time", ylab = "Cumulative Hazard")
lines(cumhaz$time, exp(0.94177)*cumhaz$hazard, col="blue", lwd = 3)

#plot survival function
bslf <- survfit(cx1)
plot(bslf, xlab = "Breeding span (years)", ylab = "Cumulative Survival")

cx2 <- coxph(brd.ob ~ MinAgeFBr, data = brd.April15)
summary(cx2)
cox.zph(cx2)
plot(cox.zph(cx2))
anova(cx2)
bslf2 <- survfit(cx2)
plot(bslf2, xlab = "Breeding span (years)", ylab = "Cumulative Survival")

#Adapted from code used in hatch year birds 
fit.age <- data.frame(predict(cx2, type="risk", se.fit=TRUE))
no.NA <- brd.April15 %>% drop_na(6)
fit.age$Age <- no.NA$MinAgeFBr

colnames(fit.age)[3] <- "Age"
predframe.age <- with(fit.age, data.frame(Age, hazard=fit, 
        lwr = fit - 1.96*se.fit, upr = fit + 1.96*se.fit))

pB <- ggplot(predframe.age, aes(Age, hazard)) +
  geom_line()
pB + ylim(0,3) +
  xlab("Age at First Breeding") +
  ylab("Predicted Relative Risk") + 
  theme_bw() +
  geom_hline(yintercept = 1.0, linetype = 2) +
  geom_ribbon(data = predframe.age, aes(ymin=lwr, ymax=upr), alpha = 0.2) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10)) 



sim3 <- coxsimLinear(cx2, b = "MinAgeFBr", Xj = c(1,12))
simGG(sim3)

cx3 <- coxph(brd.ob ~ FirstYr, data = brd.April15)
anova(cx3)
summary(cx3)
bslf3 <- survfit(cx3)
plot(bslf3, xlab = "Breeding span (years)", ylab = "Cumulative Survival")


sim4 <- coxsimLinear(cx3, b = "FY", Xj = c(1981,2016))

cx4 <- coxph(brd.ob ~ SexInd + MinAgeFBr, data=brd.April15)
anova(cx4)
anova(cx2, cx4)
anova(cx)
bslf4 <- survfit(cx4)
plot(bslf4, xlab = "Breeding span (years)", ylab = "Cumulative Survival")


cxb <- coxph(brd.ob ~ MinAgeFBr + FirstYr, data =brd.April15)
anova(cxb)
summary(cxb)

cx5 <- coxph(brd.ob ~ SexInd + MinAgeFBr + FirstYr, data=brd.April15)
anova(cx5)
cx6 <- coxph(brd.ob ~ MinAgeFBr + FirstYr + SexInd*MinAgeFBr, 
             data = brd.April15)
anova(cx6)
summary(cx6)

anova(cx2, cx5, cxb)

anova(cx2, cxb, cx5, cx6)
anova(cx2, cxb, cx6)
anova(cx2, cxb, cx5)

## AIC table
library(AICcmodavg)
Cand.models1 <- list()
Cand.models1[[1]] <- cx1
Cand.models1[[2]] <- cx2
Cand.models1[[3]] <- cx3
Cand.models1[[4]] <- cx4
Cand.models1[[5]] <- cx5
Cand.models1[[6]] <- cx6
Cand.models1[[7]] <- cxb
Modnames1 <- paste("Model", 1:length(Cand.models1), sep="")
hymodsA <- aictab(Cand.models1, Modnames1, sort = TRUE, second.ord=FALSE)
hymodsA
hymodsB <- aictab(Cand.models1, Modnames1, sort = FALSE)
hymodsB
#I think there is a "penalty" for large numbers of parameters
#The AIC scores are the exact opposite from the deviance reduction
#I suspect it's because of the # of parameters for year

Cand.models2 <- list()
Cand.models2[[1]] <- cx1
Cand.models2[[2]] <- cx3
Modnames2 <- paste("Model", 1:length(Cand.models2), sep ="")
hymodsC <- aictab(Cand.models2, Modnames2, sort =TRUE)
hymodsC

anova(cxb, cx5, cx6, test = "chisq")

extractAIC(cx1)
extractAIC(cx2)
extractAIC(cx3)
extractAIC(cx4)
extractAIC(cx5)
extractAIC(cx6)
extractAIC(cxb)

#Deviance 
anova(cx1)
anova(cx2)
anova(cx3)
anova(cx4)
anova(cx5)
anova(cx6)

anova(cx2,cx4)



#Analysis of deviance table comparing first three Cox PH models 
dev.compare <- anova(cx1, cx2, cx3, cx4, cx5, cxb, cx6, test="Chisq")
dev.compare

dev2 <- anova(cx3,cxb, cx5,cx6)
dev2

#Calculate residuals for Coxph fit 
resd <- residuals(cx1, type="deviance", collapse=brd$ID)
#Ask GF about model diagnostics 
plot(resd)

#Frailty models where cohort year is a random effect 
frail1 <- coxme(brd.ob ~ MinAgeFBr + (1|FirstYr), data = brd.April15)
summary(frail1)
anova(cx2, frail1)
frail2 <- coxme(brd.ob ~ SexInd + (1|FirstYr), data= brd.April15)
summary(frail2)
frailnull <- coxme(brd.ob ~ (1|FirstYr), data= brd.April15)

frail.AS <- coxme(brd.ob ~ MinAgeFBr + SexInd + (1|FirstYr), 
                  data = brd.April15)
summary(frail.AS)
anova(cx5)


anova(frail2, frail1)

## AFT frailty models 

AFT1 <- survreg(brd.ob ~ MinAgeFBr + frailty(FirstYr, dist='gamma'), 
                data = brd.April15, dist = "weibull")
summary(AFT1)

AFT2 <- survreg(brd.ob ~ Sex + MinAgeFBr + frailty(FirstYr, dist='gamma'),
                data=brd.April15, dist = "weibull")
summary(AFT2)

AFT3 <- survreg(brd.ob ~ MinAgeFBr + MinAgeFBr*Sex + frailty(FirstYr, dist='gamma'),
                 data=brd.April15, dist = "weibull")
summary(AFT3)
anova(AFT2, AFT3)

#################################################################
#all.brd$band[which(all.brd$JayID=="-SW_")] <- "1000-10000"
#all.brd$band[which(all.brd$JayID==6962)] <- "2000-20000"
#all.brd$band[which(all.brd$JayID=="RLC-U")] <- "3000-30000"

#all.brd$band <- as.character(all.brd$band)

#Create the final stop interval 
#newdat <- tmerge(data1=all.brd, data2=all.brd, id=JayID, tstop=Days)

#colnames(brd)[1] <- "JayID"
#colnames(brd)[2] <- "band"
#newdat <- tmerge(data1=brd, data2=all.brd, id=band, tstop=Days)

#Change helper variables to factors
#all.brd$Help <- as.factor(all.brd$Help)
#all.brd$HelpNum <- as.factor(all.brd$HelpNum)




#####################################################################
all.brd[,18] <- NULL
#Uplaod data set with time intervals 
breeders.real <- read.csv("bdtv.csv")
str(breeders.real)

breeders.real$BrdrDate <- as.Date(breeders.real$BrdrDate, format = "%m/%d/%Y")
breeders.real$LastObsDate <- as.Date(breeders.real$LastObsDate, format = "%m/%d/%Y")

breeders.real$Status[which(breeders.real$LastObsDate == "2016-04-12")] <- 0
colnames(breeders.real)[1] <- "USFWBand"
#####################################################################
########example with subsample
samp <- subset(birds, birds$Help == 0)
samp1 <- subset(birds, birds$Help == 1)

bsamp <- rbind(samp, samp1)
str(bsamp)

bsamp <- na.omit(bsamp)
bsamp$HelpNum <- as.factor(bsamp$HelpNum)
bsamp$tstart <- as.numeric(bsamp$tstart)
bsamp$tstop <- as.numeric(bsamp$tstop)

bb <- subset(bsamp, bsamp$tstop > bsamp$tstart)

test <- coxph(Surv(tstart, tstop, Status) ~ Help, data = bb)
summary(test)
test2 <- coxph(Surv(tstart, tstop, Status) ~ HelpNum, data = bb)
summary(test2)

########################################################################
brd.sam <- all.brd[1:100,]
write.csv(brd.sam, file = "brdsample.csv")

#Adding standardized date column starting from beginning of study 

#first year breeding in date format
myyears <- format(all.brd$BrdrDate, '%Y')
myyears <- as.POSIXlt(all.brd$BrdrDate)
#Breeder year in years since 1900
years <- myyears$year


#### program to calculate time intervals
newdf <- data.frame(matrix(NA, nrow=3429, ncol = 9))
tstart <- vector()
tstop <- vector()

for (i in levels(all.brd["JayID"]))
  if 
    tstart = 0
    
  

#all.brd["tstart"] <- ""
#all.brd["tstop"] <- ""

#How to add intervals 


all.brd <- all.brd[,c(1,2,3,4,5,6,13,14,15,16,22,23,
                      7,8,9,10,11,12,17,18,19,20,21)]

write.csv(all.brd, "Jan30.csv")



######################################################################

# need to find a way to get the terr info in without deleting the records
# where this info is not available 

#get band # and terryr for breeders!!! 
#bandty <- data.frame(brd.April15$JayID, brd.April15$USFWBand, brd.April15$TerrYr)
bandty <- data.frame(brd$JayID, brd$USFWBand, brd$TerrYr)
colnames(bandty) <- c("JayID", "USFWBand", "TerrYr")

#seems to be correct 
birdsAB <- merge(x = breeders.real, y = terr.info, by ="TerrYr", all.x =TRUE)
colnames(breeders.real)[1] <- "USFWBand"
birdsABC <- merge(x = birdsAB, y = bandty, by = "USFWBand", 
                                                  all.x = TRUE)
birdsABC[26] <- NULL
birdsABC <- birdsABC[, c(1,25,2,3,4,5,6,7,8,9,10,11,12,13,14,
                         15,16,17,18,19,20,21,22,23,24)]

length(unique(birdsABC$USFWBand))
length(unique(birdsABC$JayID))

# terr.info from above is the territory info (size, oak, tsf) by terryr
#########################################################################
#losing rows but I don't understand why, I think bc no breeding at terr
#Data frame with most information needed 
all <- merge(all.brd, terr.info, all.x=TRUE)
### seems to have worked.......
length(unique(all$JayID))
length(unique(all$USFWBand))

#sample size reduced 
length(unique(all$band))
fem <- subset(all, all$Sex=="F")
length(unique(fem$band))
mal <- subset(all, all$Sex=="M")
length(unique(mal$band))

#reorganize data frame 
#two columns for sex and some other unnecessary cols 
#change order to better facilitate analyses

#brd <- brd[,c(1,2,8,7,6,3,4,5)]    #to rearrange cols
names(all)
all.info <- all[,c(2,4,8,9,5,6,7,18,11,3,12,14,15,16,1,17,21,22,23,19,20)]
names(all.info)
colnames(all.info)[17] <- "OakScrub"
colnames(all.info)[19] <- "TSF"

all.info$Help <- as.factor(all.info$Help)
all.info$FirstYr <- as.factor(all.info$FirstYr)
length(unique(all.info$JayID))

#######################################################################
# read in csv files with fledge data 
# males and females in separate files 

femflg <- read.csv("Erin_Females_fledge_year.csv")
malflg<- read.csv("Erin_Males_fledge_year.csv")

str(femflg)
str(malflg)

# combine territory name with last two digits of year for TerrYr
malflg["TerrYr"] <- paste(malflg$Terr, str_sub(malflg$NestYear, start= -2),
                          sep ="")
femflg["TerrYr"] <- paste(femflg$Terr, str_sub(femflg$NestYear, start = -2),
                          sep="")

#remove failed nests, have to look at who drops out to find years where 
#birds truly had zero fledglings 

#set blanks to NA's
malflg[malflg==""] <- NA
femflg[femflg==""] <- NA

#Keep rows with a fledge date, so only includes birds with fledglings 
#need some way to know who had zero fledglings 

#malflg <- subset(malflg, !malflg$FldgDate == "NA")
#femflg <- subset(femflg, !femflg$FldgDate == "NA")

colnames(malflg)[1] <-"JayID"
colnames(femflg)[1] <- "JayID"

fldgs <- rbind(malflg, femflg)

#Sum fledge number by terryear and bird 
newdata <- ddply(fldgs, .(TerrYr, JayID), summarise, FldgNum = sum(FldgNum))
newdata[2] <- NULL

#need this step for correct merge
fldgs2<- newdata[!duplicated(newdata$TerrYr),]

colnames(birdsABC)[3] <- "TerrYr"

#works 
jays.df <- merge(x = birdsABC, y = fldgs2, by = "TerrYr", all.x = TRUE)

#But, this needs to be sorted by band number and in order because I merged
#by terryr
real.df <- arrange(jays.df, USFWBand, NestYear)

#add in id-less birds by hand
#fix(real.df)

colnames(real.df)[23] <- "OakScrub"
colnames(real.df)[25] <- "TSF"
real.df["SdOak"] <- scale(real.df$OakScrub, center = FALSE, 
                          scale = TRUE)
real.df["SdSize"] <- scale(real.df$TerrSize, center = FALSE, 
                          scale = TRUE)
real.df["SdTSF"] <- scale(real.df$TSF, center = FALSE, 
                          scale = TRUE)
names(real.df)
str(real.df)
real.df$FldgNum[is.na(real.df$FldgNum)] <- 0

length(unique(real.df$USFWBand))
length(unique(real.df$JayID))
#N = 968

#change NA's to blank in cols where imputation not necessary
real.df[,3] <- as.character(real.df[,3])
real.df[,3][is.na(real.df[,3])] <- ""
real.df$AgeKnown <- as.character(real.df$AgeKnown)
real.df$AgeKnown[is.na(real.df$AgeKnown)] <- ""
real.df[,14] <- as.character(real.df[,14])
real.df[,14][is.na(real.df[,14])] <- ""
real.df[,15] <- as.character(real.df[,15])
real.df[,15][is.na(real.df[,14])] <- ""
real.df[,16] <- as.character(real.df[,16])
real.df[,16][is.na(real.df[,16])] <- ""
real.df[,17] <- as.character(real.df[,17])
real.df[,17][is.na(real.df[,17])] <- ""
real.df[,19] <- as.character(real.df[,19])
real.df[,19][is.na(real.df[,19])] <- ""

length(unique(real.df$USFWBand))
length(unique(real.df$JayID))
# 4 15 2017 N = 968
str(real.df)
#real.df$Help <- as.factor(real.df$Help)
#real.df$HelpNum <- as.factor(real.df$HelpNum)
real.df$TerrSize <- as.numeric(real.df$TerrSize)
real.df$MinAgeFBreed <- as.numeric(real.df$MinAgeFBreed)

real.df$YrsExp <- as.numeric(real.df$YrsExp)
real.df$CurrentAge <- as.numeric(real.df$CurrentAge)
#######################################################################
## Imputation then modeling, FINALLY  
library(mice)
library(VIM)
library(mi)

countNA(real.df$Help)
242/3368
#6% missing, general rule is less than 5% and don't have to impute
countNA(real.df$OakScrub)
#50% missing
countNA(real.df$TSF)
#same value

## 2 24 2017
# MAR - missing at random 
# Info partly depends on terr year and the location of terrs in space

#pool()  and with()
msdat.df <- real.df[,c(2,3,21,22,23,24,25)]
str(msdat.df)

imput.df <- mice(msdat.df, MaxNWts = 7000)

#Example from mice package
data(nhanes)
head(nhanes)
md.pattern(nhanes)
p <- md.pairs(nhanes)

######################################################################
#2 27 2017 Time varying models with complete cases, probably can get helper
#data from the censuses for the missing years
# Also have to adjust time intervals for birds who died right away in subs.
# breeding seasons, where next interval ends less than a multiple of 365 days

names(real.df)
str(real.df)

#write.table(real.df, "JaysBrdrs.txt", sep="\t")

## 5 4 2017 - added in missing terryrs
## 5 30 2017 -  removed duplicate rows
new.realdf <- read.table("JaysBrdrs.txt", header = TRUE)

## Add breeder output info
names(jays.fec)
str(jays.fec)

#merge by terryr
test <- merge(new.realdf, jays.fec, by = "TerrYr", all.x=TRUE)

#need to add group size, not the same as helper number 
str(group.size)
group.size <- group.size[,c(1,3)]
colnames(group.size) <- c("GroupSize", "TerrYr")

new.realdf <- merge(test, group.size, all.x=TRUE)
#real.df <- merge(real.df, group.size, by = "TerrYr", all.x = TRUE)
new.realdf$TerrYr <- as.character(new.realdf$TerrYr)
new.realdf$USFWBand <- as.character(new.realdf$USFWBand)
#real.df$TerrYr <- as.character(real.df$TerrYr)
#real.df$USFWBand <- as.character(real.df$USFWBand)



#length(unique(real.df$JayID))
#length(unique(real.df$USFWBand))

# 5 4 2017
## a handful of records where year and time interval out of order
#reprows <- ddply(hybrd, .(hybrd$JayID), nrow)
reps <- ddply(new.realdf, .(new.realdf$JayID, new.realdf$NestYear.x), nrow)
       
#write.csv(real.df, file = "ErinApril20.csv", sep="", col.names = TRUE)

#t.mess <-  subset(real.df, real.df$tstart > real.df$tstop)
#t.mess$tstart <- t.mess$tstop - 1
#real.df$tstart[which(real.df$tstart >= real.df$tstop)] <- real.df$tstop 
#info <- data.frame(t.mess$band, t.mess$JayID, t.mess$tstart, t.mess$tstop)
#colnames(info) <- c("band", "JayID", "tstart", "tstop")
# adjust time intervals
#t.mess$tstart <- t.mess$tstop - 1 
#newish <- subset(real.df, real.df$tstop >= real.df$tstart)
#real.df$tstart[which(real.df$tstart >= real.df$tstop)] <- real.df$tstop - 1

### Without imputation, just with the data I have 



#Fledge number 
mod1 <- coxph(Surv(tstart, tstop, Status) ~ FldgNum, data = new.realdf)
summary(mod1)
anova(mod1)
plot(cox.zph(mod1))
extractAIC(mod1)

#Helpers 0 or 1, doesn't really make sense 
#mod2 <- coxph(Surv(tstart, tstop, Status) ~ Help, data = real.df)
#summary(mod2)

#Area of oak scrub
mod3 <- coxph(Surv(tstart, tstop, Status) ~ SdOak, data = new.realdf)
summary(mod3)
anova(mod3)
extractAIC(mod3)

#Area of 2-9 tsf
mod4 <- coxph(Surv(tstart, tstop, Status) ~ SdTSF, data = new.realdf)
summary(mod4)
anova(mod4)
extractAIC(mod4)



tsize <- na.omit(new.realdf)
mod.size <- coxph(Surv(tstart, tstop, Status) ~ MinAgeFBreed +
          SdSize + GroupSize, data = tsize)
summary(mod.size)
anova(mod.size)
extractAIC(mod.size)


fit.ter <- data.frame(predict(mod.size, type="risk", se.fit=TRUE))
fit.ter$SdSize <- tsize$SdSize 
predframe.terr <- with(fit.ter, data.frame(SdSize, hazard=fit,
            lwr = fit - 1.96*se.fit, upr = fit + 1.96*se.fit))
pTS <- ggplot(predframe.terr, aes(SdSize, hazard)) + 
  geom_line()

#Breeder experience 
mod5 <- coxph(Surv(tstart, tstop, Status) ~ YrsExp, data = new.realdf)
summary(mod5)
anova(mod5)

fit.ye <- data.frame(predict(mod5, type = "risk", se.fit=TRUE))
fit.ye$YE <- real.df$YrsExp
colnames(fit.ye)[3] <- "YrsExp"

predframe.ye <- with(fit.ye, data.frame(YrsExp, hazard=fit, 
        lwr = fit - 1.96*se.fit, upr = fit + 1.96*se.fit))
pYE <- ggplot(predframe.ye, aes(YrsExp, hazard)) +
  geom_line()
pYE + ylim(0,2) +
  xlab("Years Experienced") +
  ylab("Predicted Relative Hazard") + 
  theme_bw() +
  geom_hline(yintercept = 1.0, linetype = 2) +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  geom_ribbon(data = predframe.ye, aes(ymin=lwr, ymax=upr), alpha = 0.2) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 


mod.age <- coxph(Surv(tstart, tstop, Status) ~ CurrentAge, data = new.realdf)
summary(mod.age)
anova(mod.age)

fit.CA <- data.frame(predict(mod.age, type = "risk", se.fit=TRUE))
fit.CA$CA <- na.omit(real.df$CurrentAge)
colnames(fit.CA)[3] <- "CurrentAge"

predframe.CA <- with(fit.CA, data.frame(CurrentAge, hazard=fit, 
                                        lwr = fit - 1.96*se.fit, upr = fit + 1.96*se.fit))
pCA <- ggplot(predframe.CA, aes(CurrentAge, hazard)) +
  geom_line()
pCA + ylim(0,2) +
  xlab("Age") +
  ylab("Predicted Relative Hazard") + 
  theme_bw() +
  geom_hline(yintercept = 1.0, linetype = 2) +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  geom_ribbon(data = predframe.CA, aes(ymin=lwr, ymax=upr), alpha = 0.2) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 

#Number of helpers - this isn't good, use group size
#mod6 <- coxph(Surv(tstart, tstop, Status) ~ HelpNum, data = real.df)
#summary(mod6)
#anova(mod6)

mod7 <- coxph(Surv(tstart, tstop, Status) ~ GroupSize, data = new.realdf)
summary(mod7)
anova(mod7)

fit.gs <- data.frame(predict(mod7, type = "risk", se.fit=TRUE))
fit.gs$gs <- na.omit(real.df$GroupSize)
colnames(fit.gs)[3] <- "GroupSize"


predframe.gs <- with(fit.gs, data.frame(GroupSize, hazard=fit, 
        lwr = fit - 1.96*se.fit, upr = fit + 1.96*se.fit))
pGS <- ggplot(predframe.gs, aes(GroupSize, hazard)) +
  geom_line()
pGS + ylim(0,2) +
  xlab("Group Size") +
  ylab("Predicted Relative Hazard") + 
  theme_bw() +
  geom_hline(yintercept = 1.0, linetype = 2) +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  geom_ribbon(data = predframe.gs, aes(ymin=lwr, ymax=upr), alpha = 0.2) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 


test <-coxme(Surv(tstart, tstop, Status) ~ GroupSize + (1|FirstYr), data=new.realdf)
summary(test)
anova(mod7, test)

#terryr
test2 <- coxme(Surv(tstart, tstop, Status) ~ GroupSize + (1|TerrYr), 
               data = new.realdf)
summary(test2)
anova(mod7, test2)

test3 <- coxme(Surv(tstart, tstop, Status) ~ GroupSize + (1|NestYear),
               data = new.realdf)
summary(test3)
anova(mod7, test3)
anova(test,test2,test3)

#3 30 2017
#Just for fun

just4fun <- coxme(Surv(tstart, tstop, Status) ~ MinAgeFBreed +
                    SdOak + GroupSize + YrsExp + (1|NestYear), 
                  data = new.realdf)
summary(just4fun)

wait <- coxme(Surv(tstart, tstop, Status) ~ MinAgeFBreed + SdSize +
                     (1|NestYear), data=new.realdf)
summary(wait)
wait2 <- coxme(Surv(tstart, tstop, Status) ~ MinAgeFBreed + SdOak + 
                 (1|NestYear), data = new.realdf)
summary(wait2)

just4fun2 <- coxme(Surv(tstart, tstop, Status) ~ MinAgeFBreed + 
                     YrsExp + (1|NestYear), data = new.realdf)
summary(just4fun2)
just4fun3 <- coxph(Surv(tstart, tstop, Status) ~ MinAgeFBreed + 
                          YrsExp, data = new.realdf)
anova(just4fun3, just4fun2)

morefun <- coxme(Surv(tstart, tstop, Status) ~ YrsExp + GroupSize + 
                   (1|FirstYr) + (1|NestYear), data= new.realdf)
summary(morefun)
morefun2 <- coxme(Surv(tstart, tstop, Status) ~ GroupSize + (1|FirstYr) +
                    (1|NestYear), data = new.realdf)
summary(morefun2)
anova(morefun2, morefun)

hmm <- coxme(Surv(tstart, tstop, Status) ~ YrsExp + (1|FirstYr) + (1|NestYear), 
        data = new.realdf)
summary(hmm)

compare <- coxme(Surv(tstart, tstop, Status) ~ (1|FirstYr) + (1|NestYear), 
                 data = new.realdf)
compare2 <- coxme(Surv(tstart, tstop, Status) ~ (1|FirstYr), data = new.realdf)
compare3 <- coxme(Surv(tstart, tstop, Status) ~ GroupSize + (1|FirstYr), 
                  data = new.realdf)
anova(mod7, morefun2)
anova(compare2, compare)
anova(compare3, morefun2)

hmm2 <- coxme(Surv(tstart, tstop, Status) ~ MinAgeFBreed + (1|FirstYr) +
                (1|NestYear), data =new.realdf)
summary(hmm2) # large redution in LL (integrated to Null) - dev = 16.0

hmm3 <- coxme(Surv(tstart, tstop, Status) ~ MinAgeFBreed + GroupSize + 
                (1|FirstYr) + (1|NestYear), data = new.realdf)
summary(hmm3) #reduction from LL, dev= 16.0

###########################################################

#upload pedigree
jayped <- read.csv("Demo_Pedigree_2017.csv")
names(jayped)
jayped$SexIND[which(jayped$Sex=="F")] <- 2
jayped$SexIND[which(jayped$Sex=="M")] <- 1
jayped$SexIND[which(jayped$Sex=="")] <- 3
jayped[,11:12] <- NULL

colnames(real.df)[2] <-  "USFWBand"

#All birds with and without pedigree information
## I don't think this step is necessary because the pedigree needs to come 
# from a df with one record per subject 
#brp <- merge(x = real.df, y = jayped, intersect = "USFWBand", all.x =TRUE)
colnames(brd)[2] <- "USFWBand"
colnames(brd)[1] <- "JayID"

#N = 444, so these are all known age birds with known parentage
#Need to add the founders
#brp <- merge(x = brd, y = jayped, intersect = "USFWBand", all.x =TRUE)
#4 21 2017 - added a few birds
brp <- merge(x = brd.April15, y = jayped, intersect = "USFWBand", all.x =TRUE)


#at this point, everyone either has parent info or they don't
#If they don't, already coded correctly (at least as far as NA)
#Need to add birds in the mom and dad ID list that ARENT in the focal list 
#to the focal list
str(brp)
#At this point, the birds who have blanks for parents are founders
#They don't have parental info on the jayped from the database but they are
#in the focal ID list 

#Need the birds that are in breeder list but not focal id list
# take moms and dads, get rid of NAs and see who is not in the focal list?


brp$USFWBand <- as.character(brp$USFWBand)
brp$FBreeder <- as.character(brp$FBreeder)
brp$MBreeder <- as.character(brp$MBreeder)
brp$MUSFWBand <- as.character(brp$MUSFWBand)
brp$FUSFWBand <- as.character(brp$FUSFWBand)

ids.only <- data.frame(brp$JayID, brp$USFWBand, brp$FBreeder,
                       brp$FUSFWBand, brp$MBreeder, brp$MUSFWBand)
colnames(ids.only) <- c("JayID", "USFWBand", "MomID", "MomBand", "DadID", "DadBand")
work <- subset(ids.only, !(ids.only$MomID %in% ids.only$JayID))

fuckyou <- setdiff(ids.only$MomID, ids.only$JayID)
WHY <- data.frame(fuckyou)
colnames(WHY) <- "JayID"
mom.founders <- na.omit(WHY)
mom.founders$Sex <- "F"

ugh <- setdiff(ids.only$DadID, ids.only$JayID)
WHY2 <- data.frame(ugh)
colnames(WHY2) <- "JayID"
dad.founders <- na.omit(WHY2)
dad.founders$Sex <- "M"

founders <- rbind(mom.founders, dad.founders)
#####################################################################
#list of parents, some NAs
parents1 <- data.frame(brp$JayID, brp$FBreeder, brp$FUSFWBand, 
                       brp$MBreeder, brp$MUSFWBand)
colnames(parents1) <- c("JayID","FBreeder", "FUSFWBand", "MBreeder", "MUSFWBand")

parents1 <-parents1[complete.cases(parents1), ]
parents1$FBreeder <- as.character(parents1$FBreeder)
parents1$MBreeder <- as.character(parents1$MBreeder)
parents1$FUSFWBand <- as.character(parents1$FUSFWBand)
parents1$MUSFWBand <- as.character(parents1$MUSFWBand)
#4 21 2017 at this point, the birds with parental info (both parents known)
#are in parents 1, now I want the birds that are in the focal ID list
#marked by "JayID" who do not have parent info to be the founders
#somehow I'm losing them - see MIA dataframe

##3 15 2017 there's a problem here in that this is adding pairs
# 4 21 2017 adjustments to code to get IDS
foundersA <- subset(brp, !(brp$FUSFWBand %in% parents1$FUSFWBand))
foundsf <- subset(foundersA, foundersA$Sex == "F")
#foundsf[,3:4] <- NULL
foundsm <- subset(foundersA, foundersA$Sex == "M")
#foundsm[,1:2] <- NULL

names(foundsf)
names(foundsm)

foundersA$MUSFWBand[which(foundersA$MBreeder == "-SW_")] <- "4"
foundersA$MUSFWBand[which(foundersA$MBreeder == "6893")] <- "1"
foundersA$MUSFWBand[which(foundersA$MBreeder == "-LQ")] <- "2"
foundersA$MUSFWBand[which(foundersA$MBreeder == "RLC-U")] <- "3"  
foundersA$FUSFWBand[which(foundersA$FBreeder == "oZ-")] <- "7"
foundersA$FUSFWBand[which(foundersA$FBreeder == "K-GPC")] <- "6"

# 42 4017 - don't need much of the code in the chunk below anymore
###################################################################
#get rid of duplicates  -  a row adds for each time they parent 
#foundsf <- foundsf[!duplicated(foundsf),]
#foundsm <- foundsm[!duplicated(foundsm),]

#foundsm$MUSFWBand[which(foundsm$MBreeder == "-SW_")] <- "4"
#foundsm$MUSFWBand[which(foundsm$MBreeder == "6893")] <- "1"
#foundsm$MUSFWBand[which(foundsm$MBreeder == "-LQ")] <- "2"
#foundsm$MUSFWBand[which(foundsm$MBreeder == "RLC-U")] <- "3"                
#have to add in more founders
#there are birds in the parent list not in the focal list so those birds need
#to be added

#foundsf$FUSFWBand[which(foundsf$FBreeder == "oZ-")] <- "7"
#foundsf$FUSFWBand[which(foundsf$FBreeder == "K-GPC")] <- "6"
#add blank columns then can rbind to main df to make pedigree - brp

#colnames(foundsf) <- c("JayID", "USFWBand")
#colnames(foundsm) <- c("JayID", "USFWBand")
#names(foundsf)
#names(foundsm)
#foundsm$USFWBand[which(foundsm$JayID == "6584")] <- "8"
#foundsf <- foundsf[,c(2,3,1)]
#foundsm <- foundsm[,c(2,3,1)]

#add blank cols 
founders["USFWBand"] <- ""
founders["AgeKnown"] <- ""
founders["MinAgeFBr"] <- ""
founders["Year"] <- ""
founders["Terr"] <- ""
founders["Fbreed"] <- ""
founders["LastObsDate"] <- ""
founders["Yrs"] <- ""
founders["Censor"] <- ""
founders["FY"] <- ""
founders["TerrYr"] <- ""
founders["SexInd"] <- ""
founders["SexIND"] <- ""
founders["NatalNest"] <- ""
founders["FBreeder"] <- 0
founders["FUSFWBand"] <- 0
founders["MBreeder"] <- 0
founders["MUSFWBand"] <- 0


#brd.founders <- rbind(foundsf, foundsm)
#colnames(brd.founders)[1] <- "JayID"
#colnames(brd.founders)[2] <- "USFWBand"
#colnames(brd.founders)[1] <- "JayID"
#colnames(brd.founders)[2] <- "USFWBand"
##Need to add the missing sexes
#found.sex <- read.csv("found_sex.csv")
#str(found.sex)
#whoru3 <- subset(found.sex, !(found.sex$JayID %in% brd.founders$JayID))
#found.sex$JayID <- as.character(found.sex$JayID)
#found.sex$USFWBand <- as.character(found.sex$USFWBand)
#found.sex1 <- subset(found.sex, (found.sex$JayID %in% brd.founders$JayID))
#brd.founders <- cbind(brd.founders, found.sex1)

founders <- founders[,c(1,3,4,2,12,8,9,5,7,10,11,13,14,15,16,17,18,19,20)]

names(brp)
names(founders)
colnames(founders)[6] <- "BrdrDate"
colnames(founders)[5] <- "FirstYr"

######################################################################

founders$BrdrDate <- as.Date(founders$BrdrDate, format = "%m/%d/%Y")
founders$LastObsDate <- as.Date(founders$LastObsDate, format = "%m/%d/%Y")

brd.new <- rbind(brp, founders)
#fuck <- subset(brd.new, !(brd.new$MUSFWBand %in% brd.new$USFWBand))
#fuck <- fuck[!duplicated(fuck$JayID),]
#wtf <- subset(fuck, !(fuck$JayID %in% foundersA$JayID))

#brd.new <- rbind(brp, brd.founders)
#make sure Sex ind is correct 


###############################################
brd.new$MUSFWBand[which(brd.new$MBreeder == "RLC-U")] <- "3"
brd.new$MUSFWBand[which(brd.new$MBreeder == "-SW_")] <- "4"
brd.new$FUSFWBand[which(brd.new$FBreeder == "oZ-")] <- "7"
brd.new$FUSFWBand[which(brd.new$FBreeder == "K-GPC")] <- "6"
brd.new$MUSFWBand[which(brd.new$MBreeder == "6584")] <- "10"
brd.new$USFWBand[which(brd.new$JayID == "6584")] <- "10"
#########################################################

brd.new["SexIND"] <- NULL
brd.new$Sex <- as.character(brd.new$Sex)
brd.new$SexInd <- as.character(brd.new$SexInd)
brd.new$SexInd[which(brd.new$Sex == "F")] <- 2
brd.new$SexInd[which(brd.new$Sex == "M")] <- 1
brd.new$SexInd <- as.integer(brd.new$SexInd)
brd.new$Sex <- as.factor(brd.new$Sex)

brd.new[is.na(brd.new)] <- 0
#brd.new <- brd.new[!duplicated(brd.new$JayID),]

ids <- makefamid(brd.new$JayID, brd.new$MBreeder, brd.new$FBreeder)
famid <- cbind(ids)
colnames(famid) <- "famid"
brd.new <- cbind(brd.new, famid)

##Change NAs of founders to 0's
breeder.ped <- with(brd.new, pedigree(id = JayID, dadid=MBreeder, 
      momid = FBreeder, sex = SexInd, famid =famid, missid = 0))
ped1 <- breeder.ped[1]
print(ped1)
#plot(ped1)

bred.kin <- kinship(breeder.ped)

#Mixed effects models with pedigree

#There are some individuals in the data but not in the pedigree
#Figure out who they are 

#MIA <- subset(real.df, !(brd.new$JayID %in% real.df$JayID))
#MIA <- MIA[!duplicated(MIA$JayID), ]

#tryit <- merge(real.df, brd.new, by= c("JayID", "JayID"), all.x = TRUE)
#fuckfuck <- with(tryit, pedigree(id = JayID, dadid=MBreeder, 
                                 #momid=FBreeder, sex=SexInd, 
                                 #famid=famid, missid=0))


### June 27 2017 

#drop NA's so everything is comparable via AIC 
brd.j27 <- new.realdf %>% drop_na(17,38)
length(unique(brd.j27$USFWBand))
length(unique(brd.j27$JayID))

#territory
brd.tj27 <- new.realdf %>% drop_na(17,23,24,25,27,28,29,38)
length(unique(brd.tj27$USFWBand))
length(unique(brd.tj27$JayID))

## Survival object with intervals
j.ob <- Surv(brd.j27$tstart, brd.j27$tstop, brd.j27$Status)
j.tob <- Surv(brd.tj27$tstart, brd.tj27$tstop, brd.tj27$Status)

#June 27 new model code
j.mod1 <- coxme(j.ob ~ MinAgeFBreed + Sex + GroupSize + (1|JayID) +
            (1|NestYear.x), varlist=coxmeMlist(bred.kin), data = brd.j27)

j.mod2 <- coxme(j.ob ~ MinAgeFBreed + (1|JayID) + (1|NestYear.x), 
                varlist=coxmeMlist(bred.kin), data = brd.j27)
j.mod3 <- coxme(j.ob ~ MinAgeFBreed + GroupSize + (1|JayID) + (1|NestYear.x),
              varlist=coxmeMlist(bred.kin), data = brd.j27)

j.mod4 <- coxme(j.ob ~ MinAgeFBreed + GroupSize + CurrentAge + 
                (1|JayID) + (1|NestYear.x), varlist=coxmeMlist(bred.kin),
                data = brd.j27)
j.mod5 <- coxme(j.ob ~ MinAgeFBreed + GroupSize + FldgNum + 
              (1|JayID) + (1|NestYear.x), varlist = coxmeMlist(bred.kin),
              data = brd.j27)
summary(j.mod5)
j.mod6 <- coxme(j.ob ~ MinAgeFBreed + GroupSize + Sex + FldgNum + 
                  CurrentAge + (1|JayID) + (1|NestYear.x), 
            varlist=coxmeMlist(bred.kin), data = brd.j27)
j.mod7 <- coxme(j.ob ~ MinAgeFBreed + GroupSize + (1|NestYear.x),
                data = brd.j27)
j.mod8 <- coxme(j.ob ~ MinAgeFBreed + GroupSize + (1|JayID),
                varlist = coxmeMlist(bred.kin), data = brd.j27)
j.mod9 <- coxme(j.ob ~ MinAgeFBreed + GroupSize + Sex + (1|NestYear.x),
                data = brd.j27)
j.mod10 <- coxme(j.ob ~ MinAgeFBreed + GroupSize + CurrentAge +
                   (1|NestYear.x), data = brd.j27)
j.mod11 <- coxme(j.ob ~ MinAgeFBreed + GroupSize + FldgNum +
                  (1|NestYear.x), data = brd.j27)
j.mod12 <- coxme(j.ob ~ MinAgeFBreed + GroupSize + Sex + CurrentAge +
                   FldgNum + (1|NestYear.x), data = brd.j27)

Cand.modelsJ1 <- list()
Cand.modelsJ1[[1]] <- j.mod1
Cand.modelsJ1[[2]] <- j.mod2
Cand.modelsJ1[[3]] <- j.mod3
Cand.modelsJ1[[4]] <- j.mod4
Cand.modelsJ1[[5]] <- j.mod5
Cand.modelsJ1[[6]] <- j.mod6
Cand.modelsJ1[[7]] <- j.mod7
Cand.modelsJ1[[8]] <- j.mod8
Cand.modelsJ1[[9]] <- j.mod9
Cand.modelsJ1[[10]] <- j.mod10
Cand.modelsJ1[[11]] <- j.mod11
Cand.modelsJ1[[11]] <- j.mod12
ModnamesJ1 <- paste("Model", 1:length(Cand.modelsJ1), sep="")
tableJ1 <- aictab(Cand.modelsJ1, ModnamesJ1, sort = TRUE, second.ord=FALSE)
tableJ1

Cand.modelsJ2 <- list()
Cand.modelsJ2[[1]] <- j.mod3
Cand.modelsJ2[[2]] <- j.mod7
Cand.modelsJ2[[3]] <- j.mod8
ModnamesJ2 <- paste("Model", 1:length(Cand.modelsJ2), sep="")
tableJ2 <- aictab(Cand.modelsJ2, ModnamesJ2, sort = TRUE, second.ord=FALSE)
tableJ2


### Territory variables 
#j.tob
#brd.tj27
names(brd.tj27)

t.mod1 <- coxme(j.tob ~ MinAgeFBreed + GroupSize + SdSize + (1|JayID) +
            (1|NestYear.x), varlist = coxmeMlist(bred.kin), data = brd.tj27)

t.mod2 <- coxme(j.tob ~ (1|NestYear.x), data = brd.tj27)
summary(t.mod2)

t.mod3 <- coxme(j.tob ~ (1|JayID) + (1|NestYear.x), varlist=coxmeMlist(bred.kin),                data = brd.tj27)
summary(t.mod3)

t.mod4 <- coxme(j.tob ~ MinAgeFBreed + SdSize + (1|JayID) + (1|NestYear.x),
              varlist = coxmeMlist(bred.kin), data = brd.tj27)
summary(t.mod4)
t.mod5 <- coxme(j.tob ~ MinAgeFBreed + SdSize + (1|NestYear.x), data = brd.tj27)
summary(t.mod5)

t.mod6 <- coxme(j.tob ~ MinAgeFBreed + GroupSize + (1|NestYear.x), data = brd.tj27)
t.mod7 <- coxme(j.tob ~ MinAgeFBreed + SdSize + GroupSize + (1|NestYear.x), 
                data = brd.tj27)
t.mod8 <- coxme(j.tob ~ MinAgeFBreed + GroupSize*SdSize + (1|NestYear.x),
                data = brd.tj27)
summary(t.mod8)
t.mod9 <- coxme(j.tob ~ MinAgeFBreed + GroupSize*SdSize + (1|JayID) + 
                  (1|NestYear.x), data = brd.tj27)
t.mod10 <- coxme(j.tob ~ MinAgeFBreed + GroupSize + (1|JayID) + (1|NestYear.x),
               varlist = coxmeMlist(bred.kin), data = brd.tj27)

Cand.modelsJ3 <- list()
Cand.modelsJ3[[1]] <- t.mod1
Cand.modelsJ3[[2]] <- t.mod2
Cand.modelsJ3[[3]] <- t.mod3
Cand.modelsJ3[[4]] <- t.mod4
Cand.modelsJ3[[5]] <- t.mod5
Cand.modelsJ3[[6]] <- t.mod6
Cand.modelsJ3[[7]] <- t.mod7
Cand.modelsJ3[[8]] <- t.mod8
Cand.modelsJ3[[9]] <- t.mod9
Cand.modelsJ3[[10]] <- t.mod10
ModnamesJ3 <- paste("Model", 1:length(Cand.modelsJ3), sep="")
tableJ3<- aictab(Cand.modelsJ3, ModnamesJ3, sort = TRUE, second.ord=FALSE)
tableJ3


fit.tsize <- data.frame(predict.coxme(t.mod4, type = "risk"))
colnames(fit.tsize) <- "Fit"
fit.tsize$SdSize <- brd.tj27$SdSize

predframe.ts <- with(fit.tsize, data.frame(SdSize, hazard=Fit))
pGS <- ggplot(predframe.ts, aes(SdSize, hazard)) +
  geom_point()
pGS + ylim(0,2) +
  xlab("Group Size") +
  ylab("Predicted Relative Hazard") + 
  theme_bw() +
  geom_hline(yintercept = 1.0, linetype = 2) +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  geom_ribbon(data = predframe.gs, aes(ymin=lwr, ymax=upr), alpha = 0.2) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) 


##########################################################################

###
#Old as of June 27 2017

#yes I do need to build a new ped
ped.test1 <- coxme(Surv(tstart, tstop, Status) ~ (1|JayID), data = new.realdf, 
                   varlist = coxmeMlist(bred.kin))
summary(ped.test1)

mx.mod2 <- coxme(Surv(tstart, tstop, Status) ~ MinAgeFBreed +
          YrsExp + GroupSize + (1|JayID) + (1|NestYear.x), data = new.realdf, 
          varlist = coxmeMlist(bred.kin))
summary(mx.mod2)
mx.mod3 <- coxme(Surv(tstart, tstop, Status) ~ MinAgeFBreed +
            (1|JayID), data = real.df, varlist = coxmeMlist(bred.kin))
summary(mx.mod3)
mx.mod3B <- coxme(Surv(tstart, tstop, Status) ~ MinAgeFBreed + 
                    (1|JayID) + (1|NestYear), data = real.df, 
                  varlist = coxmeMlist(bred.kin))
summary(mx.mod3B)
anova(mx.mod3, mx.mod3B)

mx.mod4 <- coxme(Surv(tstart, tstop, Status) ~ MinAgeFBreed + GroupSize +
            (1|JayID), data = real.df, varlist=coxmeMlist(bred.kin))
summary(mx.mod4)

mx.mod5 <- coxme(Surv(tstart, tstop, Status) ~ MinAgeFBreed + GroupSize +
                    (1|JayID) + (1|NestYear), data = real.df,
                  varlist=coxmeMlist(bred.kin))

summary(mx.mod5)
anova(mx.mod4, mx.mod5) 

mx.mod6 <- coxme(Surv(tstart, tstop, Status) ~ MinAgeFBreed + GroupSize + 
                   (1|JayID) + (1|NestYear) + (1|FirstYr), data = real.df, 
                 varlist =coxmeMlist(bred.kin))
anova(mx.mod5, mx.mod6)
mx.mod7 <- coxme(Surv(tstart, tstop, Status) ~ MinAgeFBreed + GroupSize + 
                   (1|JayID) + (1|FirstYr), data = real.df, 
                 varlist = coxmeMlist(bred.kin))
summary(mx.mod7)
anova(mx.mod7, mx.mod5)

mx.mod8 <- coxme(Surv(tstart, tstop, Status) ~ MinAgeFBreed + (1|JayID) +
                   (1|NestYear), data = real.df, varlist=coxmeMlist(bred.kin))
summary(mx.mod8)
summary(mx.mod5)

##################################################################
#all.new <- merge(all.info, newdata)
#str(all.new)

#length(unique(all.new$band))

#all.new["time1"] <- ""
#all.new["time2"] <- ""

#write.csv(all.new, file = "allbirds.csv")
#Change dates back to numeric for survival object

#all.new$BrdrDate <- as.numeric(all.new$BrdrDate)
#all.new$LastObsDate <- as.numeric(all.new$LastObsDate)

#remove records with NA for lastobsdate and Jay ID 
#all.new <- subset(all.new, !all.new$JayID == "NA" & !all.new$LastObsDate == "NA")

#str(all.new)
#Have to change the censorship indicators for each interval 

#length(unique(all.new$JayID))
#length(unique(all.new$band))
#One repeat band 
#reprows <- ddply(all.new, .(all.new$JayID), nrow)

#############################################################################

#####FUCK YEAH I had the insight to save this little peice of code 
#Merge data frames to get tstart and tstop
#This data frame contains each nest year for each bird, should include
#all years lived for each bird
ints <- data.frame(birds$band, birds$NestYear, birds$tstart, birds$tstop)
colnames(ints) <- c("band", "NestYear", "tstart", "tstop")

#need to add those missing nest years
testing123 <- merge(ints, all.new, intersect = "band")
#this works but I need the missing years, as in the ones there is no data
#look into the left and right merging ideas -that could be the way 

#####################################################################

# Block 4 - Time varying Models 
# Data is in long format, I think I need a wide format as well for tmerge
cut.points <- unique(all.new$Years[all.new$Censor==1])
surv2 <- survSplit(data=all.new, cut = cut.points, end = "time",
                   start = "time0", event = "death")


test2 <- 

test <- reshape(all.new, 
    varying = c("YrsExp", "CurrentAge", "Terr","OakScrub",
                "Terrsize", "TSF", "Help", "HelpNum","FldgNum"),
    idvar = c("JayID", "band"), timevar="TerrYr",
    direction = "wide")



#####################################################################
#Example from 2014 J. Stat Software

id <- c(1,2,3,4,5,6)
time <- c(1,4,7,10,12,13)
death <- c(1,0,1,1,0,1)
age <- c(20,21,19,22,20,24)
female <- c(0,1,0,1,0,1)

SURV <- data.frame(id, time, death, age, female)
cut.points <- unique(SURV$time[SURV$death ==1])
SURV2 <- survSplit(Surv(time, death) ~ id, data = SURV , cut = cut.points, end = "time", start = "time0", event = "death")

model.1 <- coxph(Surv(time, death) ~ female + age, data = SURV)
#covariates
covs <- data.frame(age = 21, female = 0)
summary(survfit(model.1), newdata = covs, type = "aalen")

