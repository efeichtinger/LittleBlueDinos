# Hatch Year/First Year Birds - subsample to test for sex effects 

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
#all birds from 1981 on 
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

colnames(bird.df)[7] <- "Mass"

str(bird.df)

#Add censorship column, 1 if dead before 1 yr, 0 if yr > 1 or 4/12/12016
bird.df["Censor"] <- 1

#this is the line where I could change to April census + 1
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
##Subset - I kept the name to keep all model code as is
yrlg.df <- subset(yrlg.df, yrlg.df$FldgDate > "1999-01-01")
length(unique(yrlg.df$USFWBand))

#1981 on
#full.df <- yrlg.df # go to line 481

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
yrlg.ob <- Surv(yrlg.df$Days, yrlg.df$Censor, type = c('right'))
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

###########################################################################
prop <- c(51,36,30,56,23,65,55,71,28,96,39,72,113,11,41,65,48)

years <- seq(as.Date("1999/1/1"), by = "year", length.out=17)
yr <- as.numeric(format(years, "%Y"))
year <- as.vector(yr)
year <- cbind(year)

life.table <- data.frame(fit.yr$n, prop)
life.table["Year"] <- year
colnames(life.table)[1] <- "N"
colnames(life.table)[2] <- "Nsurv"
colnames(life.table)[3] <- "Year"

life.table["p"] <- life.table$Nsurv/life.table$N

#Make a nice bar graph of p, proportion of survivors each year fledglings 
ggplot(life.table, aes(year, p)) +
  geom_bar(stat="identity")
########################################################################

#Sex -  actually looks alright, shows that the CI's for males and females 
#overlapp
sex.fit <- survfit(yrlg.ob ~ yrlg.df$Sex, conf.type = "log-log")
plot.sex <- plot(sex.fit, conf.int = TRUE,col = c("maroon2","blue3"), xlab = "Time (years)", log = "y",
                   ylim = c(0.4,1), xlim=c(0,1), ylab = "Cumulative Survival")
legend("topright", c("Females","Males"), col=c("maroon2","blue3"),
       lty = c(1,1),lwd=1)

plot.sexbw <- plot(sex.fit, conf.int = FALSE, col=c("black","darkgrey"),
                  ylab = "Survival", xlab = "Time (years)", lty = c(5,1),
                  lwd = 2, xlim = c(0,1), ylim = c(0.4,1))
legend("topright", c("Females", "Males"), col=c("black","darkgrey"),
       lty = c(5, 1), lwd = 2)

# 3 7 2017 
plot.sexci <- plot(sex.fit, conf.int = TRUE,col = c("black","darkgrey"), 
          xlab = "Time (years)", log = "y",
          lty = c(1,5), lwd = c(2,2), ylim = c(0.4,1), xlim=c(0,1), 
          ylab = "Cumulative Survival")
legend("topright", c("Females","Males"), col=c("black","darkgrey"),
       lty = c(1,5),lwd=2)


#### 3 7 2017 - use this plot or at least start here
#### unless all of plots need to be made in ggplot2... 
par(mfrow = c(1,2), mai = c(0.9, 0.8, 0.1, 0.2))

all.fit <- plot(my.fit, xlab = "Time (Days)", ylab = "Cumulative Survival",conf.int=TRUE,
                lwd = 2, ylim = c(0,1),xlim=c(0,365)) 

plot.test <- plot(sex.fit, conf.in = TRUE,
          col = c("black","black","black","gray50","gray50","gray50"),
          lty = c(1,3,3,1,3,3), lwd = 2,
          xlab = "Time (Days)", ylim = c(0,1), 
          xlim = c(0,365))
legend("topright", inset = 0.0, c("Females", "Males"), col = c("black", "gray50"), 
       lty = c(1,1), lwd = 2, cex = 0.75)




library(GGally)

legend_title <- "Sex"
hy.sex <- ggsurv(sex.fit, plot.cens = FALSE, CI = TRUE, 
                 surv.col = c("green","black"),
                 xlab = "Time (Years)",
                 ylab = "Proportion Surviving") 
hy.sex + xlim(c(0,1)) + ylim(c(0,1)) 

names(sex.fit)

hy.sex2 <- ggsurv(sex.fit)
hy.sex2 + 
  xlim(c(0,1)) + 
  ylim(c(0,1)) + 
  ggplot2::scale_color_discrete(
    name = 'Sex',
    breaks = c(1,2),
    labels = c('Male', 'Female')
  )

###########################################################################
ggsurv <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
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


########################################################################
##ggsurv function from the internet Tal Galili

beep <- ggsurv(sex.fit)
beep + xlim(0, 1) 

# 3 7 2017 - make a nice plot with different line types for CI lines 
# and lines for males and females 
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
summary(cox.sex)
cox.sex
anova(cox.sex)
svvc <- survfit(cox.sex)
plot(svvc, xlim = c(0,1), xlab = "Time (years)", ylab = "Cumulative Survival")

cox.zph(cox.sex)
plot(cox.zph(cox.sex))

#Weibull
Wsex <- survreg(yrlg.ob ~ yrlg.df$Sex, dist = 'weibull')
summary(Wsex)

cumhaz <- basehaz(cox.sex)
## is this correct? Females with a higher cumulative hazard
plot(cumhaz$time, cumhaz$hazard, type = "l", lwd = 3, xlab = "Time", ylab = "Cumulative Hazard", xlim=c(0,2))
lines(cumhaz$time, exp(0.90910)*cumhaz$hazard, col="blue", lwd = 3)

yrlg.df["SexInd"] <- 0
yrlg.df$SexInd[which(yrlg.df$Sex=="M")] <-1
cox.sex2 <- coxph(yrlg.ob ~ SexInd, data=yrlg.df)
sim1 <- coxsimLinear(cox.sex, b= "SexInd", qi="Relative Hazard", Xj=c(0,1))

cumhaz1 <- basehaz(cox.sex)
## is this correct? Females with a higher cumulative hazard
plot(cumhaz1$time, cumhaz1$hazard, main = "Hazard Rates")


cox.mass <- coxph(yrlg.ob ~ Mass, data = yrlg.df)
cox.mass
anova(cox.mass)
plot(cox.zph(cox.mass))

mass.sex <- coxph(yrlg.ob ~ Sex + Mass, data = yrlg.df)
mass.sex
anova(mass.sex)
plot(cox.zph(mass.sex))

ms.sexint <- coxph(yrlg.ob ~ Sex + Mass + Sex:Mass, data = yrlg.df)
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
anova(flg.num)

hatch.num <- coxph(yrlg.ob ~ HatchNum, data=yrlg.df)
hatch.num

fld.hatch <- coxph(yrlg.ob ~ HatchNum + FldgNum, data=yrlg.df)
fld.hatch

fhmass <- coxph(yrlg.ob ~ HatchNum + FldgNum + Mass, data=yrlg.df)
fhmass

#Mass and flg num interaction
mass.fl <- coxph(yrlg.ob ~ FldgNum + Mass + FldgNum:Mass, data = yrlg.df)
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
anova(cox.year, cox3)

#### Mixed effects models
#NestYear/Cohort Year 
mm.year <- coxme(yrlg.ob ~ Mass + (1|NestYear), data = yrlg.df)
mm.year

mm.nest <- coxme(yrlg.ob ~ Mass + (1|NatalNest), data = yrlg.df)
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
              stdscr:TerrSize, data = yrlg.df)
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

#############################################################################
str(full.df)


###############################################################################

pedAll<-pedigree(id=sample.ped$id, 
                 dadid=sample.ped$father, momid=sample.ped$mother, 
                 sex=sample.ped$sex, famid=sample.ped$ped)

demog.ped <- read.csv("Erin_Ped_1999.csv")

pedjays <- pedigree(id=demog.ped$Jays, dadid=demog.ped$MBreeder,
                    momid=demog.ped$FBreeder, sex = demog.ped$Sex)
