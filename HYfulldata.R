## Jan 18 2017

# Hatch year birds with full data set 1981-2015 
# HatchYear.R is code and analyses for the subsample of 1999-2015 when birds
# were sexed using blood samples 

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
hybrd$Censor[which(hybrd$Years >= 1)]<-0

yrlg.df <- subset(hybrd, hybrd$Years > 0 & hybrd$Days > 0)
yrlg.df$Censor[which(yrlg.df$LastObsDate == "2016-04-12")] <- 0

#Jay ID not unique but the band number is 
#yrlg.df[,2] <- NULL

yrlg.df["Expr1003"] <- NULL
yrlg.df["Day11"] <- yrlg.df$MeasDate - yrlg.df$HatchDate
yrlg.df <- subset(yrlg.df, yrlg.df$Day11 == 11)

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

yrlg.df1$FldgDate <- as.numeric(yrlg.df1$FldgDate)
yrlg.df1$LastObsDate <- as.numeric(yrlg.df1$LastObsDate)

yrlg.df1$Years <- as.numeric(yrlg.df1$Years)
yrlg.df1$Days <- as.numeric(yrlg.df1$Days)

yrlg.df1 <- na.omit(yrlg.df1)

str(yrlg.df1)

# 3 3 2017 leave these variables as integers for now, can change later
# I *think* helper number is a count so it should be an integer and not
# a category
# need this data type for the simPH program to work 
#yrlg.df1$Help <- as.factor(yrlg.df1$Help)
#yrlg.df1$HelpNum <- as.factor(yrlg.df1$HelpNum)

length(unique(yrlg.df1$USFWBand))
#2323

#basic survival object and KM plots
hy.ob <- Surv(yrlg.df1$Years, yrlg.df1$Censor, type = c('right'))
hy.fit <- survfit(hy.ob ~ 1, conf.type = "log-log")
summary(hy.fit)
hy.fit
str(summary(hy.fit))
plot(hy.fit$time, hy.fit$surv, xlim = c(0,1), type = "l")
lines(hy.fit$time, hy.fit$upper, lty = 2)
lines(hy.fit$time, hy.fit$lower, lty = 2)

yr.fit <- survfit(hy.ob ~ yrlg.df1$Year, conf.type="log-log")
yr.fit
summary(yr.fit)
str(summary(yr.fit))
plot(yr.fit, log="y", xlim = c(0,1))
plot(yr.fit$time, yr.fit$sur, xlim = c(0,1), type = "l")

mean(yr.fit$n)
sd(yr.fit$n)
range(yr.fit$n)
median(yr.fit$n)
plot(yr.fit$n)

kmcurves <- plot(hy.fit, xlab="Year", log = "y", 
               ylab = "Survival", main = "HY birds",
               ylim = c(0.3,1), xlim = c(0,1))
library(GGally)
surv.plot1 <- ggsurv(hy.fit, plot.cens = FALSE, xlab = "Time (Years)",
                     ylab = "Proportion Surviving") +
  xlim(c(0,1)) +
  ylim(c(0,1))
surv.plot1

ggsurv(yr.fit) +
  xlim(c(0,1))


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

#table of deaths, p and q for birds from KM estimates 
HYsurv <- read.csv("HYlifetale.csv")
range(HYsurv$p)
mean(HYsurv$p)

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

cat("Mean nestling mass= ", mean(yrlg.df1$Mass))
cat("Standard deviation of mass", sd(yrlg.df1$Mass))
cat("Median of mass", median(yrlg.df1$Mass))

#Fledglings 
names(HYsurv)
fldgs.plot <- ggplot(HYsurv, aes(Year, N)) +
  geom_bar(width = 0.75, stat = "identity", color = "black", fill = "red") +
  theme(axis.text.x = element_text(angle = 45)) + 
  xlab("Year") +
  ylab("Count of Individuals") +
  geom_bar(aes(Year, D), stat = "identity", width = 0.75,
           color= "black", fill = "blue")
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
#plot(ped1)
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

##################################################################
## 3 27 2017 Add in parent age 
##need parent ages to add to yrlg.df1
ages <- na.omit(new.df)
ages <- ages[,c(1,2,3,4,6,21,22,23,24)]
ages.m <- ages[,c(1,2,3,4,5,6,7)]
ages.d <- ages[,c(1,2,3,4,5,8,9)]


parent.age <- all.brd[,c(1,2,6,10,11,12,18)]
mom.age <- subset(parent.age, sex == "F")
colnames(mom.age)[1] <- "FBreeder"
dad.age <- subset(parent.age, sex == "M")
colnames(dad.age)[1] <- "MBreeder"

new.d <- merge(ages.d, dad.age, by = c("TerrYr","MBreeder"), all.x=TRUE)
new.m <- merge(ages.m, mom.age, by = c("TerrYr","FBreeder"), all.x=TRUE)

new.m <- new.m[,c(3,4,5,6,1,7,9,10,11)]
colnames(new.m)[6] <- "momid"
colnames(new.m)[8] <- "momexp"
colnames(new.m)[9] <- "momage"

new.d <- new.d[,c(3,4,5,6,1,7,9,10,11)]
colnames(new.d)[6] <- "dadid"
colnames(new.d)[8] <- "dadexp"
colnames(new.d)[9] <- "dadage"

new.d <- new.d[,c(1,2,5,8,9)]
new.m <- new.m[,c(1,2,5,8,9)]

new1 <- merge(yrlg.df1, new.d, by = c("USFWBand","JayID"), all.x=TRUE)
new2 <- merge(new1, new.m, by = c("USFWBand","JayID"), all.x = TRUE)

yrlg.df1 <- new2
yrlg.df1[,20] <- NULL
yrlg.df1[,22] <- NULL
colnames(yrlg.df1)[4] <- "TerrYr"

#COunt experience level per age 
agexp.count <- count(yrlg.df1, c("momage", "momexp"))
exp.count  <- count(yrlg.df1, "momexp")
counts <- count(yrlg.df1, c("momexp", "momage"))

#how to plot...?
mexp0 <- subset(agexp.count, momexp == 0)
mexp0 <- na.omit(mexp0)
qplot(momage, freq, data = mexp0, geom="bar")
ggplot(mexp0, aes(x=momage, y=freq)) +
geom_bar(stat = "identity")


counts <- na.omit(counts)
ggplot(counts, aes(x=momage, y=freq), fill=momexp) +
geom_bar(colour="blue", stat = "identity")

a <- ggplot(counts, aes(x=momage, y=freq))
a + geom_point(aes(colour=factor(momexp), y=freq))

b <- ggplot(counts, aes(x=momexp, y=freq))
b + geom_point(aes(colour=factor(momage), y=freq))

b+ geom_bar(aes(fill = factor(momage), y=freq), stat="identity")
  

#A few duplicates because I matched by TerrYr, there were some breakups
#of breeding pairs
#reps1 <- ddply(new1, .(new1$USFWBand), nrow)
#reps2 <- ddply(new2, .(new2$USFWBand), nrow)


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

library(corrplot)
covars1 <- yrlg.df1[,c(24,25,26)]
covars1["oakfrac"] <- covars1$scrb.count/covars1$TerrSize
covars1["tsffrac"] <- covars1$tsf.count/covars1$TerrSize
colnames(covars1) <- c("OakScrub", "TerrSize", "TSF", "RelOak", "RelTSF")
corrs1 <- cor(covars1, use="complete.obs")
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

loadings <- data.frame(vars.pca$rotation, 
                       .names = row.names(vars.pca$rotation))
p + geom_text(data = loadings, mapping = aes(x = PC1, y = PC2, 
                                             label = .names, colour = .names)) + 
  coord_fixed(ratio=1) +
  labs(x = "PCI", y = "PC2")

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
flg.sub1 <- flg.sub %>% drop_na(24,25,26,27,28,29,30,31,32)


#should be 1529
flg.sub1 <- cbind(flg.sub1, pc.scores)
plot(flg.sub1$PC1, flg.sub1$stdsize)
plot(flg.sub1$PC2, flg.sub1$stdscr)
plot(flg.sub1$PC2, flg.sub1$stdtsf)

str(flg.sub1)
names(flg.sub1)
flg.sub2 <- flg.sub1 %>% drop_na(22,23)

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
yrlg.mob <- Surv(yrlg.df1$Years, yrlg.df1$Censor, type = c('right'))

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

#Year as a factor
cox.year <- coxph(yrlg.mob ~ Year, data = yrlg.df1)
cox.year
anova(cox.year)

#Compare the first 3, brood size and fledge # add nothing 
d,ev.compare1 <- anova(hatch.num, flg.num, cox.year, test="Chisq")
dev.compare1

#Parent age and experience
#just a simple plot
plot(yrlg.df1$momage, yrlg.df1$Years)
plot(yrlg.df1$Years, yrlg.df1$momage)
plot(yrlg.df1$momexp, yrlg.df1$Years)

momage <- coxph(yrlg.mob ~ momage, data = yrlg.df1)
summary(momage)
plot(cox.zph(momage))
dadage <- coxph(yrlg.mob ~ dadage, data = yrlg.df1)
plot(cox.zph(dadage))
summary(dadage)
anova(momage)
anova(dadage)

##polynomials
momage.py <- coxph(yrlg.mob ~ momage + I(momage^2), data = yrlg.df1)
summary(momage.py)
anova(momage, momage.py)

momexp.mod <- coxph(yrlg.mob ~ momexp, data = yrlg.df1)
summary(momexp.mod)
anova(momexp.mod)


dadexp.mod <- coxph(yrlg.mob ~ dadexp, data = yrlg.df1)
summary(dadexp.mod)
anova(dadexp.mod)

sim.moage <- coxsimLinear(momage, b = "momage", Xj = 1:15)
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
temp.mom <- temp.mom %>% drop_na(22,23)
temp.ob <- Surv(temp.mom$Years, temp.mom$Censor, type = c('right'))
temp.cox <- coxph(temp.ob ~ momage, data = temp.mom)
summary(temp.cox)
temp.cox1 <- coxph(temp.ob ~ momexp, data = temp.mom)
summary(temp.cox1)
temp.cox2 <- coxph(temp.ob ~ momage + momexp, data = temp.mom)
summary(temp.cox2)

temp.cox3 <- coxph(temp.ob ~ momage + momexp + momage*momexp, 
                   data = temp.mom)
summary(temp.cox3)

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


## 3 29 2017 Getting closer with the plot below 
p.mom1 <- ggplot(trts, aes(momage, Risk)) + geom_point() 
p.mom1 + facet_grid(.~momexp)

ggplot(trts, aes(x=momage,y=momexp)) +
  geom_point(aes(size=Risk))

ggplot(trts, aes(momage, momexp)) + 
  geom_point(aes(colour=Risk))

qplot(momage, Risk, data = trts, colour = momexp)
qplot(momage, momexp, data=trts, colour= Risk)
##only bottom half of plot makes sense
##email Gordon about this 

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
summary(mass.full)
anova(mass.full)
#plot of scaled Schoenfield residuals to check PH assumption 
plot(cox.zph(mass.full))
#simPH method
sim.mass <- coxsimLinear(mass.full, b = "Mass", qi = "Relative Hazard", 
                         Xj = c(15, 65), spin = TRUE)
simGG(sim.mass,type="ribbons", alpha = 0.35)
#predicted values 
preds <- predict(mass.full, type="risk")

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
helpplot2 <- simGG(sim.help2, xlab = "Helper Number", alpha = 0.4)
plotb <- helpplot2 + theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
                  scale_x_continuous(breaks =c(0,1,2,3,4,5,6),
                  labels=c("0","1","2","3","4","5","6")) +
                  theme(axis.title.y = element_blank())
plotb

#adjust scaling
multiplot(plota, plotb, cols = 2)
#This figure is close to what I want, needs some aes. adj
sim.help3 <- coxsimLinear(cx.group, b ="HelpNum",
                  qi = "Hazard Ratio", Xj = c(0,6))
simGG(sim.help3, alpha = 0.4, xlab = "Helper Number")
#try nonlinear fit
#predict
predvals1 <- predict(cx.group, type="risk")
plot(yrlg.df1$HelpNum, predvals1,
     xlab = "Natal Territory Helper Number", ylab="Predicted Risk Scores")


#Mass and Helper (0,1)
hlp.mass <- coxph(yrlg.mob ~ Mass + Help, data = yrlg.df1)
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
ggplot(data = temp, aes(x = Mass, y = Risk, group=Help, colour=Help)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("darkgrey", "black")) 
  
ggplot(data = temp, aes(Mass, Risk, group=Help, colour=Help)) +
  geom_point(size = 1) +
  scale_color_manual(values= c("darkgrey", "black"))

ggplot(data = temp, aes(Mass, Risk, group=Help, colour=Help)) +
  geom_point(size = 1) +
  scale_color_manual(values= c("blue", "darkorange"))


#Mass and Helper number
hlpn.mass <- coxph(yrlg.mob ~ Mass + HelpNum, data = yrlg.df1)
summary(hlpn.mass)
anova(hlpn.mass)
predvals3 <- predict(hlpn.mass, type = "risk")
plot(yrlg.df1$HelpNum, predvals3)


#Plot for showing the predicted values of hazard ratio from model of
#risk as a function of mass and helper number 
treats <- expand.grid(HelpNum = levels(tmp$HelpNum), 
                      Mass = tmp$Mass)
temp2 <- treats
temp2$Risk <- predict(hlpn.mass, temp2, type = "risk")
qplot(Mass, Risk, data = temp2, colour=HelpNum)


#interaction between helper number and mass
hlpmni <- coxph(yrlg.mob ~ Mass*HelpNum, data = yrlg.df1)
summary(hlpmni)
anova(hlpmni)
predvals4 <- predict(hlpmni, type="risk")
plot(yrlg.df1$Mass, predvals4)
plot(yrlg.df1$HelpNum, predvals4)

#interaction between helpers and mass, does this even make biological sense?
hlpmi <- coxph(yrlg.mob ~ Mass*Help, data=yrlg.df1)
summary(hlpmi)
anova(hlpmi)
anova(mass.full, hlpmi)

#mom age 
mha <- coxph(yrlg.mob ~ Mass + HelpNum + momage, data= yrlg.df1)
summary(mha)

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

########################################################################
#Territory variables 
#use yrlg.tob

#PC1 and 2
cx.pc <- coxph(yrlg.tob ~ PC1 + PC2, data = flg.sub1)
summary(cx.pc)
cx.pc1 <- coxph(yrlg.tob ~ PC1, data = flg.sub1)
summary(cx.pc1)

anova(cx.pc1, cx.pc)

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

#ID as random effect, no fixed effect - KEEP 3/21/2017
kin1 <- coxme(yrlg.mob ~ (1|USFWBand), data = yrlg.df1, 
              varlist = coxmeMlist(kins))
summary(kin1)
#Add mass as fixed effect
kin.mass <- coxme(yrlg.mob ~ Mass + (1|USFWBand), data = yrlg.df1,
                 varlist = coxmeMlist(kins))
summary(kin.mass)

## Add mom age
kin.mage <- coxme(yrlg.mob ~ Mass + momage + (1|USFWBand) + (1|Year), 
                  data = yrlg.df1, varlist = coxmeMlist(kins))
summary(kin.mage)

kin.mexp <- coxme(yrlg.mob ~ Mass + momexp + (1|USFWBand) + (1|Year),
                  data = yrlg.df1, varlist = coxmeMlist(kins))
summary(kin.mexp)

kin.mah <- coxme(yrlg.mob ~ Mass + momage + HelpNum + (1|USFWBand) +
                   (1|Year), data = yrlg.df1, varlist = coxmeMlist(kins))
summary(kin.mah)
anova(kin.mage, kin.mah)

kin.maep <- coxme(yrlg.mob ~ Mass + momage + momexp + (1|USFWBand) +
                    (1|Year), data = yrlg.df1, varlist = coxmeMlist(kins))
summary(kin.maep)
#asses fit of adding mass 
anova(kin1, kin.mass, test = "chisq")
#model with mass better

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