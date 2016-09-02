#### Breeders
#### July 21 2016

library(survival)
library(ggplot2)
library(car)
library(coxme)
library(kinship2)
library(SurvRegCensCov)
library(parfm)

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
#does not work - snytax to account for date 
brd <- subset(brd, brd$FY >= "1981")


#Change dates back to numeric for surv object 
brd$Fbreed <- as.numeric(brd$Fbreed)
brd$LastObsDate <- as.numeric(brd$LastObsDate)
brd$Yrs <- as.numeric(brd$Yrs)


#summary stats for breeder lifespan, i.e. time spent as a breeder
#not total lifespan 
mean(brd$Yrs)
var(brd$Yrs)
sd(brd$Yrs)

brd.ob <- Surv(brd$Yrs, brd$Censor,type= c('right'))
brd.fit <- survfit(brd.ob ~ 1, conf.type = "log-log")

# adjust x axis it's breeding span 
kmplot <- plot(brd.fit, xlab="Breeding time span (years)", log = "y", 
               ylab = "Cumulative Survival (log)", main = "Breeders",
               ylim = c(0.01,1), xlim = c(0,15))

brd.sex <- survfit(brd.ob ~ brd$Sex, conf.type="log-log")
sxplot <- plot(brd.sex, col = c("deepskyblue3","darkorange3"), 
               xlab = "Breeding time span (years)", log="y", 
               ylab = "Cumulative Survival (log)", main = "Breeders",
               ylim = c(0.01,1), xlim = c(0,15))
legend("topright", c("Females","Males"), col=c("deepskyblue3","darkorange3"),
       lty = c(1,1),lwd=1)

cx1 <- coxph(brd.ob ~ Sex, data = brd)
summary(cx1)
cox.zph(cx1)
plot(cox.zph(cx1))


####################################################################
# Block 2 - Input of Vegetation, fire, and terr size data
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


## Merge data frames by terryr - 4 cols terryr, scrub, tsf, terr cell count






#####################################################################
# Block 3 - Input of data for each year bred 
## Read in file for breeders all years, multiple records per individual
## Time varying Cox Models 


#####################################################################

# Block 4 - Models 