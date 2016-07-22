#### Breeders
#### July 21 2016

library(survival)

## Read in file of breeders - known and unknown age

brd <- read.csv("Breeders.csv")

#remove duplicate rows 
brd <- brd[!duplicated(brd),]

colnames(brd)[1] <- "ID"
colnames(brd)[2] <- "Band"
colnames(brd)[3] <- "FY"
colnames(brd)[4] <- "Fbreed"


brd <- brd[,c(1,2,8,7,6,3,4,5)]

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

#Change dates back to numeric for surv object 
brd$Fbreed <- as.numeric(brd$Fbreed)
brd$LastObsDate <- as.numeric(brd$LastObsDate)
brd$Yrs <- as.numeric(brd$Yrs)


brd.ob <- Surv(brd$Yrs, brd$Censor,type= c('right'))
brd.fit <- survfit(brd.ob ~ 1, conf.type = "log-log")
kmplot <- plot(brd.fit, xlab="Time (years)", log = "y", 
               ylab = "Cumulative Survival (log)", main = "Breeders",
               ylim = c(0.01,1), xlim = c(0,15))
