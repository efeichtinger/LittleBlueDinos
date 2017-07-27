##Survival analyses 
##Data set has multiple rows per individual if they bred more than once

#use this on laptop
setwd("C:/Users/Erin/LittleBlueDinos")
#use this on school desktop
setwd("C:/Users/efeichtinger/LittleBlueDinos")

library(survival)
library(ggplot2)

data <- read.csv("Erin_Breeders_All_Years.csv")

#remove duplicates  - for years where there was more than one nest in a year
breed <- data[!duplicated(data),]
str(breed)

colnames(breed)[1] <- "ID"
colnames(breed)[2] <- "Band"
colnames(breed)[4] <- "FirstDateBred"

#convert dates to date format
breed$FirstDateBred <- as.Date(breed$FirstDateBred, format = "%m/%d/%Y")
breed$LastObsDate <- as.Date(breed$LastObsDate, format = "%m/%d/%Y")

#subtract dates to get number of days
date.diff<- breed$LastObsDate-breed$FirstDateBred

#and survival period in years, account for leap year 
breed["Yrs"] <- as.numeric(date.diff/365.25)

breed$FirstYr <- as.factor(breed$FirstYr)

#very important piece of code for the model to work properly, remove any 
#weird entries like birds that have negative years of experience or a negative
#survival interval 
breed <- subset(breed, breed$Yrs > 0 & breed$YrsExp >= 0)

#add column for censorship status, in survival package - 0=alive, 1=dead
breed["censorship"] <- 1
#If last observed date = 10/14/2015, 0 for still alive today
breed$censorship[which(breed$LastObsDate=="2015-10-14")]<-0

str(breed)

##Plot age at first breeding vs. years experience 
plot(breed$CurrentAge, breed$YrsExp)

#Want to get an idea of the different ages at which birds first start breeding
#GF suggests bubble plot - crime in the states example 
#I think I need to get the sum of the number of birds at each age first 
#breeding occurred

sum(breed)

##Find a different way to do this, there has got to be a better way
b1 <- subset(breed, breed$YrsExp ==0 & breed$AgeFirstBreed == 1)
b2 <- subset(breed, breed$YrsExp ==0 & breed$AgeFirstBreed == 2)
b3 <- subset(breed, breed$YrsExp ==0 & breed$AgeFirstBreed ==3)
b4 <- subset(breed, breed$YrsExp ==0 & breed$AgeFirstBreed == 4)
b5 <- subset(breed, breed$YrsExp ==0 & breed$AgeFirstBreed == 5)
b6 <- subset(breed, breed$YrsExp ==0 & breed$AgeFirstBreed == 6)
b7 <- subset(breed, breed$YrsExp ==0 & breed$AgeFirstBreed == 7)
b8 <- subset(breed, breed$YrsExp ==0 & breed$AgeFirstBreed == 8)
b9 <- subset(breed, breed$YrsExp ==0 & breed$AgeFirstBreed == 9)
b10 <- subset(breed, breed$YrsExp ==0 & breed$AgeFirstBreed == 10)

number <- c(52,494,298,75,22,6,8,0,2,1)
agename <- c("Age1","Age2","Age3","Age4","Age5","Age6","Age7","Age8","Age9",
             "Age10")
a <- data.frame(agename, number)

ggplot(a, aes(x=agename, y=number, size=number)) +
geom_point(colour="white", fill ="blue", shape=21)


#change back to numeric for survival object 
breed$MinDate <- as.numeric(breed$FirstDateBred)
breed$LastObsDate <- as.numeric(breed$LastObsDate)

