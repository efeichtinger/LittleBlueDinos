#TSF EDA - exploring the time since fire data 
#September 2015

#use this on laptop
setwd("C:/Users/Erin/Dropbox/Jay survival analyses")

#use this on school desktop
setwd("C:/Users/efeichtinger/Dropbox/Jay survival analyses/Fire and veg data Fall 2015")

tsf.df <- read.csv(file="TSF_TerrYr.csv")
head(tsf.df)
str(tsf.df)
dim(tsf.df)

#Have to filter out rows where InABS =FALSE
tsf <- subset(tsf.df, InABS==TRUE, select=TERRYR:DateOfLastFire)

#find mean time since fire in years by terryr
agg1 <- aggregate(tsf.df$TSF_years, list(terryr = tsf.df$TERRYR), mean)
agg2 <- aggregate(tsf.df$TSF_years, list(terryr= tsf.df$TERRYR), sd)

agg3 <- aggregate(tsf.df$CellCount, list(terryr = tsf.df$TERRYR), mean)
agg4 <- aggregate(tsf.df$CellCount, list(terryr = tsf.df$TERRYR), sd)




