# July 5 2016
# Pedigree simulation so I can learn the ins and outs
# Mostly to correctly code the founders
# The founders fields for dadid and momid have NA in these fields 

library(kinship2)
pdsim <- read.csv("Ped_sim.csv")
ped.df <- pdsim[-(12:185), ] 

ped.df["sexcode"] <- 0
ped.df["status"] <- 1
ped.df["affected"] <- 0

sim.ped <- pedigree(ped.df$JayID, ped.df$DadID, ped.df$MomID, ped.df$sexcode,
              ped.df$affected, ped.df$status)