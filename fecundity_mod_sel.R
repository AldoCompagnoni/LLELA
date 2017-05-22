# model selection for fecundity (seeds per flower) 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(MASS)
library(dplyr)
library(glmmADMB)
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/ac79/Documents/CODE/LLELA/model_sel_results.R")

#read in data--------------------------------------------------------------
d         <- read.csv("Data/vr.csv")
fem_seeds <- read.csv("Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")

# plot level data ---------------------------------------------------------

# Only data from 2014
d14 <- subset(d, year == 2014)
d14 <- subset(d14, surv_t1 != 0)

# fecundity data ---------------------------------------------------------------------- 
fem_seeds$focalI    <- paste("f",fem_seeds$IndividualN,sep="")
fem_seeds           <- fem_seeds[,c("Plot","focalI","SeedN")] 
fem_seeds           <- fem_seeds[!is.na(fem_seeds$SeedN),]
names(fem_seeds)[1] <- "plot"

# merge data sets 
fecund_data      <- merge(d14,fem_seeds) 
fecund_data$plot <- as.factor(fecund_data$plot)

# Preliminary analyses ---------------------------------------------------------------------

# should I use individual size or not? 
isMod = list()
isMod[[1]]=glmmadmb(SeedN ~ (1 | plot), family = "nbinom2", data=fecund_data)
isMod[[2]]=glmmadmb(SeedN ~ log_l_t0 + (1 | plot), family = "nbinom2", data=fecund_data)

# answer: no. (slightly, but no!)
AICtab(isMod, weights=T)


# MODEL SELECTION --------------------------------------------------------

# formulas for candidate models 
candidate_mods <- list(
  SeedN ~ ( 1 | plot),
  SeedN ~ TotDensity + ( 1 | plot),
  SeedN ~ TotDensity + TotDensity:sr + ( 1 | plot)
)

# model specification
fit_func <- function(x){
  return( glmmadmb(x, family="nbinom2", data=fecund_data) )
}

# fit models
nsMod       <- lapply(candidate_mods, fit_func)
fec_select  <- AICtab(nsMod, weights=T)
fec_avg     <- model_avg(fec_select, nsMod)
write.csv(fec_avg, "Results/VitalRates_3/fecuntity_best.csv", row.names = F)

# Model selection table
x <- nsMod # trick to feed AICtab the models within a function
write.csv(mod_sel_res(nsMod),
          "Results/VitalRates_3/mod_sel_fecuntity.csv",row.names=F)
rm(x)