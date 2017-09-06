##Sex predictor of flowering probability: SEX RATIO (female individuals/total individuals)
setwd("C:/cloud/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle)
library(glmmADMB)
library(dplyr)
source("C:/CODE/LLELA/model_avg.R")
source("C:/CODE/LLELA/model_sel_results.R")


# read in data -------------------------------------------------------------------------------
d       <- read.csv("Data/vr.csv")
f14     <- format_flower(d)


# Model selection ----------------------------------------------------------------------------------
candidate_mods <- list(
  flowN_t1 ~ ( 1 | plot),
  flowN_t1 ~ sex + ( 1 | plot),
  flowN_t1 ~ TotDensity + ( 1 | plot),
  flowN_t1 ~ TotDensity + sex:TotDensity + ( 1 | plot),
  flowN_t1 ~ sex + TotDensity + sex:TotDensity + ( 1 | plot),
  flowN_t1 ~ TotDensity + sr:TotDensity + ( 1 | plot),
  flowN_t1 ~ sex + TotDensity + sr:TotDensity + ( 1 | plot),
  flowN_t1 ~ TotDensity + sr:TotDensity + sex:TotDensity + ( 1 | plot),
  flowN_t1 ~ sex + TotDensity + sr:TotDensity + sex:TotDensity + ( 1 | plot)
)

# fit models and model averages
nfMod         <- lapply(candidate_mods, function(x) glmmadmb(x,data=f14,family="nbinom2") )
n_flow_select <- AICtab(nfMod, weights=T)
n_flow_avg    <- model_avg(n_flow_select, nfMod)
write.csv(n_flow_avg, "Results/VitalRates_4_Cade2015/n_flowers_best.csv", row.names = F)


# Model selection table
x <- nfMod
write.csv(mod_sel_res(nfMod), 
          "Results/VitalRates_4_Cade2015/mod_sel_n_flowers.csv",row.names=F)
rm(x)
