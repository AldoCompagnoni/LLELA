# Sex competition predictor: SEX RATIO (female individuals/total individuals)
setwd("C:/cloud/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)
library(lme4)
library(dplyr)
source("C:/CODE/LLELA/model_avg.R")
source("C:/CODE/LLELA/model_sel_results.R")

# load and format data ------------------------------------------------------------------
d     <- format_growth( read.csv("Data/vr.csv") )
d14   <- subset(d, year == 2014)
d15   <- subset(d, year == 2015)

# Analysis ----------------------------------------------------------------------
candidate_mods <- list(
  log_ratio ~ ( 1 | plot),
  log_ratio ~ sex + ( 1 | plot),
  log_ratio ~ TotDensity + ( 1 | plot),
  log_ratio ~ TotDensity + sex:TotDensity + ( 1 | plot),
  log_ratio ~ sex + TotDensity + sex:TotDensity + ( 1 | plot),
  log_ratio ~ TotDensity + sr:TotDensity + ( 1 | plot),
  log_ratio ~ sex + TotDensity + sr:TotDensity + ( 1 | plot),
  log_ratio ~ TotDensity + sr:TotDensity + sex:TotDensity + ( 1 | plot),
  log_ratio ~ sex + TotDensity + sr:TotDensity + sex:TotDensity + ( 1 | plot)
)

# 2014 fit models and cal. mod weights
m14             <- lapply(candidate_mods, function(x) lmer(x, data=d14))
lr_14_mod_sel   <- AICtab(m14, weights = T)
lr_14_avg       <- model_avg(lr_14_mod_sel, m14, d14)
write.csv(lr_14_avg, "Results/VitalRates_4_Cade2015/log_ratio_14_best.csv", row.names = F)

# model selection results
x <- m14
mod_14_w        <- rgr_mod_sel(m14)
write.csv(mod_14_w,"Results/VitalRates_4_Cade2015/mod_sel_log_ratio_14.csv",row.names=F)
rm(x)

# 2015 fit models and model weights
m15            <- lapply(candidate_mods, function(x) lmer(x, data=d15))
lr_15_mod_sel  <- AICtab(m15, weights = T)
lr_15_avg      <- model_avg(lr_15_mod_sel, m15)
write.csv(lr_15_avg, "Results/VitalRates_4_Cade2015/log_ratio_15_best.csv", row.names = F)

# model selection results
x              <- m15 
mod_15_w       <- rgr_mod_sel(m15)
write.csv(mod_15_w,"Results/VitalRates_4_Cade2015/mod_sel_log_ratio_15.csv",row.names=F)
rm(x)