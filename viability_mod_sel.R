# Model selection for seed VIABILITY
# Data: tetrazolium scoring and germination data 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(lme4) ; library(bbmle) ; library(boot) ; library(testthat)
library(dplyr) ; library(glmmADMB)
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/ac79/Documents/CODE/LLELA/model_sel_results.R")

# Read in data and format -----------------------------------------------------------
d           <- read.csv("Data/vr.csv")
viabVr      <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                        stringsAsFactors = F)

# Format data -----------------------------------------------------------------------
viabVr      <- viabVr %>%
                mutate(totN = F + M) %>%
                mutate(sr = F / totN) %>%
                mutate(plot = as.factor(plot) )
viabVr      <- subset(viabVr, totFlow < 60) #exclude extreme values


# Viability model selection ---------------------------------------------------------

# Omit NAs for glmmadmb
tetr_dat  <- na.omit(dplyr::select(viabVr,Yes,fail,sr_f,totFlow,sr,totN,plot,log_l_t0))
germ_dat  <- na.omit(dplyr::select(viabVr,germTot,germFail,sr_f,totFlow,sr,totN,plot,log_l_t0))


# should I use individual size or not? ----------------------------------------------
# germination
gSize <- list()
gSize[[1]]=glmmadmb(cbind(germTot,germFail) ~ (1 | plot), family = "binomial", data=germ_dat)
gSize[[2]]=glmmadmb(cbind(germTot,germFail) ~ log_l_t0 + (1 | plot), family = "binomial", data=germ_dat)
# answer: no. (slightly, but no!)
AICtab(gSize, weights=T)

# tetrazolium
tSize <- list()
tSize[[1]]=glmmadmb(cbind(Yes,fail) ~ (1 | plot), family = "binomial", data=tetr_dat)
tSize[[2]]=glmmadmb(cbind(Yes,fail) ~ log_l_t0 + (1 | plot), family = "binomial", data=tetr_dat)
# answer: no. (slightly, but no!)
AICtab(tSize, weights=T)


# 1. Viable Seed Number: germination assays (germTot / germFail)---------------------
candidate_mods <- list(
  cbind(germTot,germFail) ~ (1 | plot),
  cbind(germTot,germFail) ~ sr_f + (1 | plot),
  cbind(germTot,germFail) ~ totFlow + (1 | plot),
  cbind(germTot,germFail) ~ sr_f + totFlow + (1 | plot),
  cbind(germTot,germFail) ~ sr_f * totFlow + (1 | plot)
)

# fit models
germ_flowN  <- lapply(candidate_mods, function(x) glmmadmb(x, family="binomial", data=germ_dat))
germ_flowN_sel  <- AICtab(germ_flowN,weights=T) 
germ_flowN_avg  <- model_avg(germ_flowN_sel, germ_flowN)
write.csv(germ_flowN_avg, "Results/VitalRates_3/germination_best.csv", row.names = F)

# Model selection table
x <- germ_flowN
write.csv(mod_sel_res(germ_flowN),
          "Results/VitalRates_3/mod_sel_germination.csv",row.names=F)
rm(x)


# 2. viable Seed Number: tetrazolium assays  (Yes / fail) ---------------------------
candidate_mods <- list(
  cbind(Yes,fail) ~ (1 | plot),
  cbind(Yes,fail) ~ sr_f + (1 | plot),
  cbind(Yes,fail) ~ totFlow + (1 | plot),
  cbind(Yes,fail) ~ sr_f + totFlow + (1 | plot),
  cbind(Yes,fail) ~ sr_f * totFlow + (1 | plot)
)

# fit models
tetr_flowN        <- lapply(candidate_mods, function(x) glmmadmb(x,family="binomial", data=tetr_dat))
tetr_flowN_select <- AICtab(tetr_flowN,weights=T)
tetr_flowN_avg    <- model_avg(tetr_flowN_select, tetr_flowN)
write.csv(tetr_flowN_avg, "Results/VitalRates_3/tetrazolium_best.csv", row.names = F)

# Model selection table
x <- tetr_flowN
write.csv(mod_sel_res(tetr_flowN),
          "mod_sel_tetrazolium.csv",row.names=F)
rm(x)
