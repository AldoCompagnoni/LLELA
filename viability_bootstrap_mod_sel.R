# Model selection for seed VIABILITY
# Data: tetrazolium scoring and germination data 
# setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
setwd("C:/cloud/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(lme4) ; library(bbmle) ; library(boot) ; library(testthat)
library(dplyr) ; library(tidyr)
library(glmmADMB)
# source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")
# source("C:/Users/ac79/Documents/CODE/LLELA/model_sel_results.R")
source("C:/CODE/LLELA/model_avg.R")
source("C:/CODE/LLELA/model_sel_results.R")


# Read in data and format -----------------------------------------------------------
d           <- read.csv("Data/vr.csv")
viabVr      <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                        stringsAsFactors = F)


# Format data -----------------------------------------------------------------------
viabVr      <- viabVr %>%
                mutate(totN = F + M) %>%
                mutate(sr = F / totN) %>%
                mutate(perm_id = paste(plot, focalI, P_ID, sep = "_") ) %>%
                mutate(plot = as.factor(plot) )
viabVr      <- subset(viabVr, totFlow < 60) #exclude extreme values 
              

# Viability data --------------------------------------------------------------------


# Omit NAs for glmmadmb
tetr_dat  <- na.omit(dplyr::select(viabVr,Yes,fail,sr_f,totFlow,sr,totN,plot,focalI,P_ID,log_l_t0, F, M, perm_id))
germ_dat  <- na.omit(dplyr::select(viabVr,germTot,germFail,sr_f,totFlow,sr,totN,plot,focalI,P_ID,log_l_t0, F, M, perm_id))
              
# randomized one-female plots
one_f_2   <- combinations(8, 2, germ_dat$perm_id)
one_f_3   <- combinations(8, 3, germ_dat$perm_id)
one_f_4   <- combinations(8, 4, germ_dat$perm_id)
   

# # Viability model selection ---------------------------------------------------------
# 
# # Omit NAs for glmmadmb
# tetr_dat  <- na.omit(dplyr::select(viabVr,Yes,fail,sr_f,totFlow,sr,totN,plot,log_l_t0, F, M))
# germ_dat  <- na.omit(dplyr::select(viabVr,germTot,germFail,sr_f,totFlow,sr,totN,plot,log_l_t0, F, M))
# 
# 
# # should I use individual size or not? ----------------------------------------------
# # germination
# gSize <- list()
# gSize[[1]]=glmmadmb(cbind(germTot,germFail) ~ (1 | plot), family = "binomial", data=germ_dat)
# gSize[[2]]=glmmadmb(cbind(germTot,germFail) ~ log_l_t0 + (1 | plot), family = "binomial", data=germ_dat)
# # answer: no. (slightly, but no!)
# AICtab(gSize, weights=T)
# 
# # tetrazolium
# tSize <- list()
# tSize[[1]]=glmmadmb(cbind(Yes,fail) ~ (1 | plot), family = "binomial", data=tetr_dat)
# tSize[[2]]=glmmadmb(cbind(Yes,fail) ~ log_l_t0 + (1 | plot), family = "binomial", data=tetr_dat)
# # answer: no. (slightly, but no!)
# AICtab(tSize, weights=T)


# 1. Viable Seed Number: germination assays (germTot / germFail)---------------------

# model structures
candidate_mods <- list(
  cbind(germTot,germFail) ~ (1 | plot),
  cbind(germTot,germFail) ~ sr_f + (1 | plot),
  cbind(germTot,germFail) ~ totFlow + (1 | plot),
  cbind(germTot,germFail) ~ sr_f + totFlow + (1 | plot),
  cbind(germTot,germFail) ~ sr_f * totFlow + (1 | plot)
)

# # fit models
# germ_flowN      <- lapply(candidate_mods, function(x) glmmadmb(x, family="binomial", data=germ_dat))
# germ_flowN_sel  <- AICtab(germ_flowN,weights=T) 
# germ_flowN_avg  <- model_avg(germ_flowN_sel, germ_flowN)


# function to perform bootstrapped 
bootstrap_mod <- function(i, one_f_2){
  
  # get one-female plots and test that we have 2 of them 
  one_f_data      <- subset(germ_dat, perm_id %in% one_f_2[i,] )
  expect_equal(nrow(one_f_data), 2)
  
  mod_data        <- subset(germ_dat, !(F == 1 & M == 0)) %>%
                        bind_rows( one_f_data )
  germ_flowN      <- lapply(candidate_mods, function(x) glmmadmb(x, family="binomial", data=mod_data))
  germ_flowN_sel  <- AICtab(germ_flowN,weights=T) 
  germ_flowN_avg  <- model_avg(germ_flowN_sel, germ_flowN)
  
  return(germ_flowN_avg)
  
}
# bootstrapped results
bootstr_res_l <- lapply(1:nrow(one_f_2), bootstrap_mod, one_f_2) 
bootstr_res   <- Reduce(function(...) rbind(...), bootstr_res_l)

boxplot(avg ~ predictor, data = bootstr_res)
abline(h=0, lty = 2)
 




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
