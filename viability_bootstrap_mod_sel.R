# Model selection for seed VIABILITY
# Data: tetrazolium scoring and germination data 
# setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
setwd("C:/cloud/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(lme4) ; library(bbmle) ; library(boot) ; library(testthat)
library(dplyr) ; library(tidyr) ; library(gtools)
library(glmmADMB) ; library(parallel)
# source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")
# source("C:/Users/ac79/Documents/CODE/LLELA/model_sel_results.R")
source("C:/CODE/LLELA/model_avg.R")
source("C:/CODE/LLELA/model_sel_results.R")
source("C:/CODE/LLELA/model_avg.R")
source("C:/CODE/LLELA/prediction.R")


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
  
detach("package:gtools")

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


# parallel processing set up ---------------------------------------------------------------

# detect cores
cores     <- detectCores()
# set up as many clusters as detected by 'detectCores()'
cluster   <- parallel::makePSOCKcluster(cores)

# attach packages that will be needed on each cluster
clusterEvalQ(cluster, list(library(lme4), library(glmmADMB), library(bbmle) , library(boot) , 
                           library(testthat), library(dplyr) , library(tidyr)) )
# attach objects that will be needed on each cluster
clusterExport(cluster, "model_avg")


# 1. Viable Seed Number: germination assays (germTot / germFail)---------------------

# model structures
candidate_mods <- list(
  cbind(germTot,germFail) ~ (1 | plot),
  cbind(germTot,germFail) ~ sr_f + (1 | plot),
  cbind(germTot,germFail) ~ totFlow + (1 | plot),
  cbind(germTot,germFail) ~ sr_f + totFlow + (1 | plot),
  cbind(germTot,germFail) ~ sr_f * totFlow + (1 | plot)
)
candidate_mods <- setNames(candidate_mods, paste0("model",1:5))

# fit BEST model
germ_flowN      <- lapply(candidate_mods, function(x) glmmadmb(x, family="binomial", data=germ_dat))
germ_flowN_sel  <- AICtab(germ_flowN,weights=T)
germ_flowN_avg  <- model_avg(germ_flowN_sel, germ_flowN)


# function to perform bootstrapped estimates
bootstrap_mod <- function(i, one_f_2, germ_dat, candidate_mods){
  
  # get one-female plots and test that we have 2 of them 
  one_f_data      <- subset(germ_dat, perm_id %in% one_f_2[i,] )
  expect_equal(nrow(one_f_data), 2)
  
  mod_data        <- subset(germ_dat, !(F == 1 & M == 0)) %>%
                        bind_rows( one_f_data )
  germ_flowN      <- lapply(candidate_mods, function(x) glmmadmb(x, family="binomial", data=mod_data) )
  germ_flowN_sel  <- AICtab(germ_flowN, mnames = names(germ_flowN), weights=T)  
  germ_flowN_avg  <- model_avg(germ_flowN_sel, germ_flowN)
  
  return(germ_flowN_avg)
  
}

boostr_par_l <- list()
for(i in 1:nrow(one_f_2)){ 
  
  # get one-female plots and test that we have 2 of them 
  one_f_data      <- subset(germ_dat, perm_id %in% one_f_2[i,] )
  expect_equal(nrow(one_f_data), 2)
  
  mod_data        <- subset(germ_dat, !(F == 1 & M == 0)) %>%
                        bind_rows( one_f_data )
  germ_flowN      <- lapply(candidate_mods, function(x) glmmadmb(x, family="binomial", data=mod_data) )
  germ_flowN_sel  <- AICtab(germ_flowN, mnames = names(germ_flowN), weights=T)  
  germ_flowN_avg  <- model_avg(germ_flowN_sel, germ_flowN)
  
  boostr_par_l[[i]] <- germ_flowN_avg

}


# run 10 processes on the # clusters defined in "cluster"
boostr_par_l  <- parLapply(cluster, 1:nrow(one_f_2), bootstrap_mod, one_f_2, germ_dat, candidate_mods)
bootstr_res_l <- Map(function(x,y) tibble::add_column(x, sample = y, .before = 1),
                     boostr_par_l, 1:length(one_f_2) )
bootstr_res   <- Reduce(function(...) rbind(...), bootstr_res_l)

boxplot(avg ~ predictor, data = bootstr_res)
abline(h=0, lty = 2)
 


# Graph -----------------------------------------------------------------------------

# model matrix for prediction
predic_seq<- bootstr_res %>%
                subset( sample==1 ) %>%
                .$predictor %>%
                as.character
mod_des   <- expand.grid("(Intercept)" = 1, 
                         totFlow = seq(1,50,1),
                         sr_f = round(seq(0,1,0.05),2)) %>%
                mutate( 'sr_f:totFlow' = totFlow * sr_f) %>%
                .[,predic_seq]

# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(gray.colors(7,start = 0, end = 1))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

# Sex ratio as dot color 
par(mfrow=c(1,1),mar=c(2.6,2.5,0.2,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,10))

# viability (germination) data ---------------------------------------------------------
plot(jitter(viabVr$totFlow,factor = 2),jitter(viabVr$germ_ratio,factor = 2),
     pch=21,ylim=c(0,1.01),xlim=c(0,48), 
     bg = cRamp(viabVr$sr_f), cex = 1.5,
     xlab="Panicle density",ylab="Seed viability rate")

# predicted values using 'pred' function
plot_bootstr <- function(i){ 
  tmp       <- subset(bootstr_res, sample == i)[,-1]
  mod_pred  <- pred(h_des, tmp, viab, inv.logit)
  lines(h_des$totFlow,mod_pred$viab,lwd=1,lty=1, col = "grey")
}
lapply(1:nrow(one_f_2), plot_bootstr)

# MEAN prediction
h_des     <- subset(mod_des,sr_f == 0.95)
mod_pred  <- pred(h_des, germ_flowN_avg, viab, inv.logit)
lines(h_des$totFlow,mod_pred$viab,lwd=1,lty=1)




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
