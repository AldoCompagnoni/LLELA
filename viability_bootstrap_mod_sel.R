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
germ_dat  <- na.omit(dplyr::select(viabVr,germTot,germFail,sr_f,totFlow,sr,totN,plot,focalI,P_ID,log_l_t0, F, M, perm_id))
one_f_df  <- subset(germ_dat, F == 1 & M == 0)

# randomized one-female plots
one_f_2   <- combinations(8, 2, one_f_df$perm_id)
one_f_3   <- combinations(8, 3, one_f_df$perm_id)
one_f_4   <- combinations(8, 4, one_f_df$perm_id)
detach("package:gtools")


# Viable Seed Number: germination assays (germTot / germFail)---------------------

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



# parallel processing model fit -------------------------------------------------------

# detect cores
cores     <- detectCores()
# set up as many clusters as detected by 'detectCores()'
cluster   <- parallel::makePSOCKcluster(cores)

# attach packages that will be needed on each cluster
clusterEvalQ(cluster, list(library(lme4), library(glmmADMB), library(bbmle) , library(boot) , 
                           library(testthat), library(dplyr) , library(tidyr)) )
# attach objects that will be needed on each cluster
clusterExport(cluster, "model_avg")


# function to perform bootstrapped estimates
bootstrap_mod <- function(i, one_f_2, germ_dat, candidate_mods){
  
  # get one-female plots and test that we have 2 of them 
  one_f_data      <- subset(germ_dat, perm_id %in% one_f_2[i,] )
  expect_equal(nrow(one_f_data), 2)
  expect_equal(unique(one_f_data$F), 1)
  
  no_one_f        <- subset(germ_dat, !(F == 1 & M == 0) )
  expect_equal( nrow(subset(no_one_f, F == 1 & M == 0)), 0 )
    
  mod_data        <- bind_rows( no_one_f, one_f_data )
  germ_flowN      <- lapply(candidate_mods, function(x) glmmadmb(x, family="binomial", data=mod_data) )
  germ_flowN_sel  <- AICtab(germ_flowN, mnames = names(germ_flowN), weights=T)  
  germ_flowN_avg  <- model_avg(germ_flowN_sel, germ_flowN)
  
  return(germ_flowN_avg)
  
}

# run processes on the # clusters defined in "cluster"
boostr_par_1_l    <- parLapply(cluster, 1:nrow(one_f_2), bootstrap_mod, one_f_2, germ_dat, candidate_mods)
bootstr_res_1_l   <- Map(function(x,y) tibble::add_column(x, sample = y, .before = 1),
                          boostr_par_1_l, 1:nrow(one_f_2) )
bootstr_res_1     <- Reduce(function(...) rbind(...), bootstr_res_1_l) %>%
                          arrange(sample, predictor)


# Graph -----------------------------------------------------------------------------

# model matrix for prediction
mod_des   <- expand.grid("(Intercept)" = 1, 
                         totFlow = seq(1,48,1),
                         sr_f = round(seq(0,1,0.05),2)) %>%
                mutate( 'sr_f:totFlow' = totFlow * sr_f) %>%
                .[,order(names(.))]

# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(gray.colors(7,start = 0, end = 1))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

# design matrices
l_des     <- subset(mod_des,sr_f == 0.05)
m_des     <- subset(mod_des,sr_f == 0.5)
h_des     <- subset(mod_des,sr_f == 0.95)

# return N. viability predictions 
pred_bootstr <- function(i, des){ 
  tmp       <- subset(bootstr_res_1, sample == i)[,-1]
  mod_pred  <- pred(des, tmp, viab, inv.logit)
  return(mod_pred$viab)
}
preds_h <- lapply(1:nrow(one_f_2), pred_bootstr, h_des)
preds_h <- Reduce(function(...) rbind(...), preds_h)
preds_m <- lapply(1:nrow(one_f_2), pred_bootstr, m_des)
preds_m <- Reduce(function(...) rbind(...), preds_m)
preds_l <- lapply(1:nrow(one_f_2), pred_bootstr, l_des)
preds_l <- Reduce(function(...) rbind(...), preds_l)

# maximum and minimum prediction 
pred_min_max <- function(x){
  
  data.frame( upr = apply(x,2,min), 
              lwr = apply(x,2,max) )
  
}
preds_h_lim <- pred_min_max(preds_h)
preds_m_lim <- pred_min_max(preds_m)
preds_l_lim <- pred_min_max(preds_l)


# the plot 
tiff("Results/VitalRates_3/Figure5_bootstrap_one_fem.tiff",
     unit="in",width=6.3,height=4,res=600,compression="lzw")

# Sex ratio as dot color 
par(mfrow=c(1,1),mar=c(2.6,2.5,0.2,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,10))

# viability (germination) data ---------------------------------------------------------
plot(jitter(viabVr$totFlow,factor = 2),jitter(viabVr$germ_ratio,factor = 2),
     pch=21,ylim=c(0,1.01),xlim=c(0,48), 
     bg = cRamp(viabVr$sr_f), cex = 1.5,
     xlab="Panicle density",ylab="Seed viability rate")

newx        <- arrange(h_des, totFlow) %>% .$totFlow
polygon(c(rev(newx), newx), c(rev(preds_h_lim$upr), preds_h_lim$lwr), 
        col = adjustcolor( "grey80", alpha.f = 0.5), border = NA)
polygon(c(rev(newx), newx), c(rev(preds_m_lim$upr), preds_m_lim$lwr), 
        col = adjustcolor( "grey80", alpha.f = 0.5), border = NA)
polygon(c(rev(newx), newx), c(rev(preds_l_lim$upr), preds_l_lim$lwr), 
        col = adjustcolor( "grey80", alpha.f = 0.5), border = NA)

polygon(c(rev(newx), newx), c(rev(preds_h_lim$upr), preds_h_lim$lwr), density = 10, angle = 45)


# MEAN predictions
mod_pred  <- pred(l_des, germ_flowN_avg, viab, inv.logit)
lines(l_des$totFlow,mod_pred$viab,lwd=2,lty=3)
mod_pred  <- pred(m_des, germ_flowN_avg, viab, inv.logit)
lines(m_des$totFlow,mod_pred$viab,lwd=2,lty=2)
mod_pred  <- pred(h_des, germ_flowN_avg, viab, inv.logit)
lines(h_des$totFlow,mod_pred$viab,lwd=2,lty=1)

# sex ratio legend
colfunc = colorRampPalette(cRamp(unique(arrange(viabVr,sr_f)$sr_f)))
legend_image <- as.raster(matrix(colfunc(21), ncol=1))
text(x=60, y = seq(0.6,1,l=3), labels = seq(1,0,l=3), xpd=NA)
rasterImage(legend_image, 51, 0.6, 58, 1, xpd=NA)
text(61, 1, "Sex ratio", pos = 4,xpd=NA)
text(61, 0.8, "(proportion", pos = 4,xpd=NA)
text(61, 0.6, "female)", pos = 4,xpd=NA)

# prediction legend
legend(50,0.55, c("5% female plot",
                  "50% female plot",
                  "95% female plot"), cex = 1, xpd = NA,
       lty = c(3,2,1), lwd=2, bty = "n")

dev.off()



# bootstrap by treatment -----------------------------------------------------------
design  <- germ_dat %>% 
              dplyr::select(F, M) %>%
              unique

# extract at least 1 sample, at most 2 samples from  
extract_sample <- function(i, germ_dat, design){
  
  treat   <- design[i,]
  tmp     <- subset(germ_dat, F == treat$F & M == treat$M) 
  
  if( nrow(tmp) == 1 ){
    out     <- tmp  
  } else{
    keep_i  <- sample(1:nrow(tmp), 2)
    out     <- tmp[keep_i,]
  }
  
  out
  
}

# do analyses
bootstrap_mod <- function(i, design, germ_dat, candidate_mods){
  
  # get at most 2 samples per treatment 
  boot_l    <- lapply(1:nrow(design), extract_sample, germ_dat, design)
  booted_df <- Reduce(function(...) rbind(...), boot_l)
  
  # fit, select, and perform model average
  germ_flowN      <- lapply(candidate_mods, function(x) glmmadmb(x, family="binomial", data=booted_df) )
  germ_flowN_sel  <- AICtab(germ_flowN, mnames = names(germ_flowN), weights=T)  
  germ_flowN_avg  <- model_avg(germ_flowN_sel, germ_flowN)
  
  return(germ_flowN_avg)
  
}

# attach objects that will be needed on each cluster
clusterExport(cluster, "extract_sample")

boostr_par_l  <- parLapply(cluster, 1:200,  bootstrap_mod, design, germ_dat, candidate_mods)
bootstr_res_l <- Map(function(x,y) tibble::add_column(x, sample = y, .before = 1),
                     boostr_par_l, 1:200 )
bootstr_res   <- Reduce(function(...) bind_rows(...), bootstr_res_l)


# Graph -----------------------------------------------------------------------------

# model matrix for prediction
predic_seq<- bootstr_res %>%
                subset( sample==1 ) %>%
                .$predictor %>%
                as.character
mod_des   <- expand.grid("(Intercept)" = 1, 
                         totFlow = seq(1,48,1),
                         sr_f = round(seq(0,1,0.05),2)) %>%
                  mutate( 'sr_f:totFlow' = totFlow * sr_f) %>%
                  .[,predic_seq]

# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(gray.colors(7,start = 0, end = 1))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

# design matrices
l_des     <- subset(mod_des,sr_f == 0.05)
m_des     <- subset(mod_des,sr_f == 0.5)
h_des     <- subset(mod_des,sr_f == 0.95)

# return N. viability predictions 
pred_bootstr <- function(i, des){ 
  tmp       <- subset(bootstr_res, sample == i)[,-1]
  mod_pred  <- pred(des, tmp, viab, inv.logit)
  return(mod_pred$viab)
}
preds_h <- lapply(1:200, pred_bootstr, h_des)
preds_h <- Reduce(function(...) rbind(...), preds_h)
preds_l <- lapply(1:200, pred_bootstr, l_des)
preds_l <- Reduce(function(...) rbind(...), preds_l)
preds_m <- lapply(1:200, pred_bootstr, m_des)
preds_m <- Reduce(function(...) rbind(...), preds_m)

# maximum and minimum prediction 
pred_min_max <- function(x){
  
  data.frame( upr = apply(x,2,min), 
              lwr = apply(x,2,max) )
  
}
preds_h_lim <- pred_min_max(preds_h)
preds_m_lim <- pred_min_max(preds_m)
preds_l_lim <- pred_min_max(preds_l)


#Graph: total number of tillers -----------------------------------------------------
tiff("Results/VitalRates_3/Figure5_bootstrapped.tiff",
     unit="in",width=6.3,height=4,res=600,compression="lzw")

# Sex ratio as dot color 
par(mfrow=c(1,1),mar=c(2.6,2.5,0.2,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,10))

# plot in itself 
plot(jitter(viabVr$totFlow,factor = 2),jitter(viabVr$germ_ratio,factor = 2),
     pch=21,ylim=c(0,1.01),xlim=c(0,48), 
     bg = cRamp(viabVr$sr_f), cex = 1.5,
     xlab="Panicle density",ylab="Seed viability rate")

newx        <- arrange(h_des, totFlow) %>% .$totFlow
polygon(c(rev(newx), newx), c(rev(preds_h_lim$upr), preds_h_lim$lwr), 
        col = adjustcolor( "grey80", alpha.f = 0.5), border = NA)
polygon(c(rev(newx), newx), c(rev(preds_m_lim$upr), preds_m_lim$lwr), 
        col = adjustcolor( "grey80", alpha.f = 0.5), border = NA)
polygon(c(rev(newx), newx), c(rev(preds_l_lim$upr), preds_l_lim$lwr), 
        col = adjustcolor( "grey80", alpha.f = 0.5), border = NA)

# MEAN predictions
mod_pred  <- pred(l_des, germ_flowN_avg, viab, inv.logit)
lines(l_des$totFlow,mod_pred$viab,lwd=2,lty=3)
mod_pred  <- pred(m_des, germ_flowN_avg, viab, inv.logit)
lines(m_des$totFlow,mod_pred$viab,lwd=2,lty=2)
mod_pred  <- pred(h_des, germ_flowN_avg, viab, inv.logit)
lines(h_des$totFlow,mod_pred$viab,lwd=2,lty=1)

# sex ratio legend
colfunc = colorRampPalette(cRamp(unique(arrange(viabVr,sr_f)$sr_f)))
legend_image <- as.raster(matrix(colfunc(21), ncol=1))
text(x=60, y = seq(0.6,1,l=3), labels = seq(1,0,l=3), xpd=NA)
rasterImage(legend_image, 51, 0.6, 58, 1, xpd=NA)
text(61, 1, "Sex ratio", pos = 4,xpd=NA)
text(61, 0.8, "(proportion", pos = 4,xpd=NA)
text(61, 0.6, "female)", pos = 4,xpd=NA)

# prediction legend
legend(50,0.55, c("5% female plot",
                  "50% female plot",
                  "95% female plot"), cex = 1, xpd = NA,
       lty = c(3,2,1), lwd=2, bty = "n")

dev.off()
