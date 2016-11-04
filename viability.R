# Model selection for seed VIABILITY
# Data: tetrazolium scoring and germination data 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(lme4) ; library(bbmle) ; library(boot) ; library(testthat)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")

# Read in data and format -----------------------------------------------------------
viabVr      <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                        stringsAsFactors = F)
viabVr$totN <- viabVr$F + viabVr$M
viabVr$sr   <- viabVr$F / viabVr$totN
viabVr$plot <- as.factor(viabVr$plot)

# Viability model selection ---------------------------------------------------------

# 1. viable Seed Number: tetrazolium assays  (Yes / fail) ---------------------------
# This is the "standard" scoring: consider as viable any seed stained by tetrazolium

# Omit NAs for glmmadmb
tetr_dat  <- na.omit(select(viabVr,yesMaybe,failMaybe,sr_f,totFlow,sr,totN,plot))

# Number of flowers and their sex ratio as predictors
tetr_flowN=list()
tetr_flowN[[1]]  <- glmmadmb(cbind(yesMaybe,failMaybe) ~ sr_f + (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN[[2]]  <- glmmadmb(cbind(yesMaybe,failMaybe) ~ totFlow + (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN[[3]]  <- glmmadmb(cbind(yesMaybe,failMaybe) ~ sr_f + totFlow + (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN[[4]]  <- glmmadmb(cbind(yesMaybe,failMaybe) ~ sr_f * totFlow + (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN_select<- AICtab(tetr_flowN,weights=T)

# Planting density and sex ratio as predictors
tetr_dens=list()
tetr_dens[[1]]    <- glmmadmb(cbind(yesMaybe,failMaybe) ~ sr + (1 | plot),family="binomial", data=tetr_dat)
tetr_dens[[2]]    <- glmmadmb(cbind(yesMaybe,failMaybe) ~ totN + (1 | plot),family="binomial", data=tetr_dat)
tetr_dens[[3]]    <- glmmadmb(cbind(yesMaybe,failMaybe) ~ sr + totN + (1 | plot),family="binomial", data=tetr_dat)
tetr_dens[[4]]    <- glmmadmb(cbind(yesMaybe,failMaybe) ~ sr * totN + (1 | plot),family="binomial", data=tetr_dat)
tetr_dens_select  <- AICtab(tetr_dens,weights=T)

# Model 4 has ~100% support
tetr_flowN_avg  <- model_avg(tetr_flowN_select, tetr_flowN)
tetr_dens_avg   <- model_avg(tetr_dens_select, tetr_dens)

write.csv(tetr_flowN_avg, "Results/VitalRates_3/tetrazolium_best.csv", row.names = F)
write.csv(tetr_dens_avg, "Results/VitalRates_3/tetrazolium_dens_best.csv", row.names = F)


# 2. Viable Seed Number: germination assays (germTot / germFail)---------------------
# Number of flowers and their sex ratio as predictors

# Omit NAs for glmmadmb
germ_dat  <- na.omit(select(viabVr,germTot,germFail,sr_f,totFlow,sr,totN,plot))

germ_flowN      <- list()
germ_flowN[[1]] <- glmmadmb(cbind(germTot,germFail) ~ sr_f + (1 | plot),family="binomial", data=germ_dat)
germ_flowN[[2]] <- glmmadmb(cbind(germTot,germFail) ~ totFlow + (1 | plot),family="binomial", data=germ_dat)
germ_flowN[[3]] <- glmmadmb(cbind(germTot,germFail) ~ sr_f + totFlow + (1 | plot),family="binomial", data=germ_dat)
germ_flowN[[4]] <- glmmadmb(cbind(germTot,germFail) ~ sr_f * totFlow + (1 | plot),family="binomial", data=germ_dat)
germ_flowN_sel  <- AICtab(germ_flowN,weights=T) 

# Planting density and sex ratio as predictors
germ_dens       <- list()
germ_dens[[1]]  <- glmmadmb(cbind(germTot,germFail) ~ sr + (1 | plot),family="binomial", data=germ_dat)
germ_dens[[2]]  <- glmmadmb(cbind(germTot,germFail) ~ totN + (1 | plot),family="binomial", data=germ_dat)
germ_dens[[3]]  <- glmmadmb(cbind(germTot,germFail) ~ sr + totN + (1 | plot),family="binomial", data=germ_dat)
germ_dens[[4]]  <- glmmadmb(cbind(germTot,germFail) ~ sr * totN + (1 | plot),family="binomial", data=germ_dat)
germ_dens_sel   <- AICtab(germ_dens,weights=T) 

# Average best two models
germ_flowN_avg  <- model_avg(germ_flowN_sel, germ_flowN)
germ_dens_avg   <- model_avg(germ_dens_sel, germ_dens)

write.csv(germ_flowN_avg, "Results/VitalRates_3/germination_best.csv", row.names = F)
write.csv(germ_dens_avg, "Results/VitalRates_3/germination_dens_best.csv", row.names = F)


# Graphs ---------------------------------------------------------------------
tiff("Results/VitalRates_3/viability.tiff",unit="in",
     width=6.3,height=4,res=600,compression="lzw")

lower=quantile(viabVr$totN,prob=c(0.1))
upper=quantile(viabVr$totN,prob=c(0.9))

par(mfrow=c(1,2),mar=c(2.5,2.5,1.1,0.1),mgp=c(1.5,0.6,0),
    oma=c(0,0,2.5,0.3),xpd=NA)
titlePlace=0.2

# Tetrazolium data
plot(jitter(viabVr$totN, factor = 2), jitter(viabVr$tetra_maybe_ratio, factor = 2), pch=1, ylim = c(0,1),
     cex = viabVr$sr * 2, xlab="Planting density",ylab="Seed viability rate")
title("Tetrazolium data", line = titlePlace)
xSeq=seq(0,48,by = 1)
beta=tetr_dens_avg$avg
yMeanLow=inv.logit(beta[1] + beta[2]*0.1 + beta[3]*xSeq + beta[4]*xSeq*0.1)
yMeanHigh=inv.logit(beta[1] + beta[2]*1 + beta[3]*xSeq + beta[4]*xSeq*1)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

legend(0,1.35, c("100% female plot","10% female plot"),
       pch = 1, pt.cex = c(1,0.1)*1.5, bty="n")

# Germination data
plot(jitter(viabVr$totN, factor = 2),jitter(viabVr$germ_ratio, factor = 2),pch=1, ylim = c(0,1),
     cex = viabVr$sr * 2, xlab="Planting density",ylab="Seed germination rate")
title("Germination data", line = titlePlace)
xSeq=seq(0,48,by = 1)
beta=germ_dens_avg$avg
yMeanLow=inv.logit(beta[1] + beta[2]*0.1 + beta[3]*xSeq + beta[4]*xSeq*0.1)
yMeanHigh=inv.logit(beta[1] + beta[2]*1 + beta[3]*xSeq + beta[4]*xSeq*1)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

legend(0,1.35, c("100% female plots","10% female plots"),
       lty = c(1,2), lwd = 2, bty="n")

dev.off()
