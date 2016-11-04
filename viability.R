# Model selection for seed VIABILITY
# Data: tetrazolium scoring and germination data 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(lme4) ; library(bbmle) ; library(boot) ; library(testthat)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")

# Read in data and format=====================================================================================
viabVr <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                   stringsAsFactors = F)

# Viability model selection ---------------------------------------------------------

viabVr$totN <- viabVr$F + viabVr$M
viabVr$sr   <- viabVr$F / viabVr$totN

# 1. viable Seed Number: tetrazolium assays  (Yes / fail) ---------------------------
# This is the "standard" scoring: consider as viable any seed stained by tetrazolium
#l_tetr_maybe=list()
#l_tetr_maybe[[1]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f,family="binomial", data=viabVr)
#l_tetr_maybe[[2]]=glm(cbind(yesMaybe,failMaybe) ~ totFlow,family="binomial", data=viabVr)
#l_tetr_maybe[[3]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f + totFlow,family="binomial", data=viabVr)
#l_tetr_maybe[[4]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f * totFlow,family="binomial", data=viabVr)
#tetr_select      <- AICtab(l_tetr_maybe,weights=T)

l_tetr_maybe=list()
l_tetr_maybe[[1]]=glm(cbind(yesMaybe,failMaybe) ~ sr,family="binomial", data=viabVr)
l_tetr_maybe[[2]]=glm(cbind(yesMaybe,failMaybe) ~ totN,family="binomial", data=viabVr)
l_tetr_maybe[[3]]=glm(cbind(yesMaybe,failMaybe) ~ sr + totN,family="binomial", data=viabVr)
l_tetr_maybe[[4]]=glm(cbind(yesMaybe,failMaybe) ~ sr * totN,family="binomial", data=viabVr)
tetr_select      <- AICtab(l_tetr_maybe,weights=T)

# Model 4 has ~100% support
tetr_avg <- model_avg(tetr_select, l_tetr_maybe)
write.csv(tetr_avg, "Results/VitalRates_3/tetrazolium_best.csv", row.names = F)


# 2. Viable Seed Number: germination assays (germTot / germFail)---------------------
#l_viab_germ=list()
#l_viab_germ[[1]]=glm(cbind(germTot,germFail) ~ sr_f,family="binomial", data=viabVr)
#l_viab_germ[[2]]=glm(cbind(germTot,germFail) ~ totFlow,family="binomial", data=viabVr)
#l_viab_germ[[3]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow,family="binomial", data=viabVr)
#l_viab_germ[[4]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow,family="binomial", data=viabVr)
#germ_select     <- AICtab(l_viab_germ,weights=T) 

l_viab_germ=list()
l_viab_germ[[1]]=glm(cbind(germTot,germFail) ~ sr,family="binomial", data=viabVr)
l_viab_germ[[2]]=glm(cbind(germTot,germFail) ~ totN,family="binomial", data=viabVr)
l_viab_germ[[3]]=glm(cbind(germTot,germFail) ~ sr + totN,family="binomial", data=viabVr)
l_viab_germ[[4]]=glm(cbind(germTot,germFail) ~ sr * totN,family="binomial", data=viabVr)
germ_select     <- AICtab(l_viab_germ,weights=T) 


# Average best two models
germ_avg <- model_avg(germ_select, l_viab_germ)
write.csv(germ_avg, "Results/VitalRates_3/germination_best.csv", row.names = F)


#Graphs ---------------------------------------------------------------------
tiff("Results/VitalRates_3/viability.tiff",unit="in",
     width=6.3,height=4,res=600,compression="lzw")

lower=quantile(viabVr$totN,prob=c(0.1))
upper=quantile(viabVr$totN,prob=c(0.9))

par(mfrow=c(1,2),mar=c(2.5,2.5,1.1,0.1),mgp=c(1.5,0.6,0),
    oma=c(0,0,2.5,0),xpd=NA)
titlePlace=0.2

# Tetrazolium data
plot(viabVr$totN, viabVr$tetra_maybe_ratio, pch=1, ylim = c(0,1),
     cex = viabVr$sr * 2, xlab="Planting density",ylab="Seed viability rate")
title("Tetrazolium data", line = titlePlace)
xSeq=seq(0,48,by = 1)
beta=tetr_avg$avg
yMeanLow=inv.logit(beta[1] + beta[2]*0.1 + beta[3]*xSeq + beta[4]*xSeq*0.1)
yMeanHigh=inv.logit(beta[1] + beta[2]*0.9 + beta[3]*xSeq + beta[4]*xSeq*0.9)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

legend(0,1.35, c("90% female plot","10% female plot"),
       pch = 1, pt.cex = c(0.9,0.1)*1.5, bty="n")

# Germination data
plot(viabVr$totN,viabVr$germ_ratio,pch=1, ylim = c(0,1),
     cex = viabVr$sr * 2, xlab="Planting density",ylab="Seed germination rate")
title("Germination data", line = titlePlace)
xSeq=seq(0,48,by = 1)
beta=germ_avg$avg
yMeanLow=inv.logit(beta[1] + beta[2]*0.1 + beta[3]*xSeq + beta[4]*xSeq*0.1)
yMeanHigh=inv.logit(beta[1] + beta[2]*0.9 + beta[3]*xSeq + beta[4]*xSeq*0.9)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

legend(0,1.35, c("90% female plots","10% female plots"),
       lty = c(1,2), lwd = 2, bty="n")

dev.off()
