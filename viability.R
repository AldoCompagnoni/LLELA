# Model selection for seed VIABILITY
# Data: tetrazolium scoring and germination data 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(lme4) ; library(bbmle) ; library(boot) ; library(testthat)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")

# Read in data and format=====================================================================================
viabVr <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                   stringsAsFactors = F)

# Viability model selection ---------------------------------------------------------

# 1. viable Seed Number: tetrazolium assays  (Yes / fail) ---------------------------
# This is the "standard" scoring: consider as viable any seed stained by tetrazolium
l_tetr_maybe=list()
l_tetr_maybe[[1]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f,family="binomial", data=viabVr)
l_tetr_maybe[[2]]=glm(cbind(yesMaybe,failMaybe) ~ totFlow,family="binomial", data=viabVr)
l_tetr_maybe[[3]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f + totFlow,family="binomial", data=viabVr)
l_tetr_maybe[[4]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f * totFlow,family="binomial", data=viabVr)
tetr_select      <- AICtab(l_tetr_maybe,weights=T)

# Model 4 has ~100% support
tetr_avg <- model_avg(tetr_select, l_tetr_maybe)
write.csv(tetr_avg, "Results/VitalRates_3/tetrazolium_best.csv", row.names = F)


#2. Viable Seed Number: germination assays (germTot / germFail)---------------------
l_viab_germ=list()
l_viab_germ[[1]]=glm(cbind(germTot,germFail) ~ sr_f,family="binomial", data=viabVr)
l_viab_germ[[2]]=glm(cbind(germTot,germFail) ~ totFlow,family="binomial", data=viabVr)
l_viab_germ[[3]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow,family="binomial", data=viabVr)
l_viab_germ[[4]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow,family="binomial", data=viabVr)
germ_select     <- AICtab(l_viab_germ,weights=T) 

# Average best two models
germ_avg <- model_avg(germ_select, l_viab_germ)
write.csv(germ_avg, "Results/VitalRates_3/germination_best.csv", row.names = F)


#Graphs ---------------------------------------------------------------------
tiff("Results/VitalRates_3/viability.tiff",unit="in",
     width=6.3,height=4,res=600,compression="lzw")

lower=quantile(viabVr$totFlow,prob=c(0.1))
upper=quantile(viabVr$totFlow,prob=c(0.9))

par(mfrow=c(1,2),mar=c(2.5,2.5,1.1,0.1),mgp=c(1.5,0.6,0),
    oma=c(0,0,3,0),xpd=NA)
titlePlace=0.2

# Tetrazolium data
plot(viabVr$sr_f,viabVr$tetra_maybe_ratio,pch=16,
     xlab="Sex ratio",ylab="Seed viability rate")
title("Tetrazolium data", line = titlePlace)
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
beta=tetr_avg$avg
yMeanLow=inv.logit(beta[1] + beta[2]*xSeq + beta[3]*lower + beta[4]*xSeq*lower)
yMeanHigh=inv.logit(beta[1] + beta[2]*xSeq + beta[3]*upper + beta[4]*xSeq*upper)
lines(xSeq,yMeanLow,lwd=2,col="blue")
lines(xSeq,yMeanHigh,lwd=2,col="red")

# Legend
text(0,1.3,"Tetrazolium model:  Sex ratio X Density",pos=4)
text(0,1.2,"Germination model: Sex ratio + Density",pos=4)

# Germination data
plot(viabVr$sr_f,viabVr$germ_ratio,pch=16,
     xlab="Sex ratio",ylab="Seed germination rate")
title("Germination data", line = titlePlace)
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
beta=germ_avg$avg
yMeanLow=inv.logit(beta[1] + beta[2]*xSeq + beta[4]*lower + beta[3]*xSeq*lower)
yMeanHigh=inv.logit(beta[1] + beta[2]*xSeq + beta[4]*upper + beta[3]*xSeq*upper)
lines(xSeq,yMeanLow,lwd=2,col="blue")
lines(xSeq,yMeanHigh,lwd=2,col="red")

legend(0.18,1.26,c("High density plot","Low density plot"),lwd=2,lty=1,
       col=c("blue","red"),bty="n")

dev.off()
