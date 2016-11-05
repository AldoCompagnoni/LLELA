# Model selection for seed germination data
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(lme4) ; library(bbmle) ; library(boot) ; library(testthat)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")


# Read in data and format------------------------------------------------------------
viabVr <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv", 
                   stringsAsFactors = F)


# Viability model selection ---------------------------------------------------------
l_viab_germ=list()
l_viab_germ[[1]]=glm(cbind(germTot,germFail) ~ sr_f,family="binomial", data=viabVr)
l_viab_germ[[2]]=glm(cbind(germTot,germFail) ~ totFlow,family="binomial", data=viabVr)
l_viab_germ[[3]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow,family="binomial", data=viabVr)
l_viab_germ[[4]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow,family="binomial", data=viabVr)
germ_select     <- AICtab(l_viab_germ,weights=T) 

# Average best models
germ_avg <- model_avg(germ_select, l_viab_germ)
write.csv(germ_avg, "Results/VitalRates_3/germination_best.csv", row.names = F)


#Graphs ---------------------------------------------------------------------
tiff("Results/VitalRates_3/germination.tiff",unit="in",
     width=4,height=4,res=600,compression="lzw")

lower=quantile(viabVr$totFlow,prob=c(0.1))
upper=quantile(viabVr$totFlow,prob=c(0.9))

par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.5,0.6,0),xpd=NA,oma=c(0,0,0,0))

# Germination data
plot(viabVr$sr_f,viabVr$germ_ratio,pch=16,
     xlab="Sex ratio",ylab="Seed germination rate")
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
beta=germ_avg$avg
yMeanLow=inv.logit(beta[1] + beta[2]*xSeq + beta[4]*lower + beta[3]*xSeq*lower)
yMeanHigh=inv.logit(beta[1] + beta[2]*xSeq + beta[4]*upper + beta[3]*xSeq*upper)
lines(xSeq,yMeanLow,lwd=2,col="blue")
lines(xSeq,yMeanHigh,lwd=2,col="red")

legend(0.08,0.97,c("High density plot","Low density plot"),lwd=2,lty=1,
       col=c("blue","red"),bty="n")

dev.off()
