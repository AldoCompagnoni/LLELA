# Model selection for seed VIABILITY
# Data: tetrazolium scoring and germination data 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(lme4) ; library(bbmle) ; library(boot) ; library(testthat)
library(dplyr) ; library(glmmADMB)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")

# Read in data and format -----------------------------------------------------------
d           <- read.csv("Data/vr.csv")
viabVr      <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                        stringsAsFactors = F)

# Format data -----------------------------------------------------------------------
viabVr$totN <- viabVr$F + viabVr$M
viabVr$sr   <- viabVr$F / viabVr$totN
viabVr$plot <- as.factor(viabVr$plot)
viabVr      <- subset(viabVr, totFlow < 60) #exclude extreme values



# Viability model selection ---------------------------------------------------------

# 1. viable Seed Number: tetrazolium assays  (Yes / fail) ---------------------------

# Omit NAs for glmmadmb
tetr_dat  <- na.omit(dplyr::select(viabVr,Yes,fail,sr_f,totFlow,sr,totN,plot,log_l_t0))

# Number of flowers and their sex ratio as predictors
tetr_flowN=list()
tetr_flowN[[1]]  <- glmmadmb(cbind(Yes,fail) ~ sr_f + log_l_t0 + (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN[[2]]  <- glmmadmb(cbind(Yes,fail) ~ totFlow + log_l_t0 + (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN[[3]]  <- glmmadmb(cbind(Yes,fail) ~ sr_f + totFlow + log_l_t0 + (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN[[4]]  <- glmmadmb(cbind(Yes,fail) ~ sr_f * totFlow + log_l_t0 + (1 | plot),family="binomial", data=tetr_dat)

# Model averagin models
tetr_flowN_select<- AICtab(tetr_flowN,weights=T)
tetr_flowN_avg  <- model_avg(tetr_flowN_select, tetr_flowN)
write.csv(tetr_flowN_avg, "Results/VitalRates_3/tetrazolium_best.csv", row.names = F)


# 2. Viable Seed Number: germination assays (germTot / germFail)---------------------

# Omit NAs for glmmadmb
germ_dat  <- na.omit(dplyr::select(viabVr,germTot,germFail,sr_f,totFlow,sr,totN,plot,log_l_t0))

germ_flowN      <- list()
germ_flowN[[1]] <- glmmadmb(cbind(germTot,germFail) ~ log_l_t0 + (1 | plot),family="binomial", data=germ_dat)
germ_flowN[[2]] <- glmmadmb(cbind(germTot,germFail) ~ sr_f + log_l_t0 + (1 | plot),family="binomial", data=germ_dat)
germ_flowN[[3]] <- glmmadmb(cbind(germTot,germFail) ~ totFlow + log_l_t0 + (1 | plot),family="binomial", data=germ_dat)
germ_flowN[[4]] <- glmmadmb(cbind(germTot,germFail) ~ sr_f + totFlow + log_l_t0 + (1 | plot),family="binomial", data=germ_dat)
germ_flowN[[5]] <- glmmadmb(cbind(germTot,germFail) ~ sr_f * totFlow + log_l_t0 + (1 | plot),family="binomial", data=germ_dat)

# Model averagin models
germ_flowN_sel  <- AICtab(germ_flowN,weights=T) 
germ_flowN_avg  <- model_avg(germ_flowN_sel, germ_flowN)
write.csv(germ_flowN_avg, "Results/VitalRates_3/germination_best.csv", row.names = F)


# Graphs ---------------------------------------------------------------------

# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(topo.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

# graph
tiff("Results/VitalRates_3/viability.tiff",unit="in",width=6.3,height=3.15,
     res=600,compression="lzw")

lower=quantile(viabVr$totN,prob=c(0.1))
upper=quantile(viabVr$totN,prob=c(0.9))

par(mfrow=c(1,2),mar=c(2.5,2.5,1.1,0.1),mgp=c(1.5,0.6,0),
    oma=c(0,0,0,0),xpd=NA)
titlePlace=0.2

# Tetrazolium data
plot(jitter(viabVr$totFlow, factor = 2), jitter(viabVr$tetra_ratio, factor = 2), 
     pch=21, ylim = c(0,1), bg = cRamp(viabVr$sr), cex = 1.5, 
     xlab="Planting density",ylab="Seed viability rate")
title("Tetrazolium data", line = titlePlace)
xSeq=seq(0,54,by = 1)
size=mean(viabVr$log_l_t0)
beta=tetr_flowN_avg$avg
yMeanLow=inv.logit(beta[1] + beta[2]*size + beta[3]*0.1 + beta[4]*xSeq*0.1 + beta[5]*xSeq)
yMeanHigh=inv.logit(beta[1] + beta[2]*size + beta[3]*1 + beta[4]*xSeq*1 + beta[5]*xSeq)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

# Germination data
plot(jitter(viabVr$totFlow, factor = 2),jitter(viabVr$germ_ratio, factor = 2),
     pch=21, ylim = c(0,1), bg = cRamp(viabVr$sr), cex = 1.5, 
     xlab="Planting density",ylab="Seed germination rate")
title("Germination data", line = titlePlace)
xSeq=seq(0,54,by = 1)
size=median(viabVr$log_l_t0)
beta=germ_flowN_avg$avg
yMeanLow=inv.logit(beta[1] + beta[2]*size + beta[3]*0.1 + beta[4]*xSeq*0.1 + beta[5]*xSeq)
yMeanHigh=inv.logit(beta[1] + beta[2]*size + beta[3]*1 + beta[4]*xSeq*1 + beta[5]*xSeq)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

dev.off()

