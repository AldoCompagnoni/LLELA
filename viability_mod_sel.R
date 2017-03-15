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
viabVr$totN <- viabVr$F + viabVr$M
viabVr$sr   <- viabVr$F / viabVr$totN
viabVr$plot <- as.factor(viabVr$plot)
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
germ_flowN      <- list()
germ_flowN[[1]] <- glmmadmb(cbind(germTot,germFail) ~ (1 | plot),family="binomial", data=germ_dat)
germ_flowN[[2]] <- glmmadmb(cbind(germTot,germFail) ~ sr_f + (1 | plot),family="binomial", data=germ_dat)
germ_flowN[[3]] <- glmmadmb(cbind(germTot,germFail) ~ totFlow + (1 | plot),family="binomial", data=germ_dat)
germ_flowN[[4]] <- glmmadmb(cbind(germTot,germFail) ~ sr_f + totFlow + (1 | plot),family="binomial", data=germ_dat)
germ_flowN[[5]] <- glmmadmb(cbind(germTot,germFail) ~ sr_f * totFlow + (1 | plot),family="binomial", data=germ_dat)

# Model averagin models
germ_flowN_sel  <- AICtab(germ_flowN,weights=T,sort=F) 
germ_flowN_avg  <- model_avg(germ_flowN_sel, germ_flowN)
write.csv(germ_flowN_avg, "Results/VitalRates_3/germination_best.csv", row.names = F)

# Model selection table
write.csv(sel_results(germ_flowN_sel, 5, "Seeds germination rate"),
          "Results/VitalRates_3/germination_mod_sel.csv",row.names=F)


# 2. viable Seed Number: tetrazolium assays  (Yes / fail) ---------------------------

# Number of flowers and their sex ratio as predictors
tetr_flowN=list()
tetr_flowN[[1]]  <- glmmadmb(cbind(Yes,fail) ~ (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN[[2]]  <- glmmadmb(cbind(Yes,fail) ~ sr_f + (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN[[3]]  <- glmmadmb(cbind(Yes,fail) ~ totFlow + (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN[[4]]  <- glmmadmb(cbind(Yes,fail) ~ sr_f + totFlow + (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN[[5]]  <- glmmadmb(cbind(Yes,fail) ~ sr_f * totFlow + (1 | plot),family="binomial", data=tetr_dat)

# Model averagin models
tetr_flowN_select <- AICtab(tetr_flowN,weights=T)
tetr_flowN_avg    <- model_avg(tetr_flowN_select, tetr_flowN)
write.csv(tetr_flowN_avg, "Results/VitalRates_3/tetrazolium_best.csv", row.names = F)

# Model selection table
write.csv(sel_results(tetr_flowN_select, 5, "Seeds viability rate"),
          "tetrazolium_mod_sel.csv",row.names=F)


# Graphs ---------------------------------------------------------------------

# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(topo.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

# graph
#tiff("Results/VitalRates_3/viability.tiff",unit="in",width=6.3,height=3.15,
#     res=600,compression="lzw")

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
#size=mean(viabVr$log_l_t0)
beta=tetr_flowN_avg$avg
yMeanLow=inv.logit(beta[1] + beta[2]*0.2 + beta[3]*xSeq*0.2 + beta[4]*xSeq)
yMeanHigh=inv.logit(beta[1] + beta[2]*1 + beta[3]*xSeq*1 + beta[4]*xSeq)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

# Germination data
plot(jitter(viabVr$totFlow, factor = 2),jitter(viabVr$germ_ratio, factor = 2),
     pch=21, ylim = c(0,1), bg = cRamp(viabVr$sr), cex = 1.5, 
     xlab="Planting density",ylab="Seed germination rate")
title("Germination data", line = titlePlace)
xSeq=seq(0,54,by = 1)
#size=median(viabVr$log_l_t0)
beta=germ_flowN_avg$avg
yMeanLow=inv.logit(beta[1] + beta[2]*0.2 + beta[3]*xSeq*0.2 + beta[4]*xSeq)
yMeanHigh=inv.logit(beta[1] + beta[2]*1 + beta[3]*xSeq*1 + beta[4]*xSeq)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

#dev.off()

