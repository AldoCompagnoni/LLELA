##Sex predictor of flowering probability: SEX RATIO (female individuals/total individuals)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)
library(boot)
library(glmmADMB)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")


# Data --------------------------------------------------------------------------------
pans      <- read.csv("Data/Spring 2014/plot_level_panicles.csv")
pans_f    <- subset(pans, F != 0) # plots containing females
pans_m    <- subset(pans, M != 0) # plots containing males


# Model selection ---------------------------------------------------------------------


# females -----------------------------------------------------------------------------
fMod        <- list()
fMod[[1]]   <- lm(F_flow_t1 ~ 1,     data=pans)
fMod[[2]]   <- lm(F_flow_t1 ~ sr,    data=pans)
fMod[[3]]   <- lm(F_flow_t1 ~ N,     data=pans)
fMod[[4]]   <- lm(F_flow_t1 ~ N + sr,data=pans)
fMod[[5]]   <- lm(F_flow_t1 ~ N * sr,data=pans)

# Model average
f_pan_sel   <- AICtab(fMod,weights=T)
f_pan_avg   <- model_avg(f_pan_sel, fMod)
write.csv(f_pan_avg, "Results/VitalRates_3/f_flowers_plot_best.csv", row.names = F)


# males -------------------------------------------------------------------------------
mMod        <- list()
mMod[[1]]   <- lm(M_flow_t1 ~ 1,     data=pans)
mMod[[2]]   <- lm(M_flow_t1 ~ sr,    data=pans)
mMod[[3]]   <- lm(M_flow_t1 ~ N,     data=pans)
mMod[[4]]   <- lm(M_flow_t1 ~ N + sr,data=pans)
mMod[[5]]   <- lm(M_flow_t1 ~  N * sr,data=pans)

# Model average
m_pan_sel   <- AICtab(mMod,weights=T)
m_pan_avg   <- model_avg(m_pan_sel, mMod)
write.csv(m_pan_avg, "Results/VitalRates_3/m_flowers_plot_best.csv", row.names = F)


# Graph -----------------------------------------------------------------------------

# Plot data
plot(jitter(pans_f$F_flow_t1) ~ pans_f$N, xlab=expression("Planting density"),
     ylab="Number of flowers", pch=16, col="blue", ylim = c(-0.5,90))
par(new=T) ; 
plot(jitter(pans_m$M_flow_t1) ~ pans_m$N,pch=17,xlab="",ylab="",col="red",
     xaxt="n", ylim = c(-0.5,90))

# lines
xSeq  <- seq(0,48,1)
betaF <- f_pan_avg$avg
betaM <- m_pan_avg$avg

yF_h  <- (betaF[1] + betaF[2]*xSeq + betaF[4]*xSeq*1 + betaF[3]*1)
yF_l  <- (betaF[1] + betaF[2]*xSeq + betaF[4]*xSeq*0.2 + betaF[3]*0.2)
yM_h  <- (betaM[1] + betaM[2]*xSeq + betaM[4]*xSeq*1 + betaM[3]*1)
yM_l  <- (betaM[1] + betaM[2]*xSeq + betaM[4]*xSeq*0.2 + betaM[3]*0.2)

lines(xSeq, yF_h, col = "blue", lwd = 2, lty = 1)
lines(xSeq, yF_l, col = "blue", lwd = 2, lty = 2)
lines(xSeq, yM_h, col = "red",  lwd = 2, lty = 1)
lines(xSeq, yM_l, col = "red",  lwd = 2, lty = 2)

