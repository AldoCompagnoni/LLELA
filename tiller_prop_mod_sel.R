# Sex competition predictor: SEX RATIO (female individuals/total individuals)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)    
library(glmmADMB) #Fit models with a Negative Binomial
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")


# Read data ------------------------------------------------------------------
x   <- read.csv("Data/vr.csv")

d   <- subset(x, !is.na(till_t1))
# plot as a factor
d   <- mutate(d, plot = as.factor(plot) ) # glmmadmb wants plot as a factor


# proportions of tillers
prop_t    <- function(x){
  
  by_plot   <- split(d, d$plot)
  tot_t     <- lapply(by_plot, function(x) x$tot_till_t1 = sum(x$till_t1) )
  tot_f     <- lapply(by_plot, function(x) sum(subset(x, sex == "f")$till_t1) )
  prop_f    <- data.frame( plot = as.factor(names(tot_t)), tot_f = unlist(tot_f), 
                           tot_t = unlist(tot_t), prop_f = unlist(tot_f)/unlist(tot_t) )
  pror_d    <- merge(unique(dplyr::select(d,plot,TotDensity,sr)), prop_f)
  
  return(pror_d)
  
}



# 2014 --------------------------------------------------------------------------------
x   <- subset(d, year == 2014)
#x <- prop_t(x)
#x <- subset(x, sr != 0)
#d14 <- subset(x, sr != 1)
d14 <- prop_t(x)

# Planting density
t14=list() # lMod stands for "leaf model" (density is quantified by N. of leaves)
t14[[1]]= glm(cbind(tot_f, tot_t-tot_f) ~ sr,data=d14,family="binomial")
t14[[2]]= glm(cbind(tot_f, tot_t-tot_f) ~ TotDensity,data=d14,family="binomial")
t14[[3]]= glm(cbind(tot_f, tot_t-tot_f) ~ TotDensity + sr,data=d14,family="binomial")
t14[[4]]= glm(cbind(tot_f, tot_t-tot_f) ~ TotDensity * sr,data=d14,family="binomial")
AICtab(t14, weights=T)


# 2015 --------------------------------------------------------------------------------
x <- subset(d, year == 2015)
#x <- prop_t(x)
#x = subset(x, sr != 0)
#d15 = subset(x, sr != 1)
d15 <- prop_t(x)

# Planting density
t15=list() # lMod stands for "leaf model" (density is quantified by N. of leaves)
t15[[1]]= glm(cbind(tot_f, tot_t-tot_f) ~ sr,data=d15,family="binomial")
t15[[2]]= glm(cbind(tot_f, tot_t-tot_f) ~ TotDensity,data=d15,family="binomial")
t15[[3]]= glm(cbind(tot_f, tot_t-tot_f) ~ TotDensity + sr,data=d15,family="binomial")
t15[[4]]= glm(cbind(tot_f, tot_t-tot_f) ~ TotDensity * sr,data=d15,family="binomial")
AICtab(t15, weights=T)


# Proportion Graphs ##########################################################
par(mfrow=c(1,2), mar = c(4,4,0.5,0.5))
cexlab = 1.1

#2014 -----------------------------------------------------------------------
scal = ((d14$TotDensity/48)+0.3)*3
plot(d14$sr, d14$prop_f, ylab = expression("Proportion of female tillers"[t+1]),
     cex=scal, cex.lab = cexlab,
     xlab = expression("Proportion of female tillers"[t]),ylim=c(0,1),xlim=c(0,1))

xSeq  <- seq(0,1,0.1) 
beta  <- coef(t14[[4]])
y_l   <- inv.logit( beta[1] + beta[2]*5 + beta[3]*xSeq + beta[4]*5*xSeq )
y_h   <- inv.logit( beta[1] + beta[2]*40 + beta[3]*xSeq + beta[4]*40*xSeq )
lines(xSeq,y_h,lwd=3)
lines(xSeq,y_l,lwd=2,lty=2)
abline(0,1,lty=2)

#2015 -----------------------------------------------------------------------
scal = ((d15$TotDensity/48)+0.3)*3
plot(d15$sr, d15$prop_f, ylab = expression("Proportion of female tillers"[t+1]),
     cex=scal,cex.lab = cexlab,
     xlab = expression("Proportion of female tillers"[t]),ylim=c(0,1),xlim=c(0,1))

xSeq  <- seq(0,1,0.1) 
beta  <- coef(t15[[4]])
y_l   <- inv.logit( beta[1] + beta[2]*5 + beta[3]*xSeq + beta[4]*5*xSeq )
y_h   <- inv.logit( beta[1] + beta[2]*40 + beta[3]*xSeq + beta[4]*40*xSeq )
lines(xSeq,y_h,lwd=3)
lines(xSeq,y_l,lwd=2,lty=2)
abline(0,1,lty=2)


# Density Graphs ##########################################################
tiff("Results/VitalRates_3/tiller_proport.tiff",unit="in",
     width=6.3,height=3.15,res=600,compression="lzw")

par(mfrow=c(1,2), mar = c(4,4,2,0.5), mgp = c(1.5,0.5,0))
cexlab = 1.1

#2014 -----------------------------------------------------------------------
scal <- ((d14$sr/1)+0.1)*2
plot(d14$TotDensity, d14$prop_f, xlab = "Planting density",
     cex=scal,cex.lab=cexlab,
     ylab = expression("Proportion of female tillers"[t+1]))

xSeq  <- seq(0,48,1) 
beta  <- coef(t14[[4]])
y_l   <- inv.logit( beta[1] + beta[2]*xSeq + beta[3]*0.1 + beta[4]*0.1*xSeq )
y_h   <- inv.logit( beta[1] + beta[2]*xSeq + beta[3]*0.9 + beta[4]*0.9*xSeq )
lines(xSeq,y_h,lwd=3)
lines(xSeq,y_l,lwd=2,lty=2)

pt = (c(0.1,0.9)+0.1)*2
legend(0, 1.32,c("10% female plots","90% female plots"),lty=c(2,1),lwd=c(1,2),
       pch=1,xpd=NA, bty = "n", pt.cex = pt) 

#2015 -----------------------------------------------------------------------
scal <- ((d15$sr/1)+0.1)*2
plot(d15$TotDensity, d15$prop_f, xlab = "Planting density",
     cex=scal,cex.lab=cexlab,
     ylab = expression("Proportion of female tillers"[t+1]))

xSeq  <- seq(0,48,1) 
beta  <- coef(t15[[4]])
y_l   <- inv.logit( beta[1] + beta[2]*xSeq + beta[3]*0.1 + beta[4]*0.1*xSeq )
y_h   <- inv.logit( beta[1] + beta[2]*xSeq + beta[3]*0.9 + beta[4]*0.9*xSeq )
lines(xSeq,y_h,lwd=3)
lines(xSeq,y_l,lwd=2,lty=2)
title("2015")

dev.off()
