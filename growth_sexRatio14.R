##Growth analyses using one-sex plots only
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) #For AIC weights, because I'm lazy!
library(glmmADMB) #Fit models with a Negative Binomial
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")

#read in data
d=read.csv("Data/vr.csv")
#remove dead individuals (this is a GROWTH model!)
d=subset(d, surv_t1 != 0)

#logtransform leaf numbers
d$mC_t0[is.na(d$mC_t0)]=0
d$fC_t0[is.na(d$fC_t0)]=0
d$plot=as.factor(d$plot) #glmmadmb wants plot as a factor


##################################################################################################
#1.Compare total growth between female only and male only plots#########
##################################################################################################
d14=subset(d, year==2014)


#2014----------------------------------------------------------------------------------------------  

# Sex ratio as original sex ratio
lMod=list() # lMod stands for "leaf model" (density is quantified by N. of leaves)
lMod[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=d14,family="nbinom2")
lMod[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=d14,family="nbinom2")
lMod[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=d14,family="nbinom2")
# Target fitness + effect = tot density 
lMod[[4]]=glmmadmb(l_t1 ~ log_l_t0 + TotDensity + (1 | plot),data=d14,family="nbinom2")
lMod[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + TotDensity + (1 | plot),data=d14,family="nbinom2")
lMod[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + TotDensity + (1 | plot),data=d14,family="nbinom2")
# Target fitness + effect = tot density + response=sex 
lMod[[7]]=glmmadmb(l_t1 ~ log_l_t0 + TotDensity + sex:TotDensity + (1 | plot),data=d14,family="nbinom2")
lMod[[8]]=glmmadmb(l_t1 ~ log_l_t0 + sex + TotDensity + TotDensity:sex + (1 | plot),data=d14,family="nbinom2")
lMod[[9]]=glmmadmb(l_t1 ~ log_l_t0 * sex + TotDensity + TotDensity:sex + (1 | plot),data=d14,family="nbinom2")
# Target fitness + effect = sex 
lMod[[10]]=glmmadmb(l_t1 ~ log_l_t0 + sr + TotDensity + (1 | plot),data=d14,family="nbinom2")
lMod[[11]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + TotDensity + (1 | plot),data=d14,family="nbinom2")
lMod[[12]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + TotDensity + (1 | plot),data=d14,family="nbinom2")
# Target fitness + effect = sex + response=sex 
lMod[[13]]=glmmadmb(l_t1 ~ log_l_t0 + sr + TotDensity + sr:sex + TotDensity:sex +(1 | plot),data=d14,family="nbinom2")
lMod[[14]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + TotDensity + sr:sex + TotDensity:sex + (1 | plot),data=d14,family="nbinom2")
lMod[[15]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + TotDensity + sr:sex + TotDensity:sex  + (1 | plot),data=d14,family="nbinom2")

gr14_mod_sel  <- AICtab(lMod, weights = T)
gr14_avg      <- model_avg(gr14_mod_sel, lMod)

select(gr14_avg,predictor,avg)


# GRAPH ----------------------------------------------------------------------

# Set up colors for plots
# 2014
d14$col=as.integer(d14$sex)
d14$col=as.character(factor(d14$col,labels=c("blue","red")))
d14$symb=as.integer(as.character(factor(d14$col,labels=c("17","16"))))



par(mar=c(2.8,2.5,0.1,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0.1))

plot(d14$TotDensity, d14$l_t1, pch=16, ylab="Target individuals: Number of leaves",
     xlab="Proportion of female individuals",col=d14$col, ylim = c(0,82))
beta = gr14_avg[,c("predictor","avg")]$avg
xSeq <- seq(0,48,by=1)
sr   <- 0.5
size <- mean(d14$log_l_t0)
y_m <- exp( beta[1] + beta[2]*size + beta[3] + beta[4]*xSeq + 
              beta[5]*xSeq + beta[6]*sr + beta[7]*size)
y_f <- exp( beta[1] + beta[2]*size + beta[4]*xSeq + beta[6]*sr )
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")


# Actual graphs
par(mar=c(2.8,2.5,0.1,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0.1))

plot(d14$sr,d14$l_t1, pch=16, ylab="Target individuals: Number of leaves",
     xlab="Proportion of female individuals",col=d14$col)
beta = gr14_avg[,c("predictor","avg")]$avg
xSeq <- seq(0,1,by=0.1)
size <- mean(d14$log_l_t0)
dens <- 1
y_m_l <- exp( beta[1] + beta[2]*size + beta[3] + beta[4]*dens + 
              beta[5]*dens + beta[6]*xSeq + beta[7]*size)
y_f_l <- exp( beta[1] + beta[2]*size + beta[4]*dens + beta[6]*xSeq)
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")
dens <- 48
y_m_h <- exp( beta[1] + beta[2]*size + beta[3] + beta[4]*dens + 
                beta[5]*dens + beta[6]*xSeq + beta[7]*size)
y_f_h <- exp( beta[1] + beta[2]*size + beta[4]*dens + beta[6]*xSeq)
lines(xSeq,y_f_l,lwd=3,lty=2,col="blue")
lines(xSeq,y_m_l,lwd=3,lty=2,col="red")

lines(xSeq,y_f_h,lwd=3,lty=1,col="blue")
lines(xSeq,y_m_h,lwd=3,lty=1,col="red")

legend(-0.05,88,c("Males","Females"),
       lty=1,lwd=3,col=c("red","blue"),bty="n")







# Actual graphs
par(mar=c(2.8,2.5,0.1,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0.1))

plot(d14$TotDensity,d14$l_t1, pch=16, ylab="Target individuals: Number of leaves",
     xlab="Proportion of female individuals",col=d14$col)
beta = gr14_avg[,c("predictor","avg")]$avg
xSeq <- seq(0,1,by=0.1)
size <- mean(d14$log_l_t0)
dens <- 1
y_m_l <- exp( beta[1] + beta[2]*size + beta[3] + beta[4]*dens + 
                beta[5]*dens + beta[6]*xSeq + beta[7]*size)
y_f_l <- exp( beta[1] + beta[2]*size + beta[4]*dens + beta[6]*xSeq)
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")
dens <- 48
y_m_h <- exp( beta[1] + beta[2]*size + beta[3] + beta[4]*dens + 
                beta[5]*dens + beta[6]*xSeq + beta[7]*size)
y_f_h <- exp( beta[1] + beta[2]*size + beta[4]*dens + beta[6]*xSeq)
lines(xSeq,y_f_l,lwd=3,lty=2,col="blue")
lines(xSeq,y_m_l,lwd=3,lty=2,col="red")

lines(xSeq,y_f_h,lwd=3,lty=1,col="blue")
lines(xSeq,y_m_h,lwd=3,lty=1,col="red")

legend(-0.05,88,c("Males","Females"),
       lty=1,lwd=3,col=c("red","blue"),bty="n")
