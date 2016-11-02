##Sex predictor of flowering probability: SEX RATIO (female individuals/total individuals)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) #For AIC weights, because I'm lazy!
library(lme4) ; library(glmmADMB)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")

# read in data
d <- read.csv("Data/vr.csv")

# make "plot" a factor to fit models 
d$plot <- as.factor(d$plot) 


# Select YEAR ONE data ------------------------------------------------------------------------------------
tmp                               <- subset(d,year==2014)
tmp[which(tmp$log_l_t1==-Inf),]   <- NA
#prepare the 'new_t1' column
tmp$new_t1[tmp$new_t1=="SKIPPED"] <- NA
tmp$new_t1[tmp$new_t1=="cnf"]     <- NA
tmp$new_t1[tmp$new_t1==""]        <- NA
tmp$new_t1                        <- as.numeric(as.character(tmp$new_t1))
f14                               <- na.omit(tmp[,c("plot","flowN_t1","log_l_t0","sex","sr","TotDensity")])

# Model selection ----------------------------------------------------------------------------------

nfMod=list()
#Target fitness
nfMod[[1]]=glmmadmb(flowN_t1 ~ log_l_t0 + (1 | plot),data=f14,family="poisson")
nfMod[[2]]=glmmadmb(flowN_t1 ~ log_l_t0 + sex + (1 | plot),data=f14,family="poisson")
nfMod[[3]]=glmmadmb(flowN_t1 ~ log_l_t0 * sex + (1 | plot),data=f14,family="poisson")
#Target fitness + effect = tot density 
nfMod[[4]]=glmmadmb(flowN_t1 ~ log_l_t0 + TotDensity + (1 | plot),data=f14,family="poisson")
nfMod[[5]]=glmmadmb(flowN_t1 ~ log_l_t0 + sex + TotDensity + (1 | plot),data=f14,family="poisson")
nfMod[[6]]=glmmadmb(flowN_t1 ~ log_l_t0 * sex + TotDensity + (1 | plot),data=f14,family="poisson")
#Target fitness + effect = tot density + response=sex 
nfMod[[7]]=glmmadmb(flowN_t1 ~ log_l_t0 + TotDensity + TotDensity:sex + (1 | plot),data=f14,family="poisson")
nfMod[[8]]=glmmadmb(flowN_t1 ~ log_l_t0 + sex*TotDensity + (1 | plot),data=f14,family="poisson")
nfMod[[9]]=glmmadmb(flowN_t1 ~ log_l_t0 * sex + TotDensity + TotDensity:sex + (1 | plot),data=f14,family="poisson")
#Target fitness + effect = sex 
nfMod[[10]]=glmmadmb(flowN_t1 ~ log_l_t0 + sr + TotDensity + (1 | plot),data=f14,family="poisson")
nfMod[[11]]=glmmadmb(flowN_t1 ~ log_l_t0 + sex + sr + TotDensity + (1 | plot),data=f14,family="poisson")
nfMod[[12]]=glmmadmb(flowN_t1 ~ log_l_t0 * sex + sr + TotDensity + (1 | plot),data=f14,family="poisson")
#Target fitness + effect = sex + response=sex 
nfMod[[13]]=glmmadmb(flowN_t1 ~ log_l_t0 + sr + TotDensity + sr:sex + TotDensity:sex +(1 | plot),data=f14,family="poisson")
nfMod[[14]]=glmmadmb(flowN_t1 ~ log_l_t0 + sex + sr + TotDensity + sr:sex + TotDensity:sex + (1 | plot),data=f14,family="poisson")
nfMod[[15]]=glmmadmb(flowN_t1 ~ log_l_t0 * sex + sr + TotDensity + sr:sex + TotDensity:sex  + (1 | plot),data=f14,family="poisson")

# Model average
n_flow_select <- AICtab(nfMod,weights=T)
n_flow_avg    <- model_avg(n_flow_select, nfMod)
write.csv(n_flow_avg, "Results/VitalRates_3/n_flowers_best.csv", row.names = F)


# GRAPHS ---------------------------------------------------------------------------------------------------------------

tiff("Results/VitalRates_3/flowering.tiff",unit="in",width=3.5,height=3.5,res=600,compression="lzw")

par(mfrow=c(1,1),mar=c(3,3,0.1,0.1),mgp=c(1.4,0.35,0),cex.lab=0.8,cex.axis=0.8,
    cex.main=0.9)

# Set up colors for plots
f14$col <- as.character(factor(as.integer(f14$sex),labels=c("blue","red")))
mal14   <- subset(f14,sex=="m")
fem14   <- subset(f14,sex=="f")

plot( fem14$flowN_t1+0.05 ~ fem14$sr, pch = 16,ylim=c(0,7),xlim=c(0,1),
      ylab = "Number of flowers", xlab="Proportion of female individuals", col = "blue")
par(new=T) ; plot( mal14$flowN_t1-0.05 ~ mal14$sr,pch = 16,ylim=c(0,7),xlim=c(0,1),
                  ylab = "Number of flowers", xlab="",col = "red")

low   <- 5
high  <- 42
size  <- mean(f14$log_l_t0)
xSeq  <- seq(0,1,by = 0.1)
beta  <- n_flow_avg[,c("predictor","avg")]$avg

y_m_l <- exp( beta[1] + beta[2]*size + beta[3]*low + 
              beta[4]*low + beta[5]*xSeq + beta[6] + 
              beta[7]*size + beta[8]*xSeq)
y_m_h <- exp( beta[1] + beta[2]*size + beta[3]*high + 
              beta[4]*high + beta[5]*xSeq + beta[6] + 
              beta[7]*size + beta[8]*xSeq)
y_f_l <- exp( beta[1] + beta[2]*size + beta[3]*low + 
              beta[5]*xSeq)
y_f_h <- exp( beta[1] + beta[2]*size + beta[3]*high + 
              beta[5]*xSeq)

lines(xSeq,y_m_l,lty=2,lwd=2,col="red")
lines(xSeq,y_m_h,lty=1,lwd=2,col="red")
lines(xSeq,y_f_l,lty=2,lwd=2,col="blue")
lines(xSeq,y_f_h,lty=1,lwd=2,col="blue")

legend(0,7.5,c("high density","low density","male","female"),
       lty=c(1,2,1,1),lwd=2,col=c("black","black","red","blue"),bty="n")

dev.off()
