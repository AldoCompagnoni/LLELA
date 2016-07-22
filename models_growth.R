##1. Selection of models of POAR growth using data from LLELA's competition experiment.  
##2. Year specific data anlayzed SEPARATELY.
##3. BEWARE:  the fields "c_t0","fC_t0", and "mC_t0" refer to FALL 2013, even when "year==2015".
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) #For AIC weights, because I'm lazy!
library(lme4)
library(nlme)

#read in data
d=read.csv("Analysis/vr.csv")

#remove died individuals (this is a GROWTH model!)
d=subset(d, surv_t1 != 0)

#Remove a clear mistake: missing data in 2013 and 2014, but 19 leaves in 2015
d=d[-which(d$plot ==36 & d$focalI =="m3"),]

#logtransform leaf numbers
d$log_l_t0=log(d$l_t0)
d$log_l_t1=log(d$l_t1)

d$mC_t0[is.na(d$mC_t0)]=0
d$fC_t0[is.na(d$fC_t0)]=0

#Quadratic terms
#d$M2=d$M^2
#d$F2=d$F^2
#d$mC_t0_2=d$mC_t0^2
#d$fC_t0_2=d$fC_t0^2


#YEAR TWO##################################################################################################
d$plot=as.factor(d$plot)
dat=subset(d,year==2015)

#NEGATIE BINOMIAL
library(glmmADMB)
l2=list()
l2[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=dat,family="nbinom2")
l2[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=dat,family="nbinom2")
l2[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=dat,family="nbinom2")
l2[[3]]=glmmadmb(l_t1 ~ log_l_t0 + log_l_t0:sex + (1 | plot),data=dat,family="nbinom2")
l2[[4]]=glmmadmb(l_t1 ~ log_l_t0 + mC_t0 + fC_t0 + (1 | plot),data=dat,family="nbinom2")
l2[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + mC_t0 + fC_t0 + (1 | plot),data=dat,family="nbinom2")
l2[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + mC_t0 + fC_t0 + (1 | plot),data=dat,family="nbinom2")
l2[[7]]=glmmadmb(l_t1 ~ log_l_t0 + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=dat,family="nbinom2")
l2[[8]]=glmmadmb(l_t1 ~ log_l_t0 + sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=dat,family="nbinom2")
l2[[9]]=glmmadmb(l_t1 ~ log_l_t0 * sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=dat,family="nbinom2")
AICtab(l2,weights=T)


tiff("Results/VitalRates_simple/growth.tiff",unit="in",width=4,height=8,res=600,compression="lzw")

dat$col=as.character(dat$sex)
dat$col[dat$col=="f"]="blue"
dat$col[dat$col=="m"]="red"

par(mfrow=c(2,1),mar=c(3,3,1,0.1),mgp=c(1.4,0.5,0))

#best models
plot(dat$mC_t0,dat$l_t1,pch=16,xlab="N. of total plot leaves 2013",
     ylab="log(N. individual leaves 2014)")
par(new=T) ;  plot(dat$l_t1 ~ dat$fC_t0,pch=16,xlab="",ylab="",xaxt="n")
title(main = "2015: best mod (53% weight)", line=0.2,cex=0.9)
legend(170,86,c("Male competition","Female competition"),
       col=c("grey50"),lty=c(1,2),lwd=3,bty="n")

xSeq <- seq(0,max(c(dat$fC_t0,dat$mC_t0)),by=1)
y_m <- exp(coef(l2[[4]])[1] + coef(l2[[4]])[2]*mean(dat$log_l_t0,na.rm=T) + 
              coef(l2[[4]])[3]*xSeq + coef(l2[[4]])[4]*mean(dat$fC_t0,na.rm=T))
y_f <- exp(coef(l2[[4]])[1] + coef(l2[[4]])[2]*mean(dat$log_l_t0,na.rm=T) + 
              coef(l2[[4]])[4]*xSeq + coef(l2[[4]])[3]*mean(dat$mC_t0,na.rm=T))
lines(xSeq,y_m,lwd=3,lty=1,col="grey50")
lines(xSeq,y_f,lwd=3,lty=2,col="grey50")


#2nd best models
plot(dat$mC_t0,dat$l_t1,pch=16,xlab="N. of total plot leaves 2013",
     ylab="log(N. individual leaves 2014)",col="red")
par(new=T) ;  plot(dat$l_t1 ~ dat$fC_t0,pch=16,xlab="",ylab="",col="blue",xaxt="n")
title(main = "2015: 2nd best (24% weight)", line=0.2,cex=0.9)
legend(80,85,c("Focal male, male comp.",
                 "Focal female, male comp.",
                 "Focal male, female comp.",
                 "Focal female, female comp."),
       col=c("red","blue"),lty=c(1,1,2,2),lwd=2,bty="n")

xSeq <- seq(0,max(c(dat$fC_t0,dat$mC_t0)),by=1)
yFem_m <- exp(fixef(l2[[5]])[1] + fixef(l2[[5]])[2]*mean(dat$log_l_t0,na.rm=T) + 
              fixef(l2[[5]])[4]*xSeq + fixef(l2[[5]])[5]*mean(dat$fC_t0,na.rm=T))
yFem_f <- exp(fixef(l2[[5]])[1] + fixef(l2[[5]])[2]*mean(dat$log_l_t0,na.rm=T) + 
              fixef(l2[[5]])[5]*xSeq + fixef(l2[[5]])[4]*mean(dat$mC_t0,na.rm=T))
yMal_m <- exp(fixef(l2[[5]])[1] + fixef(l2[[5]])[2]*mean(dat$log_l_t0,na.rm=T) + fixef(l2[[5]])[3] +
                fixef(l2[[5]])[4]*xSeq + fixef(l2[[5]])[5]*mean(dat$fC_t0,na.rm=T))
yMal_f <- exp(fixef(l2[[5]])[1] + fixef(l2[[5]])[2]*mean(dat$log_l_t0,na.rm=T) + fixef(l2[[5]])[3] +
                fixef(l2[[5]])[5]*xSeq + fixef(l2[[5]])[4]*mean(dat$mC_t0,na.rm=T))
lines(xSeq,yFem_m,lwd=2,col="blue")
lines(xSeq,yFem_f,lwd=2,lty=2,col="blue")
lines(xSeq,yMal_m,lwd=2,col="red")
lines(xSeq,yMal_f,lwd=2,lty=2,col="red")

dev.off()

