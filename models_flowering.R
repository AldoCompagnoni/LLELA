##1. Selection of models of POAR flowering using data from LLELA's competition experiment.  
##2. Year specific data anlayzed SEPARATELY.
##3. BEWARE:  the fields "c_t0","fC_t0", and "mC_t0" refer to FALL 2013, even when "year==2015".
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) #For AIC weights, because I'm lazy!
library(lme4)
library(nlme)

#read in data
d=read.csv("Data/vr.csv")

#Remove resuscitated individuals (dead in spring 2014, alive in 2015)
#NOTE: possibly these are NOT DEAD in 2014, because they all have
#few leaves: resprouted from base?
d=d[-which(d$plot ==18 & d$focalI =="m2"),]
d=d[-which(d$plot ==38 & d$focalI =="f1"),]
d=d[-which(d$plot ==46 & d$focalI =="f1"),]
d=d[-which(d$plot ==83 & d$focalI =="f5"),]
d=d[-which(d$plot ==36 & d$focalI =="m3"),]

#logtransform leaf numbers
d$log_l_t0=log(d$l_t0)
d$log_l_t1=log(d$l_t1)
d$plot=as.factor(d$plot)


#YEAR ONE##########################################################################################################
dat=subset(d,year==2014)
dat[which(dat$log_l_t1==-Inf),]=NA

library(glmmADMB)
#Density as "number of leaves"-----------------------------------------------------------
lMod=list() #lMod stands for "leaf model" (density is quantified by N. of leaves)
tmp=na.omit(dat[,c("plot","flow_t1","log_l_t1","sex","c_t0","mC_t0","fC_t0","M","F")])
#Target fitness
lMod[[1]]=glmmadmb(flow_t1 ~ log_l_t1 + (1 | plot),data=tmp,family="binomial")
lMod[[2]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + (1 | plot),data=tmp,family="binomial")
lMod[[3]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = tot density 
lMod[[4]]=glmmadmb(flow_t1 ~ log_l_t1 + c_t0 + (1 | plot),data=tmp,family="binomial")
lMod[[5]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + c_t0 + (1 | plot),data=tmp,family="binomial")
lMod[[6]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + c_t0 + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = tot density + response=sex 
lMod[[7]]=glmmadmb(flow_t1 ~ log_l_t1 + c_t0 + c_t0:sex + (1 | plot),data=tmp,family="binomial")
lMod[[8]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + c_t0 + c_t0:sex + (1 | plot),data=tmp,family="binomial")
lMod[[9]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + c_t0 + c_t0:sex + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = sex 
lMod[[10]]=glmmadmb(flow_t1 ~ log_l_t1 + mC_t0 + fC_t0 + (1 | plot),data=tmp,family="binomial")
lMod[[11]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + mC_t0 + fC_t0 + (1 | plot),data=tmp,family="binomial")
lMod[[12]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + mC_t0 + fC_t0 + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = sex + response=sex 
lMod[[13]]=glmmadmb(flow_t1 ~ log_l_t1 + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=tmp,family="binomial")
lMod[[14]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=tmp,family="binomial")
lMod[[15]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=tmp,family="binomial")
AICtab(lMod,weights=T)



#Density as "number of individuals"-----------------------------------------------------------
#RESULTS ARE ALMOST THE SAME, best two models make up 60% of weight.
dMod=list() #lMod stands for "leaf model" (density is quantified by N. of leaves)
tmp=na.omit(dat[,c("plot","flow_t1","log_l_t1","sex","TotDensity","M","F","M","F")])
#Target fitness
dMod[[1]]=glmmadmb(flow_t1 ~ log_l_t1 + (1 | plot),data=tmp,family="binomial")
dMod[[2]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + (1 | plot),data=tmp,family="binomial")
dMod[[3]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = tot density 
dMod[[4]]=glmmadmb(flow_t1 ~ log_l_t1 + TotDensity + (1 | plot),data=tmp,family="binomial")
dMod[[5]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + TotDensity + (1 | plot),data=tmp,family="binomial")
dMod[[6]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + TotDensity + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = tot density + response=sex 
dMod[[7]]=glmmadmb(flow_t1 ~ log_l_t1 + TotDensity + TotDensity:sex + (1 | plot),data=tmp,family="binomial")
dMod[[8]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + TotDensity + TotDensity:sex + (1 | plot),data=tmp,family="binomial")
dMod[[9]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + TotDensity + TotDensity:sex + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = sex 
dMod[[10]]=glmmadmb(flow_t1 ~ log_l_t1 + M + F + (1 | plot),data=tmp,family="binomial")
dMod[[11]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + M + F + (1 | plot),data=tmp,family="binomial")
dMod[[12]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + M + F + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = sex + response=sex 
dMod[[13]]=glmmadmb(flow_t1 ~ log_l_t1 + M + F + M:sex + F:sex + (1 | plot),data=tmp,family="binomial")
dMod[[14]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + M + F + M:sex + F:sex + (1 | plot),data=tmp,family="binomial")
dMod[[15]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + M + F + M:sex + F:sex + (1 | plot),data=tmp,family="binomial")
AICtab(dMod,weights=T)






#
tiff("Results/VitalRates_simple/flowering_simple.tiff",unit="in",width=4,height=8,res=600,compression="lzw")

par(mfrow=c(2,1),mar=c(3,3,1,0.1),mgp=c(1.4,0.5,0))
invlogit<-function(x){exp(x)/(1+exp(x))}

#Best model
plot(tmp$flow_t1 ~ tmp$mC_t0,pch=16,xlab="Competitor leaves 2013",
     ylab="Flowering probability, spring 2014",col="red")
par(new=T) ; plot(tmp$flow_t1 ~ tmp$fC_t0,pch=16,xlab="",ylab="",col="blue",xaxt="n")
title(main = "2015: best mod (64% weight)", line=0.2,cex=0.9)
xSeq <- seq(0,max(c(tmp$fC_t0,tmp$mC_t0)),by=1)
yFem_f=invlogit(coef(mod[[6]])[1] + coef(mod[[6]])[2]*mean(tmp$log_l_t1,na.rm=T) +
                  coef(mod[[6]])[4]*mean(tmp$mC_t0,na.rm=T) + coef(mod[[6]])[5]*xSeq)
yFem_m=invlogit(coef(mod[[6]])[1] + coef(mod[[6]])[2]*mean(tmp$log_l_t1,na.rm=T) +
                  coef(mod[[6]])[4]*xSeq + coef(mod[[6]])[5]*mean(tmp$fC_t0,na.rm=T))
yMal_f=invlogit(coef(mod[[6]])[1] + coef(mod[[6]])[2]*mean(tmp$log_l_t1,na.rm=T) + coef(mod[[6]])[3] + 
                  coef(mod[[6]])[4]*mean(tmp$mC_t0,na.rm=T) + coef(mod[[6]])[5]*xSeq + 
                  coef(mod[[6]])[6]*mean(tmp$log_l_t1,na.rm=T))
yMal_m=invlogit(coef(mod[[6]])[1] + coef(mod[[6]])[2]*mean(tmp$log_l_t1,na.rm=T) + coef(mod[[6]])[3] + 
                  coef(mod[[6]])[4]*xSeq + coef(mod[[6]])[5]*mean(tmp$fC_t0,na.rm=T) + 
                  coef(mod[[6]])[6]*mean(tmp$log_l_t1,na.rm=T))
lines(xSeq,yMal_m,lwd=2,col="red")
lines(xSeq,yMal_f,lwd=2,lty=2,col="red")
lines(xSeq,yFem_m,lwd=2,col="blue")
lines(xSeq,yFem_f,lwd=2,lty=2,col="blue")
legend(-15,0.5,c("Focal male, male competition",
                 "Focal female, male competition",
                 "Focal male, female competition",
                 "Focal female, female competition"),
       col=c("red","blue"),lty=c(1,1,2,2),lwd=2,bty="n")

#2nd best model
xSeq <- seq(0,max(c(tmp$fC_t0,tmp$mC_t0)),by=1)
yFem_f=invlogit(coef(mod[[9]])[1] + coef(mod[[9]])[2]*mean(tmp$log_l_t1,na.rm=T) +
                  coef(mod[[9]])[4]*mean(tmp$mC_t0,na.rm=T) + coef(mod[[9]])[5]*xSeq)
yFem_m=invlogit(coef(mod[[9]])[1] + coef(mod[[9]])[2]*mean(tmp$log_l_t1,na.rm=T) +
                  coef(mod[[9]])[4]*xSeq + coef(mod[[9]])[5]*mean(tmp$fC_t0,na.rm=T))

yMal_f=invlogit(coef(mod[[9]])[1] + coef(mod[[9]])[2]*mean(tmp$log_l_t1,na.rm=T) + coef(mod[[9]])[3] + 
                  coef(mod[[9]])[4]*mean(tmp$mC_t0,na.rm=T) + coef(mod[[9]])[5]*xSeq + 
                  coef(mod[[9]])[6]*mean(tmp$log_l_t1,na.rm=T))
yMal_m=invlogit(coef(mod[[9]])[1] + coef(mod[[9]])[2]*mean(tmp$log_l_t1,na.rm=T) + coef(mod[[9]])[3] + 
                  coef(mod[[9]])[4]*xSeq + coef(mod[[9]])[5]*mean(tmp$fC_t0,na.rm=T) + coef(mod[[9]])[7]*xSeq + 
                  coef(mod[[9]])[6]*mean(tmp$log_l_t1,na.rm=T))
plot(tmp$flow_t1 ~ tmp$mC_t0,pch=16,xlab="Competitor leaves 2013",
     ylab="Flowering probability, spring 2014",col="red")
par(new=T) ; plot(tmp$flow_t1 ~ tmp$fC_t0,pch=16,xlab="",ylab="",col="blue",xaxt="n")
title(main = "2015: best mod (16% weight)", line=0.2,cex=0.9)
lines(xSeq,yMal_m,lwd=2,col="red")
lines(xSeq,yMal_f,lwd=2,lty=2,col="red")
lines(xSeq,yFem_m,lwd=2,col="blue")
lines(xSeq,yFem_f,lwd=2,lty=2,col="blue")

dev.off()

