##1. Selection of models of POAR survival using data from LLELA's competition experiment.
##2. Year specific data anlayzed SEPARATELY.
##3. BEWARE:  the fields "c_t0","fC_t0", and "mC_t0" refer to FALL 2013, even when "year==2015".
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) #For AIC weights, because I'm lazy!
library(lme4)
library(glmmADMB)

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


#YEAR TWO##########################################################################################################
d$plot=as.factor(d$plot)
dat=subset(d,year==2015)
sum(dat$surv_t1==0,na.rm=T)


#Density as "number of leaves"-----------------------------------------------------------
lMod=list() #lMod stands for "leaf model" (density is quantified by N. of leaves)

#prepare the 'new_t1' column
dat$new_t1[dat$new_t1=="SKIPPED"]=NA
dat$new_t1[dat$new_t1=="cnf"]=NA
dat$new_t1[dat$new_t1==""]=NA
dat$new_t1=as.numeric(as.character(dat$new_t1))
tmp=na.omit(dat[,c("plot","surv_t1","log_l_t0","sex","c_t0","mC_t0","fC_t0","M","F","new_t1")])

############################################################################################
#Preliminary model selection: what's best?
############################################################################################

#1. plot random effect only
#2. plot + new tillers
#3. new tillers, but no plot

#1. plot random effect only--------------------------------------------------------------------------
#Target fitness
lMod[[1]]=glmmadmb(surv_t1 ~ log_l_t0 + (1 | plot),data=tmp,family="binomial")
lMod[[2]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + (1 | plot),data=tmp,family="binomial")
lMod[[3]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = tot density 
lMod[[4]]=glmmadmb(surv_t1 ~ log_l_t0 + c_t0 + (1 | plot),data=tmp,family="binomial")
lMod[[5]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + c_t0 + (1 | plot),data=tmp,family="binomial")
lMod[[6]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + c_t0 + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = tot density + response=sex 
lMod[[7]]=glmmadmb(surv_t1 ~ log_l_t0 + c_t0 + c_t0:sex + (1 | plot),data=tmp,family="binomial")
lMod[[8]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + c_t0 + c_t0:sex + (1 | plot),data=tmp,family="binomial")
lMod[[9]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + c_t0 + c_t0:sex + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = sex 
lMod[[10]]=glmmadmb(surv_t1 ~ log_l_t0 + mC_t0 + fC_t0 + (1 | plot),data=tmp,family="binomial")
lMod[[11]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + mC_t0 + fC_t0 + (1 | plot),data=tmp,family="binomial")
lMod[[12]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + mC_t0 + fC_t0 + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = sex + response=sex 
lMod[[13]]=glmmadmb(surv_t1 ~ log_l_t0 + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=tmp,family="binomial")
lMod[[14]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=tmp,family="binomial")
lMod[[15]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=tmp,family="binomial")

#2. plot + new tillers--------------------------------------------------------------------------
#Use "new tillers" instead of 'plot random effect'
lMod[[16]]=glmmadmb(surv_t1 ~ log_l_t0 + new_t1 + (1 | plot),data=tmp,family="binomial")
lMod[[17]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + new_t1 + (1 | plot),data=tmp,family="binomial")
lMod[[18]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + new_t1 + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = tot density 
lMod[[19]]=glmmadmb(surv_t1 ~ log_l_t0 + c_t0 + new_t1 + (1 | plot),data=tmp,family="binomial")
lMod[[20]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + c_t0 + new_t1 + (1 | plot),data=tmp,family="binomial")
lMod[[21]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + c_t0 + new_t1 + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = tot density + response=sex 
lMod[[22]]=glmmadmb(surv_t1 ~ log_l_t0 + c_t0 + c_t0:sex + new_t1 + (1 | plot),data=tmp,family="binomial")
lMod[[23]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + c_t0 + c_t0:sex + new_t1 + (1 | plot),data=tmp,family="binomial")
lMod[[24]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + c_t0 + c_t0:sex + new_t1 + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = sex 
lMod[[25]]=glmmadmb(surv_t1 ~ log_l_t0 + mC_t0 + fC_t0 + new_t1 + (1 | plot),data=tmp,family="binomial")
lMod[[26]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + mC_t0 + fC_t0 + new_t1 + (1 | plot),data=tmp,family="binomial")
lMod[[27]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + mC_t0 + fC_t0 + new_t1 + (1 | plot),data=tmp,family="binomial")
#Target fitness + effect = sex + response=sex 
lMod[[28]]=glmmadmb(surv_t1 ~ log_l_t0 + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + new_t1 + (1 | plot),data=tmp,family="binomial")
lMod[[29]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + new_t1 + (1 | plot),data=tmp,family="binomial")
lMod[[30]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + new_t1 + (1 | plot),data=tmp,family="binomial")

#3. new tillers, but no plot--------------------------------------------------------------------------
#Use "new tillers" instead of 'plot random effect'
lMod[[31]]=glmmadmb(surv_t1 ~ log_l_t0 + new_t1,data=tmp,family="binomial")
lMod[[32]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + new_t1,data=tmp,family="binomial")
lMod[[33]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + new_t1,data=tmp,family="binomial")
#Target fitness + effect = tot density 
lMod[[34]]=glmmadmb(surv_t1 ~ log_l_t0 + c_t0 + new_t1,data=tmp,family="binomial")
lMod[[35]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + c_t0 + new_t1,data=tmp,family="binomial")
lMod[[36]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + c_t0 + new_t1,data=tmp,family="binomial")
#Target fitness + effect = tot density + response=sex 
lMod[[37]]=glmmadmb(surv_t1 ~ log_l_t0 + c_t0 + c_t0:sex + new_t1,data=tmp,family="binomial")
lMod[[38]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + c_t0 + c_t0:sex + new_t1,data=tmp,family="binomial")
lMod[[39]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + c_t0 + c_t0:sex + new_t1,data=tmp,family="binomial")
#Target fitness + effect = sex 
lMod[[40]]=glmmadmb(surv_t1 ~ log_l_t0 + mC_t0 + fC_t0 + new_t1,data=tmp,family="binomial")
lMod[[41]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + mC_t0 + fC_t0 + new_t1,data=tmp,family="binomial")
lMod[[42]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + mC_t0 + fC_t0 + new_t1,data=tmp,family="binomial")
#Target fitness + effect = sex + response=sex 
lMod[[43]]=glmmadmb(surv_t1 ~ log_l_t0 + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + new_t1,data=tmp,family="binomial")
lMod[[44]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + new_t1,data=tmp,family="binomial")
lMod[[45]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + new_t1,data=tmp,family="binomial")
AICtab(lMod,weights=T)

#plot + new_tillers is clearly best

#coding "experimet". Erase this as soon as upi are done with it!
plot(tmp$c_t0,sqrt(tmp$new_t1),pch=16)

prov=list()
prov[[1]]=glmmadmb(surv_t1 ~ log_l_t0 + new_t1 + (1 | plot),data=tmp,family="binomial")
prov[[2]]=glmmadmb(surv_t1 ~ log_l_t0 + c_t0 + (1 | plot),data=tmp,family="binomial")
prov[[3]]=glmmadmb(surv_t1 ~ log_l_t0 + c_t0 * new_t1 + (1 | plot),data=tmp,family="binomial")
prov[[4]]=glmmadmb(surv_t1 ~ log_l_t0*new_t1 + c_t0 + (1 | plot),data=tmp,family="binomial")
prov[[5]]=glmmadmb(surv_t1 ~ log_l_t0 + c_t0 + log_l_t0:new_t1 + c_t0:new_t1 + new_t1 + (1 | plot),data=tmp,family="binomial")
AICtab(prov,weights=T)





#
tiff("Results/VitalRates_simple/survival_simple.tiff",unit="in",width=4,height=8,res=600,compression="lzw")

par(mfrow=c(2,1),mar=c(3,3,1,0.1),mgp=c(1.4,0.5,0))
invlogit<-function(x){exp(x)/(1+exp(x))}

#best models
plot(dat$surv_t1 ~ dat$mC_t0,pch=16,xlab="N. of total plot leaves 2013",
     ylab="Survival to spring 2015",col="red")
par(new=T) ; plot(dat$surv_t1 ~ dat$fC_t0,pch=16,xlab="",ylab="",col="blue", xaxt="n")
title(main = "2015: best mod (50% weight)", line=0.2,cex=0.9)
legend(-15,0.45,c("Male density","Female density"), col="black",lty=c(1,2),lwd=2,bty="n")
#text(-15,0.1,"Predictor: Density by sex",pos=4,cex=1.2)
xSeq <- seq(0,max(c(dat$fC_t0,dat$mC_t0)),by=1)
y_m=invlogit(coef(mod[[4]])[1] + coef(mod[[4]])[2]*mean(dat$log_l_t0,na.rm=T) +
              coef(mod[[4]])[3]*xSeq + coef(mod[[4]])[4]*mean(dat$fC_t0,na.rm=T))
y_f=invlogit(coef(mod[[4]])[1] + coef(mod[[4]])[2]*mean(dat$log_l_t0,na.rm=T) +
              coef(mod[[4]])[3]*mean(dat$mC_t0,na.rm=T) + coef(mod[[4]])[4]*xSeq)
lines(xSeq,y_m,lty=1,lwd=2,col="black")
lines(xSeq,y_f,lty=2,lwd=2,col="black")


#2nd best models
plot(dat$surv_t1 ~ dat$mC_t0,pch=16,xlab="N. of total plot leaves 2013",
     ylab="Survival to spring 2015",col="red")
par(new=T) ; plot(dat$surv_t1 ~ dat$fC_t0,pch=16,xlab="",ylab="",col="blue", xaxt="n")
title(main = "2015: 2nd best (20% weight)", line=0.2,cex=0.9)
legend(-15,0.5,c("Focal male, male competition",
                 "Focal female, male competition",
                 "Focal male, female competition",
                 "Focal female, female competition"),
       col=c("red","blue"),lty=c(1,1,2,2),lwd=2,bty="n")
xSeq <- seq(0,max(c(dat$fC_t0,dat$mC_t0)),by=1)
yMal_f=invlogit(coef(mod[[5]])[1] + coef(mod[[5]])[2]*mean(dat$log_l_t0,na.rm=T) + coef(mod[[5]])[3] + 
                coef(mod[[5]])[4]*mean(dat$mC_t0,na.rm=T) + coef(mod[[5]])[5]*xSeq)
yMal_m=invlogit(coef(mod[[5]])[1] + coef(mod[[5]])[2]*mean(dat$log_l_t0,na.rm=T) + coef(mod[[5]])[3] + 
                coef(mod[[5]])[4]*xSeq + coef(mod[[5]])[5]*mean(dat$fC_t0,na.rm=T))
yFem_f=invlogit(coef(mod[[5]])[1] + coef(mod[[5]])[2]*mean(dat$log_l_t0,na.rm=T) +
                  coef(mod[[5]])[4]*mean(dat$mC_t0,na.rm=T) + coef(mod[[5]])[5]*xSeq)
yFem_m=invlogit(coef(mod[[5]])[1] + coef(mod[[5]])[2]*mean(dat$log_l_t0,na.rm=T) + 
                  coef(mod[[5]])[4]*xSeq + coef(mod[[5]])[5]*mean(dat$fC_t0,na.rm=T))
lines(xSeq,yMal_m,lwd=2,col="red")
lines(xSeq,yMal_f,lwd=2,lty=2,col="red")
lines(xSeq,yFem_m,lwd=2,col="blue")
lines(xSeq,yFem_f,lwd=2,lty=2,col="blue")

dev.off()

