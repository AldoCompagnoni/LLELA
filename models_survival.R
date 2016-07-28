##1. Selection of models of POAR survival using data from LLELA's competition experiment.
##2. Year specific data anlayzed SEPARATELY.
##3. BEWARE:  the fields "TotDensity","F", and "M" refer to FALL 2013, even when "year==2015".
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
dMod=list() #dMod stands for "leaf model" (density is quantified by N. of leaves)

#prepare the 'new_t1' column
dat$new_t1[dat$new_t1=="SKIPPED"]=NA
dat$new_t1[dat$new_t1=="cnf"]=NA
dat$new_t1[dat$new_t1==""]=NA
dat$new_t1=as.numeric(as.character(dat$new_t1))
a15=na.omit(dat[,c("plot","surv_t1","log_l_t0","sex","TotDensity","M","F","M","F","new_t1")])

############################################################################################
#Preliminary model selection: what's best?
############################################################################################

#1. plot random effect only
#2. plot + new tillers
#3. NOT REPORTED, worse in model selection: new tillers + (1 + new tillers | plot) 

#1. plot random effect only--------------------------------------------------------------------------
#Target fitness
dMod[[1]]=glmmadmb(surv_t1 ~ log_l_t0 + (1 | plot),data=a15,family="binomial")
dMod[[2]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + (1 | plot),data=a15,family="binomial")
dMod[[3]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + (1 | plot),data=a15,family="binomial")
#Target fitness + effect = tot density 
dMod[[4]]=glmmadmb(surv_t1 ~ log_l_t0 + TotDensity + (1 | plot),data=a15,family="binomial")
dMod[[5]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + TotDensity + (1 | plot),data=a15,family="binomial")
dMod[[6]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + TotDensity + (1 | plot),data=a15,family="binomial")
#Target fitness + effect = tot density + response=sex 
dMod[[7]]=glmmadmb(surv_t1 ~ log_l_t0 + TotDensity + TotDensity:sex + (1 | plot),data=a15,family="binomial")
dMod[[8]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + TotDensity + TotDensity:sex + (1 | plot),data=a15,family="binomial")
dMod[[9]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + TotDensity + TotDensity:sex + (1 | plot),data=a15,family="binomial")
#Target fitness + effect = sex 
dMod[[10]]=glmmadmb(surv_t1 ~ log_l_t0 + M + F + (1 | plot),data=a15,family="binomial")
dMod[[11]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + M + F + (1 | plot),data=a15,family="binomial")
dMod[[12]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + M + F + (1 | plot),data=a15,family="binomial")
#Target fitness + effect = sex + response=sex 
dMod[[13]]=glmmadmb(surv_t1 ~ log_l_t0 + M + F + M:sex + F:sex + (1 | plot),data=a15,family="binomial")
dMod[[14]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + M + F + M:sex + F:sex + (1 | plot),data=a15,family="binomial")
dMod[[15]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + M + F + M:sex + F:sex + (1 | plot),data=a15,family="binomial")

#2. plot + new tillers--------------------------------------------------------------------------
#Use "new tillers" instead of 'plot random effect'
dMod[[16]]=glmmadmb(surv_t1 ~ log_l_t0 + new_t1 + (1 | plot),data=a15,family="binomial")
dMod[[17]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + new_t1 + (1 | plot),data=a15,family="binomial")
dMod[[18]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + new_t1 + (1 | plot),data=a15,family="binomial")
#Target fitness + effect = tot density 
dMod[[19]]=glmmadmb(surv_t1 ~ log_l_t0 + TotDensity + new_t1 + (1 | plot),data=a15,family="binomial")
dMod[[20]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + TotDensity + new_t1 + (1 | plot),data=a15,family="binomial")
dMod[[21]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + TotDensity + new_t1 + (1 | plot),data=a15,family="binomial")
#Target fitness + effect = tot density + response=sex 
dMod[[22]]=glmmadmb(surv_t1 ~ log_l_t0 + TotDensity + TotDensity:sex + new_t1 + (1 | plot),data=a15,family="binomial")
dMod[[23]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + TotDensity + TotDensity:sex + new_t1 + (1 | plot),data=a15,family="binomial")
dMod[[24]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + TotDensity + TotDensity:sex + new_t1 + (1 | plot),data=a15,family="binomial")
#Target fitness + effect = sex 
dMod[[25]]=glmmadmb(surv_t1 ~ log_l_t0 + M + F + new_t1 + (1 | plot),data=a15,family="binomial")
dMod[[26]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + M + F + new_t1 + (1 | plot),data=a15,family="binomial")
dMod[[27]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + M + F + new_t1 + (1 | plot),data=a15,family="binomial")
#Target fitness + effect = sex + response=sex 
dMod[[28]]=glmmadmb(surv_t1 ~ log_l_t0 + M + F + M:sex + F:sex + new_t1 + (1 | plot),data=a15,family="binomial")
dMod[[29]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + M + F + M:sex + F:sex + new_t1 + (1 | plot),data=a15,family="binomial")
dMod[[30]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + M + F + M:sex + F:sex + new_t1 + (1 | plot),data=a15,family="binomial")
AICtab(dMod,weights=T)


#########################################################################################################
###GRAPHS######
#########################################################################################################

tiff("Results/VitalRates_simple/survival_tillers.tiff",unit="in",width=4,height=8,res=600,compression="lzw")

sexAsInteger=as.integer(a15$sex)
a15$col=as.character(factor(sexAsInteger,labels=c("blue","red")))
a15$symb=as.integer(as.character(factor(sexAsInteger,labels=c("17","16"))))

#set up graph parameters
par(mfrow=c(2,1),mar=c(3,3,1,0.1),mgp=c(1.4,0.5,0))
invlogit<-function(x){exp(x)/(1+exp(x))}

#best models
plot(a15$surv_t1+0.015 ~ a15$M,pch=16,xlab="N. of total plot leaves 2013",ylim=c(0,1),
     ylab="Survival to spring 2015")
par(new=T) ; plot(a15$surv_t1-0.015 ~ a15$F,pch=17,xlab="",ylab="",ylim=c(0,1),col="grey50", xaxt="n")

title(main = "2015: best mod (45% weight)", line=0.2,cex=0.9)
legend(0,0.45,c("Male density","Female density"), col="black",lty=c(1,2),lwd=2,bty="n")
#text(-15,0.1,"Predictor: Density by sex",pos=4,cex=1.2)
xSeq <- seq(0,max(c(a15$F,a15$M)),by=1)
y_m=invlogit(coef(dMod[[25]])[1] + coef(dMod[[25]])[2]*mean(a15$log_l_t0,na.rm=T) +
             coef(dMod[[25]])[3]*xSeq + coef(dMod[[25]])[4]*mean(a15$F,na.rm=T) +
             coef(dMod[[25]])[5]*mean(a15$new_t1,na.rm=T))
y_f=invlogit(coef(dMod[[25]])[1] + coef(dMod[[25]])[2]*mean(a15$log_l_t0,na.rm=T) +
             coef(dMod[[25]])[3]*mean(a15$M,na.rm=T) + coef(dMod[[25]])[4]*xSeq + 
             coef(dMod[[25]])[5]*mean(a15$new_t1,na.rm=T))
lines(xSeq,y_m,lty=1,lwd=2,col="black")
lines(xSeq,y_f,lty=2,lwd=2,col="black")


#2nd best models
plot(a15$surv_t1+0.015 ~ a15$M,pch=16,xlab="N. of total plot leaves 2013",ylim=c(0,1),
     ylab="Survival to spring 2015",col=a15$col)
par(new=T) ; plot(a15$surv_t1-0.015 ~ a15$F,pch=17,ylim=c(0,1),xlab="",ylab="",col=a15$col, xaxt="n")
title(main = "2015: 2nd best (17% weight)", line=0.2,cex=0.9)
legend(0,0.5,c("Focal male, male competition",
                 "Focal female, male competition",
                 "Focal male, female competition",
                 "Focal female, female competition"),
       col=c("red","blue"),lty=c(1,1,2,2),lwd=2,bty="n")
xSeq <- seq(0,max(c(a15$F,a15$M)),by=1)
yMal_f=invlogit(coef(dMod[[26]])[1] + coef(dMod[[26]])[2]*mean(a15$log_l_t0,na.rm=T) + coef(dMod[[26]])[3] + 
                coef(dMod[[26]])[4]*mean(a15$M,na.rm=T) + coef(dMod[[26]])[5]*xSeq + 
                coef(dMod[[26]])[6]*mean(a15$new_t1,na.rm=T))
yMal_m=invlogit(coef(dMod[[26]])[1] + coef(dMod[[26]])[2]*mean(a15$log_l_t0,na.rm=T) + coef(dMod[[26]])[3] + 
                coef(dMod[[26]])[4]*xSeq + coef(dMod[[26]])[5]*mean(a15$F,na.rm=T) + 
                coef(dMod[[26]])[6]*mean(a15$new_t1,na.rm=T))
yFem_f=invlogit(coef(dMod[[26]])[1] + coef(dMod[[26]])[2]*mean(a15$log_l_t0,na.rm=T) +
                coef(dMod[[26]])[4]*mean(a15$M,na.rm=T) + coef(dMod[[26]])[5]*xSeq +
                coef(dMod[[26]])[6]*mean(a15$new_t1,na.rm=T))
yFem_m=invlogit(coef(dMod[[26]])[1] + coef(dMod[[26]])[2]*mean(a15$log_l_t0,na.rm=T) + 
                coef(dMod[[26]])[4]*xSeq + coef(dMod[[26]])[5]*mean(a15$F,na.rm=T) +
                coef(dMod[[26]])[6]*mean(a15$new_t1,na.rm=T))
lines(xSeq,yMal_m,lwd=2,col="red")
lines(xSeq,yMal_f,lwd=2,lty=2,col="red")
lines(xSeq,yFem_m,lwd=2,col="blue")
lines(xSeq,yFem_f,lwd=2,lty=2,col="blue")

dev.off()

