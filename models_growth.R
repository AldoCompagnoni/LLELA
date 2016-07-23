##1. Selection of models of POAR growth using data from LLELA's competition experiment.  
##2. Year specific data anlayzed SEPARATELY.
##3. BEWARE:  the fields "c_t0","fC_t0", and "mC_t0" refer to FALL 2013, even when "year==2015".
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) #For AIC weights, because I'm lazy!
library(lme4)
library(nlme)

#read in data
d=read.csv("Data/vr.csv")

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
d$plot=as.factor(d$plot) #plots need to be factors
#only use year two
dat=subset(d,year==2014)
dat$logRat=log(dat$l_t1/dat$l_t0)
library(glmmADMB) #Fit models with a Negative Binomial

#########################################################################################################
###RAW RESPONSE######
#########################################################################################################

#Density as "number of leaves"-----------------------------------------------------------
lMod=list() #lMod stands for "leaf model" (density is quantified by N. of leaves)
#Target fitness
lMod[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=dat,family="nbinom2")
lMod[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=dat,family="nbinom2")
lMod[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=dat,family="nbinom2")
#Target fitness + effect = tot density 
lMod[[4]]=glmmadmb(l_t1 ~ log_l_t0 + c_t0 + (1 | plot),data=dat,family="nbinom2")
lMod[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + c_t0 + (1 | plot),data=dat,family="nbinom2")
lMod[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + c_t0 + (1 | plot),data=dat,family="nbinom2")
#Target fitness + effect = tot density + response=sex 
lMod[[7]]=glmmadmb(l_t1 ~ log_l_t0 + c_t0 + c_t0:sex + (1 | plot),data=dat,family="nbinom2")
lMod[[8]]=glmmadmb(l_t1 ~ log_l_t0 + sex + c_t0 + c_t0:sex + (1 | plot),data=dat,family="nbinom2")
lMod[[9]]=glmmadmb(l_t1 ~ log_l_t0 * sex + c_t0 + c_t0:sex + (1 | plot),data=dat,family="nbinom2")
#Target fitness + effect = sex 
lMod[[10]]=glmmadmb(l_t1 ~ log_l_t0 + mC_t0 + fC_t0 + (1 | plot),data=dat,family="nbinom2")
lMod[[11]]=glmmadmb(l_t1 ~ log_l_t0 + sex + mC_t0 + fC_t0 + (1 | plot),data=dat,family="nbinom2")
lMod[[12]]=glmmadmb(l_t1 ~ log_l_t0 * sex + mC_t0 + fC_t0 + (1 | plot),data=dat,family="nbinom2")
#Target fitness + effect = sex + response=sex 
lMod[[13]]=glmmadmb(l_t1 ~ log_l_t0 + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=dat,family="nbinom2")
lMod[[14]]=glmmadmb(l_t1 ~ log_l_t0 + sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=dat,family="nbinom2")
lMod[[15]]=glmmadmb(l_t1 ~ log_l_t0 * sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=dat,family="nbinom2")
AICtab(lMod,weights=T)
sum(AICtab(lMod,weights=T)$weight[1:3]) #first 3 models make only 66%


#Density as "number of leaves"-----------------------------------------------------------
dMod=list() #dMod stands for "density model" (density is # m/f planted)
#Target fitness
dMod[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=dat,family="nbinom2")
dMod[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=dat,family="nbinom2")
dMod[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=dat,family="nbinom2")
#Target fitness + effect = tot density 
dMod[[4]]=glmmadmb(l_t1 ~ log_l_t0 + TotDensity + (1 | plot),data=dat,family="nbinom2")
dMod[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + TotDensity + (1 | plot),data=dat,family="nbinom2")
dMod[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + TotDensity + (1 | plot),data=dat,family="nbinom2")
#Target fitness + effect = tot density + response=sex 
dMod[[7]]=glmmadmb(l_t1 ~ log_l_t0 + TotDensity + TotDensity:sex + (1 | plot),data=dat,family="nbinom2")
dMod[[8]]=glmmadmb(l_t1 ~ log_l_t0 + sex + TotDensity + TotDensity:sex + (1 | plot),data=dat,family="nbinom2")
dMod[[9]]=glmmadmb(l_t1 ~ log_l_t0 * sex + TotDensity + TotDensity:sex + (1 | plot),data=dat,family="nbinom2")
#Target fitness + effect = sex 
dMod[[10]]=glmmadmb(l_t1 ~ log_l_t0 + M + F + (1 | plot),data=dat,family="nbinom2")
dMod[[11]]=glmmadmb(l_t1 ~ log_l_t0 + sex + M + F + (1 | plot),data=dat,family="nbinom2")
dMod[[12]]=glmmadmb(l_t1 ~ log_l_t0 * sex + M + F + (1 | plot),data=dat,family="nbinom2")
#Target fitness + effect = sex + response=sex 
dMod[[13]]=glmmadmb(l_t1 ~ log_l_t0 + M + F + M:sex + F:sex + (1 | plot),data=dat,family="nbinom2")
dMod[[14]]=glmmadmb(l_t1 ~ log_l_t0 + sex + M + F + M:sex + F:sex + (1 | plot),data=dat,family="nbinom2")
dMod[[15]]=glmmadmb(l_t1 ~ log_l_t0 * sex + M + F + M:sex + F:sex + (1 | plot),data=dat,family="nbinom2")
AICtab(dMod,weights=T)



#########################################################################################################
###LOG RATIO RESPONSE######
#########################################################################################################


llMod=list() #llMod stands for "leaf model", "log ratio"
#Target fitness
llMod[[1]]=lmer(logRat ~ (1 | plot),data=dat)
llMod[[2]]=lmer(logRat ~ sex + (1 | plot),data=dat)
llMod[[3]]=lmer(logRat ~ sex + (1 | plot),data=dat)
#Target fitness + effect = tot density 
llMod[[4]]=lmer(logRat ~ c_t0 + (1 | plot),data=dat)
llMod[[5]]=lmer(logRat ~ sex + c_t0 + (1 | plot),data=dat)
llMod[[6]]=lmer(logRat ~ sex + c_t0 + (1 | plot),data=dat)
#Target fitness + effect = tot density + response=sex 
llMod[[7]]=lmer(logRat ~ c_t0 + c_t0:sex + (1 | plot),data=dat)
llMod[[8]]=lmer(logRat ~ sex + c_t0 + c_t0:sex + (1 | plot),data=dat)
llMod[[9]]=lmer(logRat ~ sex + c_t0 + c_t0:sex + (1 | plot),data=dat)
#Target fitness + effect = sex 
llMod[[10]]=lmer(logRat ~ mC_t0 + fC_t0 + (1 | plot),data=dat)
llMod[[11]]=lmer(logRat ~ sex + mC_t0 + fC_t0 + (1 | plot),data=dat)
llMod[[12]]=lmer(logRat ~ sex + mC_t0 + fC_t0 + (1 | plot),data=dat)
#Target fitness + effect = sex + response=sex 
llMod[[13]]=lmer(logRat ~ mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=dat)
llMod[[14]]=lmer(logRat ~ sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=dat)
llMod[[15]]=lmer(logRat ~ sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=dat)
AICtab(llMod,weights=T)




dlMod=list() #dlMod stands for "leaf model", "log ratio"
#Target fitness
dlMod[[1]]=lmer(logRat ~ (1 | plot),data=dat)
dlMod[[2]]=lmer(logRat ~ sex + (1 | plot),data=dat)
dlMod[[3]]=lmer(logRat ~ sex + (1 | plot),data=dat)
#Target fitness + effect = tot density 
dlMod[[4]]=lmer(logRat ~ TotDensity + (1 | plot),data=dat)
dlMod[[5]]=lmer(logRat ~ sex + TotDensity + (1 | plot),data=dat)
dlMod[[6]]=lmer(logRat ~ sex + TotDensity + (1 | plot),data=dat)
#Target fitness + effect = tot density + response=sex 
dlMod[[7]]=lmer(logRat ~ TotDensity + TotDensity:sex + (1 | plot),data=dat)
dlMod[[8]]=lmer(logRat ~ sex + TotDensity + TotDensity:sex + (1 | plot),data=dat)
dlMod[[9]]=lmer(logRat ~ sex + TotDensity + TotDensity:sex + (1 | plot),data=dat)
#Target fitness + effect = sex 
dlMod[[10]]=lmer(logRat ~ M + F + (1 | plot),data=dat)
dlMod[[11]]=lmer(logRat ~ sex + M + F + (1 | plot),data=dat)
dlMod[[12]]=lmer(logRat ~ sex + M + F + (1 | plot),data=dat)
#Target fitness + effect = sex + response=sex 
dlMod[[13]]=lmer(logRat ~ M + F + M:sex + F:sex + (1 | plot),data=dat)
dlMod[[14]]=lmer(logRat ~ sex + M + F + M:sex + F:sex + (1 | plot),data=dat)
dlMod[[15]]=lmer(logRat ~ sex + M + F + M:sex + F:sex + (1 | plot),data=dat)
AICtab(dlMod,weights=T)









#########################################################################################################
###DIFFERENCE IN N. of LEAVES######
#########################################################################################################



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

