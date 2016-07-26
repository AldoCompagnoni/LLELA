##1. Selection of models of POAR growth using data from LLELA's competition experiment.  
##2. Year specific data anlayzed SEPARATELY.
##3. BEWARE:  the fields "c_t0","fC_t0", and "mC_t0" refer to FALL 2013, even when "year==2015".
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) #For AIC weights, because I'm lazy!
library(glmmADMB) #Fit models with a Negative Binomial

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

d$plot=as.factor(d$plot) #glmmadmb wants plot as a factor


#########################################################################################################
###2014######
#########################################################################################################

d14=subset(d,year==2014)
d14$logRat=log(dat$l_t1/dat$l_t0)

#Density as "number of leaves"-----------------------------------------------------------
lMod=list() #lMod stands for "leaf model" (density is quantified by N. of leaves)
#Target fitness
lMod[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=d14,family="nbinom2")
lMod[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=d14,family="nbinom2")
lMod[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=d14,family="nbinom2")
#Target fitness + effect = tot density 
lMod[[4]]=glmmadmb(l_t1 ~ log_l_t0 + c_t0 + (1 | plot),data=d14,family="nbinom2")
lMod[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + c_t0 + (1 | plot),data=d14,family="nbinom2")
lMod[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + c_t0 + (1 | plot),data=d14,family="nbinom2")
#Target fitness + effect = tot density + response=sex 
lMod[[7]]=glmmadmb(l_t1 ~ log_l_t0 + c_t0 + c_t0:sex + (1 | plot),data=d14,family="nbinom2")
lMod[[8]]=glmmadmb(l_t1 ~ log_l_t0 + sex + c_t0 + c_t0:sex + (1 | plot),data=d14,family="nbinom2")
lMod[[9]]=glmmadmb(l_t1 ~ log_l_t0 * sex + c_t0 + c_t0:sex + (1 | plot),data=d14,family="nbinom2")
#Target fitness + effect = sex 
lMod[[10]]=glmmadmb(l_t1 ~ log_l_t0 + mC_t0 + fC_t0 + (1 | plot),data=d14,family="nbinom2")
lMod[[11]]=glmmadmb(l_t1 ~ log_l_t0 + sex + mC_t0 + fC_t0 + (1 | plot),data=d14,family="nbinom2")
lMod[[12]]=glmmadmb(l_t1 ~ log_l_t0 * sex + mC_t0 + fC_t0 + (1 | plot),data=d14,family="nbinom2")
#Target fitness + effect = sex + response=sex 
lMod[[13]]=glmmadmb(l_t1 ~ log_l_t0 + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=d14,family="nbinom2")
lMod[[14]]=glmmadmb(l_t1 ~ log_l_t0 + sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=d14,family="nbinom2")
lMod[[15]]=glmmadmb(l_t1 ~ log_l_t0 * sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=d14,family="nbinom2")
AICtab(lMod,weights=T)
#sum(AICtab(lMod,weights=T)$weight[1:4]) #first 3 models make only 66%


#########################################################################################################
###2015######
#########################################################################################################

#only use year two
d15=subset(d,year==2015)
d15$logRat=log(d15$l_t1/d15$l_t0)

#Density as "number of individuals"-----------------------------------------------------------
dMod=list() #dMod stands for "density model" (density is # m/f planted)
#Target fitness
dMod[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=d15,family="nbinom2")
dMod[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=d15,family="nbinom2")
dMod[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=d15,family="nbinom2")
#Target fitness + effect = tot density 
dMod[[4]]=glmmadmb(l_t1 ~ log_l_t0 + TotDensity + (1 | plot),data=d15,family="nbinom2")
dMod[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + TotDensity + (1 | plot),data=d15,family="nbinom2")
dMod[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + TotDensity + (1 | plot),data=d15,family="nbinom2")
#Target fitness + effect = tot density + response=sex 
dMod[[7]]=glmmadmb(l_t1 ~ log_l_t0 + TotDensity + TotDensity:sex + (1 | plot),data=d15,family="nbinom2")
dMod[[8]]=glmmadmb(l_t1 ~ log_l_t0 + sex + TotDensity + TotDensity:sex + (1 | plot),data=d15,family="nbinom2")
dMod[[9]]=glmmadmb(l_t1 ~ log_l_t0 * sex + TotDensity + TotDensity:sex + (1 | plot),data=d15,family="nbinom2")
#Target fitness + effect = sex 
dMod[[10]]=glmmadmb(l_t1 ~ log_l_t0 + M + F + (1 | plot),data=d15,family="nbinom2")
dMod[[11]]=glmmadmb(l_t1 ~ log_l_t0 + sex + M + F + (1 | plot),data=d15,family="nbinom2")
dMod[[12]]=glmmadmb(l_t1 ~ log_l_t0 * sex + M + F + (1 | plot),data=d15,family="nbinom2")
#Target fitness + effect = sex + response=sex 
dMod[[13]]=glmmadmb(l_t1 ~ log_l_t0 + M + F + M:sex + F:sex + (1 | plot),data=d15,family="nbinom2")
dMod[[14]]=glmmadmb(l_t1 ~ log_l_t0 + sex + M + F + M:sex + F:sex + (1 | plot),data=d15,family="nbinom2")
dMod[[15]]=glmmadmb(l_t1 ~ log_l_t0 * sex + M + F + M:sex + F:sex + (1 | plot),data=d15,family="nbinom2")
AICtab(dMod,weights=T)


#########################################################################################################
###DIFFERENCE IN N. of LEAVES######
#########################################################################################################


tiff("Results/VitalRates_simple/growth2Best.tiff",unit="in",width=6.3,height=7,res=600,compression="lzw")

#Set up colors for plots
sexAsInteger=as.integer(d14$sex)
d14$col=as.character(factor(sexAsInteger,labels=c("blue","red")))
d14$symb=as.integer(as.character(factor(sexAsInteger,labels=c("17","16"))))

#Start plotting
par(mfcol=c(2,2),mar=c(3,3,1,0.1),mgp=c(1.4,0.5,0))

#2014----------------------------------------------------------------------------------------
#best model
plot(d14$c_t0,d14$l_t1,pch=d14$symb,xlab="N. of total plot leaves 2013",
     ylab="Target individual: Number of leaves (2014)",col=d14$col)
title(main = "2014: best mod (19% weight)", line=0.2,cex=0.9)
legend(205,105,c("Male targets","Female targets"),
       col=c("red","blue"),lty=c(1,2),lwd=3,bty="n")
xSeq <- seq(0,max(c(d14$c_t0,d14$c_t0)),by=1)
y_m <- exp(coef(lMod[[5]])[1] + coef(lMod[[5]])[2]*mean(d14$log_l_t0,na.rm=T) + 
             coef(lMod[[5]])[3] + coef(lMod[[5]])[4]*xSeq)
y_f <- exp(coef(lMod[[5]])[1] + coef(lMod[[5]])[2]*mean(d14$log_l_t0,na.rm=T) + 
             coef(lMod[[5]])[4]*xSeq)
lines(xSeq,y_m,lwd=3,lty=1,col="red")
lines(xSeq,y_f,lwd=3,lty=2,col="blue")

#2nd best models
plot(d14$c_t0,d14$l_t1,pch=16,xlab="N. of total plot leaves 2013",
     ylab="Target individual: Number of leaves (2014)")
title(main = "2014: 2nd best (18% weight)", line=0.2,cex=0.9)

xSeq <- seq(0,max(d14$c_t0),by=1)
yPred<- exp(coef(lMod[[4]])[1] + coef(lMod[[4]])[2]*mean(d14$log_l_t0,na.rm=T) + 
              coef(lMod[[4]])[3]*xSeq)
lines(xSeq,yPred,lwd=2)


#2015----------------------------------------------------------------------------------------
#best model
plot(d15$M,d15$l_t1,pch=16,xlab="Plot density of male or female individuals",
     ylab="Target individual: Number of leaves (2015)",xlim=c(-1,48))
par(new=T) ;  plot(d15$F+0.5,d15$l_t1,pch=17,xlab="",ylab="",xaxt="n",xlim=c(-1,48),col="grey50")

title(main = "2015: best mod (30% weight)", line=0.2,cex=0.9)
legend(25,85,c("Male density","Female density"),
       col=c("black","grey50"),pch=c(16,17),bty="n")

xSeq <- seq(1,48,by=1)
y_m <- exp(coef(dMod[[10]])[1] + coef(dMod[[10]])[2]*mean(d15$log_l_t0,na.rm=T) + 
              coef(dMod[[10]])[3]*xSeq + coef(dMod[[10]])[4]*mean(d15$F,na.rm=T))
y_f <- exp(coef(dMod[[10]])[1] + coef(dMod[[10]])[2]*mean(d15$log_l_t0,na.rm=T) + 
              coef(dMod[[10]])[4]*xSeq + coef(dMod[[10]])[3]*mean(d15$M,na.rm=T))
lines(xSeq,y_m,lwd=3,lty=1)
lines(xSeq,y_f,lwd=3,lty=1,col="grey50")

#2nd best models
plot(d15$TotDensity,d15$l_t1,pch=16,xlab="Total plot density",
     ylab="Target individual: Number of leaves (2015)")
title(main = "2015: 2nd best (19% weight)", line=0.2,cex=0.9)
xSeq <- seq(1,48,by=1)
ypred <- exp(coef(dMod[[4]])[1] + coef(dMod[[4]])[2]*mean(d15$log_l_t0,na.rm=T) + coef(dMod[[4]])[3]*xSeq)
lines(xSeq,ypred,lwd=2)

dev.off()
