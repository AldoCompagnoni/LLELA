##Sex predictor of flowering probability: SEX RATIO (female individuals/total individuals)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) #For AIC weights, because I'm lazy!
library(lme4) ; library(glmmADMB)

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
#make "plot" a factor to fit models 
d$plot=as.factor(d$plot) 

#Transform densities to SEX RATIO
d$sr  <- d$F / d$TotDensity
d$sr2 <- d$sr^2


#YEAR ONE##########################################################################################################
tmp=subset(d,year==2014)
tmp[which(tmp$log_l_t1==-Inf),]=NA
#prepare the 'new_t1' column
tmp$new_t1[tmp$new_t1=="SKIPPED"]=NA
tmp$new_t1[tmp$new_t1=="cnf"]=NA
tmp$new_t1[tmp$new_t1==""]=NA
tmp$new_t1=as.numeric(as.character(tmp$new_t1))
f14=na.omit(tmp[,c("plot","flow_t1","log_l_t1","sex","sr","sr2","new_t1")])


############################################################################################
#Model selection
############################################################################################

#Target fitness
fMod=list()
fMod[[1]]=glmmadmb(flow_t1 ~ log_l_t1 + (1 | plot),data=f14,family="binomial")
fMod[[2]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + (1 | plot),data=f14,family="binomial")
fMod[[3]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + (1 | plot),data=f14,family="binomial")
#Target fitness + effect = tot density 
fMod[[4]]=glmmadmb(flow_t1 ~ log_l_t1 + new_t1 + (1 | plot),data=f14,family="binomial")
fMod[[5]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + new_t1 + (1 | plot),data=f14,family="binomial")
fMod[[6]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + new_t1 + (1 | plot),data=f14,family="binomial")
#Target fitness + effect = tot density + response=sex 
fMod[[7]]=glmmadmb(flow_t1 ~ log_l_t1 + new_t1 + new_t1:sex + (1 | plot),data=f14,family="binomial")
fMod[[8]]=glmmadmb(flow_t1 ~ log_l_t1 + sex*new_t1 + (1 | plot),data=f14,family="binomial")
fMod[[9]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + new_t1 + new_t1:sex + (1 | plot),data=f14,family="binomial")
#Target fitness + effect = sex 
fMod[[10]]=glmmadmb(flow_t1 ~ log_l_t1 + sr + new_t1 + (1 | plot),data=f14,family="binomial")
fMod[[11]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + sr + new_t1 + (1 | plot),data=f14,family="binomial")
fMod[[12]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + sr + new_t1 + (1 | plot),data=f14,family="binomial")
fMod[[13]]=glmmadmb(flow_t1 ~ log_l_t1 + sr*new_t1 + (1 | plot),data=f14,family="binomial")
fMod[[14]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + sr*new_t1 + (1 | plot),data=f14,family="binomial")
fMod[[15]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + sr*new_t1 + (1 | plot),data=f14,family="binomial")
#Target fitness + effect = sex + response=sex 
fMod[[16]]=glmmadmb(flow_t1 ~ log_l_t1 + sr + new_t1 + sr:sex + new_t1:sex +(1 | plot),data=f14,family="binomial")
fMod[[17]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + sr + new_t1 + sr:sex + new_t1:sex + (1 | plot),data=f14,family="binomial")
fMod[[18]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + sr + new_t1 + sr:sex + new_t1:sex  + (1 | plot),data=f14,family="binomial")
fMod[[19]]=glmmadmb(flow_t1 ~ log_l_t1 + sr*new_t1*sex + (1 | plot),data=f14,family="binomial")
fMod[[20]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + sr*new_t1*sex + (1 | plot),data=f14,family="binomial")
fMod[[21]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + sr*new_t1*sex + (1 | plot),data=f14,family="binomial")

##QUADRATIC TERMS
fMod[[22]]=glmmadmb(flow_t1 ~ log_l_t1 + sr + sr2 + new_t1 + (1 | plot),data=f14,family="binomial")
fMod[[23]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + sr + sr2 +new_t1 + (1 | plot),data=f14,family="binomial")
fMod[[24]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + sr + sr2 +new_t1 + (1 | plot),data=f14,family="binomial")
fMod[[25]]=glmmadmb(flow_t1 ~ log_l_t1 + sr*new_t1 + sr2 +(1 | plot),data=f14,family="binomial")
fMod[[26]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + sr*new_t1 + sr2 +(1 | plot),data=f14,family="binomial")
fMod[[27]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + sr*new_t1 + sr2 +(1 | plot),data=f14,family="binomial")
#Target fitness + effect = sex + response=sex 
fMod[[28]]=glmmadmb(flow_t1 ~ log_l_t1 + sr + sr2 +new_t1 + sr:sex + new_t1:sex +(1 | plot),data=f14,family="binomial")
fMod[[29]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + sr + sr2 +new_t1 + sr:sex + new_t1:sex + (1 | plot),data=f14,family="binomial")
fMod[[30]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + sr + sr2 +new_t1 + sr:sex + new_t1:sex  + (1 | plot),data=f14,family="binomial")
fMod[[31]]=glmmadmb(flow_t1 ~ log_l_t1 + sr*new_t1*sex + sr2 +(1 | plot),data=f14,family="binomial")
fMod[[32]]=glmmadmb(flow_t1 ~ log_l_t1 + sex + sr*new_t1*sex + sr2 +(1 | plot),data=f14,family="binomial")
fMod[[33]]=glmmadmb(flow_t1 ~ log_l_t1 * sex + sr*new_t1*sex + sr2 +(1 | plot),data=f14,family="binomial")

AICtab(fMod,weights=T)



#########################################################################################################
###GRAPHS#
#########################################################################################################

tiff("Results/VitalRates_Sept16/flowering.tiff",unit="in",width=6.3,height=3.15,res=600,compression="lzw")

#Set up colors for plots
f14$col=as.character(factor(as.integer(f14$sex),labels=c("blue","red")))
mal14=subset(f14,sex=="m") ; fem14=subset(f14,sex=="f")

par(mfrow=c(1,2),mar=c(3,3,1,0.1),mgp=c(1.4,0.35,0),cex.lab=0.8,cex.axis=0.8,
    cex.main=0.9)

#Best model
plot(mal14$flow_t1+0.015 ~ mal14$log_l_t1,
     xlab=expression("Target individuals: log(leaves in 2014)"),ylim=c(0,1),
     ylab="Flowering probability, spring (2014)",pch=16,col=mal14$col)
par(new=T) ; plot(fem14$flow_t1-0.015 ~ fem14$log_l_t1,pch=17,xlab="",ylab="",col=fem14$col,xaxt="n",ylim=c(0,1))
title(main = "2015: best mod (16% weight)", line=0.2,cex=0.5)

beta=coef(fMod[[6]])
xSeq <- seq(min(f14$log_l_t1),max(f14$log_l_t1),by=0.1)
hD=quantile(f14$new_t1,probs=c(0.8)) ; lD=quantile(f14$new_t1,probs=c(0.2)) 
yf_h=boot::inv.logit(beta[1] + beta[2]*xSeq + beta[3] + beta[4]*hD) 
ym_h=boot::inv.logit(beta[1] + beta[2]*xSeq + beta[3] + beta[4]*hD + beta[5]*xSeq)
yf_l=boot::inv.logit(beta[1] + beta[2]*xSeq + beta[3] + beta[4]*lD) 
ym_l=boot::inv.logit(beta[1] + beta[2]*xSeq + beta[3] + beta[4]*lD + beta[5]*xSeq)
lines(xSeq,yf_h,lwd=2,col="blue")
lines(xSeq,ym_h,lwd=2,col="red")
lines(xSeq,yf_l,lwd=2,lty=2,col="blue")
lines(xSeq,ym_l,lwd=2,lty=2,col="red")

legend(2.38,.35,c("Females, high density","Males, high density",
                 "Females, low density","Males, low density"),cex=0.7,
                  col=c("red","blue","red","blue"),pch=c(16,17,16,17),
                  lty=c(1,1,2,2),lwd=2,bty="n")


#2nd best model
plot(mal14$flow_t1+0.015 ~ mal14$log_l_t1,
     xlab=expression("Target individuals: log(leaves in 2014)"),ylim=c(0,1),
     ylab="Flowering probability, spring (2014)",pch=16,col=mal14$col)
par(new=T) ; plot(fem14$flow_t1-0.015 ~ fem14$log_l_t1,pch=17,xlab="",ylab="",col=fem14$col,xaxt="n",ylim=c(0,1))
title(main = "2015: best mod (12% weight)", line=0.2,cex=0.7)

beta=coef(fMod[[9]])
xSeq <- seq(min(f14$log_l_t1),max(f14$log_l_t1),by=0.1)
hD=quantile(f14$new_t1,probs=c(0.8)) ; lD=quantile(f14$new_t1,probs=c(0.2)) 
yf_h=boot::inv.logit(beta[1] + beta[2]*xSeq + beta[3] + beta[4]*hD) 
ym_h=boot::inv.logit(beta[1] + beta[2]*xSeq + beta[3] + beta[4]*hD + beta[5]*xSeq + beta[6]*hD)
yf_l=boot::inv.logit(beta[1] + beta[2]*xSeq + beta[3] + beta[4]*lD) 
ym_l=boot::inv.logit(beta[1] + beta[2]*xSeq + beta[3] + beta[4]*lD + beta[5]*xSeq + beta[6]*lD)
lines(xSeq,yf_h,lwd=2,col="blue")
lines(xSeq,ym_h,lwd=2,col="red")
lines(xSeq,yf_l,lwd=2,lty=2,col="blue")
lines(xSeq,ym_l,lwd=2,lty=2,col="red")

dev.off()
