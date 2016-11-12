##Sex predictor of survival: SEX RATIO (female individuals/total individuals)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) #For AIC weights, because I'm lazy!
library(lme4) ; library(glmmADMB)

#read in data
d=read.csv("Data/vr.csv")

#Remove resuscitated individuals (dead in spring 2014, alive in 2015)
#NOTE: possibly these are NOT DEAD in 2014, because they all have
#few leaves: resprouted from base?
d=d[-which(d$plot == 18 & d$focalI == "m2"),]
d=d[-which(d$plot == 38 & d$focalI == "f1"),]
d=d[-which(d$plot == 46 & d$focalI == "f1"),]
d=d[-which(d$plot == 83 & d$focalI == "f5"),]
d=d[-which(d$plot == 36 & d$focalI == "m3"),]

#logtransform leaf numbers
d$log_l_t0=log(d$l_t0)
d$log_l_t1=log(d$l_t1)

#Transform densities to SEX RATIO
d$sr  <- d$F / d$TotDensity
d$sr2 <- d$sr^2


#YEAR TWO##########################################################################################################
d$plot=as.factor(d$plot)
tmp=subset(d,year==2015)
#prepare the 'new_t1' column
tmp$new_t1[tmp$new_t1=="SKIPPED"]=NA
tmp$new_t1[tmp$new_t1=="cnf"]=NA
tmp$new_t1[tmp$new_t1==""]=NA
tmp$new_t1=as.numeric(as.character(tmp$new_t1))
s15=na.omit(tmp[,c("plot","surv_t1","log_l_t0","sex","new_t1","sr","sr2")])


############################################################################################
#Model selection
############################################################################################

#Target fitness
sMod=list()
sMod[[1]]=glmmadmb(surv_t1 ~ log_l_t0 + (1 | plot),data=s15,family="binomial")
sMod[[2]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + (1 | plot),data=s15,family="binomial")
sMod[[3]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + (1 | plot),data=s15,family="binomial")
#Target fitness + effect = tot density 
sMod[[4]]=glmmadmb(surv_t1 ~ log_l_t0 + new_t1 + (1 | plot),data=s15,family="binomial")
sMod[[5]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + new_t1 + (1 | plot),data=s15,family="binomial")
sMod[[6]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + new_t1 + (1 | plot),data=s15,family="binomial")
#Target fitness + effect = tot density + response=sex 
sMod[[7]]=glmmadmb(surv_t1 ~ log_l_t0 + new_t1 + new_t1:sex + (1 | plot),data=s15,family="binomial")
sMod[[8]]=glmmadmb(surv_t1 ~ log_l_t0 + sex*new_t1 + (1 | plot),data=s15,family="binomial")
sMod[[9]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + new_t1 + new_t1:sex + (1 | plot),data=s15,family="binomial")
#Target fitness + effect = sex 
sMod[[10]]=glmmadmb(surv_t1 ~ log_l_t0 + sr + new_t1 + (1 | plot),data=s15,family="binomial")
sMod[[11]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + sr + new_t1 + (1 | plot),data=s15,family="binomial")
sMod[[12]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + sr + new_t1 + (1 | plot),data=s15,family="binomial")
sMod[[13]]=glmmadmb(surv_t1 ~ log_l_t0 + sr*new_t1 + (1 | plot),data=s15,family="binomial")
sMod[[14]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + sr*new_t1 + (1 | plot),data=s15,family="binomial")
sMod[[15]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + sr*new_t1 + (1 | plot),data=s15,family="binomial")
#Target fitness + effect = sex + response=sex 
sMod[[16]]=glmmadmb(surv_t1 ~ log_l_t0 + sr + new_t1 + sr:sex + new_t1:sex +(1 | plot),data=s15,family="binomial")
sMod[[17]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + sr + new_t1 + sr:sex + new_t1:sex + (1 | plot),data=s15,family="binomial")
sMod[[18]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + sr + new_t1 + sr:sex + new_t1:sex  + (1 | plot),data=s15,family="binomial")
sMod[[19]]=glmmadmb(surv_t1 ~ log_l_t0 + sr*new_t1*sex + (1 | plot),data=s15,family="binomial")
sMod[[20]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + sr*new_t1*sex + (1 | plot),data=s15,family="binomial")
sMod[[21]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + sr*new_t1*sex + (1 | plot),data=s15,family="binomial")

##QUADRATIC TERMS
sMod[[22]]=glmmadmb(surv_t1 ~ log_l_t0 + sr + sr2 + new_t1 + (1 | plot),data=s15,family="binomial")
sMod[[23]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + sr + sr2 +new_t1 + (1 | plot),data=s15,family="binomial")
sMod[[24]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + sr + sr2 +new_t1 + (1 | plot),data=s15,family="binomial")
sMod[[25]]=glmmadmb(surv_t1 ~ log_l_t0 + sr*new_t1 + sr2 +(1 | plot),data=s15,family="binomial")
sMod[[26]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + sr*new_t1 + sr2 +(1 | plot),data=s15,family="binomial")
sMod[[27]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + sr*new_t1 + sr2 +(1 | plot),data=s15,family="binomial")
#Target fitness + effect = sex + response=sex 
sMod[[28]]=glmmadmb(surv_t1 ~ log_l_t0 + sr + sr2 +new_t1 + sr:sex + new_t1:sex +(1 | plot),data=s15,family="binomial")
sMod[[29]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + sr + sr2 +new_t1 + sr:sex + new_t1:sex + (1 | plot),data=s15,family="binomial")
sMod[[30]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + sr + sr2 +new_t1 + sr:sex + new_t1:sex  + (1 | plot),data=s15,family="binomial")
sMod[[31]]=glmmadmb(surv_t1 ~ log_l_t0 + sr*new_t1*sex + sr2 +(1 | plot),data=s15,family="binomial")
sMod[[32]]=glmmadmb(surv_t1 ~ log_l_t0 + sex + sr*new_t1*sex + sr2 +(1 | plot),data=s15,family="binomial")
sMod[[33]]=glmmadmb(surv_t1 ~ log_l_t0 * sex + sr*new_t1*sex + sr2 +(1 | plot),data=s15,family="binomial")

AICtab(sMod,weights=T)



#########################################################################################################
###GRAPHS#
#########################################################################################################


tiff("Results/VitalRates_Sept16/survival_tillers.tiff",unit="in",width=6.3,height=3.15,res=600,compression="lzw")

s15$col=as.character(factor(as.integer(s15$sex),labels=c("blue","red")))
s15$symb=as.integer(as.character(factor(s15$col,labels=c("17","16"))))

#set up graph parameters
par(mfrow=c(1,2),mar=c(3,3,1,0.1),mgp=c(1.4,0.5,0))

#best model
plot(s15$surv_t1+0.015 ~ s15$sr,pch=16,xlab="Proportion of female individuals",ylim=c(0,1),
     ylab="Survival of target individuals")
title(main = "Best model (23% weight)",line=0.2,cex=0.8)
legend(-0.05,0.35,c("High density","Low density"), col="black",
       lty=c(1,2),lwd=2,bty="n") 
#text(-15,0.1,"Predictor: Density by sex",pos=4,cex=1.2)
xSeq <- seq(0,max(c(s15$sr,s15$sr)),by=0.1)
hD <- quantile(s15$new_t1,probs=0.8) ; lD <- quantile(s15$new_t1,probs=0.2) 
yH=boot::inv.logit(coef(sMod[[10]])[1] + coef(sMod[[10]])[2]*mean(s15$log_l_t0,na.rm=T) +
                    coef(sMod[[10]])[3]*xSeq + coef(sMod[[10]])[4]*hD)
yL=boot::inv.logit(coef(sMod[[10]])[1] + coef(sMod[[10]])[2]*mean(s15$log_l_t0,na.rm=T) +
                     coef(sMod[[10]])[3]*xSeq + coef(sMod[[10]])[4]*lD)
lines(xSeq,yH,lty=1,lwd=2,col="black")
lines(xSeq,yL,lty=2,lwd=2,col="black")

#2nd best models
plot(s15$surv_t1+0.015 ~ s15$sr,pch=16,xlab="Proportion of female individuals",ylim=c(0,1),
     ylab="Survival of target individuals")
title(main = "2nd best mod (11% weight)",line=0.2,cex=0.8)

xSeq <- seq(0,max(c(s15$sr,s15$sr)),by=0.1)
hD <- quantile(s15$new_t1,probs=0.8) ; lD <- quantile(s15$new_t1,probs=0.2) 
yH=boot::inv.logit(coef(sMod[[22]])[1] + coef(sMod[[22]])[2]*mean(s15$log_l_t0,na.rm=T) +
                    coef(sMod[[22]])[3]*xSeq + coef(sMod[[22]])[4]*xSeq^2 + 
                    coef(sMod[[22]])[5]*hD)
yL=boot::inv.logit(coef(sMod[[22]])[1] + coef(sMod[[22]])[2]*mean(s15$log_l_t0,na.rm=T) +
                    coef(sMod[[22]])[3]*xSeq + coef(sMod[[22]])[4]*xSeq^2 + 
                    coef(sMod[[22]])[5]*lD)
lines(xSeq,yH,lty=1,lwd=2,col="black")
lines(xSeq,yL,lty=2,lwd=2,col="black")

dev.off()
