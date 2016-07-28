## AC: 7.20.2015
##1. Selection of models of POAR survival using data from LLELA's competition experiment.
##2. Year specific data anlayzed SEPARATELY.
##3.The fields "c_t0","fC_t0", and "mC_t0" refer to FALL 2013, even when "year==2015".
## Therefore, final file IS MISLEADING. Take this into account.
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) #For AIC weights, because I'm lazy!
library(lme4)
library(nlme)

#read in data--------------------------------------------------------------
d=read.csv("Data/vr.csv")
femPanicules=read.csv("Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")
malPanicules=read.csv("Data/Spring 2014/maleCounts/malePaniculesSpring2014.csv")

#Remove resuscitated individuals (dead in spring 2014, alive in 2015)
#NOTE: probably more adequate to consider these NOT DEAD in 2014, because they all have
#few leaves: they probably resprouted from base?
d=d[-which(d$plot ==18 & d$focalI =="m2"),]
d=d[-which(d$plot ==38 & d$focalI =="f1"),]
d=d[-which(d$plot ==46 & d$focalI =="f1"),]
d=d[-which(d$plot ==83 & d$focalI =="f5"),]
d=d[-which(d$plot ==36 & d$focalI =="m3"),]

#logtransform leaf numbers
d$log_l_t0=log(d$l_t0)
d$log_l_t1=log(d$l_t1)

##format female and male panicule data--------------------------------------------------------------
femPanicules$focalI=paste("f",femPanicules$IndividualN,sep="")
names(femPanicules)[c(1,5,7)]=c("plot","panicule_Length_cm","seed_weight_mg")
tmp=as.numeric(matrix(unlist(strsplit(as.character(malPanicules$Individual),"[A-Z]")),
                      nrow(malPanicules),2,byrow=T)[,2])
malPanicules$focalI=paste("m",tmp,sep="")
names(malPanicules)[1]="plot"
malPanicules=aggregate(panicule_Length_cm ~ plot + focalI, sum, data=malPanicules)


##########################################################################################################
#ANALYSIS#######################################################
##########################################################################################################
d14=subset(d,year==2014)
d14=subset(d14,surv_t1!=0)

#Panicule length - male/female-----------------------------------------------------------------------
femPanicules=aggregate(cbind(seed_weight_mg,panicule_Length_cm) ~ plot + focalI, sum, data=femPanicules)
malPanicules$seed_weight_mg=NA
panicules=rbind(femPanicules,malPanicules)
pan14=merge(d14,panicules,all=T)

pl=list() #lMod stands for "leaf model" (density is quantified by N. of leaves)
#Target fitness
pl[[1]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 + (1 | plot),REML = F,data=pan14)
pl[[2]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 + sex + (1 | plot),REML = F,data=pan14)
pl[[3]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 * sex + (1 | plot),REML = F,data=pan14)
#Target fitness + effect = tot density 
pl[[4]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 + c_t0 + (1 | plot),REML = F,data=pan14)
pl[[5]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 + sex + c_t0 + (1 | plot),REML = F,data=pan14)
pl[[6]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 * sex + c_t0 + (1 | plot),REML = F,data=pan14)
#Target fitness + effect = tot density + response=sex 
pl[[7]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 + c_t0 + c_t0:sex + (1 | plot),REML = F,data=pan14)
pl[[8]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 + sex + c_t0 + c_t0:sex + (1 | plot),REML = F,data=pan14)
pl[[9]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 * sex + c_t0 + c_t0:sex + (1 | plot),REML = F,data=pan14)
#Target fitness + effect = sex 
pl[[10]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 + mC_t0 + fC_t0 + (1 | plot),REML = F,data=pan14)
pl[[11]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 + sex + mC_t0 + fC_t0 + (1 | plot),REML = F,data=pan14)
pl[[12]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 * sex + mC_t0 + fC_t0 + (1 | plot),REML = F,data=pan14)
#Target fitness + effect = sex + response=sex 
pl[[13]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),REML = F,data=pan14)
pl[[14]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 + sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),REML = F,data=pan14)
pl[[15]]=lmer(panicule_Length_cm^(1/4) ~ log_l_t1 * sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),REML = F,data=pan14)
pl[[16]]=lmer(panicule_Length_cm^(1/4) ~ (1 | plot),REML = F,data=pan14) #null model
AICtab(pl,weights=T) #BEST MODEL IS MODEL 1

#Graph
#sexAsInteger=as.integer(pan14$sex)
#pan14$col=as.character(factor(sexAsInteger,labels=c("blue","red")))

#plot(pan14$log_l_t1,pan14$panicule_Length_cm,pch=16,col=pan14$col,
#     ylab="log(number of target leaves)")
#xSeq=seq(min(pan14$log_l_t1,na.rm=T),max(pan14$log_l_t1,na.rm=T),length.out=100)
#betas=fixef(pl[[3]])
#ym=betas[1] + betas[2]*xSeq + betas[3] + betas[4]*xSeq
#yf=betas[1] + betas[2]*xSeq
#lines(xSeq,ym,col="red",lwd=2)  
#lines(xSeq,yf,col="blue",lwd=2)


#Female seed weight----------------------------------------------------------------------------------
fSeed=merge(d14,femPanicules)
fSeed$logSw=log(fSeed$seed_weight_mg)
fSeed$logSw[fSeed$logSw==-Inf]=NA

fs=list() #lMod stands for "leaf model" (density is quantified by N. of leaves)
fs[[1]]=lmer(logSw ~ log_l_t1 + (1 | plot),REML = F,data=fSeed)
#Effect of total density 
fs[[2]]=lmer(logSw ~ log_l_t1 + c_t0 + (1 | plot),REML = F,data=fSeed)
#Effect of male/female density
fs[[3]]=lmer(logSw ~ log_l_t1 + mC_t0 + fC_t0 + (1 | plot),REML = F,data=fSeed)
AICtab(fs,weights=T)

plot(fSeed$fC_t0,fSeed$logSw,pch=17,col="grey50",xaxt="n",xlab="")
par(new=T) ; plot(fSeed$mC_t0,fSeed$logSw,pch=16)
xSeq=seq(0,429,1)
betas=fixef(fs[[3]])
ym=betas[1] + betas[2]*mean(fSeed$log_l_t1,na.rm=T) +  betas[3]*xSeq + betas[4]*mean(fSeed$fC_t0,na.rm=T)
yf=betas[1] + betas[2]*mean(fSeed$log_l_t1,na.rm=T) +  betas[3]*mean(fSeed$fC_t0,na.rm=T) + betas[4]*xSeq
lines(xSeq,ym,lwd=2)
lines(xSeq,yf,lwd=2,col="grey50")
