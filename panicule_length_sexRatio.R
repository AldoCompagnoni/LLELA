setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) #For AIC weights, because I'm lazy!
library(lme4)
library(nlme)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")

# read in data --------------------------------------------------------------
d             <- read.csv("Data/vr.csv")
femPanicules  <- read.csv("Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")
malPanicules  <- read.csv("Data/Spring 2014/maleCounts/malePaniculesSpring2014.csv")


# format female and male panicule data--------------------------------------------------------------
femPanicules$focalI=paste("f",femPanicules$IndividualN,sep="")
# format column names
names(femPanicules)[c(1,5,7)]=c("plot","panicule_Length_cm","seed_weight_mg")
# format focal individual for male data
tmp=as.numeric(matrix(unlist(strsplit(as.character(malPanicules$Individual),"[A-Z]")),
                      nrow(malPanicules),2,byrow=T)[,2])
malPanicules$focalI=paste("m",tmp,sep="")
names(malPanicules)[1]="plot"
# Add up length of panicules for each individual 
malPanicules=aggregate(panicule_Length_cm ~ plot + focalI, sum, data=malPanicules)


# MODEL SELECTION ######################################################################

# FEMALE panicule lengths -------------------------------------------------------------------------------------
d14=subset(d,year==2014)
d14=subset(d14,surv_t1!=0)
fem_i_pan_data <- merge(d14,femPanicules)

plMod=list()
plMod[[1]]=lmer(log(panicule_Length_cm) ~ log_l_t0+ (1 | plot), data=fem_i_pan_data)
plMod[[2]]=lmer(log(panicule_Length_cm) ~ log_l_t0 + TotDensity+ (1 | plot),data=fem_i_pan_data)
plMod[[3]]=lmer(log(panicule_Length_cm) ~ log_l_t0 + sr+ (1 | plot),data=fem_i_pan_data)
plMod[[4]]=lmer(log(panicule_Length_cm) ~ log_l_t0 + sr + TotDensity + (1 | plot),data=fem_i_pan_data)
plMod[[5]]=lmer(log(panicule_Length_cm) ~ log_l_t0 + sr * TotDensity + (1 | plot),data=fem_i_pan_data)

# Model average
f_i_pl_mod_select <- AICtab(plMod,weights=T)
f_i_pl_avg        <- model_avg(f_i_pl_mod_select, plMod)


# Individual CUMULATIVE female panicule lengths -------------------------------------------------------------------------------------
femPanicules   <- aggregate(cbind(seed_weight_mg,panicule_Length_cm) ~ 
                              plot + focalI, sum, data=femPanicules)
fem_pan_data <- merge(d14,femPanicules)

plMod=list()
plMod[[1]]=lmer(log(panicule_Length_cm) ~ log_l_t0 + (1 | plot), data=fem_pan_data)
plMod[[2]]=lmer(log(panicule_Length_cm) ~ log_l_t0 + TotDensity + (1 | plot),data=fem_pan_data)
plMod[[3]]=lmer(log(panicule_Length_cm) ~ log_l_t0 + sr+ (1 | plot),data=fem_pan_data)
plMod[[4]]=lmer(log(panicule_Length_cm) ~ log_l_t0 + sr + TotDensity + (1 | plot),data=fem_pan_data)
plMod[[5]]=lmer(log(panicule_Length_cm) ~ log_l_t0 + sr * TotDensity + (1 | plot),data=fem_pan_data)

# Model average
f_pl_mod_select <- AICtab(plMod,weights=T)
f_pl_avg        <- model_avg(f_pl_mod_select, plMod)


#CUMULATIVE Panicule lengths for both males and females -----------------------------------------------------------------------
malPanicules$seed_weight_mg=NA
panicules=rbind(femPanicules,malPanicules)
pan14=merge(d14,panicules,all=T)

plMod=list()
#Target fitness
plMod[[1]]=lmer(log(panicule_Length_cm) ~log_l_t0 + (1 | plot),data=pan14)
plMod[[2]]=lmer(log(panicule_Length_cm) ~log_l_t0 + sex + (1 | plot),data=pan14)
plMod[[3]]=lmer(log(panicule_Length_cm) ~log_l_t0 * sex + (1 | plot),data=pan14)
#Target fitness + effect = tot density 
plMod[[4]]=lmer(log(panicule_Length_cm) ~log_l_t0 + TotDensity + (1 | plot),data=pan14)
plMod[[5]]=lmer(log(panicule_Length_cm) ~log_l_t0 + sex + TotDensity + (1 | plot),data=pan14)
plMod[[6]]=lmer(log(panicule_Length_cm) ~log_l_t0 * sex + TotDensity + (1 | plot),data=pan14)
#Target fitness + effect = tot density + response=sex 
plMod[[7]]=lmer(log(panicule_Length_cm) ~log_l_t0 + TotDensity + TotDensity:sex + (1 | plot),data=pan14)
plMod[[8]]=lmer(log(panicule_Length_cm) ~log_l_t0 + sex*TotDensity + (1 | plot),data=pan14)
plMod[[9]]=lmer(log(panicule_Length_cm) ~log_l_t0 * sex + TotDensity + TotDensity:sex + (1 | plot),data=pan14)
#Target fitness + effect = sex 
plMod[[10]]=lmer(log(panicule_Length_cm) ~log_l_t0 + sr + TotDensity + (1 | plot),data=pan14)
plMod[[11]]=lmer(log(panicule_Length_cm) ~log_l_t0 + sex + sr + TotDensity + (1 | plot),data=pan14)
plMod[[12]]=lmer(log(panicule_Length_cm) ~log_l_t0 * sex + sr + TotDensity + (1 | plot),data=pan14)
#Target fitness + effect = sex + response=sex 
plMod[[13]]=lmer(log(panicule_Length_cm) ~log_l_t0 + sr + TotDensity + sr:sex + TotDensity:sex +(1 | plot),data=pan14)
plMod[[14]]=lmer(log(panicule_Length_cm) ~log_l_t0 + sex + sr + TotDensity + sr:sex + TotDensity:sex + (1 | plot),data=pan14)
plMod[[15]]=lmer(log(panicule_Length_cm) ~log_l_t0 * sex + sr + TotDensity + sr:sex + TotDensity:sex  + (1 | plot),data=pan14)

# Model average
f_m_pl_mod_select <- AICtab(plMod,weights=T)
f_m_pl_avg        <- model_avg(f_m_pl_mod_select, plMod)




### OLD CODE, PLEASE IGNORE! #####################################################################

tiff("Results/VitalRates_simple/fertility.tiff",unit="in",width=3.5,height=7,res=600,compression="lzw")

#Start plotting
par(mfcol=c(2,1),mar=c(2.8,3,1,0.1),mgp=c(1.4,0.5,0))

sexAsInteger=as.integer(pan14$sex)
pan14$col=as.character(factor(sexAsInteger,labels=c("blue","red")))

#best model
plot(pan14$log_l_t1,log(pan14$panicule_Length_cm),pch=16,
     ylab="log(Panicule length)",xlab="log(size target)")
xSeq=seq(min(pan14$log_l_t1,na.rm=T),max(pan14$log_l_t1,na.rm=T),length.out=100)
yPred=fixef(pl[[1]])[1] + fixef(pl[[1]])[2]*xSeq
lines(xSeq,yPred,lwd=2)
title(main = "2014: model 1 (26% weight)", line=0.2,cex=0.9)

#2nd best model
plot(pan14$log_l_t1,log(pan14$panicule_Length_cm),pch=16,col=pan14$col,
     ylab="log(Panicule length)",xlab="log(size target)")
xSeq=seq(min(pan14$log_l_t1,na.rm=T),max(pan14$log_l_t1,na.rm=T),length.out=100)
yF=fixef(pl[[3]])[1] + fixef(pl[[3]])[2]*xSeq
yM=fixef(pl[[3]])[1] + fixef(pl[[3]])[2]*xSeq + fixef(pl[[3]])[3] + fixef(pl[[3]])[4]*xSeq
lines(xSeq,yF,lwd=2,col="blue")
lines(xSeq,yM,lwd=2,col="red")
title(main = "2014: model 2 (26% weight)", line=0.2,cex=0.9)

dev.off()



#Graph
sexAsInteger=as.integer(pan14$sex)
pan14$col=as.character(factor(sexAsInteger,labels=c("blue","red")))

plot(pan14$log_l_t1,pan14$panicule_Length_cm,pch=16,col=pan14$col,
     ylab="log(number of target leaves)")
xSeq=seq(min(pan14$log_l_t1,na.rm=T),max(pan14$log_l_t1,na.rm=T),length.out=100)
betas=fixef(pl[[3]])
ym=betas[1] + betas[2]*xSeq + betas[3] + betas[4]*xSeq
yf=betas[1] + betas[2]*xSeq
lines(xSeq,ym,col="red",lwd=2)  
lines(xSeq,yf,col="blue",lwd=2)

dev.off()



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