##1. Selection of models of POAR growth using data from LLELA's competition experiment.  
##2. Year specific data anlayzed SEPARATELY.
##3. BEWARE:  the fields "c_t0","fC_t0", and "mC_t0" refer to FALL 2013, even when "year==2015".
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) #For AIC weights, because I'm lazy!
library(lme4) #Fit models with a Negative Binomial


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
d14$logRat=log(d14$l_t1/d14$l_t0)

ll14=list() #ll14 stands for "leaf model", "log ratio"
#Target fitness
ll14[[1]]=lmer(logRat ~ (1 | plot),data=d14)
ll14[[2]]=lmer(logRat ~ sex + (1 | plot),data=d14)

#Target fitness + effect = tot density 
ll14[[3]]=lmer(logRat ~ c_t0 + (1 | plot),data=d14)
ll14[[4]]=lmer(logRat ~ sex + c_t0 + (1 | plot),data=d14)

#Target fitness + effect = tot density + response=sex 
ll14[[5]]=lmer(logRat ~ c_t0 + c_t0:sex + (1 | plot),data=d14)
ll14[[6]]=lmer(logRat ~ sex + c_t0 + c_t0:sex + (1 | plot),data=d14)

#Target fitness + effect = sex 
ll14[[7]]=lmer(logRat ~ mC_t0 + fC_t0 + (1 | plot),data=d14)
ll14[[8]]=lmer(logRat ~ sex + mC_t0 + fC_t0 + (1 | plot),data=d14)

#Target fitness + effect = sex + response=sex 
ll14[[9]]=lmer(logRat ~ mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=d14)
ll14[[10]]=lmer(logRat ~ sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=d14)

AICtab(ll14,weights=T)



ld14=list() #ld14 stands for "leaf model", "log ratio"
#Density as "number of individuals"-----------------------------------------------------------
#Target fitness
ld14[[1]]=lmer(logRat ~ (1 | plot),data=d14)
ld14[[2]]=lmer(logRat ~ sex + (1 | plot),data=d14)
#Target fitness + effect = tot density 
ld14[[3]]=lmer(logRat ~ TotDensity + (1 | plot),data=d14)
ld14[[4]]=lmer(logRat ~ sex + TotDensity + (1 | plot),data=d14)
#Target fitness + effect = tot density + response=sex 
ld14[[5]]=lmer(logRat ~ TotDensity + TotDensity:sex + (1 | plot),data=d14)
ld14[[6]]=lmer(logRat ~ sex + TotDensity + TotDensity:sex + (1 | plot),data=d14)
#Target fitness + effect = sex 
ld14[[7]]=lmer(logRat ~ M + F + (1 | plot),data=d14)
ld14[[8]]=lmer(logRat ~ sex + M + F + (1 | plot),data=d14)
#Target fitness + effect = sex + response=sex 
ld14[[9]]=lmer(logRat ~ M + F + M:sex + F:sex + (1 | plot),data=d14)
ld14[[10]]=lmer(logRat ~ sex + M + F + M:sex + F:sex + (1 | plot),data=d14)
AICtab(ld14,weights=T)




#########################################################################################################
###2015######
#########################################################################################################
d15=subset(d,year==2015)
d15$logRat=log(d15$l_t1/d15$l_t0)

ll15=list() #ll15 stands for "leaf model", "log ratio"
#Target fitness
ll15[[1]]=lmer(logRat ~ (1 | plot),data=d15)
ll15[[2]]=lmer(logRat ~ sex + (1 | plot),data=d15)
#Target fitness + effect = tot density 
ll15[[3]]=lmer(logRat ~ c_t0 + (1 | plot),data=d15)
ll15[[4]]=lmer(logRat ~ sex + c_t0 + (1 | plot),data=d15)
#Target fitness + effect = tot density + response=sex 
ll15[[5]]=lmer(logRat ~ c_t0 + c_t0:sex + (1 | plot),data=d15)
ll15[[6]]=lmer(logRat ~ sex + c_t0 + c_t0:sex + (1 | plot),data=d15)
#Target fitness + effect = sex 
ll15[[7]]=lmer(logRat ~ mC_t0 + fC_t0 + (1 | plot),data=d15)
ll15[[8]]=lmer(logRat ~ sex + mC_t0 + fC_t0 + (1 | plot),data=d15)
#Target fitness + effect = sex + response=sex 
ll15[[9]]=lmer(logRat ~ mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=d15)
ll15[[10]]=lmer(logRat ~ sex + mC_t0 + fC_t0 + mC_t0:sex + fC_t0:sex + (1 | plot),data=d15)
AICtab(ll15,weights=T)



ld15=list() #ld15 stands for "leaf model", "log ratio"
#Density as "number of individuals"-----------------------------------------------------------
#Target fitness
ld15[[1]]=lmer(logRat ~ (1 | plot),data=d15)
ld15[[2]]=lmer(logRat ~ sex + (1 | plot),data=d15)
#Target fitness + effect = tot density 
ld15[[3]]=lmer(logRat ~ TotDensity + (1 | plot),data=d15)
ld15[[4]]=lmer(logRat ~ sex + TotDensity + (1 | plot),data=d15)
#Target fitness + effect = tot density + response=sex 
ld15[[5]]=lmer(logRat ~ TotDensity + TotDensity:sex + (1 | plot),data=d15)
ld15[[6]]=lmer(logRat ~ sex + TotDensity + TotDensity:sex + (1 | plot),data=d15)
#Target fitness + effect = sex 
ld15[[7]]=lmer(logRat ~ M + F + (1 | plot),data=d15)
ld15[[8]]=lmer(logRat ~ sex + M + F + (1 | plot),data=d15)
#Target fitness + effect = sex + response=sex 
ld15[[9]]=lmer(logRat ~ M + F + M:sex + F:sex + (1 | plot),data=d15)
ld15[[10]]=lmer(logRat ~ sex + M + F + M:sex + F:sex + (1 | plot),data=d15)
AICtab(ld15,weights=T)


#########################################################################################################
###DIFFERENCE IN N. of LEAVES######
#########################################################################################################


tiff("Results/VitalRates_simple/growthLogRat.tiff",unit="in",width=6.3,height=7,res=600,compression="lzw")

#Start plotting
par(mfcol=c(2,2),mar=c(2.8,3,1,0.1),mgp=c(1.4,0.5,0))

#2014----------------------------------------------------------------------------------------
#best model
boxplot(d14$logRat,ylab="log(target leaves 2014 / target leaves 2013)")
title(main = "2014: mod 1 (86% weight)", line=0.2,cex=0.9)

#2nd best model
boxplot(d14$logRat ~ d14$sex,col=c("blue","red"),ylab="log(target leaves 2014 / target leaves 2013)",
        names=c("Females","Males"))
title(main = "2014: mod 2 (14% weight)", line=0.2,cex=0.9)


#2015----------------------------------------------------------------------------------------
#best model
boxplot(d15$logRat,ylab="log(target leaves 2015 / target leaves 2014)")
title(main = "2015: mod 1 (48% weight)", line=0.2,cex=0.9)

#2nd best model
plot(d15$logRat ~ d15$TotDensity,pch=16,xlab="Total plot density",
     ylab="log(target leaves 2015 / target leaves 2014)")
xSeq <- seq(0,48,by=1)
yPred<- fixef(ld15[[3]])[1] + fixef(ld15[[3]])[2]*xSeq
lines(xSeq,yPred,lwd=2)
title(main = "2015: mod 3 (43% weight)", line=0.2,cex=0.9)

dev.off()
