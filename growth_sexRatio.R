##Sex competition predictor: SEX RATIO (female individuals/total individuals)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) #For AIC weights, because I'm lazy!
library(glmmADMB) #Fit models with a Negative Binomial
library(dplyr)

#read in data
d=read.csv("Data/vr.csv")
#remove dead individuals (this is a GROWTH model!)
d=subset(d, surv_t1 != 0)
#Remove a clear mistake: missing data in 2013 and 2014, but 19 leaves in 2015
d=d[-which(d$plot ==36 & d$focalI =="m3"),]

#logtransform leaf numbers
d$log_l_t0=log(d$l_t0)
d$log_l_t1=log(d$l_t1)
d$mC_t0[is.na(d$mC_t0)]=0
d$fC_t0[is.na(d$fC_t0)]=0
d$plot=as.factor(d$plot) #glmmadmb wants plot as a factor

#Transform densities to SEX RATIO
d$sr  <- d$F / d$TotDensity
d$sr2 <- d$sr^2

##################################################################################################
#1.Compare total growth between female only and male only plots#########
##################################################################################################

#only use year two
tmp15=subset(d,year==2015)
tmp15$new_t1[tmp15$new_t1=="SKIPPED"]=NA
tmp15$new_t1[tmp15$new_t1=="cnf"]=NA
tmp15$new_t1[tmp15$new_t1==""]=NA
tmp15$new_t1=as.numeric(as.character(tmp15$new_t1))
d15=na.omit(tmp15[,c("l_t1","log_l_t0","plot","sex","new_t1","sr","sr2")])


#Density as "new tillers"-------------------------------------------------------------------------
#Sex ratio as original sex ratio
lMod=list() #lMod stands for "leaf model" (density is quantified by N. of leaves)
#Target fitness
lMod[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=d15,family="nbinom2")
lMod[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=d15,family="nbinom2")
lMod[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=d15,family="nbinom2")
#Target fitness + effect = tot density 
lMod[[4]]=glmmadmb(l_t1 ~ log_l_t0 + new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + new_t1 + (1 | plot),data=d15,family="nbinom2")
#Target fitness + effect = tot density + response=sex 
lMod[[7]]=glmmadmb(l_t1 ~ log_l_t0 + new_t1 + new_t1:sex + (1 | plot),data=d15,family="nbinom2")
lMod[[8]]=glmmadmb(l_t1 ~ log_l_t0 + sex*new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[9]]=glmmadmb(l_t1 ~ log_l_t0 * sex + new_t1 + new_t1:sex + (1 | plot),data=d15,family="nbinom2")
#Target fitness + effect = sex 
lMod[[10]]=glmmadmb(l_t1 ~ log_l_t0 + sr + new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[11]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[12]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[13]]=glmmadmb(l_t1 ~ log_l_t0 + sr*new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[14]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr*new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[15]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr*new_t1 + (1 | plot),data=d15,family="nbinom2")
#Target fitness + effect = sex + response=sex 
lMod[[16]]=glmmadmb(l_t1 ~ log_l_t0 + sr + new_t1 + sr:sex + new_t1:sex +(1 | plot),data=d15,family="nbinom2")
lMod[[17]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + new_t1 + sr:sex + new_t1:sex + (1 | plot),data=d15,family="nbinom2")
lMod[[18]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + new_t1 + sr:sex + new_t1:sex  + (1 | plot),data=d15,family="nbinom2")
lMod[[19]]=glmmadmb(l_t1 ~ log_l_t0 + sr*new_t1*sex + (1 | plot),data=d15,family="nbinom2")
lMod[[20]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr*new_t1*sex + (1 | plot),data=d15,family="nbinom2")
lMod[[21]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr*new_t1*sex + (1 | plot),data=d15,family="nbinom2")

##QUADRATIC TERMS
lMod[[22]]=glmmadmb(l_t1 ~ log_l_t0 + sr + sr2 + new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[23]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + sr2 +new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[24]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + sr2 +new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[25]]=glmmadmb(l_t1 ~ log_l_t0 + sr*new_t1 + sr2 +(1 | plot),data=d15,family="nbinom2")
lMod[[26]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr*new_t1 + sr2 +(1 | plot),data=d15,family="nbinom2")
lMod[[27]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr*new_t1 + sr2 +(1 | plot),data=d15,family="nbinom2")
#Target fitness + effect = sex + response=sex 
lMod[[28]]=glmmadmb(l_t1 ~ log_l_t0 + sr + sr2 +new_t1 + sr:sex + new_t1:sex +(1 | plot),data=d15,family="nbinom2")
lMod[[29]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + sr2 +new_t1 + sr:sex + new_t1:sex + (1 | plot),data=d15,family="nbinom2")
lMod[[30]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + sr2 +new_t1 + sr:sex + new_t1:sex  + (1 | plot),data=d15,family="nbinom2")
lMod[[31]]=glmmadmb(l_t1 ~ log_l_t0 + sr*new_t1*sex + sr2 +(1 | plot),data=d15,family="nbinom2")
lMod[[32]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr*new_t1*sex + sr2 +(1 | plot),data=d15,family="nbinom2")
lMod[[33]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr*new_t1*sex + sr2 +(1 | plot),data=d15,family="nbinom2")

AICtab(lMod,weights=T)



#########################################################################################################
###GRAPHS######
#########################################################################################################

#Top three models (only 2015)----------------------------------------------------------------------------
tiff("Results/VitalRates_Sept16/growth_SexRatio.tiff",unit="in",width=6.3,height=3.15,res=600,compression="lzw")

#Set up colors for plots
#2015
d15$col=as.integer(d15$sex)
d15$col=as.character(factor(d15$col,labels=c("blue","red")))
d15$symb=as.integer(as.character(factor(d15$col,labels=c("17","16"))))

#Start plotting
par(mfcol=c(1,3),mar=c(3,3,1,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0.1))

#2015----------------------------------------------------------------------------------------
#best model
boxplot(d15$l_t1,ylab="Target individuals: Number of leaves")
title(main = "Best model (17% weight)", line=0.2,cex=0.9)

#2nd best
plot(d15$sr,d15$l_t1, pch=16, ylab="Target individuals: Number of leaves",
     xlab="Proportion of female individuals",col=d15$col)
title(main = "2nd best model (13% weight)", line=0.2,cex=0.9)
legend(-0.05,85,c("Male individuals","Female individuals"),
       lty=1,lwd=3,col=c("red","blue"),bty="n")
xSeq <- seq(0,max(c(d15$sr,d15$sr)),by=0.1)
y_f <- exp(coef(lMod[[11]])[1] + coef(lMod[[11]])[2]*mean(d15$log_l_t0) + 
           coef(lMod[[11]])[4]*xSeq + coef(lMod[[11]])[5]*mean(d15$new_t1))
y_m <- exp(coef(lMod[[11]])[1] + coef(lMod[[11]])[2]*mean(d15$log_l_t0) + coef(lMod[[11]])[3] +
             coef(lMod[[11]])[4]*xSeq + coef(lMod[[11]])[5]*mean(d15$new_t1))
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")

#3nd best
plot(d15$sr,d15$l_t1, pch=16, ylab="Target individuals: Number of leaves",
     xlab="Proportion of female individuals")
title(main = "3rd best model (11% weight)", line=0.2,cex=0.9,xpd=T)
xSeq <- seq(0,max(c(d15$sr,d15$sr)),by=0.1)
y <- exp(coef(lMod[[10]])[1] + coef(lMod[[10]])[2]*mean(d15$log_l_t0) + 
             coef(lMod[[10]])[3]*xSeq + coef(lMod[[10]])[4]*mean(d15$new_t1))
lines(xSeq,y,lwd=3,lty=1)

dev.off()

