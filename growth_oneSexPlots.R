##Growth analyses using one-sex plots only
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
d$sr = d$F / d$TotDensity


##################################################################################################
#1.Compare total growth between female only and male only plots#########
##################################################################################################
d$oneSex=0
d$oneSex[d$F==0 | d$M==0]=1 #flag treatments with only one sex 
oneSexD=subset(d,oneSex==1)
oneSexD14=subset(oneSexD,year==2014)
tmp15=subset(oneSexD,year==2015)
#exclude plots for which we do not have "new tillers"
tmp15$new_t1[tmp15$new_t1=="SKIPPED"]=NA
tmp15$new_t1[tmp15$new_t1=="cnf"]=NA
tmp15$new_t1[tmp15$new_t1==""]=NA
tmp15$new_t1=as.numeric(as.character(tmp15$new_t1))
oneSexD15=na.omit(tmp15[,c("l_t1","log_l_t0","plot","sex","new_t1")])

ssm14=ssm15=list() #set up lists for models


#2014----------------------------------------------------------------------------------------------  
#Effect of sex
ssm14[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=oneSexD14,family="nbinom2")
ssm14[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=oneSexD14,family="nbinom2")
ssm14[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=oneSexD14,family="nbinom2")
#Sex + density
ssm14[[4]]=glmmadmb(l_t1 ~ log_l_t0 + TotDensity + (1 | plot),data=oneSexD14,family="nbinom2")
ssm14[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + TotDensity + (1 | plot),data=oneSexD14,family="nbinom2")
ssm14[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + TotDensity + (1 | plot),data=oneSexD14,family="nbinom2")
#Sex + density + sex response to density 
ssm14[[7]]=glmmadmb(l_t1 ~ log_l_t0 + sex + TotDensity + TotDensity:sex +(1 | plot),data=oneSexD14,family="nbinom2")
ssm14[[8]]=glmmadmb(l_t1 ~ log_l_t0 * sex + TotDensity + TotDensity:sex +(1 | plot),data=oneSexD14,family="nbinom2")
AICtab(ssm14,weights=T)


#2015----------------------------------------------------------------------------------------------  
#Effect of sex
ssm15[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=oneSexD15,family="nbinom2")
ssm15[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=oneSexD15,family="nbinom2")
ssm15[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=oneSexD15,family="nbinom2")
#Sex + density
ssm15[[4]]=glmmadmb(l_t1 ~ log_l_t0 + new_t1 + (1 | plot),data=oneSexD15,family="nbinom2")
ssm15[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + new_t1 + (1 | plot),data=oneSexD15,family="nbinom2")
ssm15[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + new_t1 + (1 | plot),data=oneSexD15,family="nbinom2")
#Sex + density + sex response to density 
ssm15[[7]]=glmmadmb(l_t1 ~ log_l_t0 + sex + new_t1 + new_t1:sex +(1 | plot),data=oneSexD15,family="nbinom2")
ssm15[[8]]=glmmadmb(l_t1 ~ log_l_t0 * sex + new_t1 + new_t1:sex +(1 | plot),data=oneSexD15,family="nbinom2")
AICtab(ssm15,weights=T)


#########################################################################################################
###GRAPHS######
#########################################################################################################

#Top two models 2014 and 2015----------------------------------------------------------------------------
tiff("Results/VitalRates_Sept16/growth_OneSexPlots.tiff",unit="in",width=6.3,height=7,res=600,compression="lzw")

#Set up colors for plots
#2015
oneSexD15$col=as.integer(oneSexD15$sex)
oneSexD15$col=as.character(factor(oneSexD15$col,labels=c("blue","red")))
oneSexD15$symb=as.integer(as.character(factor(oneSexD15$col,labels=c("17","16"))))
#2014
oneSexD14$col=as.integer(oneSexD14$sex)
oneSexD14$col=as.character(factor(oneSexD14$col,labels=c("blue","red")))
oneSexD14$symb=as.integer(as.character(factor(oneSexD14$col,labels=c("17","16"))))

#Start plotting
par(mfcol=c(2,2),mar=c(3,3,1,0.1),mgp=c(1.4,0.5,0))

#2014----------------------------------------------------------------------------------------
#best model
plot(oneSexD14$TotDensity,oneSexD14$l_t1,pch=16,xlab="N. of total plot leaves 2013",
     ylab="Target individuals: Number of leaves at time t")
title(main = "2014: Model 4 (36% weight)", line=0.2,cex=0.9)
xSeq <- seq(0,max(c(oneSexD14$TotDensity,oneSexD14$TotDensity)),by=1)
y <- exp(coef(ssm14[[4]])[1] + coef(ssm14[[4]])[2]*mean(oneSexD14$log_l_t0,na.rm=T) + 
             coef(ssm14[[4]])[3]*xSeq)
lines(xSeq,y,lwd=3,lty=1,col="black")

#2nd best models
boxplot(oneSexD14$l_t1,ylab="Target individuals: Number of leaves at time t")
title(main = "2014: Model 1 (16% weight)", line=0.2,cex=0.9)


#2015----------------------------------------------------------------------------------------
#best model
boxplot(oneSexD15$l_t1,ylab="Target individuals: Number of leaves at time t")
title(main = "2015: Model 1 (43% weight)", line=0.2,cex=0.9)

plot(oneSexD15$log_l_t0,oneSexD15$l_t1, pch=16, ylab="Target individuals: Number of leaves at time t",
     xlab="Target individuals: size at time t-1",col=oneSexD15$col)
title(main = "2015: Model 1 (17% weight)", line=0.2,cex=0.9)
legend(1.5,80,c("Male individuals","Female individuals"),
       lty=1,lwd=3,col=c("blue","red"),bty="n")
xSeq <- seq(0,max(c(oneSexD15$log_l_t0,oneSexD15$log_l_t0)),by=0.1)
y_f <- exp(coef(ssm15[[2]])[1] + coef(ssm15[[2]])[2]*xSeq)
y_m <- exp(coef(ssm15[[2]])[1] + coef(ssm15[[2]])[2]*xSeq + coef(ssm15[[2]])[3])
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")

dev.off()

