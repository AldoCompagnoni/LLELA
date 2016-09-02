##Response variable: total number of tillers 
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
#squared density
d$TotDensity2 = d$TotDensity^2


#Year two
tmp15=subset(d,year==2015)
#remove 
tmp15$new_t1[tmp15$new_t1=="SKIPPED"]=NA
tmp15$new_t1[tmp15$new_t1=="cnf"]=NA
tmp15$new_t1[tmp15$new_t1==""]=NA
tmp15$new_t1=as.numeric(as.character(tmp15$new_t1))
d15=tmp15
#d15=na.omit(tmp15[,c("l_t1","log_l_t0","plot","sex","new_t1","sr","sr2")])


##################################################################################################
#Response variable: total number of tillers 
##################################################################################################

#INITIAL STRUCTURES
#Effect of plot treatment on new tillers#################################
nt=list()
#Effect of density
nt[[1]]=lm(new_t1 ~ TotDensity,data=d15)
nt[[2]]=lm(new_t1 ~ TotDensity + TotDensity2,data=d15)
#Effect of density + sex ratio
nt[[3]]=lm(new_t1 ~ TotDensity + sr,data=d15)
nt[[4]]=lm(new_t1 ~ TotDensity + TotDensity2 + sr,data=d15)
#Interactions b/w density and sex ratio
nt[[5]]=lm(new_t1 ~ sr * TotDensity,data=d15)
nt[[6]]=lm(new_t1 ~ sr * TotDensity + TotDensity2,data=d15)
nt[[7]]=lm(new_t1 ~ sr * TotDensity + TotDensity2*sr,data=d15)
#nt[[6]]=lm(new_t1 ~ sr + sr2 + TotDensity,data=d15)
#nt[[7]]=lm(new_t1 ~ sr * TotDensity + sr2,data=d15)
#nt[[8]]=lm(new_t1 ~ sr + sr2 + TotDensity + sr:TotDensity + sr2:TotDensity,data=d15)
AICtab(nt,weights=T)


##################################################################################################
#Graph: total number of tillers 
##################################################################################################

tiff("Results/VitalRates_Sept16/growth_newTillers.tiff",unit="in",width=6.3,height=3.15,res=600,compression="lzw")


par(mfcol=c(1,2),mar=c(3,3,1,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0.1))

#Best model (61%)
plot(d15$TotDensity,d15$new_t1,pch=16,ylab="Number of new tillers in 2015",
     xlab="Plot density")
title(main = "Best model (61% weight)", line=0.2,cex=0.9)
legend(0,230,c("10% females","90% females"),lwd=2,col=c("red","blue"),bty="n")

#Calculate lines
beta=coef(nt[[7]])
xSeq=seq(0,max(d15$TotDensity),by=1)
hSr=0.9 ; lSr=0.1 #high vs low female sex ratio
yL <-  beta[1] + beta[2]*lSr + beta[3]*xSeq + beta[4]*xSeq^2 + 
       beta[5]*(lSr*xSeq) + beta[6]*(lSr*xSeq^2)
yH <-  beta[1] + beta[2]*hSr + beta[3]*xSeq + beta[4]*xSeq^2 + 
       beta[5]*(hSr*xSeq) + beta[6]*(hSr*xSeq^2)
lines(xSeq,yL,col="red",lwd=2)
lines(xSeq,yH,col="blue",lwd=2)

#2nd best model (33%)
plot(d15$TotDensity,d15$new_t1,pch=16,ylab="Number of new tillers in 2015",
     xlab="Plot density")
title(main = "2nd best model (33% weight)", line=0.2,cex=0.9)

beta=coef(nt[[6]])
xSeq=seq(0,max(d15$TotDensity),by=1)
hSr=0.9
lSr=0.1
yL <-  beta[1] + beta[2]*lSr + beta[3]*xSeq + beta[4]*xSeq^2 + 
       beta[5]*(lSr*xSeq)
yH <-  beta[1] + beta[2]*hSr + beta[3]*xSeq + beta[4]*xSeq^2 + 
       beta[5]*(hSr*xSeq)
lines(xSeq,yL,col="red",lwd=2)
lines(xSeq,yH,col="blue",lwd=2)


dev.off()

