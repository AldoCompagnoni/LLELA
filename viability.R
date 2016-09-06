##Model selection for seed VIABILITY (response variable with support [0,1])
##This script uses tetrazolium scoring and germination data 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(lme4) ; library(bbmle) ; library(testthat)
library(mgcv)

#Read in data and format###########################################################################################
germ=read.csv("Data/Spring 2014/viability/germination_data_LLELA.csv")
viab=read.csv("Data/Spring 2014/viability/viabilityData.csv")
vr=subset(read.csv("Data/vr.csv"),year==2014) #data from 2013/2014 only.

#1. format Germ
germ$germTot=apply(germ[,grep("Census",names(germ))],1,sum)
#sum up regrowth as a check - this need be less than germTot
germ$germRegrowthTot=apply(germ[,grep("Re.growth",names(germ))],1,sum) 
expect_equal(sum(germ$germRegrowthTot > germ$germTot,na.rm=T),0) 
germ$germFail=25-germ$germTot
names(germ)[1]="plot"

#2. Format tetrazolium data
## First, remove rows that contain no data 
rI=which(is.na(viab$Yes) & is.na(viab$No) & is.na(viab$Maybe) & is.na(viab$Brown))
viab=viab[-rI,]
viab$Brown[is.na(viab$Brown)]=0
viab$totS=viab$Yes + viab$No + viab$Maybe + viab$Brown  #Total seeds (for binomial model)
viab$yesMaybe=viab$Yes + viab$Maybe               #Conservative way to identify "viable seeds"
viab$fail=viab$totS-viab$Yes                      #"Fail" column for the glm() analyses.
viab$failMaybe=viab$totS-viab$yesMaybe
names(viab)[1]="plot" 

#3. Format plot data
#Test that all there is no "cnf" (could not find) data points 
expect_equal(sum(vr$F_flow_t1=="cnf",na.rm=T),0)
expect_equal(sum(vr$M_flow_t1=="cnf",na.rm=T),0)
#select and format data
plotD=unique(vr[,c("plot","F","M","F_flow_t1","M_flow_t1")])
plotD$M_flow_t1=as.numeric(as.character(plotD$M_flow_t1))
plotD$F_flow_t1=as.numeric(as.character(plotD$F_flow_t1))
#Calculate total number of flowers in each plot 
focal_f_flow=aggregate(flowN_t1 ~ plot, sum, data=subset(vr,sex=="f")) ; names(focal_f_flow)[2]="F_flow_f"
focal_m_flow=aggregate(flowN_t1 ~ plot, sum, data=subset(vr,sex=="m")) ; names(focal_m_flow)[2]="M_flow_f"
ch=merge(plotD,focal_f_flow,all=T) #Merge data 
ch=merge(ch,focal_m_flow,all=T)
ch$F_flow_f[is.na(ch$F_flow_f)]=0 #Substitute NAs with zeros
ch$M_flow_f[is.na(ch$M_flow_f)]=0
#IF #1. F or M are <=5, then #2. Substitute #_flow with #_flow_f  
ch[which(ch$F<=5),]$F_flow_t1 = ch[which(ch$F<=5),]$F_flow_f 
ch[which(ch$M<=5),]$M_flow_t1 = ch[which(ch$M<=5),]$M_flow_f 
#Calculate sex ratios
ch$totFlow=ch$M_flow_f + ch$F_flow_f #Tot number of flowers
ch$sr_f=ch$F_flow_f/ch$totFlow 
ch$sr_f2=ch$sr_f ^ 2

#Merge the three datasets
tmp=merge(viab,germ[,c("plot","F_ID","P_ID","germTot","germFail")],all=T)
tmp$focalI=matrix(unlist(strsplit(as.character(tmp$F_ID),"F-")),nrow(tmp),2,byrow=T)[,2]
tmp$focalI=paste0("f",as.numeric(tmp$focalI))
viabVr=merge(ch,tmp) #Final file

#Calculate viability/germination ratios
viabVr$tetra_ratio        <- viabVr$Yes / viabVr$totS 
viabVr$tetra_maybe_ratio  <- viabVr$yesMaybe / viabVr$totS
viabVr$germ_ratio         <- viabVr$germTot / (viabVr$germTot + viabVr$germFail)

##########################################################################################################################
#VIABILITY models. Uses absolute number of viable seeds (SeedN / fertRatio) for viability and germination data. 
##########################################################################################################################

#invlogit for graphs
invlogit<-function(x){exp(x)/(1+exp(x))}

#1. viable Seed Number according to tetrazolium assays  (Yes / fail) 
l_viab_tetr=list()
l_viab_tetr[[1]]=glm(cbind(Yes,fail) ~ sr_f,family="binomial", data=viabVr)
l_viab_tetr[[2]]=glm(cbind(Yes,fail) ~ sr_f2,family="binomial", data=viabVr)
l_viab_tetr[[3]]=glm(cbind(Yes,fail) ~ totFlow,family="binomial", data=viabVr)
l_viab_tetr[[4]]=glm(cbind(Yes,fail) ~ sr_f + sr_f2,family="binomial", data=viabVr)
l_viab_tetr[[5]]=glm(cbind(Yes,fail) ~ sr_f + totFlow,family="binomial", data=viabVr)
l_viab_tetr[[6]]=glm(cbind(Yes,fail) ~ sr_f + totFlow + sr_f2,family="binomial", data=viabVr)
l_viab_tetr[[7]]=glm(cbind(Yes,fail) ~ sr_f * totFlow,family="binomial", data=viabVr)
l_viab_tetr[[8]]=glm(cbind(Yes,fail) ~ sr_f * totFlow + sr_f2,family="binomial", data=viabVr)
AICtab(l_viab_tetr,weights=T) #mod 6, 71% support, mod 6 + 8 have 100% support


#2. viable Seed Number according to tetrazolium assays  (Yes / fail) 
# This is a **less conservative** way to identify viable seeds. 
#"Maybe" refers to seeds that were barely stained by tetrazolium.
l_tetr_maybe=list()
l_tetr_maybe[[1]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f,family="binomial", data=viabVr)
l_tetr_maybe[[2]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f2,family="binomial", data=viabVr)
l_tetr_maybe[[3]]=glm(cbind(yesMaybe,failMaybe) ~ totFlow,family="binomial", data=viabVr)
l_tetr_maybe[[4]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f + sr_f2,family="binomial", data=viabVr)
l_tetr_maybe[[5]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f + totFlow,family="binomial", data=viabVr)
l_tetr_maybe[[6]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f + totFlow + sr_f2,family="binomial", data=viabVr)
l_tetr_maybe[[7]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f * totFlow,family="binomial", data=viabVr)
l_tetr_maybe[[8]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f * totFlow + sr_f2,family="binomial", data=viabVr)
AICtab(l_tetr_maybe,weights=T) #mod 6, 63% support (mod 6 + 8 have 100%)


#3. viable Seed Number according to germination assays (germTot / germFail) 
l_viab_germ=list()
l_viab_germ[[1]]=glm(cbind(germTot,germFail) ~ sr_f,family="binomial", data=viabVr)
l_viab_germ[[2]]=glm(cbind(germTot,germFail) ~ sr_f2,family="binomial", data=viabVr)
l_viab_germ[[3]]=glm(cbind(germTot,germFail) ~ totFlow,family="binomial", data=viabVr)
l_viab_germ[[4]]=glm(cbind(germTot,germFail) ~ sr_f + sr_f2,family="binomial", data=viabVr)
l_viab_germ[[5]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow,family="binomial", data=viabVr)
l_viab_germ[[6]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow + sr_f2,family="binomial", data=viabVr)
l_viab_germ[[7]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow,family="binomial", data=viabVr)
l_viab_germ[[8]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow + sr_f2,family="binomial", data=viabVr)
AICtab(l_viab_germ,weights=T) #mod 6, 67% support (mod 6 + 8 have 100% support)


#Graph results
tiff("Results/VitalRates_Sept16/viability.tiff",unit="in",width=6.3,height=6,res=600,compression="lzw")

lower=quantile(viabVr$totFlow,prob=c(0.1))
upper=quantile(viabVr$totFlow,prob=c(0.9))

par(mfrow=c(2,2),mar=c(2.5,3,1.1,0.1),mgp=c(1.5,0.6,0))
titlePlace=0.2

plot(viabVr$sr_f,viabVr$tetra_ratio,pch=16,
     xlab="Sex ratio",ylab="Viability rate")
title("Tetrazolium data", line = titlePlace)
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
xSeq2=xSeq^2
yMeanLow=invlogit(coef(l_viab_tetr[[6]])[1] + coef(l_viab_tetr[[6]])[2]*xSeq + coef(l_viab_tetr[[6]])[4]*xSeq2 + 
                    coef(l_viab_tetr[[6]])[3]*lower)
yMeanHigh=invlogit(coef(l_viab_tetr[[6]])[1] + coef(l_viab_tetr[[6]])[2]*xSeq + coef(l_viab_tetr[[6]])[4]*xSeq2 + 
                     coef(l_viab_tetr[[6]])[3]*upper)
lines(xSeq,yMeanLow,lwd=2,col="blue")
lines(xSeq,yMeanHigh,lwd=2,col="red")

plot(viabVr$sr_f,viabVr$tetra_maybe_ratio,pch=16,
     xlab="Sex ratio",ylab="Viability rate")
title("Tetrazolium maybe", line = titlePlace)
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
xSeq2=xSeq^2
yMeanLow=invlogit(coef(l_tetr_maybe[[6]])[1] + coef(l_tetr_maybe[[6]])[2]*xSeq + coef(l_tetr_maybe[[6]])[4]*xSeq2 + 
                    coef(l_tetr_maybe[[6]])[3]*lower)
yMeanHigh=invlogit(coef(l_tetr_maybe[[6]])[1] + coef(l_tetr_maybe[[6]])[2]*xSeq + coef(l_tetr_maybe[[6]])[4]*xSeq2 + 
                     coef(l_tetr_maybe[[6]])[3]*upper)
lines(xSeq,yMeanLow,lwd=2,col="blue")
lines(xSeq,yMeanHigh,lwd=2,col="red")


plot(viabVr$sr_f,viabVr$germ_ratio,pch=16,
     xlab="Sex ratio",ylab="Germination rate")
title("Germination data", line = titlePlace)
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
xSeq2=xSeq^2
yMeanLow=invlogit(coef(l_viab_germ[[6]])[1] + coef(l_viab_germ[[6]])[2]*xSeq + coef(l_viab_germ[[6]])[4]*xSeq2 + 
                    coef(l_viab_germ[[6]])[3]*lower)
yMeanHigh=invlogit(coef(l_viab_germ[[6]])[1] + coef(l_viab_germ[[6]])[2]*xSeq + coef(l_viab_germ[[6]])[4]*xSeq2 + 
                     coef(l_viab_germ[[6]])[3]*upper)
lines(xSeq,yMeanLow,lwd=2,col="blue")
lines(xSeq,yMeanHigh,lwd=2,col="red")

plot(c(1:10),type="n",bty="n",xaxt="n",yaxt="n",ylab="",xlab="") #
text(0.3,10,"Tetrazolium data model:",pos=4)
text(0.3,9.3,expression("OSR * Density + OSR"^2),pos=4)

text(0.3,8,"Tetrazolium maybe model:",pos=4)
text(0.3,7.3,expression("OSR * Density + OSR"^2),pos=4)

text(0.3,6,"Germination model:",pos=4)
text(0.3,5.3,expression("OSR + Density + OSR"^2),pos=4)

legend(0.3,4.5,c("High density plot","Low density plot"),lwd=2,lty=1,
       col=c("blue","red"),bty="n")

dev.off()


###################################################################################
#DEBUG
###################################################################################

##Why do probabilities of fertilization decrease so much as we get close to 0 female contribution?

#Hypothesis to explain pattern: 
#1. Is the decrease in prob. of fertilization a consequence of sr_f=100% data points 
#   biasing the sr_f2 term of the model? I 'chopped' data set (viability_choppedData.tiff) to check if that's the case   
#2. Related to above: are there categortical differnces between "bins" of sex ratios (sr_f)?>
#3. Residuals of best glm() model: do they change with sr_f? (no)
#   if they did at low sr_f values, this would have suggested data points at sr_f=100% bias regression estimates 
#4. Does fitting data withy GAMs tell a different story? (yes it does, no decrease in probability towards low sr_f)

#1.ORIGINAL germination analyses--------------------------------------------------------
l_viab_germ=list()
l_viab_germ[[1]]=glm(cbind(germTot,germFail) ~ sr_f,family="binomial", data=viabVr)
l_viab_germ[[2]]=glm(cbind(germTot,germFail) ~ sr_f2,family="binomial", data=viabVr)
l_viab_germ[[3]]=glm(cbind(germTot,germFail) ~ totFlow,family="binomial", data=viabVr)
l_viab_germ[[4]]=glm(cbind(germTot,germFail) ~ sr_f + sr_f2,family="binomial", data=viabVr)
l_viab_germ[[5]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow,family="binomial", data=viabVr)
l_viab_germ[[6]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow + sr_f2,family="binomial", data=viabVr)
l_viab_germ[[7]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow,family="binomial", data=viabVr)
l_viab_germ[[8]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow + sr_f2,family="binomial", data=viabVr)
AICtab(l_viab_germ,weights=T) #mod 6, 67% support (mod 6 + 8 have 100% support)


##1. Is the decrease in prob. of fertilization a consequence of sr_f=100% data points biasing sr_f2 results?
tiff("Results/VitalRates_Sept16/viability_choppedData.tiff",unit="in",width=6.3,height=6,res=600,compression="lzw")

lower=quantile(viabVr$totFlow,prob=c(0.3))
upper=quantile(viabVr$totFlow,prob=c(0.7))

par(mfrow=c(3,2),mar=c(2.5,3,1.1,0.1),mgp=c(1.5,0.6,0))

#1.ORIGINAL germination analyses--------------------------------------------------------
l_viab_germ=list()
l_viab_germ[[1]]=glm(cbind(germTot,germFail) ~ sr_f,family="binomial", data=viabVr)
l_viab_germ[[2]]=glm(cbind(germTot,germFail) ~ sr_f2,family="binomial", data=viabVr)
l_viab_germ[[3]]=glm(cbind(germTot,germFail) ~ totFlow,family="binomial", data=viabVr)
l_viab_germ[[4]]=glm(cbind(germTot,germFail) ~ sr_f + sr_f2,family="binomial", data=viabVr)
l_viab_germ[[5]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow,family="binomial", data=viabVr)
l_viab_germ[[6]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow + sr_f2,family="binomial", data=viabVr)
l_viab_germ[[7]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow,family="binomial", data=viabVr)
l_viab_germ[[8]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow + sr_f2,family="binomial", data=viabVr)
AICtab(l_viab_germ,weights=T) #mod 6, 67% support (mod 6 + 8 have 100% support)

plotViab=viabVr
beta=coef(l_viab_germ[[6]])
plot(plotViab$sr_f,plotViab$germ_ratio,pch=16,xlab="Sex ratio",ylab="Germination rate")
title("All data", line = titlePlace)
xSeq=seq(min(plotViab$sr_f),max(plotViab$sr_f),length.out=100)
xSeq2=xSeq^2
yMeanLow=invlogit(beta[1] + beta[2]*xSeq + beta[3]*lower + 
                    beta[4]*xSeq2 + beta[5]*xSeq*lower)
yMeanHigh=invlogit(beta[1] + beta[2]*xSeq + beta[3]*upper + 
                     beta[4]*xSeq2 + beta[5]*xSeq*upper)
lines(xSeq,yMeanLow,lwd=2,col="blue")
lines(xSeq,yMeanHigh,lwd=2,col="red")


#2.Exclude less than 20% females--------------------------------------------------------
l_viab_germ=list()
l_viab_germ[[1]]=glm(cbind(germTot,germFail) ~ sr_f,family="binomial", data=subset(viabVr,sr_f>0.2))
l_viab_germ[[2]]=glm(cbind(germTot,germFail) ~ sr_f2,family="binomial", data=subset(viabVr,sr_f>0.2))
l_viab_germ[[3]]=glm(cbind(germTot,germFail) ~ totFlow,family="binomial", data=subset(viabVr,sr_f>0.2))
l_viab_germ[[4]]=glm(cbind(germTot,germFail) ~ sr_f + sr_f2,family="binomial", data=subset(viabVr,sr_f>0.2))
l_viab_germ[[5]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow,family="binomial", data=subset(viabVr,sr_f>0.2))
l_viab_germ[[6]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow + sr_f2,family="binomial", data=subset(viabVr,sr_f>0.2))
l_viab_germ[[7]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow,family="binomial", data=subset(viabVr,sr_f>0.2))
l_viab_germ[[8]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow + sr_f2,family="binomial", data=subset(viabVr,sr_f>0.2))
AICtab(l_viab_germ,weights=T) #mod 6, 67% support (mod 6 + 8 have 100% support)

#plotViab=subset(viabVr,sr_f>0.2)
plotViab=viabVr
beta=coef(l_viab_germ[[8]])
plot(plotViab$sr_f,plotViab$germ_ratio,pch=16,xlab="Sex ratio",ylab="Germination rate")
title("No <20% fem.", line = titlePlace)
xSeq=seq(min(plotViab$sr_f),max(plotViab$sr_f),length.out=100)
xSeq2=xSeq^2
yMeanLow=invlogit(beta[1] + beta[2]*xSeq + beta[3]*lower + 
                    beta[4]*xSeq2 + beta[5]*xSeq*lower)
yMeanHigh=invlogit(beta[1] + beta[2]*xSeq + beta[3]*upper + 
                     beta[4]*xSeq2 + beta[5]*xSeq*upper)
lines(xSeq,yMeanLow,lwd=2,col="blue")
lines(xSeq,yMeanHigh,lwd=2,col="red")
abline(v=0.2)

#3.Exclude less than 50% females--------------------------------------------------------
l_viab_germ=list()
l_viab_germ[[1]]=glm(cbind(germTot,germFail) ~ sr_f,family="binomial", data=subset(viabVr,sr_f>0.5))
l_viab_germ[[2]]=glm(cbind(germTot,germFail) ~ sr_f2,family="binomial", data=subset(viabVr,sr_f>0.5))
l_viab_germ[[3]]=glm(cbind(germTot,germFail) ~ totFlow,family="binomial", data=subset(viabVr,sr_f>0.5))
l_viab_germ[[4]]=glm(cbind(germTot,germFail) ~ sr_f + sr_f2,family="binomial", data=subset(viabVr,sr_f>0.5))
l_viab_germ[[5]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow,family="binomial", data=subset(viabVr,sr_f>0.5))
l_viab_germ[[6]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow + sr_f2,family="binomial", data=subset(viabVr,sr_f>0.5))
l_viab_germ[[7]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow,family="binomial", data=subset(viabVr,sr_f>0.5))
l_viab_germ[[8]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow + sr_f2,family="binomial", data=subset(viabVr,sr_f>0.5))
AICtab(l_viab_germ,weights=T) #mod 6, 67% support (mod 6 + 8 have 100% support)

#plotViab=subset(viabVr,sr_f>0.5)
plotViab=viabVr
beta=coef(l_viab_germ[[8]])
plot(plotViab$sr_f,plotViab$germ_ratio,pch=16,xlab="Sex ratio",ylab="Germination rate")
title("No <50% fem.", line = titlePlace)
xSeq=seq(min(plotViab$sr_f),max(plotViab$sr_f),length.out=100)
xSeq2=xSeq^2
yMeanLow=invlogit(beta[1] + beta[2]*xSeq + beta[3]*lower + 
                    beta[4]*xSeq2 + beta[5]*xSeq*lower)
yMeanHigh=invlogit(beta[1] + beta[2]*xSeq + beta[3]*upper + 
                     beta[4]*xSeq2 + beta[5]*xSeq*upper)
lines(xSeq,yMeanLow,lwd=2,col="blue")
lines(xSeq,yMeanHigh,lwd=2,col="red")
abline(v=0.5)


#4.Exclude less than 60% females--------------------------------------------------------


l_viab_germ=list()
l_viab_germ[[1]]=glm(cbind(germTot,germFail) ~ sr_f,family="binomial", data=subset(viabVr,sr_f>0.6))
l_viab_germ[[2]]=glm(cbind(germTot,germFail) ~ sr_f2,family="binomial", data=subset(viabVr,sr_f>0.6))
l_viab_germ[[3]]=glm(cbind(germTot,germFail) ~ totFlow,family="binomial", data=subset(viabVr,sr_f>0.6))
l_viab_germ[[4]]=glm(cbind(germTot,germFail) ~ sr_f + sr_f2,family="binomial", data=subset(viabVr,sr_f>0.6))
l_viab_germ[[5]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow,family="binomial", data=subset(viabVr,sr_f>0.6))
l_viab_germ[[6]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow + sr_f2,family="binomial", data=subset(viabVr,sr_f>0.6))
l_viab_germ[[7]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow,family="binomial", data=subset(viabVr,sr_f>0.6))
l_viab_germ[[8]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow + sr_f2,family="binomial", data=subset(viabVr,sr_f>0.6))
AICtab(l_viab_germ,weights=T) #mod 6, 67% support (mod 6 + 8 have 100% support)

#plotViab=subset(viabVr,sr_f>0.5)
plotViab=viabVr
beta=coef(l_viab_germ[[8]])
plot(plotViab$sr_f,plotViab$germ_ratio,pch=16,xlab="Sex ratio",ylab="Germination rate")
title("No <60% fem.", line = titlePlace)
xSeq=seq(min(plotViab$sr_f),max(plotViab$sr_f),length.out=100)
xSeq2=xSeq^2
yMeanLow=invlogit(beta[1] + beta[2]*xSeq + beta[3]*lower + 
                    beta[4]*xSeq2 + beta[5]*xSeq*lower)
yMeanHigh=invlogit(beta[1] + beta[2]*xSeq + beta[3]*upper + 
                     beta[4]*xSeq2 + beta[5]*xSeq*upper)
lines(xSeq,yMeanLow,lwd=2,col="blue")
lines(xSeq,yMeanHigh,lwd=2,col="red")
abline(v=0.6)


#5.ORIGINAL germination analyses--------------------------------------------------------
l_viab_germ=list()
l_viab_germ[[1]]=glm(cbind(germTot,germFail) ~ sr_f,family="binomial", data=subset(viabVr,sr_f>0.25 & sr_f<1))
l_viab_germ[[2]]=glm(cbind(germTot,germFail) ~ sr_f2,family="binomial", data=subset(viabVr,sr_f>0.25 & sr_f<1))
l_viab_germ[[3]]=glm(cbind(germTot,germFail) ~ totFlow,family="binomial", data=subset(viabVr,sr_f>0.25 & sr_f<1))
l_viab_germ[[4]]=glm(cbind(germTot,germFail) ~ sr_f + sr_f2,family="binomial", data=subset(viabVr,sr_f>0.25 & sr_f<1))
l_viab_germ[[5]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow,family="binomial", data=subset(viabVr,sr_f>0.25 & sr_f<1))
l_viab_germ[[6]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow + sr_f2,family="binomial", data=subset(viabVr,sr_f>0.25 & sr_f<1))
l_viab_germ[[7]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow,family="binomial", data=subset(viabVr,sr_f>0.25 & sr_f<1))
l_viab_germ[[8]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow + sr_f2,family="binomial", data=subset(viabVr,sr_f>0.25 & sr_f<1))
AICtab(l_viab_germ,weights=T) #mod 6, 67% support (mod 6 + 8 have 100% support)

#plotViab=subset(viabVr,sr_f>0.5)
plotViab=viabVr
beta=coef(l_viab_germ[[8]])
plot(plotViab$sr_f,plotViab$germ_ratio,pch=16,xlab="Sex ratio",ylab="Germination rate")
title("No <25 or 100% fem.", line = titlePlace)
xSeq=seq(min(plotViab$sr_f),max(plotViab$sr_f),length.out=100)
xSeq2=xSeq^2
yMeanLow=invlogit(beta[1] + beta[2]*xSeq + beta[3]*lower + 
                    beta[4]*xSeq2 + beta[5]*xSeq*lower)
yMeanHigh=invlogit(beta[1] + beta[2]*xSeq + beta[3]*upper + 
                     beta[4]*xSeq2 + beta[5]*xSeq*upper)
lines(xSeq,yMeanLow,lwd=2,col="blue")
lines(xSeq,yMeanHigh,lwd=2,col="red")
abline(v=0.25) ; abline(v=1) 


#6.Exclude female-only plots--------------------------------------------------------
l_viab_germ=list()
l_viab_germ[[1]]=glm(cbind(germTot,germFail) ~ sr_f,family="binomial", data=subset(viabVr,sr_f<1))
l_viab_germ[[2]]=glm(cbind(germTot,germFail) ~ sr_f2,family="binomial", data=subset(viabVr,sr_f<1))
l_viab_germ[[3]]=glm(cbind(germTot,germFail) ~ totFlow,family="binomial", data=subset(viabVr,sr_f<1))
l_viab_germ[[4]]=glm(cbind(germTot,germFail) ~ sr_f + sr_f2,family="binomial", data=subset(viabVr,sr_f<1))
l_viab_germ[[5]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow,family="binomial", data=subset(viabVr,sr_f<1))
l_viab_germ[[6]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow + sr_f2,family="binomial", data=subset(viabVr,sr_f<1))
l_viab_germ[[7]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow,family="binomial", data=subset(viabVr,sr_f<1))
l_viab_germ[[8]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow + sr_f2,family="binomial", data=subset(viabVr,sr_f<1))
AICtab(l_viab_germ,weights=T) #mod 6, 67% support (mod 6 + 8 have 100% support)

#plotViab=subset(viabVr,sr_f>0.5)
plotViab=viabVr
beta=coef(l_viab_germ[[8]])
plot(plotViab$sr_f,plotViab$germ_ratio,pch=16,xlab="Sex ratio",ylab="Germination rate")
title("No 100% fem.", line = titlePlace)
xSeq=seq(min(plotViab$sr_f),max(plotViab$sr_f),length.out=100)
xSeq2=xSeq^2
yMeanLow=invlogit(beta[1] + beta[2]*xSeq + beta[3]*lower + 
                    beta[4]*xSeq2 + beta[5]*xSeq*lower)
yMeanHigh=invlogit(beta[1] + beta[2]*xSeq + beta[3]*upper + 
                     beta[4]*xSeq2 + beta[5]*xSeq*upper)
lines(xSeq,yMeanLow,lwd=2,col="blue")
lines(xSeq,yMeanHigh,lwd=2,col="red")
abline(v=1)

dev.off()


#2.are there categortical differnces between "bins" of sex ratios (sr_f)?>
germChop=subset(viabVr,sr_f<=0.75)
germChop$sr_level=1 ; germChop$sr_level[germChop$sr_f>0.25 & germChop$sr_f<=0.5]=2 
germChop$sr_level[germChop$sr_f>0.5]=3 
germChop$sr_level=as.factor(germChop$sr_level)
m0=glm(cbind(germTot,germFail) ~ sr_level,family="binomial",data=germChop)


#3. Residuals of best glm() model: do they change with sr_f? Answer is: NO! 
germNaOmit=na.omit(select(viabVr,germTot,germFail,sr_f,totFlow,sr_f2))
m0=glm(cbind(germTot,germFail) ~ sr_f + totFlow + sr_f2,family="binomial", data=germNaOmit)
par(mfrow=c(1,1))
plot(germNaOmit$sr_f,resid(m0))
abline(h=0)


#4. GAM fits tell a different story: no decrease in probability of viability towards low sr_f
tiff("Results/VitalRates_Sept16/viability_GAM.tiff",unit="in",width=6.3,height=6,res=600,compression="lzw")

germGamData=na.omit(viabVr[,c('germTot','germFail',"sr_f")])
dat=data.frame(sr_f=seq(0,1,by=0.01))

par(mfrow=c(3,2),mar=c(2.5,3,1.1,0.1),mgp=c(1.5,0.6,0))
modGam=gam(cbind(germTot,germFail) ~ s(sr_f,k=3),family="binomial", data=germGamData)
plot(dat$sr_f,predict(modGam, newdata = dat,type="response"),type="l",main="Basis dimension: 3",xlab="Sex ratio (% of females)",
         ylab="P. of viability (logit scale)")
modGam=gam(cbind(germTot,germFail) ~ s(sr_f,k=4),family="binomial", data=germGamData)
plot(dat$sr_f,predict(modGam, newdata = dat,type="response"),type="l",main="Basis dimension: 4",xlab="Sex ratio (% of females)",
         ylab="P. of viability (logit scale)")
modGam=gam(cbind(germTot,germFail) ~ s(sr_f,k=5),family="binomial", data=germGamData)
plot(dat$sr_f,predict(modGam, newdata = dat,type="response"),type="l",main="Basis dimension: 5",xlab="Sex ratio (% of females)",
         ylab="P. of viability (logit scale)")
modGam=gam(cbind(germTot,germFail) ~ s(sr_f,k=6),family="binomial", data=germGamData)
plot(dat$sr_f,predict(modGam, newdata = dat,type="response"),type="l",main="Basis dimension: 6",xlab="Sex ratio (% of females)",
         ylab="P. of viability (logit scale)")
modGam=gam(cbind(germTot,germFail) ~ s(sr_f,k=7),family="binomial", data=germGamData)
plot(dat$sr_f,predict(modGam, newdata = dat,type="response"),type="l",main="Basis dimension: 7",xlab="Sex ratio (% of females)",
         ylab="P. of viability (logit scale)")
modGam=gam(cbind(germTot,germFail) ~ s(sr_f,k=8),family="binomial", data=germGamData)
plot(dat$sr_f,predict(modGam, newdata = dat,type="response"),type="l",main="Basis dimension: 8",xlab="Sex ratio (% of females)",
         ylab="P. of viability (logit scale)")

dev.off()
