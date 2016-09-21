# Model selection for seed VIABILITY
# Data: tetrazolium scoring and germination data 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(lme4) ; library(bbmle)

# Read in data and format=====================================================================================
germ=read.csv("Data/Spring 2014/viability/germination_data_LLELA.csv")
viab=read.csv("Data/Spring 2014/viability/viabilityData.csv")
vr=subset(read.csv("Data/vr.csv"),year==2014) #data from 2013/2014 only.

# 1. format Germ
germ$germTot=apply(germ[,grep("Census",names(germ))],1,sum)
# sum up regrowth as a check - this need be less than germTot
germ$germRegrowthTot=apply(germ[,grep("Re.growth",names(germ))],1,sum) 
expect_equal(sum(germ$germRegrowthTot > germ$germTot,na.rm=T),0) 
germ$germFail=25-germ$germTot
names(germ)[1]="plot"

# 2. Format tetrazolium data
# First, remove rows that contain no data 
rI=which(is.na(viab$Yes) & is.na(viab$No) & is.na(viab$Maybe) & is.na(viab$Brown))
viab=viab[-rI,]
viab$Brown[is.na(viab$Brown)]=0
viab$totS=viab$Yes + viab$No + viab$Maybe + viab$Brown  #Total seeds (for binomial model)
viab$yesMaybe=viab$Yes + viab$Maybe               #Conservative way to identify "viable seeds"
viab$fail=viab$totS-viab$Yes                      #"Fail" column for the glm() analyses.
viab$failMaybe=viab$totS-viab$yesMaybe
names(viab)[1]="plot" 

# 3. Format plot data
# Test that all there is no "cnf" (could not find) data points 
expect_equal(sum(vr$F_flow_t1=="cnf",na.rm=T),0)
expect_equal(sum(vr$M_flow_t1=="cnf",na.rm=T),0)
#select and format data
plotD=unique(vr[,c("plot","F","M","F_flow_t1","M_flow_t1")])
plotD$M_flow_t1=as.numeric(as.character(plotD$M_flow_t1))
plotD$F_flow_t1=as.numeric(as.character(plotD$F_flow_t1))
# Calculate total number of flowers in each plot 
focal_f_flow=aggregate(flowN_t1 ~ plot, sum, data=subset(vr,sex=="f")) ; names(focal_f_flow)[2]="F_flow_f"
focal_m_flow=aggregate(flowN_t1 ~ plot, sum, data=subset(vr,sex=="m")) ; names(focal_m_flow)[2]="M_flow_f"
ch=merge(plotD,focal_f_flow,all=T) #Merge data 
ch=merge(ch,focal_m_flow,all=T)
ch$F_flow_f[is.na(ch$F_flow_f)]=0 #Substitute NAs with zeros
ch$M_flow_f[is.na(ch$M_flow_f)]=0
# IF #1. F or M are <=5, then #2. Substitute #_flow with #_flow_f  
ch[which(ch$F<=5),]$F_flow_t1 = ch[which(ch$F<=5),]$F_flow_f 
ch[which(ch$M<=5),]$M_flow_t1 = ch[which(ch$M<=5),]$M_flow_f 
# Calculate sex ratios
ch$totFlow=ch$M_flow_f + ch$F_flow_f #Tot number of flowers
ch$sr_f=ch$F_flow_f/ch$totFlow 
ch$sr_f2=ch$sr_f ^ 2

# Merge the three datasets
tmp=merge(viab,germ[,c("plot","F_ID","P_ID","germTot","germFail")],all=T)
tmp$focalI=matrix(unlist(strsplit(as.character(tmp$F_ID),"F-")),nrow(tmp),2,byrow=T)[,2]
tmp$focalI=paste0("f",as.numeric(tmp$focalI))
viabVr=merge(ch,tmp) #Final file

# Calculate viability/germination ratios
viabVr$tetra_ratio        <- viabVr$Yes / viabVr$totS 
viabVr$tetra_maybe_ratio  <- viabVr$yesMaybe / viabVr$totS
viabVr$germ_ratio         <- viabVr$germTot / (viabVr$germTot + viabVr$germFail)


# Viability model selection - NO QUADRATIC TERM=============================================

# invlogit for graphs
invlogit<-function(x){exp(x)/(1+exp(x))}

# Too conservative scoring, we chose to ignore it.
# 1. viable Seed Number according to tetrazolium assays  (Yes / fail)
#l_viab_tetr=list()
#l_viab_tetr[[1]]=glm(cbind(Yes,fail) ~ sr_f,family="binomial", data=viabVr)
#l_viab_tetr[[2]]=glm(cbind(Yes,fail) ~ totFlow,family="binomial", data=viabVr)
#l_viab_tetr[[3]]=glm(cbind(Yes,fail) ~ sr_f + totFlow,family="binomial", data=viabVr)
#l_viab_tetr[[4]]=glm(cbind(Yes,fail) ~ sr_f * totFlow,family="binomial", data=viabVr)
#AICtab(l_viab_tetr,weights=T) #mod 6, 71% support, mod 6 + 8 have 100% support


# 2. viable Seed Number according to tetrazolium assays  (Yes / fail) 
# This is the "standard" scoring: consider as viable any seed stained by tetrazolium
l_tetr_maybe=list()
l_tetr_maybe[[1]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f,family="binomial", data=viabVr)
l_tetr_maybe[[2]]=glm(cbind(yesMaybe,failMaybe) ~ totFlow,family="binomial", data=viabVr)
l_tetr_maybe[[3]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f + totFlow,family="binomial", data=viabVr)
l_tetr_maybe[[4]]=glm(cbind(yesMaybe,failMaybe) ~ sr_f * totFlow,family="binomial", data=viabVr)
AICtab(l_tetr_maybe,weights=T) #mod 6, 63% support (mod 6 + 8 have 100%)


#3. viable Seed Number according to germination assays (germTot / germFail)---------------------
l_viab_germ=list()
l_viab_germ[[1]]=glm(cbind(germTot,germFail) ~ sr_f,family="binomial", data=viabVr)
l_viab_germ[[2]]=glm(cbind(germTot,germFail) ~ totFlow,family="binomial", data=viabVr)
l_viab_germ[[3]]=glm(cbind(germTot,germFail) ~ sr_f + totFlow,family="binomial", data=viabVr)
l_viab_germ[[4]]=glm(cbind(germTot,germFail) ~ sr_f * totFlow,family="binomial", data=viabVr)
AICtab(l_viab_germ,weights=T) #mod 6, 67% support (mod 6 + 8 have 100% support)


#Graph results
tiff("Results/VitalRates_Sept16/viability_noQuadratic.tiff",unit="in",
     width=6.3,height=4,res=600,compression="lzw")

lower=quantile(viabVr$totFlow,prob=c(0.1))
upper=quantile(viabVr$totFlow,prob=c(0.9))

par(mfrow=c(1,2),mar=c(2.5,2.5,1.1,0.1),mgp=c(1.5,0.6,0),
    oma=c(0,0,3,0),xpd=NA)
titlePlace=0.2

# Tetrazolium data
plot(viabVr$sr_f,viabVr$tetra_maybe_ratio,pch=16,
     xlab="Sex ratio",ylab="Seed viability rate")
title("Tetrazolium data", line = titlePlace)
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
beta=coef(l_tetr_maybe[[4]])
yMeanLow=invlogit(beta[1] + beta[2]*xSeq + beta[3]*lower + beta[4]*xSeq*lower)
yMeanHigh=invlogit(beta[1] + beta[2]*xSeq + beta[3]*upper + beta[4]*xSeq*upper)
lines(xSeq,yMeanLow,lwd=2,col="blue")
lines(xSeq,yMeanHigh,lwd=2,col="red")

# Legend
text(0,1.3,"Tetrazolium model:  Sex ratio X Density",pos=4)
text(0,1.2,"Germination model: Sex ratio + Density",pos=4)

# Germination data
plot(viabVr$sr_f,viabVr$germ_ratio,pch=16,
     xlab="Sex ratio",ylab="Seed germination rate")
title("Germination data", line = titlePlace)
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
beta=coef(l_viab_germ[[3]])
yMeanLow=invlogit(beta[1] + beta[2]*xSeq + beta[3]*lower)
yMeanHigh=invlogit(beta[1] + beta[2]*xSeq + beta[3]*upper)
lines(xSeq,yMeanLow,lwd=2,col="blue")
lines(xSeq,yMeanHigh,lwd=2,col="red")

legend(0.18,1.26,c("High density plot","Low density plot"),lwd=2,lty=1,
       col=c("blue","red"),bty="n")

dev.off()
