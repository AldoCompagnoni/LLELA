##Model selection for FERTILITY, whereby: 
##1. Fertility == "Viable Seed Number / panicule length"
##2. Viable Seed Number == Seed Number * viabiliy based on:
##   a. conservative tetrazolium scoring 
##   b. Germination assays.
##   c. A conmbination of tetrazolium scoring and germination assays
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(lme4)
library(bbmle)

#Read in data and format###########################################################################################
vr=subset(read.csv("Data/vr.csv"),year==2014) #Data from spring 2014.
viab=read.csv("Data/Spring 2014/viability/viabilityData.csv")
germ=read.csv("Data/Spring 2014/viability/germination_data_LLELA.csv")
femPanicules=read.csv("Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")


#IMPORTANT NOTE about formatting tetrazolium and germination data:
#Here, I AGGREGATE successes/failures of tetrazolium and germination data at the individual level.
#Individuals that produced more than 1 panicule can have multiple entries. 
#I aggregated these records for these analyses.

#1. format Germ
germ$germTot=apply(germ[,grep("Census",names(germ))],1,sum)
#sum up regrowth as a check - this need be less than germTot
germ$germRegrowthTot=apply(germ[,grep("Re.growth",names(germ))],1,sum)
sum(germ$germRegrowthTot > germ$germTot,na.rm=T) #this should be 0
germ$germFail=25-germ$germTot
#Produce a focalI  ("focal Individual") column to identify individuals 
germ$focalI=paste0("f",matrix(unlist(strsplit(as.character(germ$F_ID),"F-")),nrow(germ),2,byrow=T)[,2])
names(germ)[1]="plot"
germData=aggregate(cbind(germTot,germFail) ~ plot + focalI, sum, data=germ)

#2. Format tetrazolium data
rI=which(is.na(viab$Yes) & is.na(viab$No) & is.na(viab$Maybe) & is.na(viab$Brown))
viab=viab[-rI,]
viab$Brown[is.na(viab$Brown)]=0
viab$totS=viab$Yes + viab$No + viab$Maybe + viab$Brown  #Total seeds (for binomial model)
viab$yesMaybe=viab$Yes + viab$Maybe               #Conservative way to identify "viable"\
viab$fail=viab$totS-viab$Yes #"Fail" column for the glm() analyses.
viab$failMaybe=viab$totS-viab$yesMaybe
names(viab)[1]="plot" 
#"format" focal IDs
viab$focalI=paste0("f",matrix(unlist(strsplit(as.character(viab$F_ID),"F-")),nrow(viab),2,byrow=T)[,2])
#Aggregate the (two!) data points that belong to two panicules of same individual 
viabData=aggregate(cbind(Yes,No,Maybe,Brown,totS,yesMaybe,fail,failMaybe) ~ plot + focalI, 
                   sum, data=viab)

#3. Format seed weight data
names(femPanicules)[1]="plot"
femPanicules$focalI=paste("f",femPanicules$IndividualN,sep="")
names(femPanicules)[c(1,5,7)]=c("plot","panicule_Length_cm","seed_weight_mg")
seedWeights=aggregate(cbind(seed_weight_mg,SeedN,panicule_Length_cm) ~ plot + focalI, sum, data=femPanicules)
seedWeights$meanSeedWeight=seedWeights$seed_weight_mg / seedWeights$SeedN 


#4. Format plot data
plotD=unique(vr[,c("plot","F","M","F_flow_t1","M_flow_t1")])
plotD$F_flow_t1=as.numeric(as.character(plotD$F_flow_t1))
plotD$M_flow_t1=as.numeric(as.character(plotD$M_flow_t1))
#Calculate number of flowers on focal individuals 
focal_f_flow=aggregate(flowN_t1 ~ plot, sum, data=subset(vr,sex=="f")) ; names(focal_f_flow)[2]="F_flow_f"
focal_m_flow=aggregate(flowN_t1 ~ plot, sum, data=subset(vr,sex=="m")) ; names(focal_m_flow)[2]="M_flow_f"
ch=merge(plotD,focal_f_flow,all=T) #Merge data 
ch=merge(ch,focal_m_flow,all=T)
ch$F_flow_f[is.na(ch$F_flow_f)]=0 #Substitute NAs with zeros
ch$M_flow_f[is.na(ch$M_flow_f)]=0
#IF: 1. the absolute number of planted females ("F") or males ("M")
#are <=5, then 
#2. Substitute F_flow/M_flow with F_flow_f/M_flow_f.
#I do this because there is mismatch between number of flowers counted at the plot level
#And numbers counted at the level of "census individuals". Counts for census individuals are more precise 
ch[which(ch$F<=5),]$F_flow_t1 = ch[which(ch$F<=5),]$F_flow_f 
ch[which(ch$M<=5),]$M_flow_t1 = ch[which(ch$M<=5),]$M_flow_f 
#Calculate sex ratios
ch$totFlow=ch$M_flow_t1 + ch$F_flow_t1 #Tot number of flowers
ch$sr_f=ch$F_flow_t1/ch$totFlow 
ch$sr_f2=ch$sr_f ^ 2

#Merge three data sets#############################################################
viabFert=merge(seedWeights,viabData[,-4],all=T)
tmp=merge(ch,viabFert,all.y=T)
viabVr=merge(tmp,germ[,c("plot","focalI","germFail","germTot")])

#Calculate total fertility
viabVr$fertRatio_tetra  <- viabVr$Yes / viabVr$totS 
viabVr$fertRatio_germ   <- viabVr$germTot / (viabVr$germTot + viabVr$germFail)
viabVr$fertRatio_comb   <- (viabVr$fertRatio_tetra + viabVr$fertRatio_germ) / 2

viabVr$fertSeedN_tetra  <- viabVr$SeedN * viabVr$fertRatio_tetra
viabVr$fertSeedN_germ   <- viabVr$SeedN * viabVr$fertRatio_germ
viabVr$fertSeedN_comb   <- viabVr$SeedN * viabVr$fertRatio_comb

viabVr$fert_tetra   <- viabVr$fertSeedN_tetra / viabVr$panicule_Length_cm
viabVr$fert_germ    <- viabVr$fertSeedN_germ / viabVr$panicule_Length_cm
viabVr$fert_comb    <- viabVr$fertSeedN_comb / viabVr$panicule_Length_cm


##########################################################################################################################
#Models#############################################################################################
##########################################################################################################################

#Fertility models. Uses absolute number of viable seeds (SeedN / fertRatio) for viability and germination data. 
#1. viableSeedN_tetra   
#2. viableSeedN_germ
#3. viableSeedN_comb

#1. viable Seed Number/panicule length according to tetrazolium assays
l_tetr=list()
l_tetr[[1]]=lm(fert_tetra ~ sr_f,data=viabVr)
l_tetr[[2]]=lm(fert_tetra ~ sr_f2,data=viabVr)
l_tetr[[3]]=lm(fert_tetra ~ totFlow,data=viabVr)
l_tetr[[4]]=lm(fert_tetra ~ sr_f + sr_f2,data=viabVr)
l_tetr[[5]]=lm(fert_tetra ~ sr_f + totFlow,data=viabVr)
l_tetr[[6]]=lm(fert_tetra ~ sr_f + totFlow + sr_f2,data=viabVr)
l_tetr[[7]]=lm(fert_tetra ~ sr_f * totFlow,data=viabVr)
l_tetr[[8]]=lm(fert_tetra ~ sr_f * totFlow + sr_f2,data=viabVr)

#best is model 4 + 2 have 75% support
AICtab(l_tetr,weights=T)


#2. viable Seed Number/panicule length according to germination assays
l_germ=list()
l_germ[[1]]=lm(fert_germ ~ sr_f,data=viabVr)
l_germ[[2]]=lm(fert_germ ~ sr_f2,data=viabVr)
l_germ[[3]]=lm(fert_germ ~ totFlow,data=viabVr)
l_germ[[4]]=lm(fert_germ ~ sr_f + sr_f2,data=viabVr)
l_germ[[5]]=lm(fert_germ ~ sr_f + totFlow,data=viabVr)
l_germ[[6]]=lm(fert_germ ~ sr_f + totFlow + sr_f2,data=viabVr)
l_germ[[7]]=lm(fert_germ ~ sr_f * totFlow,data=viabVr)
l_germ[[8]]=lm(fert_germ ~ sr_f * totFlow + sr_f2,data=viabVr)

#model 2 + 4 have 73% support
AICtab(l_germ,weights=T)


#3. viable Seed Number/panicule length according to a combination of
##  tetrazolium and germination assays
l_comb=list()
l_comb[[1]]=lm(fert_comb ~ sr_f,data=viabVr)
l_comb[[2]]=lm(fert_comb ~ sr_f2,data=viabVr)
l_comb[[3]]=lm(fert_comb ~ totFlow,data=viabVr)
l_comb[[4]]=lm(fert_comb ~ sr_f + sr_f2,data=viabVr)
l_comb[[5]]=lm(fert_comb ~ sr_f + totFlow,data=viabVr)
l_comb[[6]]=lm(fert_comb ~ sr_f + totFlow + sr_f2,data=viabVr)
l_comb[[7]]=lm(fert_comb ~ sr_f * totFlow,data=viabVr)
l_comb[[8]]=lm(fert_comb ~ sr_f * totFlow + sr_f2,data=viabVr)

#model 4 + 2 have 75% support
AICtab(l_comb,weights=T)


#Prduce graphs of linear fits
tiff("Results/fertility/fertility.tiff",unit="in",width=6.3,height=5,res=600,compression="lzw")

par(mfrow=c(2,2),mar=c(2.5,3,1.1,0.1),mgp=c(1.5,0.6,0))
titlePlace=0.2

#Tetrazolium data
plot(viabVr$sr_f,viabVr$fert_tetra,pch=16,ylim=c(0,70),
     xlab="Sex ratio",ylab="viable seed N / panicule length")
title("Tetrazolium data", line = titlePlace)
text(0.7,68,expression("Sex ratio + Sex ratio"^2))
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
xSeq2=xSeq^2
yMean=coef(l_tetr[[4]])[1] + coef(l_tetr[[4]])[2]*xSeq + coef(l_tetr[[4]])[3]*xSeq2
lines(xSeq,yMean,lwd=2)

#Germination data
plot(viabVr$sr_f,viabVr$fert_germ,pch=16,ylim=c(0,70),
     xlab="Sex ratio",ylab="viable seed N / panicule length")
title("Germination data", line = titlePlace)
text(0.7,68,expression("Sex ratio"^2))
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
xSeq2=xSeq^2
yMean=coef(l_germ[[2]])[1] + coef(l_germ[[2]])[2]*xSeq2
lines(xSeq,yMean,lwd=2)

#Combined data
plot(viabVr$sr_f,viabVr$fert_comb,pch=16,ylim=c(0,70),
     xlab="Sex ratio",ylab="viable seed N / panicule length")
title("Combined data", line = titlePlace)
text(0.7,68,expression("Sex ratio + Sex ratio"^2))
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
xSeq2=xSeq^2
yMean=coef(l_comb[[4]])[1] + coef(l_comb[[4]])[2]*xSeq + coef(l_comb[[4]])[3]*xSeq2
lines(xSeq,yMean,lwd=2)

dev.off()
