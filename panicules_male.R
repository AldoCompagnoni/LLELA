setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) #For AIC weights, because I'm lazy!
library(lme4)
library(dplyr)

source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")

# read in data --------------------------------------------------------------
d             <- read.csv("Data/vr.csv")
malPanicules  <- read.csv("Data/Spring 2014/maleCounts/malePaniculesSpring2014.csv")


# format focal individual for male data
tmp=as.numeric(matrix(unlist(strsplit(as.character(malPanicules$Individual),"[A-Z]")),
                      nrow(malPanicules),2,byrow=T)[,2])
malPanicules$focalI=paste("m",tmp,sep="")
names(malPanicules)[1]="plot"


# MALE panicule lengths -------------------------------------------------------------------------------------
d14=subset(d,year==2014)
d14=subset(d14,surv_t1!=0)
mal_i_pan_data <- merge(d14,malPanicules)
mal_i_pan_data$pan_area <- mal_i_pan_data$panicule_Length_cm *  mal_i_pan_data$panicule_Width_cm * pi

plMod=list()
plMod[[1]]=lmer(pan_area ~ log_l_t0 + (1 | plot), data=mal_i_pan_data)
plMod[[2]]=lmer(pan_area ~ log_l_t0 + TotDensity + (1 | plot),data=mal_i_pan_data)
plMod[[3]]=lmer(pan_area ~ log_l_t0 + sr+ (1 | plot),data=mal_i_pan_data)
plMod[[4]]=lmer(pan_area ~ log_l_t0 + sr + TotDensity + (1 | plot),data=mal_i_pan_data)
plMod[[5]]=lmer(pan_area ~ log_l_t0 + sr * TotDensity + (1 | plot),data=mal_i_pan_data)

# Model average
m_i_pl_mod_select <- AICtab(plMod,weights=T)
m_i_pl_avg        <- model_avg(m_i_pl_mod_select, plMod)
write.csv(m_i_pl_avg, "Results/VitalRates_3/male_panicules.csv", row.names = F)

# Graph ---------------------------------------------------------------------------

tiff("Results/VitalRates_3/male_repr_alloc.tiff",unit="in",width=3.5,height=3.5,res=600,compression="lzw")

#Start plotting
par(mfcol=c(1,1),mar=c(2.8,3,1,0.1),mgp=c(1.4,0.5,0))

#best model
plot(mal_i_pan_data$TotDensity,mal_i_pan_data$pan_area,pch=1, 
     cex = mal_i_pan_data$sr * 1.5,
     ylab="Panicule area",xlab="Planting density")
size  = mean(mal_i_pan_data$log_l_t0)
beta  = m_i_pl_avg$avg
xSeq  = seq(1,48,by =1)
y_l = beta[1] + beta[2]*size + beta[3]*0.1 + beta[4]*xSeq + beta[5]*xSeq*0.1
y_h = beta[1] + beta[2]*size + beta[3]*0.9 + beta[4]*xSeq + beta[5]*xSeq*0.9

lines(xSeq,y_l,lwd=2,lty=2)
lines(xSeq,y_h,lwd=2,lty=1)

dev.off()

