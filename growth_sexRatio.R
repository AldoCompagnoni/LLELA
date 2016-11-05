# Sex competition predictor: SEX RATIO (female individuals/total individuals)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)    
library(glmmADMB) #Fit models with a Negative Binomial
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")


# read in data
d=read.csv("Data/vr.csv")
# remove dead individuals (this is a GROWTH model!)
d=subset(d, surv_t1 != 0)
# logtransform leaf numbers
d$plot=as.factor(d$plot) # glmmadmb wants plot as a factor

# only use year two
tmp15=subset(d,year==2015)
tmp15$new_t1[tmp15$new_t1=="SKIPPED"]=NA
tmp15$new_t1[tmp15$new_t1=="cnf"]=NA
tmp15$new_t1[tmp15$new_t1==""]=NA
tmp15$new_t1=as.numeric(as.character(tmp15$new_t1))
d15=na.omit(tmp15[,c("l_t1","log_l_t0","plot","sex","new_t1","sr","TotDensity")])


# Density as "new tillers"-------------------------------------------------------------------------

# Sex ratio as original sex ratio
lMod=list() # lMod stands for "leaf model" (density is quantified by N. of leaves)
lMod[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=d15,family="nbinom2")
lMod[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=d15,family="nbinom2")
lMod[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=d15,family="nbinom2")
# Target fitness + effect = tot density 
lMod[[4]]=glmmadmb(l_t1 ~ log_l_t0 + new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + new_t1 + (1 | plot),data=d15,family="nbinom2")
# Target fitness + effect = tot density + response=sex 
lMod[[7]]=glmmadmb(l_t1 ~ log_l_t0 + new_t1 + sex:new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[8]]=glmmadmb(l_t1 ~ log_l_t0 + sex + new_t1 + new_t1:sex + (1 | plot),data=d15,family="nbinom2")
lMod[[9]]=glmmadmb(l_t1 ~ log_l_t0 * sex + new_t1 + new_t1:sex + (1 | plot),data=d15,family="nbinom2")
# Target fitness + effect = sex 
lMod[[10]]=glmmadmb(l_t1 ~ log_l_t0 + sr + new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[11]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + new_t1 + (1 | plot),data=d15,family="nbinom2")
lMod[[12]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + new_t1 + (1 | plot),data=d15,family="nbinom2")
# Target fitness + effect = sex + response=sex 
lMod[[13]]=glmmadmb(l_t1 ~ log_l_t0 + sr + new_t1 + sr:sex + new_t1:sex +(1 | plot),data=d15,family="nbinom2")
lMod[[14]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + new_t1 + sr:sex + new_t1:sex + (1 | plot),data=d15,family="nbinom2")
lMod[[15]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + new_t1 + sr:sex + new_t1:sex  + (1 | plot),data=d15,family="nbinom2")

# Model average
grow_select <- AICtab(lMod, weights = T)
grow_avg    <- model_avg(grow_select, lMod)
write.csv(grow_avg, "Results/VitalRates_3/growth_best.csv", row.names = F)


# Planting density -------------------------------------------------------------------------

# Sex ratio as original sex ratio
lMod=list() # lMod stands for "leaf model" (density is quantified by N. of leaves)
lMod[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=d15,family="nbinom2")
lMod[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=d15,family="nbinom2")
lMod[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=d15,family="nbinom2")
# Target fitness + effect = tot density 
lMod[[4]]=glmmadmb(l_t1 ~ log_l_t0 + TotDensity + (1 | plot),data=d15,family="nbinom2")
lMod[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + TotDensity + (1 | plot),data=d15,family="nbinom2")
lMod[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + TotDensity + (1 | plot),data=d15,family="nbinom2")
# Target fitness + effect = tot density + response=sex 
lMod[[7]]=glmmadmb(l_t1 ~ log_l_t0 + TotDensity + sex:TotDensity + (1 | plot),data=d15,family="nbinom2")
lMod[[8]]=glmmadmb(l_t1 ~ log_l_t0 + sex + TotDensity + TotDensity:sex + (1 | plot),data=d15,family="nbinom2")
lMod[[9]]=glmmadmb(l_t1 ~ log_l_t0 * sex + TotDensity + TotDensity:sex + (1 | plot),data=d15,family="nbinom2")
# Target fitness + effect = sex 
lMod[[10]]=glmmadmb(l_t1 ~ log_l_t0 + sr + TotDensity + (1 | plot),data=d15,family="nbinom2")
lMod[[11]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + TotDensity + (1 | plot),data=d15,family="nbinom2")
lMod[[12]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + TotDensity + (1 | plot),data=d15,family="nbinom2")
# Target fitness + effect = sex + response=sex 
lMod[[13]]=glmmadmb(l_t1 ~ log_l_t0 + sr + TotDensity + sr:sex + TotDensity:sex +(1 | plot),data=d15,family="nbinom2")
lMod[[14]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + TotDensity + sr:sex + TotDensity:sex + (1 | plot),data=d15,family="nbinom2")
lMod[[15]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + TotDensity + sr:sex + TotDensity:sex  + (1 | plot),data=d15,family="nbinom2")

# Model average
grow_select <- AICtab(lMod, weights = T)
grow_avgN   <- model_avg(grow_select, lMod)
write.csv(grow_avgN, "Results/VitalRates_3/growthN_best.csv", row.names = F)



# GRAPH ------------------------------------------------------------------------------------------------

# Set up colors for plots
# 2015
d15$col=as.integer(d15$sex)
d15$col=as.character(factor(d15$col,labels=c("blue","red")))
d15$symb=as.integer(as.character(factor(d15$col,labels=c("17","16"))))

# 2015 ---------------------------------------------------------------------------------------

tiff("Results/VitalRates_3/growth_raw.tiff",unit="in",width=6.3,height=3.15,res=600,compression="lzw")

par(mfrow=c(1,2),mar=c(2.8,2.5,0.7,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0))

plot(d15$new_t1,d15$l_t1, pch=16, ylab="Individuals number of leaves (2015)",
     xlab="Total number of new tillers",col=d15$col, ylim = c(0,82))
beta <- grow_avg[,c("predictor","avg")]$avg
xSeq <- seq(0,224,by=1 )
size <- mean( d15$log_l_t0 )
dens <- 1
y_m <- exp( beta[1] + beta[2]*size + beta[3]*xSeq + beta[4] + 
            beta[5]*0.5 + beta[6]*size + beta[7]*xSeq + beta[8]*0.5)
y_f <- exp( beta[1] + beta[2]*size + beta[3]*xSeq + beta[5]*0.5 )
lines(xSeq,y_f,lwd=1.5,lty=1,col="blue")
lines(xSeq,y_m,lwd=1.5,lty=1,col="red")

legend(-10,88,c("Males","Females"),
       lty=1,lwd=1.5,col=c("red","blue"),bty="n")

# Density
plot(d15$TotDensity,d15$l_t1, pch=16, ylab="Individuals number of leaves (2015)",
     xlab="Planting density",col=d15$col, ylim = c(0,82))
beta <- grow_avgN[,c("predictor","avg")]$avg
xSeq <- seq(0,48,by=1)
size <- mean(d15$log_l_t0)
y_m <- exp( beta[1] + beta[2]*size + beta[3] + beta[4]*0.5 + 
              beta[5]*xSeq + beta[6]*size + beta[7]*0.5 + beta[8]*xSeq )
y_f <- exp( beta[1] + beta[2]*size + beta[3] + beta[4]*0.5 + beta[5]*xSeq )
lines(xSeq,y_f,lwd=1.5,lty=1,col="blue")
lines(xSeq,y_m,lwd=1.5,lty=1,col="red")

dev.off()

