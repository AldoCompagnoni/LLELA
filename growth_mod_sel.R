# Sex competition predictor: SEX RATIO (female individuals/total individuals)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)    
library(glmmADMB) #Fit models with a Negative Binomial
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")


# load and format data ------------------------------------------------------------------
d=read.csv("Data/vr.csv")
# remove dead individuals (this is a GROWTH model!)
d=subset(d, surv_t1 != 0)
# logtransform leaf numbers
d$plot=as.factor(d$plot) # glmmadmb wants plot as a factor

# Year one
d14 <- subset(d, year==2014)

# Year two
tmp15=subset(d,year==2015)
tmp15$new_t1[tmp15$new_t1=="SKIPPED"]=NA
tmp15$new_t1[tmp15$new_t1=="cnf"]=NA
tmp15$new_t1[tmp15$new_t1==""]=NA
tmp15$new_t1=as.numeric(as.character(tmp15$new_t1))
d15=na.omit(tmp15[,c("l_t1","log_l_t0","plot","focalI","sex","new_t1","sr","TotDensity")])

#Merge
tmp     <- select(d14, plot, new_t1)
tmp     <- mutate(tmp,new_t1 = as.character(as.numeric(new_t1)))   
tmp     <- mutate(tmp, new_t0 = new_t1, new_t1 = NULL)
d14_15  <- merge(tmp, d15)


# 2014 --------------------------------------------------------------------------------

# Planting density
m14=list() # lMod stands for "leaf model" (density is quantified by N. of leaves)
m14[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=d14,family="nbinom2")
m14[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=d14,family="nbinom2")
m14[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=d14,family="nbinom2")
# Target fitness + effect = tot density 
m14[[4]]=glmmadmb(l_t1 ~ log_l_t0 + TotDensity + (1 | plot),data=d14,family="nbinom2")
m14[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + TotDensity + (1 | plot),data=d14,family="nbinom2")
m14[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + TotDensity + (1 | plot),data=d14,family="nbinom2")
# Target fitness + effect = tot density + response=sex 
m14[[7]]=glmmadmb(l_t1 ~ log_l_t0 + TotDensity + sex:TotDensity + (1 | plot),data=d14,family="nbinom2")
m14[[8]]=glmmadmb(l_t1 ~ log_l_t0 + sex + TotDensity + TotDensity:sex + (1 | plot),data=d14,family="nbinom2")
m14[[9]]=glmmadmb(l_t1 ~ log_l_t0 * sex + TotDensity + TotDensity:sex + (1 | plot),data=d14,family="nbinom2")
# Target fitness + effect = sex 
m14[[10]]=glmmadmb(l_t1 ~ log_l_t0 + sr + TotDensity + (1 | plot),data=d14,family="nbinom2")
m14[[11]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + TotDensity + (1 | plot),data=d14,family="nbinom2")
m14[[12]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + TotDensity + (1 | plot),data=d14,family="nbinom2")
# Target fitness + effect = sex + response=sex 
m14[[13]]=glmmadmb(l_t1 ~ log_l_t0 + sr + TotDensity + sr:sex + TotDensity:sex +(1 | plot),data=d14,family="nbinom2")
m14[[14]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + TotDensity + sr:sex + TotDensity:sex + (1 | plot),data=d14,family="nbinom2")
m14[[15]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + TotDensity + sr:sex + TotDensity:sex  + (1 | plot),data=d14,family="nbinom2")

gr14_mod_sel  <- AICtab(m14, weights = T)
gr14_avg      <- model_avg(gr14_mod_sel, m14)
write.csv(gr14_avg, "Results/VitalRates_3/growth_14N_best.csv", row.names = F)


# N Leaves densities
m14L=list() # lMod stands for "leaf model" (density is quantified by N. of leaves)
m14L[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=d14,family="nbinom2")
m14L[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=d14,family="nbinom2")
m14L[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=d14,family="nbinom2")
# Target fitness + effect = tot density 
m14L[[4]]=glmmadmb(l_t1 ~ log_l_t0 + c_t0 + (1 | plot),data=d14,family="nbinom2")
m14L[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + c_t0 + (1 | plot),data=d14,family="nbinom2")
m14L[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + c_t0 + (1 | plot),data=d14,family="nbinom2")
# Target fitness + effect = tot density + response=sex 
m14L[[7]]=glmmadmb(l_t1 ~ log_l_t0 + c_t0 + sex:c_t0 + (1 | plot),data=d14,family="nbinom2")
m14L[[8]]=glmmadmb(l_t1 ~ log_l_t0 + sex + c_t0 + c_t0:sex + (1 | plot),data=d14,family="nbinom2")
m14L[[9]]=glmmadmb(l_t1 ~ log_l_t0 * sex + c_t0 + c_t0:sex + (1 | plot),data=d14,family="nbinom2")
# Target fitness + effect = sex 
m14L[[10]]=glmmadmb(l_t1 ~ log_l_t0 + sr + c_t0 + (1 | plot),data=d14,family="nbinom2")
m14L[[11]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + c_t0 + (1 | plot),data=d14,family="nbinom2")
m14L[[12]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + c_t0 + (1 | plot),data=d14,family="nbinom2")
# Target fitness + effect = sex + response=sex 
m14L[[13]]=glmmadmb(l_t1 ~ log_l_t0 + sr + c_t0 + sr:sex + c_t0:sex +(1 | plot),data=d14,family="nbinom2")
m14L[[14]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + c_t0 + sr:sex + c_t0:sex + (1 | plot),data=d14,family="nbinom2")
m14L[[15]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + c_t0 + sr:sex + c_t0:sex  + (1 | plot),data=d14,family="nbinom2")

# model selection
gr14_L_mod_sel  <- AICtab(m14L, weights = T)
gr14_L_avg      <- model_avg(gr14_L_mod_sel, m14L)
write.csv(gr14_L_avg, "Results/VitalRates_3/growth_14L_best.csv", row.names = F)


# 2015 --------------------------------------------------------------------------------

# Density as "new tillers"
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


# Planting density 
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

# Best models
gr14_avg   <- read.csv("Results/VitalRates_3/growth_14N_best.csv")
gr14_L_avg <- read.csv("Results/VitalRates_3/growth_14L_best.csv")
grow_avg   <- read.csv("Results/VitalRates_3/growth_best.csv")
grow_avgN  <- read.csv("Results/VitalRates_3/growthN_best.csv")


# Set up colors for plots
# 2015
d15$col=as.integer(d15$sex)
d15$col=as.character(factor(d15$col,labels=c("blue","red")))
# 2014
d14$col=as.integer(d14$sex)
d14$col=as.character(factor(d14$col,labels=c("blue","red")))

# 2015 ---------------------------------------------------------------------------------------

tiff("Results/VitalRates_3/growth_raw.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

par(mfrow=c(2,2),mar=c(3,2.5,1,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0))

# 2015 --------------------------------------------------------------------------
plot(d15$new_t1,d15$l_t1, pch=16, ylab="Target individuals: number  leaves in 2015",
     xlab="Total number of new tillers in plot (2015)",col=d15$col, ylim = c(0,99))
beta <- grow_avg[,c("predictor","avg")]$avg
xSeq <- seq(0,224,by=1 )
size <- mean( d15$log_l_t0 )
y_m <- exp( beta[1] + beta[2]*size + beta[3]*xSeq + beta[4] + 
            beta[5]*0.5 + beta[6]*size + beta[7]*xSeq + beta[8]*0.5)
y_f <- exp( beta[1] + beta[2]*size + beta[3]*xSeq + beta[5]*0.5 )
lines(xSeq,y_f,lwd=1.5,lty=1,col="blue")
lines(xSeq,y_m,lwd=1.5,lty=1,col="red")

legend(-10,99,c("Males","Females"),
       lty=1,lwd=1.5,col=c("red","blue"),bty="n")
title("2015")

# Density
plot(d15$TotDensity,d15$l_t1, pch=16, ylab="Target individuals: number  leaves in 2015",
     xlab="Planting density",col=d15$col, ylim = c(0,99))
beta <- grow_avgN[,c("predictor","avg")]$avg
xSeq <- seq(0,48,by=1)
size <- mean(d15$log_l_t0)
y_m <- exp( beta[1] + beta[2]*size + beta[3] + beta[4]*0.5 + 
              beta[5]*xSeq + beta[6]*size + beta[7]*0.5 + beta[8]*xSeq )
y_f <- exp( beta[1] + beta[2]*size + beta[4]*0.5 + beta[5]*xSeq )
lines(xSeq,y_f,lwd=1.5,lty=1,col="blue")
lines(xSeq,y_m,lwd=1.5,lty=1,col="red")
title("2015")


# 2014 -------------------------------------------------------------------------

# Total number of leaves
plot(d14$c_t0, d14$l_t1, pch=16, ylab="Target individuals: number  leaves in 2014",
     xlab="Total number of plot leaves in fall 2013",col=d14$col, ylim = c(0,99))
beta = gr14_L_avg[,c("predictor","avg")]$avg
xSeq <- seq(0,max(d14$c_t0),by=1)
sr   <- 0.5
size <- mean(d14$log_l_t0)
y_m <- exp( beta[1] + beta[2]*xSeq + beta[3]*size + beta[4] + 
              beta[5]*xSeq + beta[6]*0.5 + beta[7]*size + beta[8]*0.5)
y_f <- exp( beta[1] + beta[2]*xSeq + beta[3]*size + beta[6]*0.5)
lines(xSeq,y_f,lwd=1.5,lty=1,col="blue")
lines(xSeq,y_m,lwd=1.5,lty=1,col="red")
title("2014")

# Planting density
plot(d14$TotDensity, d14$l_t1, pch=16, ylab="Target individuals: number  leaves in 2014",
     xlab="Planting density",col=d14$col, ylim = c(0,99))
beta <- gr14_avg[,c("predictor","avg")]$avg
xSeq <- seq(0,48,by=1)
sr   <- 0.5
size <- mean(d14$log_l_t0)
y_m <- exp( beta[1] + beta[2]*size + beta[3] + beta[4]*xSeq + 
              beta[5]*xSeq + beta[6]*sr + beta[7]*size)
y_f <- exp( beta[1] + beta[2]*size + beta[4]*xSeq + beta[6]*sr )
lines(xSeq,y_f,lwd=1.5,lty=1,col="blue")
lines(xSeq,y_m,lwd=1.5,lty=1,col="red")
title("2014")

dev.off()
