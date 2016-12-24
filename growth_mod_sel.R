# Sex competition predictor: SEX RATIO (female individuals/total individuals)
#setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
setwd("C:/Users/Aldo/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)    
library(glmmADMB) #Fit models with a Negative Binomial
library(dplyr)
#source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")
source("C:/Users/Aldo/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/Aldo/Documents/CODE/LLELA/model_avg_format.R")
source("C:/Users/Aldo/Documents/CODE/LLELA/model_sel_results.R")


# Read data ------------------------------------------------------------------
d       <- read.csv("Data/vr.csv")

# remove dead individuals (this is a GROWTH model!)
d       <- subset(d, surv_t1 != 0)
# logtransform leaf numbers
d       <- mutate(d, plot = as.factor(plot) ) # glmmadmb wants plot as a factor

# Year one
d14     <- subset(d, year == 2014)
# Year two
tmp15   <- subset(d, year == 2015)
# Missing new tillers data 
tmp15   <- mutate(tmp15, new_t1 = replace(new_t1, new_t1=="SKIPPED", NA)) 
tmp15   <- mutate(tmp15, new_t1 = replace(new_t1, new_t1=="cnf", NA))
tmp15   <- mutate(tmp15, new_t1 = replace(new_t1, new_t1=="", NA))
tmp15   <- mutate(tmp15, new_t1 = as.numeric(as.character(new_t1)))
d15     <- na.omit(dplyr::select(tmp15,l_t1,log_l_t0,plot,focalI,sex,new_t1,sr,TotDensity))


# 2014 --------------------------------------------------------------------------------

# Planting density
m14=list() # lMod stands for "leaf model" (density is quantified by N. of leaves)
m14[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=d14,family="nbinom2")
m14[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=d14,family="nbinom2")
m14[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=d14,family="nbinom2")
# Target fitness: effect + tot density 
m14[[4]]=glmmadmb(l_t1 ~ log_l_t0 + TotDensity + (1 | plot),data=d14,family="nbinom2")
m14[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + TotDensity + (1 | plot),data=d14,family="nbinom2")
m14[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + TotDensity + (1 | plot),data=d14,family="nbinom2")
# Target fitness: effect + sex ratio
m14[[7]]=glmmadmb(l_t1 ~ log_l_t0 + sr + (1 | plot),data=d14,family="nbinom2")
m14[[8]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + (1 | plot),data=d14,family="nbinom2")
m14[[9]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + (1 | plot),data=d14,family="nbinom2")
# Target fitness: tot density + sex ratio
m14[[10]]=glmmadmb(l_t1 ~ log_l_t0 + sr + TotDensity + (1 | plot),data=d14,family="nbinom2")
m14[[11]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + TotDensity + (1 | plot),data=d14,family="nbinom2")
m14[[12]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + TotDensity + (1 | plot),data=d14,family="nbinom2")
# Target fitness: tot density X sex ratio
m14[[13]]=glmmadmb(l_t1 ~ log_l_t0 + sr * TotDensity + (1 | plot),data=d14,family="nbinom2")
m14[[14]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr * TotDensity + (1 | plot),data=d14,family="nbinom2")
m14[[15]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr * TotDensity + (1 | plot),data=d14,family="nbinom2")

# Model average
gr14_mod_sel  <- AICtab(m14, weights = T)
gr14_avg      <- model_avg(gr14_mod_sel, m14)
write.csv(gr14_avg, "Results/VitalRates_3/growth_14N_best.csv", row.names = F)

# Model average summary table
gr14_avg_sum <- model_avg_format(gr14_avg)
write.csv(gr14_avg_sum, "Results/VitalRates_3/growth_14N_sum.csv", row.names = F)

# model selection result table
sel_res14     <- sel_results(gr14_mod_sel, 15, "sizet")
write.csv(sel_res14, "Results/VitalRates_3/growth14N_mod_sel.csv", row.names = F)


# 2015 --------------------------------------------------------------------------------

# Planting density
m15=list() # lMod stands for "leaf model" (density is quantified by N. of leaves)
m15[[1]]=glmmadmb(l_t1 ~ log_l_t0 + (1 | plot),data=d15,family="nbinom2")
m15[[2]]=glmmadmb(l_t1 ~ log_l_t0 + sex + (1 | plot),data=d15,family="nbinom2")
m15[[3]]=glmmadmb(l_t1 ~ log_l_t0 * sex + (1 | plot),data=d15,family="nbinom2")
# Target fitness: effect + tot density 
m15[[4]]=glmmadmb(l_t1 ~ log_l_t0 + TotDensity + (1 | plot),data=d15,family="nbinom2")
m15[[5]]=glmmadmb(l_t1 ~ log_l_t0 + sex + TotDensity + (1 | plot),data=d15,family="nbinom2")
m15[[6]]=glmmadmb(l_t1 ~ log_l_t0 * sex + TotDensity + (1 | plot),data=d15,family="nbinom2")
# Target fitness: effect + sex ratio
m15[[7]]=glmmadmb(l_t1 ~ log_l_t0 + sr + (1 | plot),data=d15,family="nbinom2")
m15[[8]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + (1 | plot),data=d15,family="nbinom2")
m15[[9]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + (1 | plot),data=d15,family="nbinom2")
# Target fitness: tot density + sex ratio
m15[[10]]=glmmadmb(l_t1 ~ log_l_t0 + sr + TotDensity + (1 | plot),data=d15,family="nbinom2")
m15[[11]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr + TotDensity + (1 | plot),data=d15,family="nbinom2")
m15[[12]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr + TotDensity + (1 | plot),data=d15,family="nbinom2")
# Target fitness: tot density X sex ratio
m15[[13]]=glmmadmb(l_t1 ~ log_l_t0 + sr * TotDensity + (1 | plot),data=d15,family="nbinom2")
m15[[14]]=glmmadmb(l_t1 ~ log_l_t0 + sex + sr * TotDensity + (1 | plot),data=d15,family="nbinom2")
m15[[15]]=glmmadmb(l_t1 ~ log_l_t0 * sex + sr * TotDensity + (1 | plot),data=d15,family="nbinom2")

# Model average
gr15_mod_sel <- AICtab(m15, weights = T)
gr15_avg     <- model_avg(gr15_mod_sel, m15)
write.csv(gr15_avg, "Results/VitalRates_3/growth_best.csv", row.names = F)

# Model average summary table
gr15_avg_sum <- model_avg_format(gr15_avg)
write.csv(gr15_avg_sum, "Results/VitalRates_3/growth_15N_sum.csv", row.names = F)

# model selection result table
sel_res15     <- sel_results(gr15_mod_sel, 15, "sizet")
write.csv(sel_res15, "Results/VitalRates_3/growth15N_mod_sel.csv", row.names = F)


# GRAPH ------------------------------------------------------------------------------------------------

# Best models
gr14_avg   <- read.csv("Results/VitalRates_3/growth_14N_best.csv")
gr15_avg   <- read.csv("Results/VitalRates_3/growth_best.csv")


# Set up colors for plots
# 2015
d15$col=as.integer(d15$sex)
d15$col=as.character(factor(d15$col,labels=c("blue","red")))
# 2014
d14$col=as.integer(d14$sex)
d14$col=as.character(factor(d14$col,labels=c("blue","red")))

# 2015 ---------------------------------------------------------------------------------------

#tiff("Results/VitalRates_3/growth_raw.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

par(mfrow=c(2,1),mar=c(3,2.5,1,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0))

# 2015 --------------------------------------------------------------------------
plot(d15$TotDensity,d15$l_t1, pch=16, ylab="Target individuals: number  leaves in 2015",
     xlab="Planting density",col=d15$col, ylim = c(0,99))
beta <- gr15_avg[,c("predictor","avg")]$avg
xSeq <- seq(0,48,by=1)
size <- mean(d15$log_l_t0)
y_m <- exp( beta[1] + beta[2]*size + beta[3]*size + beta[4] + 
              beta[5]*0.5 + beta[6]*xSeq*0.5 + beta[7]*xSeq )
y_f <- exp( beta[1] + beta[2]*size + beta[5]*0.5 + beta[6]*xSeq*0.5 + beta[7]*xSeq )
lines(xSeq,y_f,lwd=1.5,lty=1,col="blue")
lines(xSeq,y_m,lwd=1.5,lty=1,col="red")
title("2015")

# 2014 -------------------------------------------------------------------------
plot(d14$TotDensity, d14$l_t1, pch=16, ylab="Target individuals: number  leaves in 2014",
     xlab="Planting density",col=d14$col, ylim = c(0,99))
beta <- gr14_avg[,c("predictor","avg")]$avg
xSeq <- seq(0,48,by=1)
sr   <- 0.5
size <- mean(d14$log_l_t0)
y_m <- exp( beta[1] + beta[2]*size + beta[3] + beta[4]*xSeq + 
              beta[5]*0.5 + beta[6]*size + beta[7]*0.5*xSeq )
y_f <- exp( beta[1] + beta[2]*size + beta[4]*xSeq + beta[5]*0.5 + beta[7]*0.5*xSeq)
lines(xSeq,y_f,lwd=1.5,lty=1,col="blue")
lines(xSeq,y_m,lwd=1.5,lty=1,col="red")
title("2014")

#dev.off()
