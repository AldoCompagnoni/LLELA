#setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
setwd("C:/Users/Aldo/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle)
library(lme4)
library(dplyr)
library(glmmADMB)
#source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")
source("C:/Users/Aldo/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/Aldo/Documents/CODE/LLELA/model_avg_format.R")
source("C:/Users/Aldo/Documents/CODE/LLELA/model_sel_results.R")


# load and format data --------------------------------------------------------------
d             <- read.csv("Data/vr.csv")
malPanicules  <- read.csv("Data/Spring 2014/maleCounts/malePaniculesSpring2014.csv")


# MALE panicule lengths --------------------------------------------------------------
d14                     <- subset(d, year==2014)
d14                     <- subset(d14, surv_t1!=0)
m_alloc                 <- merge(d14, malPanicules)
m_alloc$pan_area        <- m_alloc$panicule_Length_cm *  m_alloc$panicule_Width_cm * pi
m_alloc$plot            <- as.factor(m_alloc$plot)


# model selection ----------------------------------------------------------------------
spike_data              <- na.omit(dplyr::select(m_alloc, CountSpikelets, 
                                                 log_l_t0, sr, TotDensity, plot))

plMod=list()
plMod[[1]]=glmmadmb(CountSpikelets ~ log_l_t0 + (1 | plot), family = "nbinom2", data=spike_data)
plMod[[2]]=glmmadmb(CountSpikelets ~ log_l_t0 + TotDensity + (1 | plot), family = "nbinom2", data=spike_data)
plMod[[3]]=glmmadmb(CountSpikelets ~ log_l_t0 + sr + (1 | plot), family = "nbinom2", data=spike_data)
plMod[[4]]=glmmadmb(CountSpikelets ~ log_l_t0 + sr + TotDensity + (1 | plot), family = "nbinom2", data=spike_data)
plMod[[5]]=glmmadmb(CountSpikelets ~ log_l_t0 + sr * TotDensity + (1 | plot), family = "nbinom2", data=spike_data)

# Model average
m_alloc_select    <- AICtab(plMod,weights=T)
m_alloc_avg       <- model_avg(m_alloc_select, plMod)
write.csv(m_alloc_avg, "Results/VitalRates_3/male_spikelets.csv", row.names = F)

# Model average summary table
m_alloc_avg_sum   <- model_avg_format(m_alloc_avg)
write.csv(m_alloc_avg_sum, "Results/VitalRates_3/male_spikelets_sum.csv", row.names = F)

# model selection result table
sel_male_panic    <- sel_results(m_alloc_select, 5, "spikelets")
write.csv(sel_male_panic, "Results/VitalRates_3/male_spikelets_mod_sel.csv", row.names = F)


# Graph ---------------------------------------------------------------------------

#tiff("Results/VitalRates_3/male_repr_alloc.tiff",unit="in",width=3.5,height=3.5,res=600,compression="lzw")

#Start plotting
par(mfcol=c(1,1),mar=c(2.8,3,1,0.1),mgp=c(1.4,0.5,0))

#best model
plot(spike_data$TotDensity,spike_data$CountSpikelets,pch=1, 
     cex = spike_data$sr * 1.5,
     ylab="Spikelet number",xlab="Planting density")
size <- mean(spike_data$log_l_t0)
beta <- m_alloc_avg$avg
xSeq <- seq(1,48,by =1)
y_l  <- exp(beta[1] + beta[2]*size + beta[3]*xSeq + beta[4]*0.1 + beta[5]*xSeq*0.1)
y_h  <- exp(beta[1] + beta[2]*size + beta[3]*xSeq + beta[4]*0.9 + beta[5]*xSeq*0.9)

lines(xSeq,y_l,lwd=2,lty=2)
lines(xSeq,y_h,lwd=2,lty=1)

#dev.off()
