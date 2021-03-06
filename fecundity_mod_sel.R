# model selection for fecundity (seeds per flower) 
#setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
setwd("C:/Users/Aldo/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(MASS)
library(dplyr)
#source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")
source("C:/Users/Aldo/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/Aldo/Documents/CODE/LLELA/model_avg_format.R")
source("C:/Users/Aldo/Documents/CODE/LLELA/model_sel_results.R")

#read in data--------------------------------------------------------------
d         <- read.csv("Data/vr.csv")
fem_seeds <- read.csv("Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")

# plot level data ---------------------------------------------------------

# Only data from 2014
d14 <- subset(d, year == 2014)
d14 <- subset(d14, surv_t1 != 0)

# fecundity data ---------------------------------------------------------------------- 
fem_seeds$focalI    <- paste("f",fem_seeds$IndividualN,sep="")
fem_seeds           <- fem_seeds[,c("Plot","focalI","SeedN")] 
fem_seeds           <- fem_seeds[!is.na(fem_seeds$SeedN),]
names(fem_seeds)[1] <- "plot"

# merge data sets 
fecund_data         <- merge(d14,fem_seeds) 
fecund_data$plot    <- as.factor(fecund_data$plot)

# ANALYSIS ---------------------------------------------------------------------

nsMod=list()
nsMod[[1]]= glmmadmb(SeedN ~ log_l_t0 + ( 1 | plot), family="nbinom2", data=fecund_data)
nsMod[[2]]= glmmadmb(SeedN ~ log_l_t0 + TotDensity + ( 1 | plot), family="nbinom2", data=fecund_data)
nsMod[[3]]= glmmadmb(SeedN ~ log_l_t0 + sr + ( 1 | plot), family="nbinom2", data=fecund_data)
nsMod[[4]]= glmmadmb(SeedN ~ log_l_t0 + sr + TotDensity + ( 1 | plot), family="nbinom2", data=fecund_data)
nsMod[[5]]= glmmadmb(SeedN ~ log_l_t0 + sr * TotDensity + ( 1 | plot), family="nbinom2", data=fecund_data)

# Model average
fec_select <- AICtab(nsMod,weights=T)
fec_avg    <- model_avg(fec_select, nsMod)
write.csv(fec_avg, "Results/VitalRates_3/fecuntity_best.csv", row.names = F)

# Model average summary table
fec_avg_sum   <- model_avg_format(fec_avg)
write.csv(fec_avg_sum, "Results/VitalRates_3/fecuntity_sum.csv", row.names = F)

# model selection result table
sel_fecund <- sel_results(fec_select, 5, "seeds")
write.csv(sel_fecund, "Results/VitalRates_3/fecundity_mod_sel.csv", row.names = F)


# GRAPH -------------------------------------------------------------------------

#tiff("Results/VitalRates_3/fecundity.tiff",unit="in",width=3.5,height=3.5,res=600,compression="lzw")

#Start plotting
par(mfcol=c(1,1),mar=c(2.8,3,0.5,0.2),mgp=c(1.4,0.5,0))

plot(fecund_data$sr, fecund_data$SeedN, pch = 16,
     xlab = "Sex Ratio", ylab = "Panicule seed number")
xSeq   <- seq(0,1,by = 0.1)
mSize  <- mean(fecund_data$log_l_t0)
beta   <- fec_avg$avg

y_low  <- exp(beta[1] + beta[2]*mSize + beta[3]*xSeq + 
              beta[4]*5 + beta[5]*xSeq*5)
y_high <- exp(beta[1] + beta[2]*mSize + beta[3]*xSeq + 
                beta[4]*42 + beta[5]*xSeq*42)
lines(xSeq, y_low, lwd = 2, lty = 2)
lines(xSeq, y_high, lwd = 2)

legend(-0.05,890,c("low density","high density"),
       lwd=2,lty=c(2,1),bty = "n", cex = 1.5)

#dev.off()

