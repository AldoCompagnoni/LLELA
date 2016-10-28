# model selection for fecundity (total seed per plant) 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) #For AIC weights, because I'm lazy!
library(MASS)
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")

#read in data--------------------------------------------------------------
d         <- read.csv("Data/vr.csv")
fem_seeds <- read.csv("Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")

# plot level data ---------------------------------------------------------
#Remove resuscitated individuals (dead in spring 2014, alive in 2015)
#NOTE: probably more adequate to consider these NOT DEAD in 2014, because they all have
#few leaves: they probably resprouted from base?
d <- d[-which(d$plot == 18 & d$focalI == "m2"),]
d <- d[-which(d$plot == 38 & d$focalI == "f1"),]
d <- d[-which(d$plot == 46 & d$focalI == "f1"),]
d <- d[-which(d$plot == 83 & d$focalI == "f5"),]
d <- d[-which(d$plot == 36 & d$focalI == "m3"),]

# logtransform leaf numbers
d$log_l_t0 <- log(d$l_t0)
d$log_l_t1 <- log(d$l_t1)
# Sex ratios
d$sr <- d$F / d$TotDensity  

# Only data from 2014
d14 <- subset(d, year == 2014)
d14 <- subset(d14, surv_t1 != 0)

# fecundity data ---------------------------------------------------------------------- 
fem_seeds$focalI    <- paste("f",fem_seeds$IndividualN,sep="")
fem_seeds           <- fem_seeds[,c("Plot","focalI","SeedN")] 
fem_seeds           <- fem_seeds[!is.na(fem_seeds$SeedN),]
names(fem_seeds)[1] <- "plot"
fem_seeds           <- fem_seeds %>% 
                        group_by(plot,focalI) %>% 
                          summarise(seed_tot = sum(SeedN))
# merge data sets 
fecund_data      <- merge(d14,fem_seeds) 


# ANALYSIS ---------------------------------------------------------------------
nsMod=list()
nsMod[[1]]=glm.nb(seed_tot ~ log_l_t0, data=fecund_data)
nsMod[[2]]=glm.nb(seed_tot ~ log_l_t0 + TotDensity,data=fecund_data)
nsMod[[3]]=glm.nb(seed_tot ~ log_l_t0 + sr,data=fecund_data)
nsMod[[4]]=glm.nb(seed_tot ~ log_l_t0 + sr + TotDensity,data=fecund_data)
nsMod[[5]]=glm.nb(seed_tot ~ log_l_t0 + sr * TotDensity,data=fecund_data)

# Model average
fec_select <- AICtab(nsMod,weights=T)
fec_avg    <- model_avg(fec_select, nsMod)
write.csv(fec_avg, "Results/VitalRates_3/fecuntity_best.csv", row.names = F)


# GRAPH -------------------------------------------------------------------------

tiff("Results/VitalRates_3/fecundity.tiff",unit="in",width=3.5,height=3.5,res=600,compression="lzw")

#Start plotting
par(mfcol=c(1,1),mar=c(2.8,3,0.5,0.2),mgp=c(1.4,0.5,0))

plot(fecund_data$sr, fecund_data$seed_tot, pch = 16,
     xlab = "Sex Ratio", ylab = "Panicule seed number")
xSeq   <- seq(0,1,by = 0.1)
mSize  <- mean(fecund_data$log_l_t0)
beta   <- fec_avg$avg
y_low  <- exp(beta[1] + beta[2]*mSize + beta[3]*xSeq + 
              beta[4]*5 + beta[5]*xSeq*5)


y_low  <- exp(beta[1] + beta[2]*mSize + beta[3]*xSeq + 
                beta[4]*5 + beta[5]*xSeq*5)

predictors[1,]

y_high <- exp(beta[1] + beta[2]*mSize + beta[3]*xSeq + 
                beta[4]*42 + beta[5]*xSeq*42)
lines(xSeq, y_low, lwd = 2, lty = 2)
lines(xSeq, y_high, lwd = 2)

legend(0,1240,c("low density","high density"),
       lwd=2,lty=c(2,1),bty = "n", cex = 1.5)

dev.off()

