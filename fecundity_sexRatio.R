## AC: 7.20.2015
##1. Selection of models of POAR survival using data from LLELA's competition experiment.
##2. Year specific data anlayzed SEPARATELY.
##3.The fields "c_t0","fC_t0", and "mC_t0" refer to FALL 2013, even when "year==2015".
## Therefore, final file IS MISLEADING. Take this into account.
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) #For AIC weights, because I'm lazy!
library(MASS)

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

# merge data sets 
fecund_data      <- merge(d14,fem_seeds) 


# ANALYSIS ---------------------------------------------------------------------

# Model selection. Use negative binomial because poisson is overdispersed.
nsMod=list()
nsMod[[1]]=glm.nb(SeedN ~ log_l_t0, data=fecund_data)
nsMod[[2]]=glm.nb(SeedN ~ log_l_t0 + TotDensity,data=fecund_data)
nsMod[[3]]=glm.nb(SeedN ~ log_l_t0 + sr,data=fecund_data)
nsMod[[4]]=glm.nb(SeedN ~ log_l_t0 + sr + TotDensity,data=fecund_data)
nsMod[[5]]=glm.nb(SeedN ~ log_l_t0 + sr * TotDensity,data=fecund_data)
mod_select = AICtab(nsMod,weights=T)


# Model average (Models that make up more than 95% of weight)
betaList      <- list()
betaList[[1]] <- data.frame(predictor = names(coef(nsMod[[4]])),
                            parameter_1 = coef(nsMod[[4]]))
betaList[[2]] <- data.frame(predictor = names(coef(nsMod[[5]])),
                            parameter_2 = coef(nsMod[[5]]))
betaList[[3]] <- data.frame(predictor = names(coef(nsMod[[2]])),
                            parameter_3 = coef(nsMod[[2]]))

# Model averages
beta_avg  <- Reduce(function(...) merge(...,all=T), betaList)
beta_avg[is.na(beta_avg)] <- 0
beta_avg$avg <- (beta_avg$parameter_1 * mod_select$weight[1] + 
                 beta_avg$parameter_2 * mod_select$weight[2] + 
                 beta_avg$parameter_3 * mod_select$weight[3]) / sum(mod_select$weight[1:3]) 

write.csv(beta_avg, "Results/VitalRates_3/fecuntity_best.csv", row.names = F)

# GRAPH -------------------------------------------------------------------------

tiff("Results/VitalRates_3/fecundity.tiff",unit="in",width=3.5,height=3.5,res=600,compression="lzw")

#Start plotting
par(mfcol=c(1,1),mar=c(2.8,3,0.5,0.2),mgp=c(1.4,0.5,0))

plot(fecund_data$sr, fecund_data$SeedN, pch = 16,
     xlab = "Sex Ratio", ylab = "Panicule seed number")
xSeq   <- seq(0,1,by = 0.1)
mSize  <- mean(fecund_data$log_l_t0)
beta   <- beta_avg$avg
y_low  <- exp(beta[1] + beta[2]*mSize + beta[3]*xSeq + 
              beta[4]*11 + beta[5]*xSeq*11)
y_high <- exp(beta[1] + beta[2]*mSize + beta[3]*xSeq + 
                beta[4]*42 + beta[5]*xSeq*42)
lines(xSeq, y_low, lwd = 2, lty = 2)
lines(xSeq, y_high, lwd = 2)

legend(0,860,c("low density","high density"),
       lwd=2,lty=c(2,1),bty = "n", cex = 1.5)

dev.off()

