# Code to produce figure1 (panels a,b,c)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) #For AIC weights, because I'm lazy!
library(MASS)
library(dplyr)
library(boot)

#Read in data--------------------------------------------------------------------
d           <- read.csv("Data/vr.csv")
fem_seeds   <- read.csv("Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")
viabVr      <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                        stringsAsFactors = F)
# best models
n_flow_beta <- read.csv("Results/VitalRates_3/n_flowers_best.csv")
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv")
tetr_beta   <- read.csv("Results/VitalRates_3/germination_best.csv")

# FORMAT DATA -------------------------------------------------------------------

# Only data from 2014
d14         <- subset(d, year == 2014)
d14         <- subset(d14, surv_t1 != 0)
d14[which(d14$log_l_t1==-Inf),]   <- NA
d14$new_t1  <- as.numeric(as.character(d14$new_t1))
f14         <- na.omit(d14[,c("plot","flow_t1","log_l_t1","log_l_t0","flowN_t1",
                              "sex","sr","new_t1","TotDensity")])

# fecundity data ---------------------------------------------------------------------- 
fem_seeds$focalI    <- paste("f",fem_seeds$IndividualN,sep="")
fem_seeds           <- fem_seeds[,c("Plot","focalI","SeedN")] 
fem_seeds           <- fem_seeds[!is.na(fem_seeds$SeedN),]
names(fem_seeds)[1] <- "plot"
fecund_data         <- merge(d14,fem_seeds) 


# FIGURE 1 -----------------------------------------------------------------------------

tiff("Results/VitalRates_3/figure1.tiff",unit="in",width=6.3,height=2.5,res=600,compression="lzw")

par(mfrow=c(1,3),mar=c(2.5,2.5,0.1,0.5),mgp=c(1.4,0.35,0),cex.lab=1.1,cex.axis=0.8,
    cex.main=0.9, oma=c(0,0,0.2,0))

# Flowering ------------------------------------------------------------------------------
# Set up colors for plots
f14$col <- as.character(factor(as.integer(f14$sex),labels=c("blue","red")))
mal14   <- subset(f14,sex=="m")
fem14   <- subset(f14,sex=="f")

# means
meanF   <- fem14 %>% group_by(TotDensity) %>% 
            summarise( meanF = mean(flowN_t1), sdF = sd(flowN_t1) )
meanM   <- mal14 %>% group_by(TotDensity) %>% 
            summarise( meanM = mean(flowN_t1), sdM = sd(flowN_t1) )

# plots
plot(meanF$TotDensity , meanF$meanF, pch = 16, col = "blue", ylim = c(0,2),
     ylab = "Number of flowers", xlab="Planting density")
arrows(meanF$TotDensity, meanF$meanF - meanF$sdF*0.1, col = "blue",
       meanF$TotDensity, meanF$meanF + meanF$sdF*0.1, length=0.02, angle=90, code=3)
par(new = T)
plot(meanM$meanM ~ meanM$TotDensity , pch = 16, col = "red", ylim = c(0,2),
     ylab = "", xlab="")
arrows(meanM$TotDensity, meanM$meanM - meanM$sdM*0.1, col = "red",  
       meanM$TotDensity, meanM$meanM + meanM$sdM*0.1, length=0.02, angle=90, code=3)

size  <- mean(f14$log_l_t0)
xSeq  <- seq(0,48,by = 1)
beta  <- n_flow_beta[,c("predictor","avg")]$avg

y_m <- exp( beta[1] + beta[2]*size + beta[3]*xSeq + 
            beta[4]*xSeq + beta[5]*0.5 + beta[6] + 
            beta[7]*size + beta[8]*0.5)
y_f <- exp( beta[1] + beta[2]*size + beta[3]*xSeq + beta[5]*0.5)

lines(xSeq,y_m,lty=1,lwd=2,col="blue")
lines(xSeq,y_f,lty=1,lwd=2,col="red")

legend(1,2.15,c("male","female"), lty=1, lwd=2, 
       col=c("red","blue"),bty="n", pch = 16 )
text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.145,
     par("usr")[4]*0.96,"(a)", cex = 1.2, xpd = T)


# fecundity ----------------------------------------------------------------------------
plot(fecund_data$TotDensity, fecund_data$SeedN, pch = 1, 
     cex = fecund_data$sr * 1.5 , ylim = c(0, 1000),
     xlab = "Planting density", ylab = "Seeds per flower")
xSeq   <- seq(0,48,by = 1)
mSize  <- mean(fecund_data$log_l_t0)
beta   <- fec_beta$avg
y_low  <- exp(beta[1] + beta[2]*mSize + beta[3]*0.1 + 
                beta[4]*xSeq + beta[5]*xSeq*0.1)
y_high <- exp(beta[1] + beta[2]*mSize + beta[3]*0.9 + 
                beta[4]*xSeq + beta[5]*xSeq*0.9)

lines(xSeq, y_low, lwd = 2, lty = 2)
lines(xSeq, y_high, lwd = 2)

legend(-1,1075,c("90% female plot","10% female plot"),
       pch = 1, pt.cex = c(0.9,0.1)*1.5, bty="n")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.145,
     par("usr")[4]*0.96,"(b)", cex = 1.2, xpd = T)


# viability ----------------------------------------------------------------------------
lower=quantile(viabVr$totFlow,prob=c(0.1))
upper=quantile(viabVr$totFlow,prob=c(0.9))

# Tetrazolium data
plot(viabVr$totFlow,viabVr$germ_ratio,pch=1,ylim=c(0,1.01), 
     cex = viabVr$sr_f * 1.5 ,
     xlab="Number of flowers",ylab="Seed viability rate")

xSeq=seq(min(viabVr$totFlow),max(viabVr$totFlow),length.out=100)
beta=tetr_beta$avg
yMeanLow=inv.logit(beta[1] + beta[2]*0.1 + beta[3]*xSeq + beta[4]*xSeq*0.1)
yMeanHigh=inv.logit(beta[1] + beta[2]*0.9 + beta[3]*xSeq + beta[4]*xSeq*0.9)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)
text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.145,
     par("usr")[4]*0.96,"(c)", cex = 1.2, xpd = T)

dev.off()
