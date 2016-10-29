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

tiff("Results/VitalRates_3/figure1.tiff",unit="in",width=6.3,height=2.1,res=600,compression="lzw")

par(mfrow=c(1,3),mar=c(2.5,2.5,0.1,0.5),mgp=c(1.4,0.35,0),cex.lab=1.1,cex.axis=0.8,
    cex.main=0.9)

# Flowering ------------------------------------------------------------------------------
# Set up colors for plots
f14$col <- as.character(factor(as.integer(f14$sex),labels=c("blue","red")))
mal14   <- subset(f14,sex=="m")
fem14   <- subset(f14,sex=="f")

plot( fem14$flowN_t1+0.05 ~ fem14$sr, pch = 16,ylim=c(0,7),xlim=c(0,1),
      ylab = "Number of flowers", xlab="Proportion of female individuals", col = "blue")
par(new=T) ; plot( mal14$flowN_t1-0.05 ~ mal14$sr,pch = 16,ylim=c(0,7),xlim=c(0,1),
                   ylab = "Number of flowers", xlab="",col = "red")

low   <- 5
high  <- 42
size  <- mean(f14$log_l_t0)
xSeq  <- seq(0,1,by = 0.1)
beta  <- n_flow_beta[,c("predictor","avg")]$avg

y_m_l <- exp( beta[1] + beta[2]*size + beta[3]*low + 
                beta[4]*low + beta[5]*xSeq + beta[6] + 
                beta[7]*size + beta[8]*xSeq)
y_m_h <- exp( beta[1] + beta[2]*size + beta[3]*high + 
                beta[4]*high + beta[5]*xSeq + beta[6] + 
                beta[7]*size + beta[8]*xSeq)
y_f_l <- exp( beta[1] + beta[2]*size + beta[3]*low + 
                beta[5]*xSeq)
y_f_h <- exp( beta[1] + beta[2]*size + beta[3]*high + 
                beta[5]*xSeq)

lines(xSeq,y_m_l,lty=2,lwd=2,col="red")
lines(xSeq,y_m_h,lty=1,lwd=2,col="red")
lines(xSeq,y_f_l,lty=2,lwd=2,col="blue")
lines(xSeq,y_f_h,lty=1,lwd=2,col="blue")

legend(0.15,7.5,c("high density","low density","male","female"),
       lty=c(1,2,1,1),lwd=2,col=c("black","black","red","blue"),bty="n")
text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.145,
     par("usr")[4]*0.96,"(a)", cex = 1.2, xpd = T)


# fecundity ----------------------------------------------------------------------------
plot(fecund_data$sr, fecund_data$SeedN, pch = 16, ylim=c(0,950),
     xlab = "Proportion female individuals", ylab = "Seeds per flower")
xSeq   <- seq(0,1,by = 0.1)
mSize  <- mean(fecund_data$log_l_t0)
beta   <- fec_beta$avg
y_low  <- exp(beta[1] + beta[2]*mSize + beta[3]*xSeq + 
                beta[4]*5 + beta[5]*xSeq*5)
y_high <- exp(beta[1] + beta[2]*mSize + beta[3]*xSeq + 
                beta[4]*42 + beta[5]*xSeq*42)
lines(xSeq, y_low, lwd = 2, lty = 2)
lines(xSeq, y_high, lwd = 2)
text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.145,
     par("usr")[4]*0.96,"(b)", cex = 1.2, xpd = T)

# viability ----------------------------------------------------------------------------
lower=quantile(viabVr$totFlow,prob=c(0.1))
upper=quantile(viabVr$totFlow,prob=c(0.9))

# Tetrazolium data
plot(viabVr$sr_f,viabVr$germ_ratio,pch=16,ylim=c(0,1.01),
     xlab="Proportion of female flowers",ylab="Seed viability rate")
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
beta=tetr_beta$avg
yMeanLow=inv.logit(beta[1] + beta[2]*xSeq + beta[3]*lower + beta[4]*xSeq*lower)
yMeanHigh=inv.logit(beta[1] + beta[2]*xSeq + beta[3]*upper + beta[4]*xSeq*upper)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)
text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.145,
     par("usr")[4]*0.96,"(c)", cex = 1.2, xpd = T)

dev.off()

