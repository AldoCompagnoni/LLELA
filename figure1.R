# Code to produce figure1 (panels a-d)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(dplyr)
library(boot)

#Read data--------------------------------------------------------------------
d           <- read.csv("Data/vr.csv")
pans        <- read.csv("Data/Spring 2014/plot_level_panicles.csv")
fem_seeds   <- read.csv("Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")
m_spiklet   <- read.csv("Data/Spring 2014/maleCounts/malePaniculesSpring2014.csv")
viabVr      <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                        stringsAsFactors = F)

# best models
f_pan_avg   <- read.csv("Results/VitalRates_3/f_flowers_plot_best.csv")
m_pan_avg   <- read.csv("Results/VitalRates_3/m_flowers_plot_best.csv")
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv")
m_spik_beta <- read.csv("Results/VitalRates_3/male_spikelets.csv")
germ_beta   <- read.csv("Results/VitalRates_3/germination_best.csv")


# FORMAT DATA -------------------------------------------------------------------

# plot level data
design    <- dplyr::select(pans,plot,N,sr)

# Number of flowers in 2014
pans_f    <- subset(pans, F != 0) # plots containing females
pans_m    <- subset(pans, M != 0) # plots containing males

# fecundity data  
fem_seeds   <- mutate(fem_seeds, focalI = paste("f",IndividualN,sep=""))
fem_seeds   <- dplyr::select(fem_seeds,Plot,focalI,SeedN)
fem_seeds   <- filter(fem_seeds, !is.na(SeedN) )
fem_seeds   <- rename(fem_seeds, plot = Plot)
fecund_data <- merge(design,fem_seeds) 

# male allocation (spikelets) data  
m_alloc     <- merge(design, m_spiklet)

# remove three plots with more than 60 flowers 
viabVr      <- subset(viabVr, totFlow < 60)


# FIGURE 1 -----------------------------------------------------------------------------

# Sex ratio as dot color ------
# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(topo.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

# Graph
tiff("Results/VitalRates_3/figure1.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

par(mfrow=c(2,2),mar=c(3,2.5,0.1,0.1),mgp=c(1.4,0.35,0),cex.lab=1.1,cex.axis=0.8,
    cex.main=0.9, oma=c(0,0,0.2,0))

# Flowering ------------------------------------------------------------------------------

# Plot data
plot(jitter(pans_f$F_flow_t1) ~ pans_f$N, xlab=expression("Planting density"),
     ylab="Number of flowers", pch=16, col="blue", ylim = c(-0.5,90))
par(new=T) ; 
plot(jitter(pans_m$M_flow_t1) ~ pans_m$N,pch=17,xlab="",ylab="",col="red",
     xaxt="n", ylim = c(-0.5,90))

# lines
xSeq  <- seq(1,48,1)
betaF <- f_pan_avg$avg
betaM <- m_pan_avg$avg
yF_h  <- (betaF[1] + betaF[2]*xSeq + betaF[4]*xSeq*1 + betaF[3]*1)
yF_l  <- (betaF[1] + betaF[2]*xSeq + betaF[4]*xSeq*0.2 + betaF[3]*0.2)
yM_h  <- (betaM[1] + betaM[2]*xSeq + betaM[4]*xSeq*1 + betaM[3]*1)
yM_l  <- (betaM[1] + betaM[2]*xSeq + betaM[4]*xSeq*0.2 + betaM[3]*0.2)

lines(xSeq, yF_h, col = "blue", lwd = 2, lty = 1)
lines(xSeq, yF_l, col = "blue", lwd = 2, lty = 2)
lines(xSeq, yM_h, col = "red",  lwd = 2, lty = 1)
lines(xSeq, yM_l, col = "red",  lwd = 2, lty = 2)

# legends
legend(0,95,c("100% female plot","  20% female plot"),lty = c(1,2), lwd = 2, bty="n")
legend(1.5,82,c("male","female"), lty=1, lwd=2, col=c("red","blue"), bty="n", pch = c(17,16) )

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(a)", cex = 1.2, xpd = T)


# fecundity ----------------------------------------------------------------------------
plot(fecund_data$N, fecund_data$SeedN, pch = 21, 
     bg = cRamp(fecund_data$sr) , ylim = c(0, 1000), cex = 1.5,
     xlab = "Planting density", ylab = "Seeds per flower")
xSeq   <- seq(0,48,by = 1)
beta   <- fec_beta$avg
y_low  <- exp(beta[1] + beta[2]*0.1 + beta[3]*xSeq + beta[4]*xSeq*0.1)
y_high <- exp(beta[1] + beta[2]*1 + beta[3]*xSeq + beta[4]*xSeq*1)

lines(xSeq, y_low, lwd = 2, lty = 2)
lines(xSeq, y_high, lwd = 2)

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(b)", cex = 1.2, xpd = T)

colfunc = colorRampPalette(cRamp(unique(arrange(fecund_data,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(19), ncol=1))
text(x=24.5, y = seq(800,950,l=3), labels = seq(0,1,l=3))
rasterImage(legend_image, 17, 800, 22, 950)
text(25, 950, "Percent of", pos = 4)
text(25, 875, "males in", pos = 4)
text(25, 800, "plot", pos = 4)


# male allocation ----------------------------------------------------------------------
plot(m_alloc$N,m_alloc$CountSpikelets,pch=21, 
     bg = cRamp(m_alloc$sr), cex = 1.5,
     ylab="Spikelet number per male flower",xlab="Planting density")
beta <- m_spik_beta$avg
xSeq <- seq(1,48,by =1)
y_l  <- exp(beta[1] + beta[2]*xSeq + beta[3]*0.2 + beta[4]*xSeq*0.2)
y_h  <- exp(beta[1] + beta[2]*xSeq + beta[3]*1 + beta[4]*xSeq*1)

lines(xSeq,y_l,lwd=2,lty=2)
lines(xSeq,y_h,lwd=2,lty=1)

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(c)", cex = 1.2, xpd = T)


# viability (germination) data ---------------------------------------------------------
plot(jitter(viabVr$totFlow,factor = 2),jitter(viabVr$germ_ratio,factor = 2),pch=21,ylim=c(0,1.01), 
     bg = cRamp(viabVr$sr_f), cex = 1.5,
     xlab="Number of flowers",ylab="Seed viability rate")
xSeq=seq(min(viabVr$totFlow),max(viabVr$totFlow),length.out=100)
beta=germ_beta$avg
yMeanLow=inv.logit(beta[1] + beta[2]*0.2 + beta[3]*xSeq*0.2 + beta[4]*xSeq)
yMeanHigh=inv.logit(beta[1] + beta[2]*1  + beta[3]*xSeq*1   + beta[4]*xSeq)

lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(d)", cex = 1.2, xpd = T)

dev.off()

