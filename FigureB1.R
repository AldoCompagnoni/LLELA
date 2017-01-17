# Code to produce figure B1
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
#setwd("C:/Users/Aldo/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(dplyr)
library(boot)


#Read data--------------------------------------------------------------------
viabVr      <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                        stringsAsFactors = F)

# best models
germ_beta   <- read.csv("Results/VitalRates_3/germination_best.csv")
tetr_beta   <- read.csv("Results/VitalRates_3/tetrazolium_best.csv")

# remove three plots with more than 60 flowers 
viabVr      <- subset(viabVr, totFlow < 60)


# Figure B1 ------------------------------------------------------------

# Sex ratio as dot color ------
# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(topo.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}


# Graph
tiff("Results/VitalRates_3/figureB1.tiff",unit="in",width=6.3,height=2.7,res=600,compression="lzw")

par(mfrow=c(1,2), mar=c(2.5,2.5,1,0.1),mgp=c(1.4,0.35,0),cex.lab=1,cex.axis=0.8,
    cex.main=0.9, oma=c(0,0,0,8))

# tetrazolium data
plot(jitter(viabVr$totFlow,factor = 2),jitter(viabVr$tetra_ratio,factor = 2),pch=21,ylim=c(0,1), 
     bg = cRamp(viabVr$sr_f), cex = 1.1, main = "Tetrazolium data",
     xlab="Number of flowers",ylab="Seed viability rate")
xSeq=seq(min(viabVr$totFlow),max(viabVr$totFlow),length.out=100)
beta=tetr_beta$avg
yMeanLow=inv.logit(beta[1] + beta[2]*0.2 + beta[3]*xSeq*0.2 + beta[4]*xSeq)
yMeanHigh=inv.logit(beta[1] + beta[2]*1 + beta[3]*xSeq*1 + beta[4]*xSeq)

lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

# germination data
plot(jitter(viabVr$totFlow,factor = 2),jitter(viabVr$germ_ratio,factor = 2),pch=21,ylim=c(0,1), 
     bg = cRamp(viabVr$sr_f), cex = 1.1, main = "Germination data",
     xlab="Number of flowers",ylab="Seed viability rate")
xSeq=seq(min(viabVr$totFlow),max(viabVr$totFlow),length.out=100)
beta=germ_beta$avg
yMeanLow=inv.logit(beta[1] + beta[2]*0.2 + beta[3]*xSeq*0.2 + beta[4]*xSeq)
yMeanHigh=inv.logit(beta[1] + beta[2]*1 + beta[3]*xSeq*1 + beta[4]*xSeq)

lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

# legend
legend(55,1.05,c("100% female","  20% female"),lty = c(1,2), lwd = 2, bty="n", xpd=NA)

colfunc = colorRampPalette(cRamp(unique(arrange(viabVr,sr_f)$sr_f)))
legend_image <- as.raster(matrix(colfunc(19), ncol=1))
text(x=62, y = seq(0.45,0.65,l=3), labels = seq(0,1,l=3), xpd=NA, cex = 0.8)
rasterImage(legend_image, 67, 0.45, 75, 0.65,xpd=NA)
text(75, 0.65, "Percent of", pos = 4,xpd=NA)
text(75, 0.55, "males in", pos = 4,  xpd=NA)
text(75, 0.45, "plot", pos = 4,      xpd=NA)

dev.off()