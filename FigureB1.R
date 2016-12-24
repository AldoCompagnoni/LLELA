# Code to produce figure B1
#setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
setwd("C:/Users/Aldo/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(dplyr)
library(boot)


#Read data--------------------------------------------------------------------
viabVr      <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                        stringsAsFactors = F)

# best models
germ_beta   <- read.csv("Results/VitalRates_3/germination_best.csv")

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
tiff("Results/VitalRates_3/figureB1.tiff",unit="in",width=6.3,height=3.15,res=600,compression="lzw")

par(mfrow=c(1,2), mar=c(2.5,2.5,0.1,0.1),mgp=c(1.4,0.35,0),cex.lab=1,cex.axis=0.8,
    cex.main=0.9, oma=c(0,0,0.2,0))

plot(jitter(viabVr$totFlow,factor = 2),jitter(viabVr$tetra_ratio,factor = 2),pch=21,ylim=c(0,1.01), 
     bg = cRamp(viabVr$sr_f), cex = 1.5,
     xlab="Number of flowers",ylab="Seed viability rate")
xSeq=seq(min(viabVr$totFlow),max(viabVr$totFlow),length.out=100)
beta=germ_beta$avg
size=mean(viabVr$log_l_t0)
yMeanLow=inv.logit(beta[1] + beta[2]*size + beta[3]*0.2 + beta[4]*xSeq*0.2 + beta[5]*xSeq)
yMeanHigh=inv.logit(beta[1] + beta[2]*size + beta[3]*1 + beta[4]*xSeq*1 + beta[5]*xSeq)

lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

# legend
plot(viabVr$totFlow,viabVr$tetra_ratio, type = "n", xaxt="n", yaxt="n",
     xlab="",ylab="", bty = "n")

legend(-5,0.6,c("100% female plot","  20% female plot"),lty = c(1,2), lwd = 2, bty="n")

colfunc = colorRampPalette(cRamp(unique(arrange(viabVr,sr_f)$sr_f)))
legend_image <- as.raster(matrix(colfunc(19), ncol=1))
text(x=6.5, y = seq(0.65,0.85,l=3), labels = seq(0,1,l=3))
rasterImage(legend_image, -6, 0.65, 3, 0.85)
text(7.5, 0.85, "Percent of", pos = 4)
text(7.5, 0.75, "males in", pos = 4)
text(7.5, 0.65, "plot", pos = 4)

dev.off()
