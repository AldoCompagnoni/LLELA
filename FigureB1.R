# Code to produce figure B1
setwd("C:/cloud/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
#setwd("C:/Users/Aldo/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(dplyr)
library(boot)
source("C:/CODE/LLELA/model_avg.R")
source("C:/CODE/LLELA/prediction.R")


#Read data--------------------------------------------------------------------
viabVr      <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                        stringsAsFactors = F)
# remove three plots with more than 60 flowers 
viabVr      <- subset(viabVr, totFlow < 60)

# best models
germ_beta   <- read.csv("Results/VitalRates_4_Cade2015/germination_best.csv")
tetr_beta   <- read.csv("Results/VitalRates_4_Cade2015/tetrazolium_best.csv")


# Figure B1 ------------------------------------------------------------

# Sex ratio as dot color ------
# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(topo.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}


# Graph
tiff("Results/VitalRates_4_Cade2015/figureB1.tiff",unit="in",width=6.3,height=2.7,res=600,compression="lzw")

par(mfrow=c(1,2), mar=c(2.5,2.3,1,0.1),mgp=c(1.4,0.3,0),cex.lab=1,cex.axis=0.8,
    cex.main=0.9, oma=c(0,0,0,9), tcl = -0.3)

# tetrazolium data ----------------------------------------------------------------
plot(jitter(viabVr$totFlow,factor = 2),jitter(viabVr$tetra_ratio,factor = 2),pch=21,ylim=c(0,1), 
     bg = cRamp(viabVr$sr_f), cex = 1.1, main = "Tetrazolium data",
     xlab="",ylab="Seed viability rate")

# three lines: three model matrices
l_des <- subset(tetr_beta, sr_f == 0.05)
m_des <- subset(tetr_beta, sr_f == 0.5)
h_des <- subset(tetr_beta, sr_f == 0.95)

# plot lines
lines(l_des$totFlow,l_des$pred,lwd=2,lty=3)
lines(m_des$totFlow,m_des$pred,lwd=2,lty=2)
lines(h_des$totFlow,h_des$pred,lwd=2,lty=1)


# germination data -----------------------------------------------------------------
plot(jitter(viabVr$totFlow,factor = 2),jitter(viabVr$germ_ratio,factor = 2),pch=21,ylim=c(0,1), 
     bg = cRamp(viabVr$sr_f), cex = 1.1, main = "Germination data",
     xlab="",ylab="")

# three lines: three model matrices
l_des <- subset(germ_beta, sr_f == 0.05)
m_des <- subset(germ_beta, sr_f == 0.5)
h_des <- subset(germ_beta, sr_f == 0.95)

# plot lines
lines(l_des$totFlow,l_des$pred,lwd=2,lty=3)
lines(m_des$totFlow,m_des$pred,lwd=2,lty=2)
lines(h_des$totFlow,h_des$pred,lwd=2,lty=1)

mtext("Number of panicles", 1, line = 1.2, at = -10)

# legends -----------------------------------------------------------------
# prediction legend
legend(54,1.1, c("5% female plot",
               "50% female plot",
               "95% female plot"), cex = 1, xpd = NA,
       lty = c(3,2,1), lwd=2, bty = "n")

# color legend
colfunc = colorRampPalette(cRamp(unique(arrange(viabVr,sr_f)$sr_f)))
legend_image <- as.raster(matrix(colfunc(19), ncol=1))
text(x=62, y = seq(0.25,0.6,l=3), labels = seq(1,0,l=3), xpd=NA, cex = 0.8)
rasterImage(legend_image, 67, 0.25, 75, 0.6,xpd=NA)
text(75, 0.6, "Percent of", pos = 4,xpd=NA)
text(75, 0.48, "female", pos = 4,  xpd=NA)
text(75, 0.4, "panicles", pos = 4,  xpd=NA)
text(75, 0.25, "in plot", pos = 4,      xpd=NA)

dev.off()