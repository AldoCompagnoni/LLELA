##Response variable: total number of tillers 
setwd("C:/cloud/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) 
library(glmmADMB) # Fit models with a Negative Binomial
library(dplyr)
source("C:/CODE/LLELA/model_avg.R")

# load and format data -----------------------------------------------------
x       <- read.csv("Data/vr.csv")
d       <- format_growth(x)
d14     <- subset(d, year == 2014)
d15     <- subset(d, year == 2015)
d15     <- subset(d15, plot != 149) # Remove outlier

# Model averages
avg15   <- read.csv("Results/VitalRates_4_Cade2015/new_t_bh15_best.csv")


# Graph -----------------------------------------------------------------------------

# Sex ratio as dot color ------
# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(topo.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

#Graph: total number of tillers -----------------------------------------------------
tiff("Results/VitalRates_4_Cade2015/Figure_new_t_15.tiff",unit="in",width=4,height=4,res=600,compression="lzw")

par(mfrow=c(1,1),mar=c(2.6,2.5,0.2,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0.1))

plot(d15$TotDensity,d15$new_t1,pch=21,ylab="Number of new tillers",
     xlab="Planting density", bg = cRamp(d14$sr), cex = 1.5)

# select
l_des <- subset(avg15, sr == 0.05)
h_des <- subset(avg15, sr == 0.95)

lines(h_des$TotN, h_des$pred, col="#ABABAB",lwd=2, lty = 1)
lines(l_des$TotN, l_des$pred, col="#141414",lwd=2, lty = 3)

# sex ratio legend
colfunc = colorRampPalette(cRamp(unique(arrange(d15,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(21), ncol=1))
text(x=2, y = seq(185,225,l=3), labels = seq(1,0,l=3))
rasterImage(legend_image, 4, 185, 9.4, 225)
text(8.8, 225, "Percent of", pos = 4)
text(8.8, 205, "females in", pos = 4)
text(8.8, 185, "plot", pos = 4)

# prediction legend
# legend(22, par("usr")[4], c("5% female plot", "95% female plot"), cex = 0.9,
#        lty = c(1,2), lwd=2, col=c("#636363","#DCDCDC"), bty = "n")
legend(0, 180, c("95% female plot", "5% female plot"), cex = 1,
       lty = c(1,3), lwd=2, col=c("#ABABAB","#141414"), bty = "n")

dev.off()

