##Response variable: total number of tillers 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) 
library(glmmADMB) # Fit models with a Negative Binomial
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")

# load and format data -----------------------------------------------------
x       <- read.csv("Data/vr.csv")
d       <- format_growth(x)
d14     <- subset(d, year == 2014)
d15     <- subset(d, year == 2015)
d15     <- subset(d15, plot != 149) # Remove outlier

# Model averages
avg14   <- read.csv("Results/VitalRates_3/new_t_bh14_best.csv")
avg15   <- read.csv("Results/VitalRates_3/new_t_bh15_best.csv")


# Graph -----------------------------------------------------------------------------

# Sex ratio as dot color ------
# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(gray.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

#Graph: total number of tillers -----------------------------------------------------
tiff("Results/VitalRates_3/Figure4.tiff",unit="in",width=4,height=4,res=600,compression="lzw")

par(mfrow=c(1,1),mar=c(2.6,2.5,0.2,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0.1))

# 2014 ----------------------------------------------------------------------------
# total number of tillers
plot(d14$TotDensity,d14$new_t1,pch=21,ylab="Number of new tillers",
     xlab="Planting density", bg = cRamp(d14$sr), cex = 1.5, 
     ylim = c(0,100))
N    <- seq(0,48,1)
beta <- as.data.frame(avg14)
fem  <- N*0.9
mal  <- N*0.1
y_h  <- (beta[,"lam.f"]*fem) / (1 + beta[,"b.f"] * (fem + beta[,"a.m"]*mal)) + 
  (beta[,"lam.m"]*mal) / (1 + beta[,"b.m"]* (beta[,"a.f"]*fem + mal))
fem  <- N*0.1
mal  <- N*0.9
y_l  <- (beta[,"lam.f"]*fem) / (1 + beta[,"b.f"] * (fem + beta[,"a.m"]*mal)) + 
  (beta[,"lam.m"]*mal) / (1 + beta[,"b.m"]* (beta[,"a.f"]*fem + mal))
lines(N,y_l,col="#636363",lwd=2)
lines(N,y_h,col="#DCDCDC",lwd=2, lty = 2)

# sex ratio legend
colfunc = colorRampPalette(cRamp(unique(arrange(d15,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(21), ncol=1))
text(x=2, y = seq(60,100,l=3), labels = seq(1,0,l=3))
rasterImage(legend_image, 4, 60, 9.5, 100)
text(8.8, 100, "Percent of", pos = 4)
text(8.8, 80, "females in", pos = 4)
text(8.8, 60, "plot", pos = 4)

# prediction legend
legend("topright", c("10% female plot", "90% female plot"), cex = 1,
       lty = c(1,2), lwd=2, col=c("#636363","#DCDCDC"), bty = "n")

dev.off()

