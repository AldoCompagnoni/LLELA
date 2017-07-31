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
avg14   <- read.csv("Results/VitalRates_3/new_t_bh14_best.csv")
avg15   <- read.csv("Results/VitalRates_3/new_t_bh15_best.csv")


# Graph -----------------------------------------------------------------------------

# Sex ratio as dot color ------
# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(gray.colors(7, start = 0, end = 1))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

#Graph: total number of tillers -----------------------------------------------------
tiff("Results/VitalRates_3/FigureE1.tiff",unit="in",width=4,height=4,res=600,compression="lzw")

par(mfrow=c(1,1),mar=c(2.6,2.5,0.2,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0.1))

# 2015 ----------------------------------------------------------------------------
# total number of tillers
plot(d15$TotDensity,d15$new_t1,pch=21,ylab="Number of new tillers",
     xlab="Planting density", bg = cRamp(d15$sr), cex = 1.5)
N    <- seq(0,48,1)
beta <- as.data.frame(avg15)
fem  <- N*0.95
mal  <- N*0.05
y_h  <- (beta[,"lam.f"]*fem) / (1 + beta[,"b.f"] * (fem + beta[,"a."]*mal)) + 
  (beta[,"lam.m"]*mal) / (1 + beta[,"b.m"] * (fem + beta[,"a."]*mal))
fem  <- N*0.05
mal  <- N*0.95
y_l  <- (beta[,"lam.f"]*fem) / (1 + beta[,"b.f"] * (fem + beta[,"a."]*mal)) + 
  (beta[,"lam.m"]*mal) / (1 + beta[,"b.m"] * (fem + beta[,"a."]*mal))
lines(N,y_h,lwd=2, lty = 1) #col="#DCDCDC",
lines(N,y_l,lwd=2, lty = 3) #col="#636363",

# sex ratio legend
colfunc = colorRampPalette(cRamp(unique(arrange(d15,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(21), ncol=1))
text(x=6.5, y = seq(180,220,l=3), labels = seq(1,0,l=3))
rasterImage(legend_image, 0, 180, 5, 220)
text(7, 220, "Percent of", pos = 4)
text(7, 200, "females in", pos = 4)
text(7, 180, "plot", pos = 4)

# prediction legend
legend(20,220, c("95% female plot", "5% female plot"), cex = 1,
       lty = c(1,3), lwd=2,  bty = "n")

dev.off()



