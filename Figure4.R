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
grow_new_t14_avg  <- read.csv("Results/VitalRates_3/growth_new_t_14.csv")
grow_new_t15_avg  <- read.csv("Results/VitalRates_3/growth_new_t_15.csv")


# Graph -----------------------------------------------------------------------------

# Sex ratio as dot color ------
# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(topo.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  


#Graph: total number of tillers -----------------------------------------------------
tiff("Results/VitalRates_3/Figure4.tiff",unit="in",width=6.3,height=3.15,res=600,compression="lzw")

par(mfrow=c(1,2),mar=c(2.6,2.5,1.3,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0.1))

# 2014 ----------------------------------------------------------------------------
# total number of tillers
plot(d14$TotDensity,d14$new_t1,pch=21,ylab="Number of new tillers",
     xlab="Planting density", bg = cRamp(d14$sr), cex = 1.5, ylim = c(0,max(d15$new_t1,na.rm=T)))
beta  <-  grow_new_t14_avg[,c("predictor","avg")]$avg
xSeq  <-  seq(0,48,by=1)
yL    <-  exp( beta[1] + beta[2]*xSeq + beta[3]*xSeq^2 + beta[4]*0.1 + beta[5]*xSeq*0.1 )
yH    <-  exp( beta[1] + beta[2]*xSeq + beta[3]*xSeq^2 + beta[4]*0.9 + beta[5]*xSeq*0.9 )
lines(xSeq,yL,col="blue3",lwd=2)
lines(xSeq,yH,col="khaki3",lwd=2, lty = 2)

# sex ratio legend
colfunc = colorRampPalette(cRamp(unique(arrange(d15,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(21), ncol=1))
text(x=4, y = seq(180,220,l=3), labels = seq(0,1,l=3))
rasterImage(legend_image, 7, 180, 12.5, 220)
text(11.5, 220, "Percent of", pos = 4)
text(11.5, 200, "males in", pos = 4)
text(11.5, 180, "plot", pos = 4)

# prediction legend
legend(1, 170, c("90% male plot", "10% male plot"), cex = 1,
       lty = c(1,2), lwd=2, col=c("blue3","khaki3"), bty = "n")

# panel ID
text(-2,245,"a) Plot level data, year 2014",pos=4, xpd = NA)


# 2015 ----------------------------------------------------------------------------
# total number of tillers
plot(d15$TotDensity,d15$new_t1,pch=21,ylab="Number of new tillers",
     xlab="Planting density", bg = cRamp(d15$sr), cex = 1.5)
beta  <-  grow_new_t15_avg[,c("predictor","avg")]$avg
xSeq  <-  seq(0,48,by=1)
yL    <-  exp( beta[1] + beta[2]*xSeq + beta[3]*xSeq^2 + beta[4]*0.1 + beta[5]*xSeq*0.1 )
yH    <-  exp( beta[1] + beta[2]*xSeq + beta[3]*xSeq^2 + beta[4]*0.9 + beta[5]*xSeq*0.9 )
lines(xSeq,yL,col="blue3",lwd=2)
lines(xSeq,yH,col="khaki3",lwd=2, lty = 2)

# panel ID
text(-2,245,"b) Plot level data, year 2015",pos=4, xpd = NA)

dev.off()
