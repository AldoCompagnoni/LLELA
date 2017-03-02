# Code to produce figure1 (panels a,b,c)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(dplyr)


#Read and format data--------------------------------------------------------------------
d       <- read.csv("Data/vr.csv")
d       <- format_growth(d)
# separate years
d14     <- filter(d, year == 2014)
d15     <- filter(d, year == 2015)


# FIGURE 3 ----------------------------------------------------------------------

# Best models
lr_14_avg   <- read.csv("Results/VitalRates_3/log_ratio_14_best.csv")
lr_15_avg   <- read.csv("Results/VitalRates_3/log_ratio_15_best.csv")

# service functions
# color palette
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(topo.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

# Symbols for the sexes 
sex_symbol <- function(x){
  out <- x %>% mutate(symb = as.integer( as.character( 
    factor( as.integer(sex),labels=c("21","24")))) )
  return(out)
}
d14 <- sex_symbol(d14)
d15 <- sex_symbol(d15)


# plot
tiff("Results/VitalRates_3/Figure3.tiff",unit="in",width=6.3,height=3.15,res=600,compression="lzw")

par(mfrow=c(1,2),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0))

# 2014 
plot(jitter(d14$TotDensity), d14$log_ratio, pch = d14$symb, ylab="Log ratio (2014)", 
     xlab="Planting density",bg=cRamp(d14$sr), ylim = c(-3,2.5), 
     lwd = 1)
beta  <- lr_14_avg[,c("predictor","avg")]$avg
y_m   <- beta[1] + beta[2] + beta[3]*0.5
y_f   <- beta[1] + beta[3]*0.5
abline(h = y_f,col="khaki3",lwd=2)
abline(h = y_m,col="blue3",lwd=2, lty = 2)


colfunc = colorRampPalette(cRamp(unique(arrange(d15,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(21), ncol=1))
text(x=24, y = seq(-3,-1.5,l=3), labels = seq(0,1,l=3))
rasterImage(legend_image, 27, -3, 31, -1.5)
text(30, -1.5, "Percent of", pos = 4)
text(30, -2.25, "males in", pos = 4)
text(30, -3.0, "plot", pos = 4)

legend(0.05,-2,c("Males","Females"), cex = 0.8, pch = c(24,21),
       lty=c(2,1),lwd=1.5,col=c("blue3","khaki3"),bty="n")


# 2015
plot(jitter(d15$TotDensity), d15$log_ratio, pch = d15$symb, ylab="Log ratio (2015)", 
     xlab="Planting density",bg=cRamp(d15$sr),ylim = c(-3,2.5), 
     lwd = 1)
xSeq  <- seq(0,48,1)
beta  <- lr_15_avg[,c("predictor","avg")]$avg
y_m   <- beta[1] + beta[2]*xSeq + beta[3]*0.5 + beta[4]
y_f   <- beta[1] + beta[2]*xSeq + beta[4]*0.5
lines(xSeq,y_f,lwd=1.5,lty=1,col="khaki3")
lines(xSeq,y_m,lwd=1.5,lty=2,col="blue3")

dev.off()
