# Code to produce figure1 (panels a,b,c)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(MASS)
library(dplyr)
library(boot)

#Read in data--------------------------------------------------------------------

# load and format data ------------------------------------------------------------------
d=read.csv("Data/vr.csv")
# remove dead individuals (this is a GROWTH model!)
d=subset(d, surv_t1 != 0)
# Year one
d14 <- subset(d, year==2014)
# Year two
tmp15=subset(d,year==2015)
tmp15$new_t1[tmp15$new_t1=="SKIPPED"]=NA
tmp15$new_t1[tmp15$new_t1=="cnf"]=NA
tmp15$new_t1[tmp15$new_t1==""]=NA
tmp15$new_t1=as.numeric(as.character(tmp15$new_t1))
d15=na.omit(tmp15[,c("l_t1","log_l_t0","plot","sex","new_t1","sr","TotDensity")])

# Best models
gr14_avg   <- read.csv("Results/VitalRates_3/growth_14N_best.csv")
grow_avgN  <- read.csv("Results/VitalRates_3/growthN_best.csv")


# FIGURE 3 ----------------------------------------------------------------------


# Sex ratio as dot color ------------------------------------------------------------

# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(topo.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

# Graph
tiff("Results/VitalRates_3/Figure3.tiff",unit="in",width=6.3,height=3.5,res=600,compression="lzw")

par(mfrow=c(1,2),mar=c(2.5,2.5,1.3,0.5),mgp=c(1.4,0.35,0),cex.lab=1,cex.axis=0.8,
    cex.main=0.9,xpd=NA)

# 2014
plot(d14$TotDensity, d14$l_t1, pch=21, ylab="Number of leaves in 2014", 
     xlab="Planting density",bg=cRamp(d14$sr), ylim = c(0,99), 
     lwd = 1)
beta <- gr14_avg[,c("predictor","avg")]$avg
xSeq <- seq(0,48,by=1)
sr   <- 0.5
size <- mean(d14$log_l_t0)
y_m <- exp( beta[1] + beta[2]*size + beta[3] + beta[4]*xSeq + 
              beta[5]*xSeq + beta[6]*sr + beta[7]*size)
y_f <- exp( beta[1] + beta[2]*size + beta[4]*xSeq + beta[6]*sr )
lines(xSeq,y_f,lwd=1.5,lty=1,col="blue")
lines(xSeq,y_m,lwd=1.5,lty=1,col="red")
text(-2,107.5,"a) Target individuals, year 2014",pos=4)

legend(29,105,c("Males","Females"), cex = 0.8,
       lty=1,lwd=1.5,col=c("red","blue"),bty="n")

# 2015
plot(d15$TotDensity,d15$l_t1, pch=21, ylab="Number of leaves in 2015",
     xlab="Planting density",bg=cRamp(d15$sr), ylim = c(0,99),
     lwd = 1)
beta <- grow_avgN[,c("predictor","avg")]$avg
xSeq <- seq(0,48,by=1)
size <- mean(d15$log_l_t0)
y_m <- exp( beta[1] + beta[2]*size + beta[3] + beta[4]*0.5 + 
              beta[5]*xSeq + beta[6]*size + beta[7]*0.5 + beta[8]*xSeq )
y_f <- exp( beta[1] + beta[2]*size + beta[4]*0.5 + beta[5]*xSeq )
lines(xSeq,y_f,lwd=1.5,lty=1,col="blue")
lines(xSeq,y_m,lwd=1.5,lty=1,col="red")
text(-2,107.5,"b) Target individuals, year 2015", pos=4)

colfunc = colorRampPalette(cRamp(unique(arrange(d15,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(21), ncol=1))
text(x=23, y = seq(78,98,l=3), labels = seq(0,1,l=3))
rasterImage(legend_image, 15, 78, 20, 98)
text(25, 98, "Percent of", pos = 4)
text(25, 88, "males in", pos = 4)
text(25, 78, "plot", pos = 4)

dev.off()

