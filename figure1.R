# Code to produce figure1 (panels a-d)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(MASS)
library(dplyr)
library(boot)

#Read data--------------------------------------------------------------------
d           <- read.csv("Data/vr.csv")
fem_seeds   <- read.csv("Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")
m_spiklet   <- read.csv("Data/Spring 2014/maleCounts/malePaniculesSpring2014.csv")
viabVr      <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                        stringsAsFactors = F)

# best models
n_flow_beta <- read.csv("Results/VitalRates_3/n_flowers_best.csv")
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv")
m_spik_beta <- read.csv("Results/VitalRates_3/male_spikelets.csv")
tetr_beta   <- read.csv("Results/VitalRates_3/germination_best.csv")


# FORMAT DATA -------------------------------------------------------------------

# Number of flowers in 2014
d14         <- subset(d, year == 2014)
d14         <- subset(d14, surv_t1 != 0)
d14         <- d14 %>% mutate(new_t1 = as.numeric( as.character(new_t1)) )
f14         <- na.omit(select(d14,plot,flow_t1,log_l_t0,flowN_t1,
                              sex,sr,new_t1,TotDensity))


# fecundity data  
fem_seeds   <- mutate(fem_seeds, focalI = paste("f",IndividualN,sep=""))
fem_seeds   <- select(fem_seeds,Plot,focalI,SeedN)
fem_seeds   <- filter(fem_seeds, !is.na(SeedN) )
fem_seeds   <- rename(fem_seeds, plot = Plot)
fecund_data <- merge(d14,fem_seeds) 

# male allocation (spikelets) data  
m_alloc     <- merge(d14, m_spiklet)

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
# Set up colors for plots
f14$col <- as.character(factor(as.integer(f14$sex),labels=c("blue","red")))
mal14   <- subset(f14,sex=="m")
fem14   <- subset(f14,sex=="f")

# means
meanF   <- fem14 %>% group_by(TotDensity) %>% 
  summarise( meanF = mean(flowN_t1), sdF = sd(flowN_t1) )
meanM   <- mal14 %>% group_by(TotDensity) %>% 
  summarise( meanM = mean(flowN_t1), sdM = sd(flowN_t1) )

# plots
plot(meanF$TotDensity , meanF$meanF, pch = 16, col = "blue", ylim = c(0.3,2),
     ylab = "Number of flowers per capita", xlab="Planting density")
arrows(meanF$TotDensity, meanF$meanF - meanF$sdF*0.1, col = "blue",
       meanF$TotDensity, meanF$meanF + meanF$sdF*0.1, length=0.02, angle=90, code=3)
par(new = T)
plot(meanM$meanM ~ meanM$TotDensity , pch = 16, col = "red", ylim = c(0.3,2),
     ylab = "", xlab="")
arrows(meanM$TotDensity, meanM$meanM - meanM$sdM*0.1, col = "red",  
       meanM$TotDensity, meanM$meanM + meanM$sdM*0.1, length=0.02, angle=90, code=3)

size  <- mean(f14$log_l_t0)
xSeq  <- seq(0,48, by = 1)
beta  <- n_flow_beta[,c("predictor","avg")]$avg

y_m <- exp( beta[1] + beta[2]*size + beta[3]*xSeq + 
              beta[4]*xSeq + beta[5]*0.5 + beta[6] + 
              beta[7]*size + beta[8]*0.5)
y_f <- exp( beta[1] + beta[2]*size + beta[3]*xSeq + beta[5]*0.5)

lines(xSeq,y_m,lty=1,lwd=2,col="blue")
lines(xSeq,y_f,lty=1,lwd=2,col="red")

legend(1,2.1,c("male","female"), lty=1, lwd=2, 
       col=c("red","blue"), bty="n", pch = 16 )

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(a)", cex = 1.2, xpd = T)


# fecundity ----------------------------------------------------------------------------
plot(fecund_data$TotDensity, fecund_data$SeedN, pch = 21, 
     bg = cRamp(fecund_data$sr) , ylim = c(0, 1000), cex = 1.5,
     xlab = "Planting density", ylab = "Seeds per flower")
xSeq   <- seq(0,48,by = 1)
mSize  <- mean(fecund_data$log_l_t0)
beta   <- fec_beta$avg
y_low  <- exp(beta[1] + beta[2]*mSize + beta[3]*0.1 + 
                beta[4]*xSeq + beta[5]*xSeq*0.1)
y_high <- exp(beta[1] + beta[2]*mSize + beta[3]*1 + 
                beta[4]*xSeq + beta[5]*xSeq*1)

lines(xSeq, y_low, lwd = 2, lty = 2)
lines(xSeq, y_high, lwd = 2)

legend(14.5,830,c("100% female plot","  20% female plot"),lty = c(1,2), lwd = 2, bty="n")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(b)", cex = 1.2, xpd = T)

colfunc = colorRampPalette(cRamp(unique(arrange(fecund_data,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(19), ncol=1))
text(x=24.5, y = seq(850,1000,l=3), labels = seq(0,1,l=3))
rasterImage(legend_image, 17, 850, 22, 1000)
text(25, 1000, "Percent of", pos = 4)
text(25, 925, "males in", pos = 4)
text(25, 850, "plot", pos = 4)


# male allocation ----------------------------------------------------------------------
plot(m_alloc$TotDensity,m_alloc$CountSpikelets,pch=21, 
     bg = cRamp(m_alloc$sr), cex = 1.5,
     ylab="Spikelet number per male flower",xlab="Planting density")
size <- mean(m_alloc$log_l_t0)
beta <- m_spik_beta$avg
xSeq <- seq(1,48,by =1)
y_l  <- exp(beta[1] + beta[2]*size + beta[3]*xSeq + beta[4]*0.1 + beta[5]*xSeq*0.1)
y_h  <- exp(beta[1] + beta[2]*size + beta[3]*xSeq + beta[4]*0.9 + beta[5]*xSeq*0.9)

lines(xSeq,y_l,lwd=2,lty=2)
lines(xSeq,y_h,lwd=2,lty=1)

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(c)", cex = 1.2, xpd = T)


# viability ----------------------------------------------------------------------------
lower=quantile(viabVr$totFlow,prob=c(0.1))
upper=quantile(viabVr$totFlow,prob=c(0.9))

# Viability (germination) data
plot(jitter(viabVr$totFlow,factor = 2),jitter(viabVr$germ_ratio,factor = 2),pch=21,ylim=c(0,1.01), 
     bg = cRamp(viabVr$sr_f), cex = 1.5,
     xlab="Number of flowers",ylab="Seed viability rate")

xSeq=seq(min(viabVr$totFlow),max(viabVr$totFlow),length.out=100)
beta=tetr_beta$avg
yMeanLow=inv.logit(beta[1] + beta[2]*0.2 + beta[3]*xSeq + beta[4]*xSeq*0.2)
yMeanHigh=inv.logit(beta[1] + beta[2]*1 + beta[3]*xSeq + beta[4]*xSeq*1)

lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(d)", cex = 1.2, xpd = T)

dev.off()

