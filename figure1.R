# Code to produce figure1 (panels a-d)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(dplyr)
library(boot)
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/ac79/Documents/CODE/LLELA/unload_mass.R")
unload_mass()

#Read data--------------------------------------------------------------------
d           <- read.csv("Data/vr.csv", stringsAsFactors = F)
fem_seeds   <- read.csv("Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")
m_spiklet   <- read.csv("Data/Spring 2014/maleCounts/malePaniculesSpring2014.csv")

# best models
n_flow_beta <- read.csv("Results/VitalRates_3/n_flowers_best.csv")
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv")
m_spik_beta <- read.csv("Results/VitalRates_3/male_spikelets.csv")
avg14       <- read.csv("Results/VitalRates_3/new_t_bh14_best.csv")


# FORMAT DATA -------------------------------------------------------------------

# Number of flowers in 2014
f14   <- format_flower(d)
f14   <- f14 %>% mutate( plot = as.numeric(as.character(plot)) )

# fecundity data  
fem_seeds   <- mutate(fem_seeds, focalI = paste("f",IndividualN,sep=""))
fem_seeds   <- dplyr::select(fem_seeds,Plot,focalI,SeedN)
fem_seeds   <- subset(fem_seeds, !is.na(SeedN) )
fem_seeds   <- mutate(fem_seeds, plot = Plot)
fecund_data <- merge(fem_seeds, select(f14,plot,sr, TotDensity))
fecund_data <- fecund_data %>% unique()

# male allocation (spikelets) data  
m_alloc     <- merge(select(f14,plot,sr,TotDensity), m_spiklet)
m_alloc     <- m_alloc %>% unique()

# tiller production data
d       <- format_growth(d)
d14     <- subset(d, year == 2014)


# FIGURE 1 -----------------------------------------------------------------------------

# Sex ratio as dot color ------
# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(gray.colors(7,start=0,end=1))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

# Graph
tiff("Results/VitalRates_3/figure1.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

par(mfrow=c(2,2),mar=c(3,2.5,0.1,0.1),mgp=c(1.4,0.35,0),cex.lab=1.1,cex.axis=0.8,
    cex.main=0.9, oma=c(0,0,0.2,0))

# Flowering ------------------------------------------------------------------------------
# Set up colors for plots
f14$col <- as.character(factor(f14$sex,labels=c("blue","red")))
mal14   <- subset(f14,sex=="m")
fem14   <- subset(f14,sex=="f")

# means
meanF   <- fem14 %>% 
            group_by(TotDensity) %>% 
            summarise( meanF = mean(flowN_t1), sdF = sd(flowN_t1) )
meanM   <- mal14 %>% 
            group_by(TotDensity) %>% 
            summarise( meanM = mean(flowN_t1), sdM = sd(flowN_t1) )

# plots
plot(meanF$TotDensity , meanF$meanF, pch = 16, col = "#ABABAB", ylim = c(0.3,2),
     ylab = "Number of flowers per capita", xlab="Planting density")
arrows(meanF$TotDensity, meanF$meanF - meanF$sdF*0.1, col = "#ABABAB",
       meanF$TotDensity, meanF$meanF + meanF$sdF*0.1, length=0.02, angle=90, code=3)
par(new = T)
plot(meanM$meanM ~ meanM$TotDensity , pch = 16, col = "black", ylim = c(0.3,2),
     ylab = "", xlab="")
arrows(meanM$TotDensity, meanM$meanM - meanM$sdM*0.1, col = "black",  
       meanM$TotDensity, meanM$meanM + meanM$sdM*0.1, length=0.02, angle=90, code=3)

xSeq  <- seq(0,48, by = 1)
beta  <- n_flow_beta[,c("predictor","avg")]$avg
y_m <- exp(beta[1] + beta[2] + beta[3]*0.5 + beta[4]*xSeq + 
           beta[5]*0.5 + beta[6]*0.5*xSeq + beta[7]*xSeq)
y_f <- exp(beta[1] + beta[3]*0.5 + beta[4]*xSeq + beta[6]*0.5*xSeq)
lines(xSeq,y_m,lty=1,lwd=2,col="black")
lines(xSeq,y_f,lty=2,lwd=2,col="#ABABAB")

legend(1,2.1,c("male","female"), lty=c(1,2), lwd=2, 
       col=c("black","#ABABAB"), bty="n", pch = 16 )

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(a)", cex = 1.2, xpd = T)


# fecundity ----------------------------------------------------------------------------
plot(fecund_data$TotDensity, fecund_data$SeedN, pch = 21, 
     bg = cRamp(fecund_data$sr) , ylim = c(0, 1000), cex = 1.5,
     xlab = "Planting density", ylab = "Seeds per flower")
xSeq   <- seq(0,48,by = 1)
beta   <- fec_beta$avg
y_low  <- exp(beta[1] + beta[3]*0.05 + beta[2]*xSeq + beta[4]*xSeq*0.05)
y_high <- exp(beta[1] + beta[3]*0.95 + beta[2]*xSeq + beta[4]*xSeq*0.95)

lines(xSeq, y_low, lwd = 2, lty = 2)
lines(xSeq, y_high, lwd = 2)

legend(14.5,830,c("95% female plot","  5% female plot"),lty = c(1,2), lwd = 2, bty="n")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(b)", cex = 1.2, xpd = T)

colfunc = colorRampPalette(cRamp(unique(arrange(fecund_data,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(19), ncol=1))
text(x=24.5, y = seq(850,1000,l=3), labels = seq(1,0,l=3))
rasterImage(legend_image, 17, 850, 22, 1000)
text(25, 1000, "Percent of", pos = 4)
text(25, 925, "females in", pos = 4)
text(25, 850, "plot", pos = 4)


# male allocation ----------------------------------------------------------------------
plot(m_alloc$TotDensity,m_alloc$CountSpikelets,pch=21, 
     bg = cRamp(m_alloc$sr), cex = 1.5,
     ylab="Spikelet number per male flower",xlab="Planting density")
beta <- m_spik_beta$avg
xSeq <- seq(1,48,by =1)
y_l  <- exp(beta[1] + beta[2]*xSeq + beta[3]*0.05 + beta[4]*xSeq*0.2)
y_h  <- exp(beta[1] + beta[2]*xSeq + beta[3]*0.95 + beta[4]*xSeq*0.95)

lines(xSeq,y_l,lwd=2,lty=2)
lines(xSeq,y_h,lwd=2,lty=1)

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(c)", cex = 1.2, xpd = T)


# total number of tillers 2014 --------------------------------------------------------
plot(d14$TotDensity,d14$new_t1,pch=21,ylab="Number of new tillers",
     xlab="Planting density", bg = cRamp(d14$sr), cex = 1.5)
N    <- seq(0,48,1)
beta <- as.data.frame(avg14)
fem  <- N*0.95
mal  <- N*0.05
y_h  <- (beta[,"lam.f"]*fem) / (1 + beta[,"b.f"] * (fem + beta[,"a."]*mal)) + 
        (beta[,"lam.m"]*mal) / (1 + beta[,"b.m"] * (fem + beta[,"a."]*mal))
fem  <- N*0.05
mal  <- N*0.95
y_l  <- (beta[,"lam.f"]*fem) / (1 + beta[,"b.f"] * (fem + beta[,"a."]*mal)) + 
        (beta[,"lam.m"]*mal) / (1 + beta[,"b.m"] * (fem + beta[,"a."]*mal))
lines(N,y_h,lwd=2, lty = 1) #col="#DCDCDC",
lines(N,y_l,lwd=2, lty = 2) #col="#636363",

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(d)", cex = 1.2, xpd = T)


dev.off()
