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
lr_14_avg   <- read.csv("Results/VitalRates_3/log_ratio_14_best.csv")
n_flow_beta <- read.csv("Results/VitalRates_3/n_flowers_best.csv")
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv")
m_spik_beta <- read.csv("Results/VitalRates_3/male_spikelets.csv")
avg14       <- read.csv("Results/VitalRates_3/new_t_bh14_best.csv")

# FORMAT DATA -------------------------------------------------------------------

# Relative growth rate (log ratio)
rgr         <- format_growth(d)
rgr14       <- filter(rgr, year == 2014)

# Number of flowers in 2014
f14         <- format_flower(d)
f14         <- f14 %>% mutate( plot = as.numeric(as.character(plot)) )

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
till        <- format_new_tillers(d)
till14      <- subset(till, year == 2014)


# FIGURE 1 -----------------------------------------------------------------------------

# Sex ratio as dot color ------
# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(gray.colors(7,start=0,end=1))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

# Symbols for the sexes 
sex_symbol <- function(x){
  out <- x %>% 
    mutate(symb = factor(sex, labels=c("21","24")) ) %>%
    mutate(symb = as.character(symb) ) %>%
    mutate(symb = as.integer(symb) )
  return(out)
}


# Graph -----------------------------------------------------------------------------------

# parameters
cex_dots <- 1.5
cex_pan  <- 1.3
cex_leg  <- 1.3


tiff("Results/VitalRates_3/figure1.tiff",unit="in", width=6.3,height=8.5,res=600,compression="lzw")

par(mfrow=c(3,2),mar=c(3,3,0.1,0.4),mgp=c(1.6,0.5,0),cex.lab=1.4,cex.axis=1,
    oma=c(0,0,0.2,0))

# relative growth rate 
rgr14  <- sex_symbol(rgr14)
  
plot(jitter(rgr14$TotDensity, factor = 2), rgr14$log_ratio, pch = rgr14$symb, 
     ylab="Relative growth rate", cex = cex_dots, xlim = c(0,48.5),
     xlab="Planting density",bg=cRamp(rgr14$sr), ylim = c(-3,2.5), 
     lwd = 1)
beta  <- lr_14_avg[,c("predictor","avg")]$avg

xSeq  <- seq(1,48,1)
y_h   <- beta[1] + beta[2] + beta[3]*0.95
y_l   <- beta[1] + beta[2] + beta[3]*0.05
abline(h= y_l, lwd = 2, lty = 2)
abline(h= y_h, lwd = 2)


colfunc      <- colorRampPalette(cRamp(unique(arrange(rgr14,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(21), ncol=1))
text(x=30, y = seq(-2.75,-1.25,l=3), labels = seq(1,0,l=3), cex = 1.3)
rasterImage(legend_image, 32.5, -2.75, 36.5, -1.25)
# text(36, -1.25, "Percent of", pos = 4, cex = 1.3)
# text(36, -2, "females in", pos = 4, cex = 1.3)
# text(36, -2.75, "plot", pos = 4, cex = 1.3)
text(36, -1.25, "Sex ratio", pos = 4, cex = 1.2)
text(35.9, -2, "(proportion", pos = 4, cex = 1.2)
text(36, -2.75, "female)", pos = 4, cex = 1.2)

legend("bottomleft", c("95% female plot","  5% female plot"),
       col = c("#ABABAB", "#141414"),
       lty = c(1,2), lwd = 2, cex = cex_leg, bty="n")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.04,
     par("usr")[4]*0.94,"(a)", cex = cex_pan, xpd = T)


# Panicles per capita 
f14  <- sex_symbol(f14)

plot(jitter(f14$TotDensity), jitter(f14$flowN_t1, factor = 3), 
     pch = f14$symb, bg=cRamp(f14$sr), cex = cex_dots, xlim = c(0,48.5),
     ylab = "Number of panicles per individual", xlab="Planting density")

xSeq  <- seq(0,48, by = 1)
beta  <- n_flow_beta[,c("predictor","avg")]$avg

y_h <- exp(beta[1] + beta[2] + beta[3]*0.95 + beta[4]*xSeq + 
             beta[5]*0.95 + beta[6]*0.95*xSeq + beta[7]*xSeq)
y_l <- exp(beta[1] + beta[2] + beta[3]*0.05 + beta[4]*xSeq + 
             beta[5]*0.05 + beta[6]*0.05*xSeq + beta[7]*xSeq)
lines(xSeq,y_h,lty=1,lwd=2,col="#ABABAB")
lines(xSeq,y_l,lty=2,lwd=2,col="#141414")


text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.04,
     par("usr")[4]*0.97,"(b)", cex = cex_pan, xpd = T)


# fecundity 
plot(fecund_data$TotDensity, fecund_data$SeedN, pch = 21, xlim = c(0,48.5),
     bg = cRamp(fecund_data$sr) , ylim = c(0, 1000), cex = cex_dots,
     xlab = "Planting density", ylab = "Seeds per female panicle")
xSeq   <- seq(0,48,by = 1)
beta   <- fec_beta$avg
y_low  <- exp(beta[1] + beta[3]*0.05 + beta[2]*xSeq + beta[4]*xSeq*0.05)
y_high <- exp(beta[1] + beta[3]*0.95 + beta[2]*xSeq + beta[4]*xSeq*0.95)

lines(xSeq, y_low, lwd = 2, lty = 2, col="#141414")
lines(xSeq, y_high, lwd = 2, col="#ABABAB")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(c)", cex = cex_pan, xpd = T)


# male allocation 
plot(m_alloc$TotDensity,m_alloc$CountSpikelets,pch=21, 
     bg = cRamp(m_alloc$sr), cex = cex_dots, xlim = c(0,48.5),
     ylab="Spikelets per male panicle",xlab="Planting density")
beta <- m_spik_beta$avg
xSeq <- seq(1,48,by =1)
y_l  <- exp(beta[1] + beta[2]*xSeq + beta[3]*0.05 + beta[4]*xSeq*0.2)
y_h  <- exp(beta[1] + beta[2]*xSeq + beta[3]*0.95 + beta[4]*xSeq*0.95)

lines(xSeq,y_l,lwd=2,lty=2, col="#141414")
lines(xSeq,y_h,lwd=2,lty=1, col="#ABABAB")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(d)", cex = cex_pan, xpd = T)


# total number of tillers 2014 
plot(till14$TotDensity,till14$new_t1,pch=21,ylab="Asexual recruitment",
     xlab="Planting density", bg = cRamp(till14$sr), cex = cex_dots,
     xlim = c(0,48.5))
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
lines(N,y_h,lwd=2, lty = 1, col="#ABABAB") #col="#DCDCDC",
lines(N,y_l,lwd=2, lty = 2, col="#141414") #col="#63 6363",

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(e)", cex = cex_pan, xpd = T)


dev.off()
