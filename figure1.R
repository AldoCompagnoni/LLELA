# Code to produce figure1 (panels a-d)
setwd("C:/cloud/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(dplyr)
library(boot)
source("C:/CODE/LLELA/model_avg.R")
source("C:/CODE/LLELA/unload_mass.R")
source("C:/CODE/LLELA/prediction.R")
unload_mass()

#Read data--------------------------------------------------------------------
d           <- read.csv("Data/vr.csv", stringsAsFactors = F)
fem_seeds   <- read.csv("Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")
m_spiklet   <- read.csv("Data/Spring 2014/maleCounts/malePaniculesSpring2014.csv")
viabVr      <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                        stringsAsFactors = F)

# best models
germ_beta   <- read.csv("Results/VitalRates_4_Cade2015/germination_best.csv")
lr_14_avg   <- read.csv("Results/VitalRates_4_Cade2015/log_ratio_14_best.csv")
n_flow_beta <- read.csv("Results/VitalRates_4_Cade2015/n_flowers_best.csv")
fec_beta    <- read.csv("Results/VitalRates_4_Cade2015/fecuntity_best.csv")
m_spik_beta <- read.csv("Results/VitalRates_4_Cade2015/male_spikelets.csv")
avg14       <- read.csv("Results/VitalRates_4_Cade2015/new_t_bh14_best.csv")

# FORMAT DATA -------------------------------------------------------------------

# remove three plots with more than 60 flowers
viabVr      <- subset(viabVr, totFlow < 60)   

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

# 6-panel figure
tiff("Results/VitalRates_4_Cade2015/figure1.tiff",unit="in", width=6.3,height=8.5,res=600,compression="lzw")

par(mfrow=c(3,2),mar=c(3,3,0.1,0.4),mgp=c(1.6,0.5,0),cex.lab=1.4,cex.axis=1,
    oma=c(0,0,0.2,0))

# relative growth rate
rgr14  <- sex_symbol(rgr14)
plot(jitter(rgr14$TotDensity, factor = 2), rgr14$log_ratio, pch = rgr14$symb,
     ylab="Growth rate", cex = cex_dots, xlim = c(0,48.5),
     xlab="Planting density",bg=cRamp(rgr14$sr), ylim = c(-3,2.5),
     lwd = 1)

# predictions across male/females
lr_14 <- lr_14_avg %>%
            group_by(TotDensity, sr) %>%
            summarise( pred = mean(pred) ) %>%
            ungroup
rgr_h <- subset(lr_14, sr == 0.95)
rgr_l <- subset(lr_14, sr == 0.05)
lines(rgr_h$TotDensity, rgr_h$pred, col = "#ABABAB", lty = 1, lwd = 2 )
lines(rgr_l$TotDensity, rgr_l$pred, col = "#141414", lty = 3, lwd = 2 )

legend(-3,-2, c("95% female plot","  5% female plot"),
       col = c("#ABABAB", "#141414"),
       lty = c(1,3), lwd = 2, cex = cex_leg, bty="n")

colfunc      <- colorRampPalette(cRamp(unique(arrange(rgr14,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(21), ncol=1))
text(x=30, y = seq(-2.6,-1.4,l=3), labels = seq(1,0,l=3), cex = 1.3)
rasterImage(legend_image, 32.5, -2.6, 36.5, -1.4)
text(36, -1.4, "Sex ratio", pos = 4, cex = 1.2)
text(35.9, -2, "(proportion", pos = 4, cex = 1.2)
text(36, -2.6, "female)", pos = 4, cex = 1.2)

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.04,
     par("usr")[4]*0.94,"(a)", cex = cex_pan, xpd = T)


# Panicles per capita 
f14  <- sex_symbol(f14)
plot(jitter(f14$TotDensity), jitter(f14$flowN_t1, factor = 3), 
     pch = f14$symb, bg=cRamp(f14$sr), cex = cex_dots, xlim = c(0,48.5),
     ylab = "Number of panicles per individual", xlab="Planting density")

# predictions across male/females
pan <- n_flow_beta %>%
            group_by(TotDensity, sr) %>%
            summarise( pred = mean(pred) ) %>%
            ungroup

pan_h <- subset(pan, sr == 0.95)
pan_l <- subset(pan, sr == 0.05)
lines(pan_h$TotDensity,pan_h$pred,lty=1,lwd=2,col="#ABABAB")
lines(pan_l$TotDensity,pan_l$pred,lty=3,lwd=2,col="#141414")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.04,
     par("usr")[4]*0.97,"(b)", cex = cex_pan, xpd = T)


# fecundity 
plot(jitter(fecund_data$TotDensity), fecund_data$SeedN, pch = 21, xlim = c(0,48.5),
     bg = cRamp(fecund_data$sr) , ylim = c(0, 1000), cex = cex_dots,
     xlab = "Planting density", ylab = "Seeds per female panicle")

fec_beta_f <- subset(fec_beta, sex == "f") %>% select(-sex)
fec_beta_m <- subset(fec_beta, sex == "m") %>% select(-sex)

# predictions across male/females (again, same predictions)
fec_beta <- fec_beta %>%
                group_by(TotDensity, sr) %>%
                summarise( pred = mean(pred) ) %>%
                ungroup

fec_h <- subset(fec_beta, sr == 0.95)
fec_l <- subset(fec_beta, sr == 0.05)
lines(fec_h$TotDensity,fec_h$pred,lty=1,lwd=2,col="#ABABAB")
lines(fec_l$TotDensity,fec_l$pred,lty=3,lwd=2,col="#141414")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(c)", cex = cex_pan, xpd = T)


# male allocation
plot(jitter(m_alloc$TotDensity),m_alloc$CountSpikelets,pch=21,
     bg = cRamp(m_alloc$sr), cex = cex_dots, xlim = c(0,48.5),
     ylab="Spikelets per male panicle",xlab="Planting density")

spik_l  <- subset(m_spik_beta, sr == 0.05) 
spik_h  <- subset(m_spik_beta, sr == 0.95)

lines(spik_h$TotDensity, spik_h$pred, lwd=2, lty=1, col="#ABABAB")
lines(spik_l$TotDensity, spik_l$pred, lwd=2, lty=3, col="#141414")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(d)", cex = cex_pan, xpd = T)


# total number of tillers 2014 
plot(jitter(till14$TotDensity),till14$new_t1,pch=21,ylab="Asexual recruitment",
     xlab="Planting density", bg = cRamp(till14$sr), cex = cex_dots,
     xlim = c(0,48.5))

till_h <- select(avg14, TotN, sr, pred) %>% subset(sr == 0.95)
till_l <- select(avg14, TotN, sr, pred) %>% subset(sr == 0.05)
lines(till_h$TotN,till_h$pred,lwd=2, lty = 1, col="#ABABAB") #col="#DCDCDC",
lines(till_l$TotN,till_l$pred,lwd=2, lty = 3, col="#141414") #col="#63 6363",

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(e)", cex = cex_pan, xpd = T)

dev.off()
