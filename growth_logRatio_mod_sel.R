# Sex competition predictor: SEX RATIO (female individuals/total individuals)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)
library(lme4)
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/ac79/Documents/CODE/LLELA/model_sel_results.R")

# load and format data ------------------------------------------------------------------
d     <- format_growth( read.csv("Data/vr.csv") )
d14   <- subset(d, year == 2014)
d15   <- subset(d, year == 2015)

# NOTE: Planting density as density predictor
# 2014 --------------------------------------------------------------------------------
m14=list() # lMod stands for "leaf model" (density is quantified by N. of leaves)
m14[[1]]=lmer(log_ratio ~ (1 | plot),data=d14)
# single factors
m14[[2]]=lmer(log_ratio ~ sex + (1 | plot),data=d14)
m14[[3]]=lmer(log_ratio ~ sr + (1 | plot),data=d14)
m14[[4]]=lmer(log_ratio ~ TotDensity + (1 | plot),data=d14)
# Additive effects
m14[[5]]=lmer(log_ratio ~ sex + sr + (1 | plot),data=d14)
m14[[6]]=lmer(log_ratio ~ sr + TotDensity + (1 | plot),data=d14)
m14[[7]]=lmer(log_ratio ~ TotDensity + sex + (1 | plot),data=d14)
# two-way interactions
m14[[8]]=lmer(log_ratio ~ sex * sr  + (1 | plot),data=d14)
m14[[9]]=lmer(log_ratio ~ sr * TotDensity + (1 | plot),data=d14)
m14[[10]]=lmer(log_ratio ~ TotDensity * sex + (1 | plot),data=d14)
# two-way interaction, + remaining factor
m14[[11]]=lmer(log_ratio ~ sex * sr + TotDensity + (1 | plot),data=d14)
m14[[12]]=lmer(log_ratio ~ sr * TotDensity + sex + (1 | plot),data=d14)
m14[[13]]=lmer(log_ratio ~ TotDensity * sex + sr + (1 | plot),data=d14)
# three way interaction
m14[[14]]=lmer(log_ratio ~ sex + sr + TotDensity + (1 | plot),data=d14)
m14[[15]]=lmer(log_ratio ~ sr * sex * TotDensity + (1 | plot),data=d14)

lr_14_mod_sel   <- AICtab(m14, weights = T)
lr_14_avg       <- model_avg(lr_14_mod_sel, m14)
write.csv(lr_14_avg, "Results/VitalRates_3/log_ratio_14_best.csv", row.names = F)

# model selection results
mod_14_w        <- rgr_mod_sel(m14)
write.csv(mod_14_w,"Results/VitalRates_3/log_ratio_14_mod_sel.csv",row.names=F)


# 2015 --------------------------------------------------------------------------------
m15=list() # lMod stands for "leaf model" (density is quantified by N. of leaves)
m15[[1]]=lmer(log_ratio ~ (1 | plot),data=d15)
# single factors
m15[[2]]=lmer(log_ratio ~ sex + (1 | plot),data=d15)
m15[[3]]=lmer(log_ratio ~ sr + (1 | plot),data=d15)
m15[[4]]=lmer(log_ratio ~ TotDensity + (1 | plot),data=d15)
# Additive effects
m15[[5]]=lmer(log_ratio ~ sex + sr + (1 | plot),data=d15)
m15[[6]]=lmer(log_ratio ~ sr + TotDensity + (1 | plot),data=d15)
m15[[7]]=lmer(log_ratio ~ TotDensity + sex + (1 | plot),data=d15)
# two-way interactions
m15[[8]]=lmer(log_ratio ~ sex * sr  + (1 | plot),data=d15)
m15[[9]]=lmer(log_ratio ~ sr * TotDensity + (1 | plot),data=d15)
m15[[10]]=lmer(log_ratio ~ TotDensity * sex + (1 | plot),data=d15)
# two-way interaction, + remaining factor
m15[[11]]=lmer(log_ratio ~ sex * sr + TotDensity + (1 | plot),data=d15)
m15[[12]]=lmer(log_ratio ~ sr * TotDensity + sex + (1 | plot),data=d15)
m15[[13]]=lmer(log_ratio ~ TotDensity * sex + sr + (1 | plot),data=d15)
# three way interaction
m15[[14]]=lmer(log_ratio ~ sex + sr + TotDensity + (1 | plot),data=d15)
m15[[15]]=lmer(log_ratio ~ sr * sex * TotDensity + (1 | plot),data=d15)

lr_15_mod_sel  <- AICtab(m15, weights = T)
lr_15_avg      <- model_avg(lr_15_mod_sel, m15)
write.csv(lr_15_avg, "Results/VitalRates_3/log_ratio_15_best.csv", row.names = F)

# model selection results
mod_15_w       <- rgr_mod_sel(m15)
write.csv(mod_15_w,"Results/VitalRates_3/log_ratio_15_mod_sel.csv",row.names=F)


# GRAPH ------------------------------------------------------------------------------------------------

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
tiff("Results/VitalRates_3/log_ratio_growth_15.tiff",
     unit="in",width=6.3,height=4.5,res=600,compression="lzw")

par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,8.7))

# 2015
plot(jitter(d15$TotDensity), d15$log_ratio, pch = d15$symb, ylab="Relative growth rate (2015)", 
     xlab="Planting density",bg=cRamp(d15$sr),ylim=c(-3,2.5), lwd = 1)
xSeq  <- seq(0,48,1)
beta  <- lr_15_avg[,c("predictor","avg")]$avg
y_m   <- beta[1] + beta[2]*xSeq + beta[3]*0.5 + beta[4]
y_f   <- beta[1] + beta[2]*xSeq + beta[3]*0.5
lines(xSeq,y_f,lwd=1.5,lty=1,col="blue")
lines(xSeq,y_m,lwd=1.5,lty=2,col="red")

# legends
legend(53,2.8,c("Males","Females"), cex = 1, pch = c(24,21),
       lty=c(2,1),lwd=1.5,col=c("red","blue"),bty="n", xpd=NA)

colfunc <- colorRampPalette(cRamp(unique(arrange(d15,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(21), ncol=1))
text(x=53, y = seq(1.5,-1,l=3), labels = seq(1,0,l=3), xpd=NA)
rasterImage(legend_image, 55, 1.5, 60, -1, xpd=NA)
text(60, 1.5, "Percent of", pos = 4, xpd=NA)
text(60, 0.25, "females in", pos = 4, xpd=NA)
text(60, -1, "plot", pos = 4, xpd=NA)

dev.off()
