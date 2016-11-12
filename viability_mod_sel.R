# Model selection for seed VIABILITY
# Data: tetrazolium scoring and germination data 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(lme4) ; library(bbmle) ; library(boot) ; library(testthat)
library(dplyr) ; library(glmmADMB)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")

# Read in data and format -----------------------------------------------------------
viabVr      <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                        stringsAsFactors = F)
viabVr$totN <- viabVr$F + viabVr$M
viabVr$sr   <- viabVr$F / viabVr$totN
viabVr$plot <- as.factor(viabVr$plot)
viabVr      <- subset(viabVr, totFlow < 60) #exclude extreme values

# Viability model selection ---------------------------------------------------------

# 1. viable Seed Number: tetrazolium assays  (Yes / fail) ---------------------------
# This is the "standard" scoring: consider as viable any seed stained by tetrazolium

# Omit NAs for glmmadmb
tetr_dat  <- na.omit(dplyr::select(viabVr,Yes,fail,sr_f,totFlow,sr,totN,plot))

# Number of flowers and their sex ratio as predictors
tetr_flowN=list()
tetr_flowN[[1]]  <- glmmadmb(cbind(Yes,fail) ~ sr_f + (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN[[2]]  <- glmmadmb(cbind(Yes,fail) ~ totFlow + (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN[[3]]  <- glmmadmb(cbind(Yes,fail) ~ sr_f + totFlow + (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN[[4]]  <- glmmadmb(cbind(Yes,fail) ~ sr_f * totFlow + (1 | plot),family="binomial", data=tetr_dat)
tetr_flowN_select<- AICtab(tetr_flowN,weights=T)

# Planting density and sex ratio as predictors
tetr_dens=list()
tetr_dens[[1]]    <- glmmadmb(cbind(Yes,fail) ~ sr + (1 | plot),family="binomial", data=tetr_dat)
tetr_dens[[2]]    <- glmmadmb(cbind(Yes,fail) ~ totN + (1 | plot),family="binomial", data=tetr_dat)
tetr_dens[[3]]    <- glmmadmb(cbind(Yes,fail) ~ sr + totN + (1 | plot),family="binomial", data=tetr_dat)
tetr_dens[[4]]    <- glmmadmb(cbind(Yes,fail) ~ sr * totN + (1 | plot),family="binomial", data=tetr_dat)
tetr_dens_select  <- AICtab(tetr_dens,weights=T)

# Model 4 has ~100% support
tetr_flowN_avg  <- model_avg(tetr_flowN_select, tetr_flowN)
tetr_dens_avg   <- model_avg(tetr_dens_select, tetr_dens)

write.csv(tetr_flowN_avg, "Results/VitalRates_3/tetrazolium_best.csv", row.names = F)
write.csv(tetr_dens_avg, "Results/VitalRates_3/tetrazolium_dens_best.csv", row.names = F)


# 2. Viable Seed Number: germination assays (germTot / germFail)---------------------
# Number of flowers and their sex ratio as predictors

# Omit NAs for glmmadmb
germ_dat  <- na.omit(dplyr::select(viabVr,germTot,germFail,sr_f,totFlow,sr,totN,plot))


germ_flowN      <- list()
germ_flowN[[1]] <- glmmadmb(cbind(germTot,germFail) ~ sr_f + (1 | plot),family="binomial", data=germ_dat)
germ_flowN[[2]] <- glmmadmb(cbind(germTot,germFail) ~ totFlow + (1 | plot),family="binomial", data=germ_dat)
germ_flowN[[3]] <- glmmadmb(cbind(germTot,germFail) ~ sr_f + totFlow + (1 | plot),family="binomial", data=germ_dat)
germ_flowN[[4]] <- glmmadmb(cbind(germTot,germFail) ~ sr_f * totFlow + (1 | plot),family="binomial", data=germ_dat)
germ_flowN_sel  <- AICtab(germ_flowN,weights=T) 

# Planting density and sex ratio as predictors
germ_dens       <- list()
germ_dens[[1]]  <- glmmadmb(cbind(germTot,germFail) ~ sr + (1 | plot),family="binomial", data=germ_dat)
germ_dens[[2]]  <- glmmadmb(cbind(germTot,germFail) ~ totN + (1 | plot),family="binomial", data=germ_dat)
germ_dens[[3]]  <- glmmadmb(cbind(germTot,germFail) ~ sr + totN + (1 | plot),family="binomial", data=germ_dat)
germ_dens[[4]]  <- glmmadmb(cbind(germTot,germFail) ~ sr * totN + (1 | plot),family="binomial", data=germ_dat)
germ_dens_sel   <- AICtab(germ_dens,weights=T) 

# Average best two models
germ_flowN_avg  <- model_avg(germ_flowN_sel, germ_flowN)
germ_dens_avg   <- model_avg(germ_dens_sel, germ_dens)

write.csv(germ_flowN_avg, "Results/VitalRates_3/germination_best.csv", row.names = F)
write.csv(germ_dens_avg, "Results/VitalRates_3/germination_dens_best.csv", row.names = F)


# Germination + tetrazolium ----------------------------------------------------
# Format tetrazolium data
tetr_dat <- dplyr::select(viabVr,Yes,fail,sr_f,totFlow,sr,totN,plot)
tetr_dat <- tetr_dat %>% mutate(yes = Yes, no = fail, type = 1) %>%
              dplyr::select(sr_f, totFlow, plot, yes, no, type)
# Format germination data
germ_dat <- dplyr::select(viabVr,germTot,germFail,sr_f,totFlow,sr,totN,plot)
germ_dat <- germ_dat %>% mutate(yes = germTot, no = germFail, type = 2) %>%
              dplyr::select(sr_f, totFlow, plot, yes, no, type)
# Combine tetrazolium and germination data
all_data <- na.omit( rbind(tetr_dat, germ_dat) )
all_data <- mutate( all_data, type = as.factor(type) , germ_ratio = yes/(yes + no) )

## Fit models
all_dens       <- list()
all_dens[[1]]  <- glmmadmb(cbind(yes, no) ~ sr_f + type + (1 | plot),family="binomial", data=all_data)
all_dens[[2]]  <- glmmadmb(cbind(yes, no) ~ totFlow + type + (1 | plot),family="binomial", data=all_data)
all_dens[[3]]  <- glmmadmb(cbind(yes, no) ~ sr_f + totFlow + type + (1 | plot),family="binomial", data=all_data)
all_dens[[4]]  <- glmmadmb(cbind(yes, no) ~ sr_f * totFlow + type + (1 | plot),family="binomial", data=all_data)

# Average best two models
all_dens_sel   <- AICtab(all_dens, weights=T) 
all_avg  <- model_avg(all_dens_sel, all_dens)
write.csv(all_avg, "Results/VitalRates_3/all_viability_best.csv", row.names = F)


# Graphs ---------------------------------------------------------------------

# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(topo.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

# graph
tiff("Results/VitalRates_3/viability.tiff",unit="in",width=6.3,height=6.3,
     res=600,compression="lzw")

lower=quantile(viabVr$totN,prob=c(0.1))
upper=quantile(viabVr$totN,prob=c(0.9))

par(mfrow=c(2,2),mar=c(2.5,2.5,1.1,0.1),mgp=c(1.5,0.6,0),
    oma=c(0,0,0,0),xpd=NA)
titlePlace=0.2

# Tetrazolium data
plot(jitter(viabVr$totFlow, factor = 2), jitter(viabVr$tetra_ratio, factor = 2), 
     pch=21, ylim = c(0,1), bg = cRamp(viabVr$sr), cex = 1.5, 
     xlab="Planting density",ylab="Seed viability rate")
title("Tetrazolium data", line = titlePlace)
xSeq=seq(0,54,by = 1)
beta=tetr_flowN_avg$avg
yMeanLow=inv.logit(beta[1] + beta[2]*0.1 + beta[3]*xSeq + beta[4]*xSeq*0.1)
yMeanHigh=inv.logit(beta[1] + beta[2]*1 + beta[3]*xSeq + beta[4]*xSeq*1)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

# Germination data
plot(jitter(viabVr$totFlow, factor = 2),jitter(viabVr$germ_ratio, factor = 2),
     pch=21, ylim = c(0,1), bg = cRamp(viabVr$sr), cex = 1.5, 
     xlab="Planting density",ylab="Seed germination rate")
title("Germination data", line = titlePlace)
xSeq=seq(0,54,by = 1)
beta=germ_flowN_avg$avg
yMeanLow=inv.logit(beta[1] + beta[2]*0.1 + beta[3]*xSeq + beta[4]*xSeq*0.1)
yMeanHigh=inv.logit(beta[1] + beta[2]*1 + beta[3]*xSeq + beta[4]*xSeq*1)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

# All data
plot(jitter(all_data$totFlow, factor = 2),jitter(all_data$germ_ratio, factor = 2),
     pch=21, ylim = c(0,1), bg = cRamp(all_data$sr), cex = 1.5, 
     xlab="Planting density",ylab="Seed viability/germination rate")
title("All data", line = titlePlace)
xSeq=seq(0,54,by = 1)
beta=all_avg$avg
yMeanLow=inv.logit(beta[1] + beta[2]*0.1 + beta[3]*xSeq*0.1 + beta[4]*xSeq)
yMeanHigh=inv.logit(beta[1] + beta[2]*1 + beta[3]*xSeq*1 + beta[4]*xSeq)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)

# plot for legends
plot(jitter(all_data$totFlow, factor = 2),jitter(all_data$germ_ratio, factor = 2),
     type ="n", xaxt = "n", yaxt = "n", ylab="", xlab="", bty = 'n')

colfunc       <- colorRampPalette(cRamp(unique(arrange(viabVr,sr)$sr)))
legend_image  <- as.raster(matrix(colfunc(19), ncol=1))
text(x=8.5, y = seq(0.7,0.9,l=3), labels = seq(0,1,l=3))
rasterImage(legend_image, 0, 0.7, 5, 0.9)
text(10, 0.9, "Percent of", pos = 4)
text(10, 0.8, "males in", pos = 4)
text(10, 0.7, "plot", pos = 4)

legend(-2,0.65, c("100% female plots","10% female plots"),
       lty = c(1,2), lwd = 2, bty="n")

dev.off()

