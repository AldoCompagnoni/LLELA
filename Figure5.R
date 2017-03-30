##Response variable: ,mating (proportion of fertilized seeds)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) 
library(glmmADMB) # Fit models with a Negative Binomial
library(dplyr)
library(boot)
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/ac79/Documents/CODE/LLELA/prediction.R")

# load and format data -----------------------------------------------------
viabVr    <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                        stringsAsFactors = F)
viabVr    <- subset(viabVr, totFlow < 60) # remove three plots with more than 60 flowers 
# best model
germ_beta <- read.csv("Results/VitalRates_3/germination_best.csv")
# model matrix for prediction 
mod_des   <- expand.grid("(Intercept)" = 1, 
                         totFlow = seq(1,50,1),
                         sr_f = round(seq(0,1,0.05),2)) 
mod_des$'sr_f:totFlow' <- mod_des$totFlow * mod_des$sr_f


# Graph -----------------------------------------------------------------------------

# Sex ratio as dot color ------
# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(gray.colors(7,start = 0, end = 1))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

#Graph: total number of tillers -----------------------------------------------------
tiff("Results/VitalRates_3/Figure5.tiff",unit="in",width=6.3,height=4,res=600,compression="lzw")

par(mfrow=c(1,1),mar=c(2.6,2.5,0.2,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,10))

# viability (germination) data ---------------------------------------------------------
plot(jitter(viabVr$totFlow,factor = 2),jitter(viabVr$germ_ratio,factor = 2),
     pch=21,ylim=c(0,1.01),xlim=c(0,48), 
     bg = cRamp(viabVr$sr_f), cex = 1.5,
     xlab="Number of flowers",ylab="Seed viability rate")

# limit model matrix
l_des <- subset(mod_des,sr_f == 0.05)
h_des <- subset(mod_des,sr_f == 0.95)
# predicted values using 'pred' function
l_pred <- pred(l_des, germ_beta, viab, inv.logit)
h_pred <- pred(h_des, germ_beta, viab, inv.logit)
# plot lines
lines(l_des$totFlow,l_pred$viab,lwd=2,lty=2)
lines(h_des$totFlow,h_pred$viab,lwd=2,lty=1)

# sex ratio legend
colfunc = colorRampPalette(cRamp(unique(arrange(viabVr,sr_f)$sr_f)))
legend_image <- as.raster(matrix(colfunc(21), ncol=1))
text(x=60, y = seq(0.6,1,l=3), labels = seq(1,0,l=3), xpd=NA)
rasterImage(legend_image, 51, 0.6, 58, 1, xpd=NA)
text(61, 1, "Percent of", pos = 4,xpd=NA)
text(61, 0.8, "females in", pos = 4,xpd=NA)
text(61, 0.6, "plot", pos = 4,xpd=NA)

# prediction legend
legend(50,0.55, c("5% female plot", "95% female plot"), cex = 1, xpd = NA,
       lty = c(1,2), lwd=2, bty = "n")

dev.off()

