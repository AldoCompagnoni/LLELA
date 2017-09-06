# Code to produce figure1 (panels a,b,c)
setwd("C:/cloud/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
source("C:/CODE/LLELA/model_avg.R")
library(dplyr)


#Read and format data--------------------------------------------------------------------
d       <- read.csv("Data/vr.csv")
d       <- format_growth(d)
# separate years
d15     <- filter(d, year == 2015)


# FIGURE ----------------------------------------------------------------------

# Best model
lr_15_avg   <- read.csv("Results/VitalRates_4_Cade2015/log_ratio_15_best.csv")

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
d15 <- sex_symbol(d15)


# plot
tiff("Results/VitalRates_4_Cade2015/Figure_log_ratio_2015.tiff",unit="in",
     width=6.3,height=4.5,res=600,compression="lzw")

par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,9), xpd = NA)

# 2015
plot(jitter(d15$TotDensity), d15$log_ratio, pch = d15$symb, ylab="Growth rate (2015)", 
     xlab="Planting density",bg=cRamp(d15$sr),lwd = 1)

# predictions across male/females
lr_15 <- lr_15_avg %>%
            group_by(sr, TotDensity) %>%
            summarise( pred = mean(pred) ) %>%
            ungroup

rgr_h <- subset(lr_15, sr == 0.05 )
rgr_l <- subset(lr_15, sr == 0.95 )
lines(rgr_h$TotDensity, rgr_h$pred, col = "#ABABAB", lty = 1, lwd = 2 )
lines(rgr_l$TotDensity, rgr_l$pred, col = "#141414", lty = 3, lwd = 2 )

colfunc = colorRampPalette(cRamp(unique(arrange(d15,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(21), ncol=1))
text(x=52, y = seq(-0.5,1.5,l=3), labels = seq(0,1,l=3))
rasterImage(legend_image, 54, -0.5, 59, 1.5)
text(59, 1.5, "Percent of", pos = 4)
text(59, 0.5, "females in", pos = 4)
text(59, -0.5, "plot", pos = 4)

legend("topleft", c("95% female plot", "5% female plot"), cex = 1,
       lty = c(1,3), lwd=2, col=c("#ABABAB","#141414"), bty = "n")

legend(52,-0.7, c("Females","Males"), cex = 1,
       pch = c(21,24), bty = "n")

dev.off()
