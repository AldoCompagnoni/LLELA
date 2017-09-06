# Code to produce figure1 (panels a,b,c)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
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
plot(jitter(d15$TotDensity), d15$log_ratio, pch = d15$symb, ylab="Log ratio (2015)", 
     xlab="Planting density",bg=cRamp(d15$sr),lwd = 1)

# predictions across male/females
lr_15 <- lr_15_avg %>%
            group_by(sex, TotDensity) %>%
            summarise( pred = mean(pred) ) %>%
            ungroup

rgr_f <- subset(lr_15, sex == "f" )
rgr_m <- subset(lr_15, sex == "m" )
lines(rgr_f$TotDensity, rgr_h$pred, col = "blue", lty = 1, lwd = 2 )
lines(rgr_m$TotDensity, rgr_l$pred, col = "red", lty = 2, lwd = 2 )

colfunc = colorRampPalette(cRamp(unique(arrange(d15,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(21), ncol=1))
text(x=52, y = seq(-1.5,0.5,l=3), labels = seq(0,1,l=3))
rasterImage(legend_image, 54, -1.5, 59, 0.5)
text(59, 0.5, "Percent of", pos = 4)
text(59, -0.5, "females in", pos = 4)
text(59, -1.5, "plot", pos = 4)

legend(52,par("usr")[4],c("Males","Females"), cex = 1, pch = c(24,21),
       lty=c(2,1),lwd=1.5,col=c("red","blue"),bty="n")

dev.off()
