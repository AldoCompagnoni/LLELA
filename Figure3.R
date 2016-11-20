# Code to produce figure1 (panels a,b,c)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(dplyr)

#Read data--------------------------------------------------------------------
d       <- read.csv("Data/vr.csv")

# Format data ------------------------------------------------------------------

# remove dead individuals (this is a GROWTH model!)
d       <- filter(d, surv_t1 != 0)
# separate years
d14     <- filter(d, year == 2014)
d15     <- filter(d, year == 2015)

# Best models
gr14_avg  <- read.csv("Results/VitalRates_3/growth_14N_best.csv")
gr15_avg  <- read.csv("Results/VitalRates_3/growthN_best.csv")



# FIGURE 3 ----------------------------------------------------------------------

# Sex ratio as dot color ---------------

# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(topo.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

# Set up symbols for the sexes ---------
sex_symbol <- function(x){
  out <- x %>% mutate(symb = as.integer( as.character( 
                                factor( as.integer(sex),labels=c("21","24")))) )
  return(out)
}
d14 <- sex_symbol(d14)
d15 <- sex_symbol(d15)


# Graph
tiff("Results/VitalRates_3/Figure3.tiff",unit="in",width=6.3,height=3.5,res=600,compression="lzw")

par(mfrow=c(1,2),mar=c(2.5,2.5,1.3,0.5),mgp=c(1.4,0.35,0),cex.lab=1,cex.axis=0.8,
    cex.main=0.9,xpd=NA)

# 2014
plot(jitter(d14$TotDensity), d14$l_t1, pch = d14$symb, ylab="Number of leaves in 2014", 
     xlab="Planting density",bg=cRamp(d14$sr), ylim = c(0,99), 
     lwd = 1)
beta <- gr14_avg[,c("predictor","avg")]$avg
xSeq <- seq(0,48,by=1)
sr   <- 0.5
size <- mean(d14$log_l_t0)
y_m <- exp( beta[1] + beta[2]*size + beta[3] + beta[4]*xSeq + 
              beta[5]*sr + beta[6]*size + beta[7]*xSeq*sr )
y_f <- exp( beta[1] + beta[2]*size + beta[4]*xSeq + beta[5]*sr + beta[7]*xSeq*sr)
lines(xSeq,y_f,lwd=1.5,lty=1,col="blue")
lines(xSeq,y_m,lwd=1.5,lty=1,col="red")
text(-2,107.5,"a) Target individuals, year 2014",pos=4)

legend(29,105,c("Males","Females"), cex = 0.8, pch = c(24,21),
       lty=1,lwd=1.5,col=c("red","blue"),bty="n")

# 2015
plot(jitter(d15$TotDensity),d15$l_t1, pch = d15$symb, ylab="Number of leaves in 2015",
     xlab="Planting density",bg=cRamp(d15$sr), ylim = c(0,99),
     lwd = 1)
beta <- gr15_avg[,c("predictor","avg")]$avg
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

