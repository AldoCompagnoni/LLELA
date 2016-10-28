# Code to produce figure1 (panels a,b,c)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) #For AIC weights, because I'm lazy!
library(MASS)
library(dplyr)

#Read in data--------------------------------------------------------------------
# raw data
d         <- read.csv("Data/vr.csv")
fem_seeds <- read.csv("Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")
viabVr    <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                      stringsAsFactors = F)

# best models
flow_beta   <- read.csv("Results/VitalRates_3/flowering_best.csv")
n_flow_beta <- read.csv("Results/VitalRates_3/n_flowers_best.csv")
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv")
tetr_beta   <- read.csv("Results/VitalRates_3/tetrazolium_best.csv")


# FORMAT DATA -------------------------------------------------------------------
# plot level data 
#Remove resuscitated individuals (dead in spring 2014, alive in 2015)
#NOTE: probably more adequate to consider these NOT DEAD in 2014, because they all have
#few leaves: they probably resprouted from base?
d <- d[-which(d$plot == 18 & d$focalI == "m2"),]
d <- d[-which(d$plot == 38 & d$focalI == "f1"),]
d <- d[-which(d$plot == 46 & d$focalI == "f1"),]
d <- d[-which(d$plot == 83 & d$focalI == "f5"),]
d <- d[-which(d$plot == 36 & d$focalI == "m3"),]

# logtransform leaf numbers
d$log_l_t0 <- log(d$l_t0)
d$log_l_t1 <- log(d$l_t1)
# Sex ratios
d$sr <- d$F / d$TotDensity  

# Only data from 2014
d14         <- subset(d, year == 2014)
d14         <- subset(d14, surv_t1 != 0)
d14[which(d14$log_l_t1==-Inf),]   <- NA
d14$new_t1  <- as.numeric(as.character(d14$new_t1))
f14         <- na.omit(d14[,c("plot","flow_t1","log_l_t1","log_l_t0",
                              "sex","sr","new_t1","TotDensity")])

# fecundity data ---------------------------------------------------------------------- 
fem_seeds$focalI    <- paste("f",fem_seeds$IndividualN,sep="")
fem_seeds           <- fem_seeds[,c("Plot","focalI","SeedN")] 
fem_seeds           <- fem_seeds[!is.na(fem_seeds$SeedN),]
names(fem_seeds)[1] <- "plot"
fem_seeds           <- fem_seeds %>% 
                        group_by(plot,focalI) %>% 
                          summarise(seed_tot = sum(SeedN))
# merge data sets 
fecund_data         <- merge(d14,fem_seeds) 


# FIGURE 1 -----------------------------------------------------------------------------

tiff("Results/VitalRates_3/figure1.tiff",unit="in",width=6.3,height=2.1,res=600,compression="lzw")

par(mfrow=c(1,3),mar=c(2.5,2.5,0.1,0.5),mgp=c(1.4,0.35,0),cex.lab=1.1,cex.axis=0.8,
    cex.main=0.9)

# Flowering ------------------------------------------------------------------------------
# Set up colors for plots
f14$col <- as.character(factor(as.integer(f14$sex),labels=c("blue","red")))
mal14   <- subset(f14,sex=="m")
fem14   <- subset(f14,sex=="f")

# Plot data
plot(mal14$flow_t1+0.015 ~ mal14$sr,
     xlab="Proportion female individuals",ylim=c(0,1.1),
     ylab="Probability of flowering",pch=16,col=mal14$col)
par(new=T) ; plot(fem14$flow_t1-0.015 ~ fem14$sr,pch=17,xlab="",
                  ylab="",col=fem14$col,xaxt="n",ylim=c(0,1.1))
text(par("usr")[1] + (par("usr")[2] - par("usr")[1])*0.05,
     par("usr")[4]*0.95,"a)", cex = 1.2)

# Compute model predictions
beta  <- flow_beta$avg
low   <- 5
high  <- 42
size  <- mean(f14$log_l_t0)
xSeq  <- seq(0,1,by = 0.1)

y_m_l <- inv.logit( beta[1] + beta[2]*size + beta[3]*low + 
                      beta[4] + beta[5]*xSeq + beta[6]*low +
                      beta[7]*size + beta[8]*xSeq)
y_m_h <- inv.logit( beta[1] + beta[2]*size + beta[3]*high + 
                      beta[4] + beta[5]*xSeq + beta[6]*high +
                      beta[7]*size + beta[8]*xSeq)
y_f_l <- inv.logit( beta[1] + beta[2]*size + beta[3]*low + 
                      beta[5]*xSeq)
y_f_h <- inv.logit( beta[1] + beta[2]*size + beta[3]*high + 
                      beta[5]*xSeq)

lines(xSeq,y_m_h,lty=1,lwd=2,col="red")
lines(xSeq,y_m_l,lty=2,lwd=2,col="red")

lines(xSeq,y_f_h,lty=1,lwd=2,col="blue")
lines(xSeq,y_f_l,lty=2,lwd=2,col="blue")

legend(-0.03,0.4,c("high density","low density"),
       lty=c(1,2),lwd=2,bty="n")
legend(0.57,0.4,c("male","female"),
       lty=1,lwd=2,col=c("red","blue"),bty="n")

# fecundity ----------------------------------------------------------------------------
plot(fecund_data$sr, fecund_data$seed_tot, pch = 16,
     xlab = "Proportion female individuals", ylab = "Seeds per flower")
xSeq   <- seq(0,1,by = 0.1)
mSize  <- mean(fecund_data$log_l_t0)
beta   <- fec_beta$avg
y_low  <- exp(beta[1] + beta[2]*mSize + beta[3]*xSeq + 
                beta[4]*5 + beta[5]*xSeq*5)
y_high <- exp(beta[1] + beta[2]*mSize + beta[3]*xSeq + 
                beta[4]*42 + beta[5]*xSeq*42)
lines(xSeq, y_low, lwd = 2, lty = 2)
lines(xSeq, y_high, lwd = 2)
text(par("usr")[1] + (par("usr")[2] - par("usr")[1])*0.05,
     par("usr")[4]*0.95,"b)", cex = 1.2)

# viability ----------------------------------------------------------------------------
lower=quantile(viabVr$totFlow,prob=c(0.1))
upper=quantile(viabVr$totFlow,prob=c(0.9))

# Tetrazolium data
plot(viabVr$sr_f,viabVr$tetra_maybe_ratio,pch=16,ylim=c(0,1.01),
     xlab="Proportion of female flowers",ylab="Seed viability rate")
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
beta=tetr_beta$avg
yMeanLow=inv.logit(beta[1] + beta[2]*xSeq + beta[3]*lower + beta[4]*xSeq*lower)
yMeanHigh=inv.logit(beta[1] + beta[2]*xSeq + beta[3]*upper + beta[4]*xSeq*upper)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)
text(par("usr")[1] + (par("usr")[2] - par("usr")[1])*0.05,
     par("usr")[4]*0.95,"c)", cex = 1.2)

dev.off()


# FIGURE 2 -----------------------------------------------------------------------------
flow_beta <- flow_beta[c(1:3,5,4,6:8),]
fec_beta
tetr_beta



# flowers per plot ----------------------------------------------------------------
predict_flower <- expand.grid("(Intercept)" = 1, 
                              log_l_t0 = mean(f14$log_l_t0),
                              TotDensity = seq(1,48,1),
                              sr = seq(0,1,0.1),
                              sex = c(1,0))
predict_flower$sexmXTotDensity = predict_flower$sex * predict_flower$TotDensity
predict_flower$log_l_t0Xsex = predict_flower$log_l_t0 * predict_flower$sex
predict_flower$srXsexm = predict_flower$sr * predict_flower$sex

# Predictions
predict_flower$pred_flow  <- inv.logit( as.matrix(predict_flower) %*% flow_beta$avg )

# Plot level number of flowers
females             <- subset(predict_flower, sex == 0)
females             <- females[order(females$TotDensity,females$sr),]
females$n_fem       <- females$sr * females$TotDensity
females$n_flowers_f <- females$n_fem * females$pred_flow

males               <- subset(predict_flower, sex == 1)
males               <- males[order(males$TotDensity,males$sr),]
males$n_mal         <- (1-males$sr) * males$TotDensity
males$n_flowers_m   <- males$n_mal * males$pred_flow

n_plot_flowers      <- merge(select(females,sr,TotDensity,n_flowers_f),
                             select(males,  sr,TotDensity,n_flowers_m))
n_plot_flowers$TotFlowers <- n_plot_flowers$n_flowers_f + n_plot_flowers$n_flowers_m




# Plot level number of seeds (fecundity) ----------------------------------------------------------------
predict_fec <- expand.grid("(Intercept)" = 1, 
                           log_l_t0 = mean( f14$log_l_t0 ),
                           sr = seq(0,1,0.1),
                           TotDensity = seq(1,48,1))
predict_fec$TotDensityXsr <- predict_fec$TotDensity * predict_fec$sr
predict_fec$pred_fec      <- exp( as.matrix(predict_fec) %*% fec_beta$avg )
predict_fec               <- predict_fec[order(predict_fec$TotDensity,
                                               predict_fec$sr),]



# Examples to plot the surface
x<-unique(females$sr)
y<-unique(females$TotDensity)
z<-matrix(females$pred, nrow=length(unique(females$sr)),
          ncol=length(unique(females$TotDensity)))

par(mfrow=c(1,1))
contour(x,y,z)
filled.contour(x,y,z)
persp(x,y,z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")

