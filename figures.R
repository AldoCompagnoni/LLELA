# Code to produce figure1 (panels a,b,c)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) #For AIC weights, because I'm lazy!
library(MASS)
library(dplyr)
library(boot)

#Read in data--------------------------------------------------------------------
# raw data
d         <- read.csv("Data/vr.csv")
fem_seeds <- read.csv("Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")
viabVr    <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                      stringsAsFactors = F)
one_sex_tiller <- read.csv("Data/one_sex_plots.csv",
                           stringsAsFactors = F)

# best models
#flow_beta   <- read.csv("Results/VitalRates_3/flowering_best.csv")
n_flow_beta <- read.csv("Results/VitalRates_3/n_flowers_best.csv")
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv")
#tetr_beta   <- read.csv("Results/VitalRates_3/tetrazolium_best.csv")
tetr_beta   <- read.csv("Results/VitalRates_3/germination_best.csv")

# FORMAT DATA -------------------------------------------------------------------

# Only data from 2014
d14         <- subset(d, year == 2014)
d14         <- subset(d14, surv_t1 != 0)
d14[which(d14$log_l_t1==-Inf),]   <- NA
d14$new_t1  <- as.numeric(as.character(d14$new_t1))
f14         <- na.omit(d14[,c("plot","flow_t1","log_l_t1","log_l_t0","flowN_t1",
                              "sex","sr","new_t1","TotDensity")])

# fecundity data ---------------------------------------------------------------------- 
fem_seeds$focalI    <- paste("f",fem_seeds$IndividualN,sep="")
fem_seeds           <- fem_seeds[,c("Plot","focalI","SeedN")] 
fem_seeds           <- fem_seeds[!is.na(fem_seeds$SeedN),]
names(fem_seeds)[1] <- "plot"
#fem_seeds           <- fem_seeds %>% 
#                        group_by(plot,focalI) %>% 
#                          summarise(seed_tot = sum(SeedN))
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


plot( fem14$flowN_t1+0.05 ~ fem14$sr, pch = 16,ylim=c(0,7),xlim=c(0,1),
      ylab = "Number of flowers", xlab="Proportion of female individuals", col = "blue")
par(new=T) ; plot( mal14$flowN_t1-0.05 ~ mal14$sr,pch = 16,ylim=c(0,7),xlim=c(0,1),
                   ylab = "Number of flowers", xlab="",col = "red")

low   <- 5
high  <- 42
size  <- mean(f14$log_l_t0)
xSeq  <- seq(0,1,by = 0.1)
beta  <- n_flow_beta[,c("predictor","avg")]$avg

y_m_l <- exp( beta[1] + beta[2]*size + beta[3]*low + 
                beta[4]*low + beta[5]*xSeq + beta[6] + 
                beta[7]*size + beta[8]*xSeq)
y_m_h <- exp( beta[1] + beta[2]*size + beta[3]*high + 
                beta[4]*high + beta[5]*xSeq + beta[6] + 
                beta[7]*size + beta[8]*xSeq)
y_f_l <- exp( beta[1] + beta[2]*size + beta[3]*low + 
                beta[5]*xSeq)
y_f_h <- exp( beta[1] + beta[2]*size + beta[3]*high + 
                beta[5]*xSeq)

lines(xSeq,y_m_l,lty=2,lwd=2,col="red")
lines(xSeq,y_m_h,lty=1,lwd=2,col="red")
lines(xSeq,y_f_l,lty=2,lwd=2,col="blue")
lines(xSeq,y_f_h,lty=1,lwd=2,col="blue")

legend(0.15,7.5,c("high density","low density","male","female"),
       lty=c(1,2,1,1),lwd=2,col=c("black","black","red","blue"),bty="n")
text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.145,
     par("usr")[4]*0.96,"(a)", cex = 1.2, xpd = T)


# fecundity ----------------------------------------------------------------------------
plot(fecund_data$sr, fecund_data$SeedN, pch = 16, ylim=c(0,950),
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
text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.145,
     par("usr")[4]*0.96,"(b)", cex = 1.2, xpd = T)

# viability ----------------------------------------------------------------------------
lower=quantile(viabVr$totFlow,prob=c(0.1))
upper=quantile(viabVr$totFlow,prob=c(0.9))

# Tetrazolium data
plot(viabVr$sr_f,viabVr$germ_ratio,pch=16,ylim=c(0,1.01),
     xlab="Proportion of female flowers",ylab="Seed viability rate")
xSeq=seq(min(viabVr$sr_f),max(viabVr$sr_f),length.out=100)
beta=tetr_beta$avg
yMeanLow=inv.logit(beta[1] + beta[2]*xSeq + beta[3]*lower + beta[4]*xSeq*lower)
yMeanHigh=inv.logit(beta[1] + beta[2]*xSeq + beta[3]*upper + beta[4]*xSeq*upper)
lines(xSeq,yMeanLow,lwd=2,lty=2)
lines(xSeq,yMeanHigh,lwd=2,lty=1)
text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.145,
     par("usr")[4]*0.96,"(c)", cex = 1.2, xpd = T)

dev.off()


# FIGURE 2 -----------------------------------------------------------------------------

# flowers per plot ----------------------------------------------------------------
predict_flower <- expand.grid("(Intercept)" = 1, 
                              log_l_t0 = mean(f14$log_l_t0),
                              TotDensity = seq(1,48,1),
                              sr = seq(0,1,0.1),
                              sex = c(1,0))
predict_flower$sexmXTotDensity = predict_flower$sex * predict_flower$TotDensity
predict_flower$log_l_t0Xsex = predict_flower$log_l_t0 * predict_flower$sex
predict_flower$srXsexm = predict_flower$sr * predict_flower$sex
# change order of predictors
flow_beta = n_flow_beta[c(1,2,3,5,6,4,7,8),c("predictor","avg")]
# Predictions
predict_flower$pred_n_flow  <- exp( as.matrix(predict_flower) %*% flow_beta$avg )

# Plot level number of flowers
females             <- subset(predict_flower, sex == 0)
females             <- females[order(females$TotDensity,females$sr),]
females$n_fem       <- females$sr * females$TotDensity
females$n_flowers_f <- females$n_fem * females$pred_n_flow

males               <- subset(predict_flower, sex == 1)
males               <- males[order(males$TotDensity,males$sr),]
males$n_mal         <- (1-males$sr) * males$TotDensity
males$n_flowers_m   <- males$n_mal * males$pred_n_flow

n_plot_flowers      <- merge(select(females,sr,TotDensity,n_flowers_f),
                             select(males,  sr,TotDensity,n_flowers_m))
n_plot_flowers$totFlowers <- n_plot_flowers$n_flowers_f + n_plot_flowers$n_flowers_m
n_plot_flowers$sr_flowers <- n_plot_flowers$n_flowers_f / n_plot_flowers$totFlowers



# individual seed numbers  (fecundity) ----------------------------------------------------------------
predict_fec <- expand.grid("(Intercept)" = 1, 
                           log_l_t0 = mean( f14$log_l_t0 ),
                           sr = seq(0,1,0.1),
                           TotDensity = seq(1,48,1))
predict_fec$TotDensityXsr <- predict_fec$TotDensity * predict_fec$sr
predict_fec$pred_i_seeds  <- exp( as.matrix(predict_fec) %*% fec_beta$avg )
predict_fec               <- predict_fec[order(predict_fec$TotDensity,
                                               predict_fec$sr),]

# seed viability (fecundity) ----------------------------------------------------------------
predict_viab            <- n_plot_flowers
predict_viab$intersect  <- 1
predict_viab$srXtotFlow <- predict_viab$totFlowers * predict_viab$sr_flowers
pred_viab_mat           <- select(predict_viab,intersect,sr_flowers,totFlowers,srXtotFlow)
predict_viab$pred_viab  <- inv.logit(as.matrix(pred_viab_mat) %*% tetr_beta$avg)

final           <- merge(select(predict_viab,sr,TotDensity,n_flowers_f,pred_viab),
                         select(predict_fec,sr,TotDensity,pred_i_seeds))
final$fertility <- final$n_flowers_f * final$pred_viab * final$pred_i_seeds
final           <- final[order(final$TotDensity,final$sr),]



# Examples to plot the surface
x<-unique(final$sr)
y<-unique(final$TotDensity)
z<-matrix(final$fertility, nrow=length(unique(final$sr)),
          ncol=length(unique(females$TotDensity)))


tiff("Results/VitalRates_3/figure2.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

persp(x,y,z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",
      zlab = "Viable seeds",
      main = "Fertility")

dev.off()

filled.contour(x,y,z, color.palette=heat.colors,
               xlab = "Proportion of female individuals", 
               ylab = "Density",
               main = "Viable seeds")

#contour(x,y,z,
#        xlab = "Planting Sex ratio", 
#        ylab = "Density")


# FIGURE 3 -----------------------------------------------------------------------------

#remove dead individuals (this is a GROWTH model!)
d=subset(d, surv_t1 != 0)

#only use year two
tmp15=subset(d, year==2015)
tmp15$new_t1[tmp15$new_t1=="SKIPPED"]=NA
tmp15$new_t1[tmp15$new_t1=="cnf"]=NA
tmp15$new_t1[tmp15$new_t1==""]=NA
tmp15$new_t1=as.numeric(as.character(tmp15$new_t1))
d15=na.omit(tmp15[,c("l_t1","log_l_t0","plot","sex","new_t1","sr")])


tiff("Results/VitalRates_3/Figure3.tiff",unit="in",width=6.3,height=3.5,res=600,compression="lzw")

par(mfrow=c(1,2),mar=c(2.5,2.5,1.3,0.5),mgp=c(1.4,0.35,0),cex.lab=1.1,cex.axis=0.8,
    cex.main=0.9)

# Panel a -  Color code sexes --------------------------------------------------------------
tillers$col=as.integer(tillers$sex)
tillers$col=as.character(factor(tillers$col,labels=c("blue","red")))

plot(tillers$tot_till_t0,tillers$tot_till_t1, pch=16, 
     ylab="Tillers per plot (spring 2015)",
     xlab="Tillers per plot (spring 2014)",col=tillers$col)
legend(0,235,c("Male individuals","Female individuals"),
       lty=1,pch=16,lwd=3,col=c("red","blue"),bty="n")

xSeq <- seq(0,max(c(tillers$tot_till_t0,tillers$tot_till_t0)),by=0.1)
beta=one_sex_grow_avg[,c("predictor","avg")]$avg
y_m    <- exp( betas[1] + betas[2] * xSeq + betas[3] + betas[4]*xSeq)
y_f    <- exp( betas[1] + betas[2] * xSeq )

lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.03,
     par("usr")[4]*1.05,"a) One sex plots", cex = 1.2, xpd = T, pos = 4)


# Panel b - Set up colors for plots -------------------------------------------------------
d15$col=as.integer(d15$sex)
d15$col=as.character(factor(d15$col,labels=c("blue","red")))
d15$symb=as.integer(as.character(factor(d15$col,labels=c("17","16"))))

plot(d15$sr,d15$l_t1, pch=16, ylab="Number of leaves per individual",
     xlab="Proportion of female individuals",col=d15$col, ylim = c(0,82))
beta = grow_avg[,c("predictor","avg")]$avg
xSeq <- seq(0,1,by=0.1)
size <- mean(d15$log_l_t0)
dens <- 1
y_m <- exp( beta[1] + beta[2]*size + beta[3]*dens + beta[4] + 
              beta[5]*xSeq + beta[6]*size + beta[7]*dens + beta[8]*xSeq)
y_f <- exp( beta[1] + beta[2]*size + beta[3]*dens + beta[5]*xSeq )
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")
dens <- 48
y_m <- exp( beta[1] + beta[2]*size + beta[3]*dens + beta[4] + 
              beta[5]*xSeq + beta[6]*size + beta[7]*dens + beta[8]*xSeq)
y_f <- exp( beta[1] + beta[2]*size + beta[3]*dens + beta[5]*xSeq )
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.03,
     par("usr")[4]*1.05,"b) Intersexual competition", cex = 1.2, xpd = T, pos = 4)

dev.off()
