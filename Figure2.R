# Code to produce figure1 (panels a,b,c)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(MASS)
library(dplyr)
library(boot)

# Read in data--------------------------------------------------------------------
# raw data
d         <- read.csv("Data/vr.csv")

# best models
n_flow_beta <- read.csv("Results/VitalRates_3/n_flowers_best.csv")
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv")
tetr_beta   <- read.csv("Results/VitalRates_3/germination_best.csv")

# FORMAT DATA -------------------------------------------------------------------------

# plot level data, including number of flowers per individual -------------------------
d14         <- subset(d, year == 2014)
d14         <- subset(d14, surv_t1 != 0)
f14         <- na.omit(d14[,c("plot","flow_t1","log_l_t1","log_l_t0","flowN_t1",
                              "sex","sr","TotDensity")])


# MODEL PREDICTIONS -----------------------------------------------------------------------------

# flowers per individual (flowering model, Figure 1a) ---------------------------------------------
predict_flower <- expand.grid("(Intercept)" = 1, 
                              log_l_t0 = mean(f14$log_l_t0),
                              TotDensity = seq(1,48,1),
                              sr = seq(0,1,0.1),
                              sex = c(1,0))
predict_flower$sexmXTotDensity = predict_flower$sex * predict_flower$TotDensity
predict_flower$log_l_t0Xsex = predict_flower$log_l_t0 * predict_flower$sex
predict_flower$srXsexm = predict_flower$sr * predict_flower$sex
# change order of predictors to match prediction matrix
flow_beta = n_flow_beta[c(1,2,3,5,6,4,7,8),c("predictor","avg")]
# Predictions
predict_flower$pred_n_flow  <- exp( as.matrix(predict_flower) %*% flow_beta$avg )


# flowers per plot ---------------------------------------------------------
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


# seeds per flower (fecundity model, Figure 1b) ---------------------------------------------
predict_fec <- expand.grid("(Intercept)" = 1, 
                           log_l_t0 = mean( f14$log_l_t0 ),
                           sr = seq(0,1,0.1),
                           TotDensity = seq(1,48,1))
predict_fec$TotDensityXsr <- predict_fec$TotDensity * predict_fec$sr
predict_fec$pred_i_seeds  <- exp( as.matrix(predict_fec) %*% fec_beta$avg )


# seed viability (viability model, Figure 1c) -----------------------------------------------
predict_viab            <- n_plot_flowers
predict_viab$intersect  <- 1
predict_viab$srXtotFlow <- predict_viab$totFlowers * predict_viab$sr_flowers
pred_viab_mat           <- select(predict_viab,intersect,sr_flowers,totFlowers,srXtotFlow)
predict_viab$pred_viab  <- inv.logit( as.matrix(pred_viab_mat) %*% tetr_beta$avg)


# Final data frame --------------------------------------------------------------------------
final           <- merge(select(predict_viab,sr,TotDensity,n_flowers_f,pred_viab),
                         select(predict_fec,sr,TotDensity,pred_i_seeds))
final$fertility <- final$n_flowers_f * final$pred_viab * final$pred_i_seeds
final           <- final[order(final$TotDensity,final$sr),]



# FIGURE 2 -----------------------------------------------------------------------------------

# Prepare data
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

# Alternatives to persp() function
tiff("Results/VitalRates_3/figure2_contour.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

par(mar=c(3,3,1,0.5))
filled.contour(x,y,z, color.palette=heat.colors,
               xlab = "Proportion of female individuals", 
               ylab = "Density",
               main = "Viable seeds")
dev.off()


contour(x,y,z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",
      main = "Viable seeds")
