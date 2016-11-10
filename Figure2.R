# Code to produce Figure 2 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(MASS)
library(dplyr)
library(boot)

# Read in data--------------------------------------------------------------------
# raw data
d           <- read.csv("Data/vr.csv")

# best models
n_flow_beta <- read.csv("Results/VitalRates_3/n_flowers_best.csv")
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv")
germ_beta   <- read.csv("Results/VitalRates_3/germination_dens_best.csv")


# FORMAT DATA ---------------------------------------------------------------------

# plot level data, including number of flowers per individual ---------------------
d14         <- subset(d, year == 2014)
d14         <- subset(d14, surv_t1 != 0)
f14         <- na.omit(d14[,c("plot","log_l_t0","flowN_t1","sex","sr","TotDensity")])


# MODEL PREDICTIONS ---------------------------------------------------------------

# flowers per individual (flowering model, Figure 1a) -----------------------------
predict_flower <- expand.grid("Intercept" = 1, 
                              log_l_t0 = mean(f14$log_l_t0),
                              TotDensity = seq(1,48,1),
                              sr = seq(0,1,0.1),
                              sex = c(1,0))
predict_flower$sexmXTotDensity  <- predict_flower$sex * predict_flower$TotDensity
predict_flower$log_l_t0Xsex     <- predict_flower$log_l_t0 * predict_flower$sex
predict_flower$srXsexm          <- predict_flower$sr * predict_flower$sex
# change order of predictors to match prediction matrix
flow_beta     <- n_flow_beta[c(1,2,3,5,6,4,7,8),c("predictor","avg")]
# Per capita number of flowers
predict_flower$pc_flowers      <- exp( as.matrix(predict_flower) %*% flow_beta$avg )


# Number of female flowers per plot -----------------------------------------------
females                   <- subset(predict_flower, sex == 0)
females                   <- females[order(females$TotDensity,females$sr),]
females$n_fem             <- females$sr * females$TotDensity
females$n_f_flowers       <- females$n_fem * females$pc_flowers


# Seeds per flower (fecundity model, Figure 1b) -----------------------------------
predict_fec <- expand.grid("Intercept" = 1, 
                           log_l_t0 = mean( f14$log_l_t0 ),
                           sr = seq(0,1,0.1),
                           TotDensity = seq(1,48,1) )
predict_fec$TotDensityXsr <- predict_fec$TotDensity * predict_fec$sr
predict_fec$per_f_seeds   <- exp( as.matrix(predict_fec) %*% fec_beta$avg )


# Seed viability (viability model, Figure 1c) -------------------------------------
predict_viab              <- select(predict_fec,Intercept,sr,TotDensity,TotDensityXsr)
predict_viab$pred_viab    <- inv.logit( as.matrix(predict_viab) %*% germ_beta$avg)


# Seeds per plot ----------------------------------------------------------------
seed_x_plot               <- merge(select(females,    sr, TotDensity, n_f_flowers),
                                   select(predict_fec,sr, TotDensity, per_f_seeds))
seed_x_plot$n_seeds       <-  seed_x_plot$n_f_flowers * seed_x_plot$per_f_seeds

# Final data frame ----------------------------------------------------------------
fert                 <- merge(select(seed_x_plot, sr,TotDensity,n_seeds),
                              select(predict_viab,sr,TotDensity,pred_viab))
fert$viable_seeds    <- fert$pred_viab * fert$n_seeds
fert                 <- fert[order(fert$TotDensity,fert$sr),]

# FIGURE 2 ------------------------------------------------------------------------

# Prepare data
x<-unique(fert$sr)
y<-unique(fert$TotDensity)
z<-matrix(fert$viable_seeds, nrow=length(unique(fert$sr)),
          ncol=length(unique(females$TotDensity)))
#z<-matrix(fert$ciao, nrow=length(unique(fert$sr)),
#          ncol=length(unique(females$TotDensity)))



tiff("Results/VitalRates_3/figure2_a_tetra.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

persp(x,y,z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",
      zlab = "Viable seeds",
      main = "Fertility")

dev.off()

# Alternatives to persp() function
tiff("Results/VitalRates_3/figure2_a_tetra_contour.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

par(mar=c(4,4,1,0.5),mgp=c(2,0.7,0))
filled.contour(x,y,z, color.palette=heat.colors, cex.lab = 1.4,
               xlab = "Proportion of female individuals", 
               ylab = "Density", 
               main = "Viable seeds")

dev.off()


contour(x,y,z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",
      main = "Viable seeds")


# Per capita fertility ---------------------------------------------------------

z<-matrix(fert$ciao, nrow=length(unique(fert$sr)),
          ncol=length(unique(females$TotDensity)))


tiff("Results/VitalRates_3/figure2_a.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

persp(x,y,z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",
      zlab = "Viable seeds per flower",
      main = "Fertility")

dev.off()

# Alternatives to persp() function
tiff("Results/VitalRates_3/figure2_a_contour.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

par(mar=c(4,4,1,0.5),mgp=c(2,0.7,0))
filled.contour(x,y,z, color.palette=heat.colors, cex.lab = 1.4,
               xlab = "Proportion of female individuals", 
               ylab = "Density", 
               main = "Viable seeds per flower")

dev.off()
