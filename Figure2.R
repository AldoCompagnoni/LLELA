# Code to produce Figure 2 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(dplyr)
library(boot)
library(testthat) # to test "order/correspondence of predictors" 

# Read in data--------------------------------------------------------------------
# raw data
d           <- read.csv("Data/vr.csv")

# best models
n_flow_beta <- read.csv("Results/VitalRates_3/n_flowers_best.csv", stringsAsFactors = F)
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv", stringsAsFactors = F)
germ_beta   <- read.csv("Results/VitalRates_3/germination_best.csv", stringsAsFactors = F) 

# FORMAT DATA ---------------------------------------------------------------------

# plot level data, including number of flowers per individual ---------------------
d14         <- subset(d, year == 2014)
d14         <- subset(d14, surv_t1 != 0)
f14         <- na.omit(dplyr::select(d14,plot,log_l_t0,flowN_t1,sex,sr,TotDensity))


# MODEL PREDICTIONS ---------------------------------------------------------------

# Flowers per individual (flowering model, Figure 1a) -----------------------------
predict_flower <- expand.grid("(Intercept)" = 1, 
                              log_l_t0 = mean(f14$log_l_t0),
                              TotDensity = seq(1,48,1),
                              sr = seq(0,1,0.01),
                              sexm = c(1,0))
predict_flower <- mutate(predict_flower, "sr:TotDensity" = sr * TotDensity)
predict_flower <- mutate(predict_flower, "log_l_t0:sexm" = log_l_t0 * sexm)

# change order of predictors to match prediction matrix
flow_beta      <- n_flow_beta[c(1,2,3,4,6,5,7), c("predictor","avg")]
# TEST correspondence of predictors
expect_equal( all(names(predict_flower) == flow_beta$predictor), TRUE )
# Per capita number of flowers
predict_flower <- mutate(predict_flower, 
                         pc_n_flowers = as.vector(exp( as.matrix(predict_flower) %*% flow_beta$avg )))


# Seeds per flower (fecundity model, Figure 1b) -----------------------------------

# Number of female flowers (panicles) per plot 
females        <- subset(predict_flower, sexm == 0)
females        <- females[order(females$TotDensity,females$sr),]
females        <- mutate(females, n_fem = sr * TotDensity) # n. of female flowers in plot
females        <- mutate(females, n_f_flowers = n_fem * pc_n_flowers) #n. of fem. panicles

# Seeds per flower (panicle)
predict_fec <- expand.grid("(Intercept)" = 1, 
                           sr = seq(0,1,0.01),
                           TotDensity = seq(1,48,1) )
predict_fec <- mutate(predict_fec, "sr:TotDensity" = TotDensity * sr)
# TEST correspondence of predictors
expect_equal( all(names(predict_fec) == fec_beta$predictor), TRUE )
predict_fec <- mutate(predict_fec, seeds_per_f = 
                        as.vector(exp( as.matrix(predict_fec) %*% fec_beta$avg )) )


# Seed viability (viability model, Figure 1d) -------------------------------------
predict_viab  <- dplyr::select(predict_fec,-seeds_per_f)
predict_viab  <- predict_viab[,order(names(predict_viab))] 
germ_beta     <- mutate(germ_beta, predictor = gsub("_f","",predictor) )
germ_beta     <- mutate(germ_beta, predictor = gsub("totFlow","TotDensity",predictor) )
expect_equal( all(names(predict_viab) == germ_beta$predictor), TRUE )
predict_viab  <- mutate(predict_viab, pred_viab  = 
                          as.vector( inv.logit( as.matrix(predict_viab) %*% germ_beta$avg)) )


# Final data frame ----------------------------------------------------------------

# Seeds per plot
seed_x_plot   <- merge(dplyr::select(females,    sr, TotDensity, n_f_flowers),
                       dplyr::select(predict_fec,sr, TotDensity, seeds_per_f))
seed_x_plot   <- mutate(seed_x_plot, n_seeds = n_f_flowers * seeds_per_f)

# fertility (viable seeds)
fert          <- merge(dplyr::select(seed_x_plot, sr,TotDensity,n_seeds),
                       dplyr::select(predict_viab,sr,TotDensity,pred_viab))
fert          <- mutate(fert, viable_seeds = pred_viab * n_seeds )
fert          <- fert[order(fert$TotDensity,fert$sr),] #order requires by persp & contour


# FIGURE 2 ------------------------------------------------------------------------

# Prepare data
x<-unique(fert$sr)
y<-unique(fert$TotDensity)
z<-matrix(fert$viable_seeds, nrow=length(unique(fert$sr)),
          ncol=length(unique(females$TotDensity)))


tiff("Results/VitalRates_3/figure2_target.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

persp(x,y,z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",
      zlab = "Viable seeds",
      main = "Fertility")

dev.off()

# Alternatives to persp() function
tiff("Results/VitalRates_3/figure2_contour_target.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

par(mar=c(4,4,1,0.5),mgp=c(2,0.7,0))
filled.contour(x,y,z, color.palette=heat.colors, cex.lab = 1.4,
               xlab = "Proportion of female individuals", 
               ylab = "Planting density", 
               main = "Number of viable seeds")

dev.off()