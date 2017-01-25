# Code to produce Figure 2 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(dplyr)
library(boot)
library(testthat) # to test "order/correspondence of predictors" 
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")

# Read in data--------------------------------------------------------------------
# raw data
d           <- read.csv("Data/vr.csv")

# best models
n_flow_beta <- read.csv("Results/VitalRates_3/n_flowers_best.csv", stringsAsFactors = F)
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv", stringsAsFactors = F)
germ_beta   <- read.csv("Results/VitalRates_3/germination_best.csv", stringsAsFactors = F) 

# FORMAT DATA ---------------------------------------------------------------------

# plot level data, including number of flowers per individual ---------------------
f14   <- format_flower(d)

# MODEL PREDICTIONS ---------------------------------------------------------------

# Flowers per individual (flowering model, Figure 1a) -----------------------------
predict_flower <- expand.grid("(Intercept)" = 1, 
                              TotDensity = seq(1,48,1),
                              sexm = c(1,0),
                              sr = seq(0,1,0.01))
predict_flower <- mutate(predict_flower, "sr:TotDensity" = sr * TotDensity)
predict_flower <- mutate(predict_flower, "sexm:TotDensity" = sexm * TotDensity)
predict_flower <- mutate(predict_flower, "sexm:sr" = sr * sexm)

# TEST correspondence of predictors
expect_equal( all(names(predict_flower) == n_flow_beta$predictor), TRUE )
# Per capita number of flowers
predict_flower <- mutate(predict_flower, 
                         pc_n_flowers = as.vector(exp( as.matrix(predict_flower) %*% 
                                                       n_flow_beta$avg ))
                         )

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
                        as.vector(exp( as.matrix(predict_fec) %*% fec_beta$avg )) 
                      )

# Seed viability (viability model, Figure 1d) -------------------------------------
predict_viab  <- dplyr::select(predict_fec,-seeds_per_f)
predict_viab  <- predict_viab[,order(names(predict_viab))] 
germ_beta     <- mutate(germ_beta, predictor = gsub("_f","",predictor) )
germ_beta     <- mutate(germ_beta, predictor = gsub("totFlow","TotDensity",predictor) )
expect_equal( all(names(predict_viab) == germ_beta$predictor), TRUE )
predict_viab  <- mutate(predict_viab, pred_viab  = 
                          as.vector( inv.logit( as.matrix(predict_viab) %*% germ_beta$avg)) )


# COMBINE MODELS #####################################################################

# Per capita and total seed production ------------------------------------------------
seed_prod <- merge(select(subset(predict_flower, sexm == 0), sr,TotDensity,pc_n_flowers), 
                   select(predict_fec, sr,TotDensity,seeds_per_f) )
seed_prod <- mutate(seed_prod, pc_seed_prod = seeds_per_f * pc_n_flowers)
seed_prod <- mutate(seed_prod, seed_prod    = pc_seed_prod * sr * TotDensity)

# fertility (viable seeds)
fert          <- merge(dplyr::select(seed_prod, sr,TotDensity,seed_prod),
                       dplyr::select(predict_viab,sr,TotDensity,pred_viab))
fert          <- mutate(fert, viable_seeds = pred_viab * seed_prod )


# FIGURE 2 ########################################################################

# data for 3d graph
seed_3d <- form_3d_surf(seed_prod,seed_prod)
viab_3d <- form_3d_surf(predict_viab,pred_viab)
vs_3d   <- form_3d_surf(fert,viable_seeds)


tiff("Results/VitalRates_3/figure2_decomp.tiff",unit="in",width=6.3,height=2.1,res=600,compression="lzw")

par(mfrow=c(1,3),mar=c(1,2,1,0),cex=0.5)

persp(seed_3d$x,seed_3d$y,seed_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",
      zlab = "Seeds",
      main = "Seed production")
persp(viab_3d$x,viab_3d$y,viab_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",
      zlab = "Proportion viable",
      main = "Seed viability")
persp(vs_3d$x,vs_3d$y,vs_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",
      zlab = "Viable seeds",
      main = "Fertility")

dev.off()


# Alternatives to persp() function
#tiff("Results/VitalRates_3/figure2_contour_decomp.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

par(mfrow=c(2,2),mar=c(4,4,1,0.5),mgp=c(2,0.7,0))
filled.contour(seed_3d$x,seed_3d$y,seed_3d$z, color.palette=heat.colors, cex.lab = 1.4,
               xlab = "Proportion of female individuals", 
               ylab = "Planting density", 
               main = "Number of viable seeds")
filled.contour(viab_3d$x,viab_3d$y,viab_3d$z, color.palette=heat.colors, cex.lab = 1.4,
               xlab = "Proportion of female individuals", 
               ylab = "Planting density", 
               main = "Number of viable seeds")
filled.contour(vs_3d$x,vs_3d$y,vs_3d$z, color.palette=heat.colors, cex.lab = 1.4,
               xlab = "Proportion of female individuals", 
               ylab = "Planting density", 
               main = "Number of viable seeds")

#dev.off()



#tiff("Results/VitalRates_3/figure2_target.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

#persp(x,y,z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
#      xlab = "Proportion of female individuals", 
#      ylab = "Density",ticktype = "detailed",
#      zlab = "Viable seeds",
#      main = "Fertility")

#dev.off()

# Alternatives to persp() function
#tiff("Results/VitalRates_3/figure2_contour_target.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

#par(mar=c(4,4,1,0.5),mgp=c(2,0.7,0))
#filled.contour(x,y,z, color.palette=heat.colors, cex.lab = 1.4,
#               xlab = "Proportion of female individuals", 
#               ylab = "Planting density", 
#               main = "Number of viable seeds")

#dev.off()
