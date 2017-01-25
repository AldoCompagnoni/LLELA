# Code to produce Figure 2 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(dplyr)
library(boot)
library(testthat) # to test "order/correspondence of predictors" 
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/ac79/Documents/CODE/LLELA/prediction.R")

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
des_flower    <- new_design(n_flow_beta)
pred_flower   <- pred(des_flower, n_flow_beta, pc_n_flowers)
pred_flower   <- mutate(pred_flower, pc_n_flowers = exp(pc_n_flowers) )

# Seeds per flower (panicle) Figure 1b) 
des_seed      <- new_design(fec_beta)
pred_seed     <- pred(des_seed, fec_beta, seeds_per_f)
pred_seed     <- mutate(pred_seed, seeds_per_f = exp(seeds_per_f) )

# Seed viability (viability model, Figure 1d)
germ_beta     <- mutate(germ_beta, predictor = gsub("_f","",predictor) )
germ_beta     <- mutate(germ_beta, predictor = gsub("totFlow","TotDensity",predictor) )
des_viab      <- new_design(germ_beta)
pred_viab     <- pred(des_viab, germ_beta, pred_viab)
pred_viab     <- mutate(pred_viab, pred_viab = inv.logit(pred_viab) )


# COMBINE MODELS #####################################################################

# PANEL 1: Per capita/total seed production ------------------------------------------------
seed_prod <- merge(dplyr::select(subset(pred_flower, sexm == 0), sr,TotDensity,pc_n_flowers), 
                   dplyr::select(pred_seed, sr,TotDensity,seeds_per_f) )
seed_prod <- mutate(seed_prod, pc_seed_prod = seeds_per_f * pc_n_flowers)
seed_prod <- mutate(seed_prod, seed_prod    = pc_seed_prod * sr * TotDensity)

# PANEL 2: viability (simply "pred_viab") -------------------------------------

# PANEL 3: Viable Seed production ----------------------------------------------------

# Number of female flowers (panicles) per plot 
females       <- subset(pred_flower, sexm == 0)
females       <- mutate(females, n_fem       = sr * TotDensity) # n. of female flowers in plot
females       <- mutate(females, n_f_flowers = n_fem * pc_n_flowers) #n. of fem. panicles
# Number of male flowers (panicles) per plot 
males         <- subset(pred_flower, sexm == 1)
males         <- mutate(males, n_mal       = (1-sr) * TotDensity) # n. of female flowers in plot
males         <- mutate(males, n_m_flowers = n_mal * pc_n_flowers) #n. of fem. panicles


germ_beta   <- read.csv("Results/VitalRates_3/germination_best.csv", stringsAsFactors = F) 
des_sr_f    <- merge(dplyr::select(females,TotDensity,sr,n_fem,n_f_flowers), 
                        dplyr::select(males,TotDensity,sr,n_mal,n_m_flowers))
des_sr_f    <- des_sr_f %>% mutate(
                              Intercept = 1,
                              totFlow = n_m_flowers + n_f_flowers,
                              sr_f = n_f_flowers / (n_m_flowers + n_f_flowers)
                              )
des_sr_f    <- des_sr_f %>% mutate(sr_f_x_totFlow = totFlow * sr_f)
pre_sr_f    <- dplyr::select(des_sr_f, Intercept, sr_f, sr_f_x_totFlow, totFlow)
pre_sr_f    <- setNames(pre_sr_f, c("(Intercept)","sr_f","sr_f:totFlow","totFlow") )  
expect_equal( all(names(pre_sr_f) == germ_beta$predictor), TRUE )

pred_val <- as.vector( as.matrix(pre_sr_f) %*% germ_beta$avg )
pre_sr_f <- mutate(pre_sr_f, viab = pred_val)
names(des_sr_f)[c(7,10)] = c("(Intercept)","sr_f:totFlow")

pred_viab_seed = merge(des_sr_f,pre_sr_f)

fert = merge(pred_viab_seed, seed_prod)
fert <- mutate(fert, viable_seeds = seed_prod * viab)

# fertility (viable seeds)
fert          <- merge(dplyr::select(seed_prod, sr,TotDensity,seed_prod),
                       dplyr::select(predict_viab,sr,TotDensity,pred_viab))
fert          <- mutate(fert, viable_seeds = pred_viab * seed_prod )


# FIGURE 2 ########################################################################

# data for 3d graph
seed_3d <- form_3d_surf(seed_prod,pc_seed_prod)
viab_3d <- form_3d_surf(pred_viab,pred_viab)
vs_3d   <- form_3d_surf(fert,viable_seeds)


tiff("Results/VitalRates_3/figure2_decomp.tiff",unit="in",width=6.3,height=2.1,res=600,compression="lzw")

par(mfrow=c(1,3),mar=c(1,2,1,0),cex=0.5)

persp(seed_3d$x,seed_3d$y,seed_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",
      zlab = "Seeds per individual",
      main = "Per-capita seed production")
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
