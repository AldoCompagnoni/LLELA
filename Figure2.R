# Code to produce Figure 2 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(dplyr)
library(boot)
library(testthat) # to test "order/correspondence of predictors" 
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/ac79/Documents/CODE/LLELA/prediction.R")

# Read in best models--------------------------------------------------------------------
n_flow_beta <- read.csv("Results/VitalRates_3/n_flowers_best.csv", stringsAsFactors = F)
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv", stringsAsFactors = F)
germ_beta   <- read.csv("Results/VitalRates_3/germination_best.csv", stringsAsFactors = F) 
newt_beta   <- read.csv("Results/VitalRates_3/new_t_best.csv", stringsAsFactors = F)

# MODEL PREDICTIONS ---------------------------------------------------------------

# Flowers per individual (flowering model, Figure 1a) -----------------------------
des_flower    <- new_design(n_flow_beta)
pred_flower   <- pred(des_flower, n_flow_beta, pc_n_flowers, exp)

# Seeds per flower (panicle) Figure 1b) 
des_seed      <- new_design(fec_beta)
pred_seed     <- pred(des_seed, fec_beta, seeds_per_f, exp)

# Seed viability (viability model, Figure 1d)
germ_beta     <- mutate(germ_beta, predictor = gsub("_f","",predictor) )
germ_beta     <- mutate(germ_beta, predictor = gsub("totFlow","TotDensity",predictor) )
des_viab      <- new_design(germ_beta)
pred_viab     <- pred(des_viab, germ_beta, pred_viab, inv.logit)

# new tillers
des_till      <- data.frame(TotDensity = seq(1,48,1))
des_till      <- des_till %>% mutate(pred_new_t = (newt_beta[,1]*TotDensity)/(1+newt_beta[,2]*TotDensity))


# COMBINE MODELS #####################################################################

# PANEL 1: Per capita/total seed production ------------------------------------------------
seed_prod <- merge(dplyr::select(subset(pred_flower, sexm == 0), sr,TotDensity,pc_n_flowers), 
                   dplyr::select(pred_seed, sr,TotDensity,seeds_per_f) )
seed_prod <- mutate(seed_prod, pc_seed_prod = seeds_per_f * pc_n_flowers)
seed_prod <- mutate(seed_prod, seed_prod    = pc_seed_prod * sr * TotDensity)

# PANEL 3: Viable Seed production ----------------------------------------------------

# Calculate absolute number of flowers in plots
# Number of female flowers (panicles) per plot 
females       <- subset(pred_flower, sexm == 0)
females       <- mutate(females, n_fem       = sr * TotDensity) # n. of female flowers in plot
females       <- mutate(females, n_f_flowers = n_fem * pc_n_flowers) #n. of fem. panicles
# Number of male flowers (panicles) per plot 
males         <- subset(pred_flower, sexm == 1)
males         <- mutate(males, n_mal       = (1-sr) * TotDensity) # n. of female flowers in plot
males         <- mutate(males, n_m_flowers = n_mal * pc_n_flowers) #n. of fem. panicles
# combine the data sets
flow_num      <- merge(dplyr::select(females,TotDensity,sr,n_fem,n_f_flowers), 
                       dplyr::select(males,  TotDensity,sr,n_mal,n_m_flowers))


#design matrix to calculate NUMBER of viable seeds 
des_v_s    <- flow_num %>% mutate(Intercept = 1,
                                  totFlow = n_m_flowers + n_f_flowers,
                                  sr_f = n_f_flowers / (n_m_flowers + n_f_flowers)
                                  )
des_v_s    <- des_v_s %>% mutate(sr_f_x_totFlow = totFlow * sr_f)

# Predicted NUMBER of viable seeds  
pred_v_seed   <- dplyr::select(des_v_s, Intercept, sr_f, sr_f_x_totFlow, totFlow)
pred_v_seed   <- setNames(pred_v_seed, c("(Intercept)","sr_f","sr_f:totFlow","totFlow") )  
germ_beta     <- read.csv("Results/VitalRates_3/germination_best.csv", stringsAsFactors = F) 
expect_equal( all(names(pred_v_seed) == germ_beta$predictor), TRUE )
tmp_viab      <- as.vector( as.matrix(pred_v_seed) %*% germ_beta$avg )
pred_v_seed   <- mutate(pred_v_seed, viab = inv.logit(tmp_viab) )

# test correspondence des_sr_f
names(des_v_s)[c(7,10)] = c("(Intercept)","sr_f:totFlow")
pred_viab_s   <- merge(des_v_s,pred_v_seed)
fert          <- merge(pred_viab_s, seed_prod)
fert          <- mutate(fert, viable_seeds     = seed_prod * viab,
                              viable_seeds_pc  = (seed_prod * viab)/n_fem)

# reproduction
reprod    <- merge(fert,des_till)
reprod    <- mutate(reprod, reprod = viable_seeds + pred_new_t)

# FIGURE 2 ------------------------------------------------------------------------------

# data for 3d graph
seed_3d <- form_3d_surf(seed_prod,pc_seed_prod)
viab_3d <- form_3d_surf(pred_viab,pred_viab)
vs_3d   <- form_3d_surf(fert,viable_seeds)
#vs_3d   <- form_3d_surf(fert,viable_seeds_pc)
#vs_3d   <- form_3d_surf(reprod,reprod)

tiff("Results/VitalRates_3/figure2_decomp.tiff",unit="in",width=6.3,height=2.1,res=600,compression="lzw")
#tiff("Results/VitalRates_3/figure2_decomp_pc.tiff",unit="in",width=6.3,height=2.1,res=600,compression="lzw")
#tiff("Results/VitalRates_3/figure2_decomp_reprod.tiff",unit="in",width=6.3,height=2.1,res=600,compression="lzw")

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
