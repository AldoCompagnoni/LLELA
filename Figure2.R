# Code to produce Figure 2 
setwd("C:/Users/Aldo/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(dplyr)
library(boot)
library(testthat) # to test "order/correspondence of predictors" 
source("C:/Users/Aldo/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/Aldo/Documents/CODE/LLELA/prediction.R")
source("C:/Users/Aldo/Documents/CODE/LLELA/unload_mass.R")
unload_mass()

# Read in best models--------------------------------------------------------------------
n_flow_beta <- read.csv("Results/VitalRates_3/n_flowers_best.csv", stringsAsFactors = F)
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv", stringsAsFactors = F)
germ_beta   <- read.csv("Results/VitalRates_3/germination_best.csv", stringsAsFactors = F) 
newt_beta   <- read.csv("Results/VitalRates_3/new_t_bh14_best.csv", stringsAsFactors = F)

# MODEL PREDICTIONS ---------------------------------------------------------------

# Flowers per individual (flowering model, Figure 1a) -----------------------------
des_flower    <- new_design(n_flow_beta)
pred_flower   <- pred(des_flower, n_flow_beta, n_flowers_pc, exp)

# Seeds per flower (panicle) Figure 1b) 
des_seed      <- new_design(fec_beta)
pred_seed     <- pred(des_seed, fec_beta, seeds_per_f, exp)

# new tillers
des_till      <- expand.grid(TotDensity = seq(1,48,1), sr = seq(0,1,0.05))
des_till      <- mutate(des_till, F = TotDensity*sr, 
                                  M = TotDensity*(1-sr))
pred_till     <- des_till %>% mutate(pred_new_t = 
                                    (newt_beta[,"lam.f"]*F) / (1+newt_beta[,"b.f"]*(F + newt_beta[,"a.m"]*M)) + 
                                    (newt_beta[,"lam.m"]*M) / (1+newt_beta[,"b.m"]*(newt_beta[,"a.f"]*F + M))
                                     )
pred_till     <- mutate(pred_till, pred_new_t_pc = pred_new_t/TotDensity)
pred_till     <- mutate(pred_till, pred_new_t_pf = pred_new_t/(TotDensity*sr) )
pred_till     <- mutate(pred_till, pred_new_t_pf = replace(pred_new_t_pf, pred_new_t_pf==Inf, NA))
                        
# COMBINE MODELS #####################################################################

# PANEL 1: Per capita/total seed production ------------------------------------------------
seed_prod <- merge(select(subset(pred_flower, sexm == 0), sr,TotDensity,n_flowers_pc), 
                   select(pred_seed, sr,TotDensity,seeds_per_f) )
seed_prod <- mutate(seed_prod, seed_prod    = seeds_per_f * n_flowers_pc * sr * TotDensity)
seed_prod <- mutate(seed_prod, seed_prod_pf = seed_prod / (sr * TotDensity))
seed_prod <- mutate(seed_prod, seed_prod_pc = seed_prod / TotDensity)

# PANEL 2: Per capita/total seed production ------------------------------------------------

# Number of female flowers (panicles) per plot 
females       <- subset(pred_flower, sexm == 0)
females       <- mutate(females, n_fem       = sr * TotDensity) # n. of female flowers in plot
females       <- mutate(females, n_f_flowers = n_fem * n_flowers_pc) #n. of fem. panicles
# Number of male flowers (panicles) per plot 
males         <- subset(pred_flower, sexm == 1)
males         <- mutate(males, n_mal       = (1-sr) * TotDensity) # n. of female flowers in plot
males         <- mutate(males, n_m_flowers = n_mal * n_flowers_pc) #n. of fem. panicles

# combine the data sets
flow_num      <- merge(dplyr::select(females,TotDensity,sr,n_fem,n_f_flowers), 
                       dplyr::select(males,  TotDensity,sr,n_mal,n_m_flowers))

#design matrix to calculate seed viability (proportion) 
des_v     <- flow_num %>% mutate(Intercept = 1,
                                  totFlow = n_m_flowers + n_f_flowers,
                                  sr_f = n_f_flowers / (n_m_flowers + n_f_flowers)
                                 )
des_v     <- des_v %>% mutate(sr_f_x_totFlow = totFlow * sr_f)
names(des_v)[c(7,10)] <- c("(Intercept)","sr_f:totFlow")
pred_viab  <- pred(des_v[,c("(Intercept)","sr_f","sr_f:totFlow","totFlow")], 
                   germ_beta, viab, inv.logit)
all.equal(select(pred_viab,sr_f,totFlow),select(des_v,sr_f,totFlow))
pred_viab  <- cbind(select(des_v,TotDensity,sr,
                           n_fem,n_f_flowers,
                           n_mal,n_m_flowers),
                    pred_viab)


# PANEL 3: Viable Seed production ----------------------------------------------------

# Predicted NUMBER of viable seeds  
#pred_v_seed   <- des_v[,c("(Intercept)","sr_f","sr_f:totFlow","totFlow")]
#germ_beta     <- read.csv("Results/VitalRates_3/germination_best.csv", stringsAsFactors = F) 
#expect_equal( all(names(pred_v_seed) == germ_beta$predictor), TRUE )
#tmp_viab      <- as.vector( as.matrix(pred_v_seed) %*% germ_beta$avg )
#pred_v_seed   <- mutate(pred_v_seed, viab = inv.logit(tmp_viab) )

# test correspondence des_sr_f
#names(des_v_s)[c(7,10)] = c("(Intercept)","sr_f:totFlow")
#pred_viab_s   <- merge(des_v_s,pred_v_seed)
#fert          <- merge(pred_viab_s, seed_prod)

# n. of flowers + viability
fert          <- merge(pred_viab, seed_prod)
fert          <- mutate(fert, viable_seeds     = seed_prod * viab,
                              viable_seeds_pf  = (seed_prod * viab)/n_fem)

# reproduction
reprod    <- merge(fert,pred_till)
reprod    <- mutate(reprod, reprod    = viable_seeds + pred_new_t)
reprod    <- mutate(reprod, reprod_pf = (viable_seeds + pred_new_t)/n_fem)
reprod    <- mutate(reprod, reprod_pc = (viable_seeds + pred_new_t)/TotDensity)


# FIGURE 2 ------------------------------------------------------------------------------

# data for PER CAPITA data

seed_3d <- form_3d_surf(seed_prod,seed_prod_pc)
viab_3d <- form_3d_surf(pred_viab,viab)
till_3d <- form_3d_surf(pred_till,pred_new_t_pc)
vs_3d   <- form_3d_surf(reprod,reprod_pc)

tiff("Results/VitalRates_3/figure2.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

m_line <- -2
par(mfrow=c(2,2),mar=c(0,2,0,0),mgp=c(3,1,0), oma=c(0,0,0,0),cex=0.7,
    cex.axis = 0.7)

persp(seed_3d$x,seed_3d$y,seed_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals",line = m_line,
      ylab = "Density",ticktype = "detailed",
      zlab = "Seed produced",mgp=c(2,1,0),
      main = "Seed production")
persp(viab_3d$x,viab_3d$y,viab_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", line = m_line,
      ylab = "Density",ticktype = "detailed",
      zlab = "Seed viability",
      main = "Seed viability")
persp(till_3d$x,till_3d$y,till_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", line = m_line,
      ylab = "Density",ticktype = "detailed",
      zlab = "New Tillers",
      main = "Production of new tillers")
persp(vs_3d$x,vs_3d$y,vs_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",line = m_line,
      zlab = "Viable recruits",
      main = "Reproduction")

dev.off()



# data for PER FEMALE data
seed_3d <- form_3d_surf(seed_prod,seed_prod_pf)
viab_3d <- form_3d_surf(pred_viab,pred_viab)
till_3d <- form_3d_surf(pred_till,pred_new_t_pf)
vs_3d   <- form_3d_surf(reprod,reprod_pf)

tiff("Results/VitalRates_3/figure2_pf.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

m_line <- -2
par(mfrow=c(2,2),mar=c(0,2,0,0),mgp=c(3,1,0), oma=c(0,0,0,0),cex=0.7,
    cex.axis = 0.7)

persp(seed_3d$x,seed_3d$y,seed_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals",line = m_line,
      ylab = "Density",ticktype = "detailed",
      zlab = "Seed produced_PF",mgp=c(2,1,0),
      main = "Seed production_PF")
persp(viab_3d$x,viab_3d$y,viab_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", line = m_line,
      ylab = "Density",ticktype = "detailed",
      zlab = "Seed viability",
      main = "Seed viability")
persp(till_3d$x,till_3d$y,till_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", line = m_line,
      ylab = "Density",ticktype = "detailed",
      zlab = "New Tillers_PF",
      main = "Production of new tillers_PF")
persp(vs_3d$x,vs_3d$y,vs_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",line = m_line,
      zlab = "Viable recruits_PF",
      main = "Reproduction_PF")

dev.off()



# data for ABSOLUTE data
seed_3d <- form_3d_surf(seed_prod,seed_prod)
viab_3d <- form_3d_surf(pred_viab,pred_viab)
till_3d <- form_3d_surf(pred_till,pred_new_t)
vs_3d   <- form_3d_surf(reprod,reprod)

tiff("Results/VitalRates_3/figure2_abs.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

m_line <- -2
par(mfrow=c(2,2),mar=c(0,2,0,0),mgp=c(3,1,0), oma=c(0,0,0,0),cex=0.7,
    cex.axis = 0.7)

persp(seed_3d$x,seed_3d$y,seed_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals",line = m_line,
      ylab = "Density",ticktype = "detailed",
      zlab = "Seed produced",mgp=c(2,1,0),
      main = "Plot-level seed production")
persp(viab_3d$x,viab_3d$y,viab_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", line = m_line,
      ylab = "Density",ticktype = "detailed",
      zlab = "Seed viability",
      main = "Seed viability")
persp(till_3d$x,till_3d$y,till_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", line = m_line,
      ylab = "Density",ticktype = "detailed",
      zlab = "New Tillers",
      main = "Production of new tillers")
persp(vs_3d$x,vs_3d$y,vs_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",line = m_line,
      zlab = "Viable recruits",
      main = "Reproduction")

dev.off()



par(mar=c(4,4,1,0.5),mgp=c(2,0.7,0))
filled.contour(vs_3d$x,vs_3d$y,vs_3d$z, color.palette=heat.colors, cex.lab = 1.4,
               xlab = "Proportion of female individuals", 
               ylab = "Planting density", 
               main = "Number of viable seeds")
