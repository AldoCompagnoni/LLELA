# Code to produce Figure 2 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(dplyr)
library(boot)
library(testthat) # to test "order/correspondence of predictors" 
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/ac79/Documents/CODE/LLELA/prediction.R")
source("C:/Users/ac79/Documents/CODE/LLELA/unload_mass.R")
unload_mass()
source("C:/Users/ac79/Documents/CODE/LLELA/filled.contour3.R")
source("http://wiki.cbr.washington.edu/qerm/sites/qerm/images/2/25/Filled.legend.R")


# Read in best models--------------------------------------------------------------------
n_flow_beta <- read.csv("Results/VitalRates_3/n_flowers_best.csv", stringsAsFactors = F)
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv", stringsAsFactors = F)
germ_beta   <- read.csv("Results/VitalRates_3/germination_best.csv", stringsAsFactors = F) 
newt_beta   <- read.csv("Results/VitalRates_3/new_t_bh14_best.csv", stringsAsFactors = F)

# "remember" graphical parameters
opar  <- par()

# MODEL PREDICTIONS ---------------------------------------------------------------

# Flowers per individual (flowering model, Figure 1a) -----------------------------
des_flower  <- new_design(n_flow_beta)
pred_flower <- pred(des_flower, n_flow_beta, n_flowers_pc, exp)

# Seeds per flower (panicle) Figure 1b) 
des_seed    <- new_design(fec_beta)
pred_seed   <- pred(des_seed, fec_beta, seeds_per_f, exp)

# new tillers
des_till    <- expand.grid(TotDensity = seq(1,48,1), sr = seq(0,1,0.05))
des_till    <- mutate(des_till, F = TotDensity*sr, 
                                  M = TotDensity*(1-sr))
pred_till   <- des_till %>% mutate(pred_new_t = 
                                    (newt_beta[,"lam.f"]*F) / (1+ newt_beta[,"b.f"]*(F + newt_beta[,"a."]*M)) + 
                                    (newt_beta[,"lam.m"]*M) / (1+ newt_beta[,"b.m"]*(F + newt_beta[,"a."]*M))
                                     )
pred_till   <- mutate(pred_till, pred_new_t_pc = pred_new_t/TotDensity)
pred_till   <- mutate(pred_till, pred_new_t_pf = pred_new_t/(TotDensity*sr) )
pred_till   <- mutate(pred_till, pred_new_t_pf = replace(pred_new_t_pf, pred_new_t_pf==Inf, NA))

                        
# COMBINE MODELS #####################################################################

# PANEL 1: Per capita/total seed production ------------------------------------------------
seed_prod <- merge(select(subset(pred_flower, sexm == 0), sr,TotDensity,n_flowers_pc), 
                   select(pred_seed, sr,TotDensity,seeds_per_f) )
seed_prod <- mutate(seed_prod, seed_prod    = seeds_per_f * n_flowers_pc * sr * TotDensity)
seed_prod <- mutate(seed_prod, seed_prod_pf = seed_prod / (sr * TotDensity))
seed_prod <- mutate(seed_prod, seed_prod_pc = seed_prod / TotDensity)

# PANEL 2: Per capita/total seed production ------------------------------------------------

# Number of female flowers (panicles) per plot 
females     <- subset(pred_flower, sexm == 0)
females     <- mutate(females, n_fem       = sr * TotDensity) # n. of female flowers in plot
females     <- mutate(females, n_f_flowers = n_fem * n_flowers_pc) #n. of fem. panicles
# Number of male flowers (panicles) per plot 
males       <- subset(pred_flower, sexm == 1)
males       <- mutate(males, n_mal       = (1-sr) * TotDensity) # n. of female flowers in plot
males       <- mutate(males, n_m_flowers = n_mal * n_flowers_pc) #n. of fem. panicles

# combine the data sets
flow_num    <- merge(dplyr::select(females,TotDensity,sr,n_fem,n_f_flowers), 
                     dplyr::select(males,  TotDensity,sr,n_mal,n_m_flowers))


#design matrix to calculate seed viability (proportion) 
des_v     <- flow_num %>% mutate(Intercept = 1,
                                  totFlow = n_m_flowers + n_f_flowers,
                                  sr_f = n_f_flowers / (n_m_flowers + n_f_flowers)
                                 )
des_v     <- des_v %>% 
              mutate(sr_f_x_totFlow = totFlow * sr_f)
names(des_v)[c(7,10)] <- c("(Intercept)","sr_f:totFlow")
viab_raw <- pred(des_v[,c("(Intercept)","sr_f","sr_f:totFlow","totFlow")], 
                  germ_beta, viab, inv.logit)
all.equal(select(viab_raw,sr_f,totFlow),select(des_v,sr_f,totFlow))
pred_viab  <- des_v %>%
                select(TotDensity,sr, n_fem,n_f_flowers, n_mal,n_m_flowers) %>%
                cbind(viab_raw)
  
# PANEL 3: Viable Seed production ----------------------------------------------------

# n. of flowers + viability
fert      <- merge(pred_viab, seed_prod)
fert      <- mutate(fert, viable_seeds     = seed_prod * viab,
                          viable_seeds_pf  = (seed_prod * viab)/n_fem)

# reproduction
s_d       <- 1 #per capita seedling death
reprod    <- merge(fert,pred_till)
reprod    <- mutate(reprod, reprod    = (viable_seeds * s_d) + pred_new_t)
reprod    <- mutate(reprod, reprod_pf = reprod /n_fem )
reprod    <- mutate(reprod, reprod_pc = reprod /TotDensity)

# sex ratios at which reproduction maximizes 
low_d     <- subset(reprod, TotDensity == 1)
max_low_d <- which(low_d$reprod_pc == max(low_d$reprod_pc) )
high_d    <- subset(reprod, TotDensity == 48)
max_high_d<- which(high_d$reprod_pc == max(high_d$reprod_pc) )
c(low_d$sr[max_low_d], high_d$sr[max_high_d])

# FIGURE 2 ------------------------------------------------------------------------------

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
      zlab = expression(Psi),mgp=c(2,1,0), sub = "a",
      main = expression("a)   Reproductive potential ("*Psi*")"))
persp(viab_3d$x,viab_3d$y,viab_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", line = m_line,
      ylab = "Density",ticktype = "detailed",
      zlab = expression(omega),
      main = expression("b)     Seed viability("*omega*")"))
persp(till_3d$x,till_3d$y,till_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", line = m_line,
      ylab = "Density",ticktype = "detailed",
      zlab = expression(gamma),
      main = expression("c)    Production of new tillers("*gamma*")"))
persp(vs_3d$x,vs_3d$y,vs_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",line = m_line,
      zlab = "R",
      main = expression("d)   Regeneration(R)"))

dev.off()


# countour
tiff("Results/VitalRates_3/figure2_contour.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

par(mfrow=c(2,2),mar=c(2,2,1,0.2),mgp=c(2,0.5,0), oma=c(1,1,0,0),cex=0.7,
    cex.axis = 0.7)

contour(seed_3d$x,seed_3d$y,seed_3d$z, 
        main = expression("a) Reproductive potential ("*Psi*")"))
contour(viab_3d$x,viab_3d$y,viab_3d$z,
        main = expression("b) Seed viability ("*omega*")"))
contour(till_3d$x,till_3d$y,till_3d$z,
        main = expression("c) Production of new tillers ("*gamma*")"))
contour(vs_3d$x,vs_3d$y,vs_3d$z,
        main = expression("d) Regeneration (R)"))

mtext("Density",side=2,line=-0.2,outer=T,cex=1)
mtext("Proportion of female individuals",at=0.5,line=-0.2,side=1,outer=T,cex=1)

dev.off()


# FIGURE 2 COUNTOUR ##################################################################################
# NOTE: you need to export this file MANUALLY using R studio.
par(opar)
dev.off()
# gplots has the function colorpanel, which is handy for making gray-scale contour plots
library(gplots)

# Data 
seed_3d <- form_3d_surf(seed_prod,seed_prod_pc)
viab_3d <- form_3d_surf(pred_viab,viab)
till_3d <- form_3d_surf(pred_till,pred_new_t_pc)
vs_3d   <- form_3d_surf(reprod,reprod_pc)

# coordinates
y_rng_1 <- c(0.56,0.95)
y_rng_2 <- c(0.08,0.47)
x_rng_1 <- c(0.08,0.4)
x_rng_2 <- c(0.56,0.9)

# other graphic parameters
off_m   <- 0.05
mgps    <- c(1,0.3,0)

#plot.new() is necessary if using the modified versions of filled.contour
plot.new()


# Reproductive potential --------------------------------------------
xcoords <- unique(seed_3d$x)
ycoords <- unique(seed_3d$y)
surface.matrix <- seed_3d$z


#I am organizing where the plots appear on the page using the "plt" argument in "par()"
par(new = "TRUE", plt = c(x_rng_1,y_rng_1),
    las = 1, cex.axis = 0.8, tck = -0.02, mgp = mgps) 

gray.colors_me <- function(n, start = 0, end =1, gamma = 2.2, alpha = NULL){
  gray.colors(n = n, start = start, end = end, gamma = gamma, alpha = alpha)
}

# Top left plot:
filled.contour3(xcoords,
                ycoords,
                surface.matrix,
                color=gray.colors_me,
                xlab = "",        # suppress x-axis annotation
                ylab = "",        # suppress y-axis annotation
                xlim = c(min(xcoords),max(xcoords)),
                ylim = c(min(ycoords),max(ycoords)),
                zlim = c(min(surface.matrix),max(surface.matrix))
)

text_coord <- par("usr")
text_y_d   <- (text_coord[4]-text_coord[3])*0.05
text(x=text_coord[1],y=text_coord[4]+text_y_d,xpd = NA,
     expression("a) Initiated seeds ("*italic(s)*")"),
     cex = 1.2,font = 2, pos = 4, offset = off_m)

# Y-axis for the whole graph (mtext does not function here)
par(new = "TRUE", plt = c(0,1,0,1)) 
text(x=-0.17,y=-5,"Density",cex=1.5,xpd = NA, srt=90)

#Add a legend:
par(new = "TRUE", mgp = mgps, 
    plt = c(x_rng_1[2]+0.02,x_rng_1[2]+0.06,
            y_rng_1),   # define plot region for legend
    las = 1, cex.axis = 0.8)

#
filled.legend(
  xcoords,
  ycoords,
  surface.matrix,
  color = gray.colors_me,
  xlab = "",
  ylab = "",
  xlim = c(min(xintercepts),max(xintercepts)),
  ylim = c(min(slopes),max(slopes)),
  zlim = range(surface.matrix))

# SECOND --------------------------------------------------------------
xcoords <- unique(viab_3d$x)
ycoords <- unique(viab_3d$y)
surface.matrix <- viab_3d$z

#I am organizing where the plots appear on the page using the "plt" argument in "par()"
par(new = "TRUE", plt = c(x_rng_2,y_rng_1), mgp = mgps,   
    las = 1, cex.axis = 0.8, tck = -0.02 )             

# Top left plot:
filled.contour3(xcoords,
                ycoords,
                surface.matrix,
                color=gray.colors_me,
                xlab = "",        # suppress x-axis annotation
                ylab = "",        # suppress y-axis annotation
                xlim = c(min(xcoords),max(xcoords)),
                ylim = c(min(ycoords),max(ycoords)),
                zlim = c(min(surface.matrix),max(surface.matrix)))

text_coord <- par("usr")
text_y_d   <- (text_coord[4]-text_coord[3])*0.05
text(x=text_coord[1],y=text_coord[4]+text_y_d,xpd = NA,
     expression("b) Seed viability ("*italic(z)*")"),
     cex = 1.2,font = 2, pos = 4, offset = off_m)


#Add a legend:
par(new = "TRUE", mgp = mgps, 
    plt = c(0.91,0.95,y_rng_1),   # define plot region for legend
    las = 1,
    cex.axis = 0.8)
#
filled.legend(
  xcoords,
  ycoords,
  surface.matrix,
  color = gray.colors_me,
  xlab = "",
  ylab = "",
  xlim = c(min(xintercepts),max(xintercepts)),
  ylim = c(min(slopes),max(slopes)),
  zlim = range(surface.matrix))


# THIRD --------------------------------------------------------------
xcoords <- unique(till_3d$x)
ycoords <- unique(till_3d$y)
surface.matrix <- till_3d$z

#I am organizing where the plots appear on the page using the "plt" argument in "par()"
par(new = "TRUE", plt = c(x_rng_1,y_rng_2), mgp = mgps,  
    las = 1, cex.axis = 0.8, tck = -0.02 )                 

# Top left plot:
filled.contour3(xcoords,
                ycoords,
                surface.matrix,
                color=gray.colors_me,
                xlab = "",        # suppress x-axis annotation
                ylab = "",        # suppress y-axis annotation
                xlim = c(min(xcoords),max(xcoords)),
                ylim = c(min(ycoords),max(ycoords)),
                zlim = c(min(surface.matrix),max(surface.matrix)))


text_coord <- par("usr")
text_y_d   <- (text_coord[4]-text_coord[3])*0.05
text(x=text_coord[1],y=text_coord[4]+text_y_d,xpd = NA,
     expression("c) Asexual recruits ("*italic(a)*")"),
     cex = 1.2,font = 2, pos = 4, offset = off_m)


#Add a legend:
par(new = "TRUE", mgp = mgps, 
    plt = c(x_rng_1[2]+0.02,x_rng_1[2]+0.06,
            y_rng_2),   # define plot region for legend
    las = 1, cex.axis = 0.8)

#
filled.legend(
  xcoords,
  ycoords,
  surface.matrix,
  color = gray.colors_me,
  xlab = "",
  ylab = "",
  xlim = c(min(xintercepts),max(xintercepts)),
  ylim = c(min(slopes),max(slopes)),
  zlim = range(surface.matrix))


# FOUR --------------------------------------------------------------
xcoords <- unique(vs_3d$x)
ycoords <- unique(vs_3d$y)
surface.matrix <- vs_3d$z

#I am organizing where the plots appear on the page using the "plt" argument in "par()"
par(new = "TRUE", plt = c(x_rng_2,y_rng_2), mgp = mgps, 
    las = 1,cex.axis = 0.8, tck = -0.02 ) 

# Top left plot:
filled.contour3(xcoords,
                ycoords,
                surface.matrix,
                color=gray.colors_me,
                xlab = "",        # suppress x-axis annotation
                ylab = "",        # suppress y-axis annotation
                xlim = c(min(xcoords),max(xcoords)),
                ylim = c(min(ycoords),max(ycoords)),
                zlim = c(min(surface.matrix),max(surface.matrix)))


text_coord <- par("usr")
text_y_d   <- (text_coord[4]-text_coord[3])*0.05
text(x=text_coord[1],y=text_coord[4]+text_y_d,xpd = NA,
     expression("d) Total recruitment ("*italic(sz)*" + "*italic(a)*")"),cex = 1.2,pos = 4, offset = off_m)


#Add a legend:
par(new = "TRUE", mgp = mgps, 
    plt = c(0.91,0.95,y_rng_2),   # define plot region for legend
    las = 1, cex.axis = 0.8)
#
filled.legend(
  xcoords,
  ycoords,
  surface.matrix,
  color = gray.colors_me,
  xlab = "",
  ylab = "",
  xlim = c(min(xintercepts),max(xintercepts)),
  ylim = c(min(slopes),max(slopes)),
  zlim = range(surface.matrix))

mtext("Proportion of female individuals",side=1, 
      line=-1, outer=T, adj = c(0.43) ,cex=1.5)





# Varying Per capita seedling survival ---------------------------------
s_d       <- 1 #per capita seedling death
reprod    <- merge(fert,pred_till)
reprod    <- mutate(reprod, reprod    = (viable_seeds * s_d) + pred_new_t)
reprod    <- mutate(reprod, reprod_pf = reprod /n_fem )
reprod    <- mutate(reprod, reprod_pc = reprod /TotDensity)
vs_3d_1   <- form_3d_surf(reprod,reprod_pc)

s_d       <- 0.7 #per capita seedling death
reprod    <- merge(fert,pred_till)
reprod    <- mutate(reprod, reprod    = (viable_seeds * s_d) + pred_new_t)
reprod    <- mutate(reprod, reprod_pf = reprod /n_fem )
reprod    <- mutate(reprod, reprod_pc = reprod /TotDensity)
vs_3d_2   <- form_3d_surf(reprod,reprod_pc)

s_d       <- 0.46 #per capita seedling death
reprod    <- merge(fert,pred_till)
reprod    <- mutate(reprod, reprod    = (viable_seeds * s_d) + pred_new_t)
reprod    <- mutate(reprod, reprod_pf = reprod /n_fem )
reprod    <- mutate(reprod, reprod_pc = reprod /TotDensity)
vs_3d_3   <- form_3d_surf(reprod,reprod_pc)

s_d       <- 0.31 #per capita seedling death
reprod    <- merge(fert,pred_till)
reprod    <- mutate(reprod, reprod    = (viable_seeds * s_d) + pred_new_t)
reprod    <- mutate(reprod, reprod_pf = reprod /n_fem )
reprod    <- mutate(reprod, reprod_pc = reprod /TotDensity)
vs_3d_4   <- form_3d_surf(reprod,reprod_pc)


tiff("Results/VitalRates_3/figure2_seedling_surv_estimates.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

m_line <- -2
par(mfrow=c(2,2),mar=c(0,2,0,0),mgp=c(3,1,0), oma=c(0,0,0,0),cex=0.7,
    cex.axis = 0.7)

persp(vs_3d_1$x,vs_3d_1$y,vs_3d_1$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals",line = m_line,
      ylab = "Density",ticktype = "detailed",
      zlab = "Regeneration",mgp=c(2,1,0), sub = "a",
      main = expression("a)   100% seedling surv"))
persp(vs_3d_2$x,vs_3d_2$y,vs_3d_2$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", line = m_line,
      ylab = "Density",ticktype = "detailed",
      zlab = "Regeneration",
      main = expression("b)   70% seedling surv"))
persp(vs_3d_3$x,vs_3d_3$y,vs_3d_3$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", line = m_line,
      ylab = "Density",ticktype = "detailed",
      zlab = "Regeneration",
      main = expression("c)   20% seedling surv"))
persp(vs_3d_4$x,vs_3d_4$y,vs_3d_4$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",line = m_line,
      zlab = "Regeneration",
      main = expression("d)   5% seedling surv"))

dev.off()



par(mfrow=c(2,2),mar=c(2,2,1,0.2),mgp=c(2,0.5,0), oma=c(1,1,0,0),cex=0.7,
    cex.axis = 0.7)

contour(vs_3d_1$x,vs_3d_1$y,vs_3d_1$z, 
        main = expression("a) Regeneration: seedling survival 100%"))
contour(vs_3d_2$x,vs_3d_2$y,vs_3d_2$z,
        main = expression("b) Regeneration: seedling survival 70%"))
contour(vs_3d_3$x,vs_3d_3$y,vs_3d_3$z,
        main = expression("c) Regeneration: seedling survival 20%"))
contour(vs_3d_4$x,vs_3d_4$y,vs_3d_4$z,
        main = expression("d) Regeneration: seedling survival 5%"))

mtext("Density",side=2,line=-0.2,outer=T,cex=1)
mtext("Proportion of female individuals",at=0.5,line=-0.2,side=1,outer=T,cex=1)




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


# contour
par(mar=c(4,4,1,0.5),mgp=c(2,0.7,0))
filled.contour(vs_3d$x,vs_3d$y,vs_3d$z, color.palette=heat.colors, cex.lab = 1.4,
               xlab = "Proportion of female individuals", 
               ylab = "Planting density", 
               main = "Number of viable seeds")

