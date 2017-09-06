# Code to produce Figure 2 
setwd("C:/cloud/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(dplyr)
library(boot)
library(testthat) # to test "order/correspondence of predictors" 
library(glmmADMB)
source("C:/CODE/LLELA/model_avg.R")
source("C:/CODE/LLELA/prediction.R")
source("C:/CODE/LLELA/unload_mass.R")
source("C:/CODE/LLELA/filled.contour3.R")
source("http://wiki.cbr.washington.edu/qerm/sites/qerm/images/2/25/Filled.legend.R")


# Read in best models / new data(for viability)----------------------------------------------------------
n_flow_beta <- read.csv("Results/VitalRates_4_Cade2015/n_flowers_best.csv", stringsAsFactors = F)
fec_beta    <- read.csv("Results/VitalRates_4_Cade2015/fecuntity_best.csv", stringsAsFactors = F)
germ_beta   <- read.csv("Results/VitalRates_4_Cade2015/germination_best.csv", stringsAsFactors = F) 
newt_beta   <- read.csv("Results/VitalRates_4_Cade2015/new_t_bh14_best.csv", stringsAsFactors = F)


# "remember" graphical parameters
opar  <- par()

# change prediction names
n_flow_beta <- rename(n_flow_beta, panicles_pc  = pred )
fec_beta    <- rename(fec_beta,    s_per_panic  = pred )
germ_beta   <- rename(germ_beta,   pred_viab    = pred )
newt_beta   <- rename(newt_beta,   pred_new_t   = pred,
                                   TotDensity   = TotN )


# PANEL 1: Per capita/total seed production ------------------------------------------------
seed_prod   <- full_join( subset(n_flow_beta, sex == "f"), fec_beta ) %>%
                  mutate( F = round(TotDensity * sr, 4) ) %>%
                  mutate( seed_prod = F * panicles_pc * s_per_panic ) %>%
                  mutate( seed_prod_pc = seed_prod / TotDensity,
                          seed_prod_pf = seed_prod / F ) 


# PANEL 2: Seed viability ------------------------------------------------------------------

# Number of female flowers (panicles) per plot 
females     <- subset(n_flow_beta, sex == "f" ) %>%
                  mutate(F = round(sr * TotDensity,4) ) %>%
                  mutate(n_f_panic = F * panicles_pc)
# Number of male flowers (panicles) per plot 
males       <- subset(n_flow_beta, sex == "m") %>%
                  mutate( M = round((1-sr) * TotDensity,4) ) %>%
                  mutate( n_m_panic = M * panicles_pc)
# combine the data sets
flow_num    <- merge(dplyr::select(females,TotDensity,sr,F,n_f_panic), 
                     dplyr::select(males,  TotDensity,sr,M,n_m_panic))


# get viability estimates "ex novo"
d           <- read.csv("Data/vr.csv")
viabVr      <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv", stringsAsFactors = F)
# Format data 
viabVr      <- viabVr %>%
                  mutate(totN = F + M) %>%
                  mutate(sr = F / totN) %>%
                  mutate(perm_id = paste(plot, focalI, P_ID, sep = "_") ) %>%
                  mutate(plot = as.factor(plot) ) %>%
                  subset(totFlow < 60) #exclude extreme values 
germ_dat    <- na.omit(dplyr::select(viabVr,germTot,germFail,sr_f,totFlow,sr,totN,plot,focalI,P_ID,log_l_t0, F, M, perm_id))

# model structures
candidate_mods <- list(
  cbind(germTot,germFail) ~ (1 | plot),
  cbind(germTot,germFail) ~ sr_f + (1 | plot),
  cbind(germTot,germFail) ~ totFlow + (1 | plot),
  cbind(germTot,germFail) ~ sr_f + totFlow + (1 | plot),
  cbind(germTot,germFail) ~ sr_f * totFlow + (1 | plot)
)
candidate_mods <- setNames(candidate_mods, paste0("model",1:5))
# fit BEST model, no bootstrap
germ_mod_l      <- lapply(candidate_mods, function(x) glmmadmb(x, family="binomial", data=germ_dat))
germ_mod_sel    <- AICtab(germ_mod_l,weights=T)
# extract model rank
mod_rank <- gsub("model", "", attributes(germ_mod_sel)$row.names) %>% as.numeric
# Models that make up more than 95% of weight
k <- 0 ; sumW <- 0 
while(sumW < 0.95){ 
  k <- k + 1
  sumW <- sum(germ_mod_sel$weight[1:k])
}
# Prepare "newdata" 
new_data    <- flow_num %>% 
                  mutate(totFlow = n_m_panic + n_f_panic,
                         sr_f    = n_f_panic / (n_m_panic + n_f_panic) )
# create prediction data frame
pred_l    <- lapply(germ_mod_l[mod_rank[1:k]],
                  function(x) predict(x, newdata = new_data, re.form = NA, type = "response") )
preds     <- ( do.call(cbind, pred_l) %*% germ_mod_sel$weight[1:k] ) / sum( germ_mod_sel$weight[1:k] )
pred_viab <- mutate(new_data, viab = as.numeric(preds) )
unload_mass()


 
# PANEL 3: new tillers ---------------------------------------------------------------------

# round F and M (to allow merging)
newt_beta   <- mutate(newt_beta, 
                      F = round(F, 4),
                      M = round(M, 4))
pred_till   <- mutate( newt_beta, 
                       pred_new_t_pc = pred_new_t / TotDensity,
                       pred_new_t_pf = pred_new_t/ (TotDensity*sr) )
pred_till   <- mutate(pred_till, pred_new_t_pf = replace(pred_new_t_pf, pred_new_t_pf==Inf, NA)) %>%
                      dplyr::select(TotDensity:pred_new_t_pf )



# PANEL 4: Viable Seed production ----------------------------------------------------------

# n. of flowers + viability
fert      <- merge(pred_viab, seed_prod)
fert      <- mutate(fert, viable_seeds     = seed_prod * viab,
                          viable_seeds_pf  = (seed_prod * viab)/ F)

# reproduction
s_d       <- 1 #per capita seedling death
reprod    <- merge(fert, pred_till)
reprod    <- mutate(reprod, reprod    = (viable_seeds * s_d) + pred_new_t)
reprod    <- mutate(reprod, reprod_pf = reprod / F)
reprod    <- mutate(reprod, reprod_pc = reprod / TotDensity)

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

# set up graphic parameters
par(opar)
dev.off()
# gplots has the function colorpanel, which is handy for making gray-scale contour plots
library(gplots)

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



#I am organizing where the plots appear on the page using the "plt" argument in "par()"
tiff("Results/VitalRates_4_Cade2015/Figure2_grey_contour.tiff",
     unit="in",width=6.3,height=6.3,res=600,compression="lzw")


par(new = "TRUE", plt = c(x_rng_1,y_rng_1),
    las = 1, cex.axis = 0.8, tck = -0.02, mgp = mgps) 

gray.colors_me <- function(n, start = 0, end =1, gamma = 2.2, alpha = NULL){
  gray.colors(n = n, start = start, end = end, gamma = gamma, alpha = alpha)
}

# Reproductive potential --------------------------------------------
xcoords <- unique(seed_3d$x)
ycoords <- unique(seed_3d$y)
surface.matrix <- seed_3d$z


# Top left plot:
filled.contour3(xcoords,
                ycoords,
                surface.matrix,
                color=gray.colors_me,
                xlab = "",        # suppress x-axis annotation
                ylab = "",        # suppress y-axis annotation
                xlim = c(min(xcoords),max(xcoords)),
                ylim = c(min(ycoords),max(ycoords)),
                zlim = c(min(surface.matrix),max(surface.matrix)) )

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

mtext("Sex ratio (proportion female)",side=1, 
      line=-1, outer=T, adj = c(0.43) ,cex=1.5)

dev.off()



# Varying Per capita seedling survival ---------------------------------
s_d       <- 1 #per capita seedling death
reprod    <- merge(fert,pred_till)
reprod    <- mutate(reprod, reprod    = (viable_seeds * s_d) + pred_new_t)
reprod    <- mutate(reprod, reprod_pf = reprod / F )
reprod    <- mutate(reprod, reprod_pc = reprod /TotDensity)
vs_3d_1   <- form_3d_surf(reprod,reprod_pc)

s_d       <- 0.7 #per capita seedling death
reprod    <- merge(fert,pred_till)
reprod    <- mutate(reprod, reprod    = (viable_seeds * s_d) + pred_new_t)
reprod    <- mutate(reprod, reprod_pf = reprod / F )
reprod    <- mutate(reprod, reprod_pc = reprod /TotDensity)
vs_3d_2   <- form_3d_surf(reprod,reprod_pc)

s_d       <- 0.46 #per capita seedling death
reprod    <- merge(fert,pred_till)
reprod    <- mutate(reprod, reprod    = (viable_seeds * s_d) + pred_new_t)
reprod    <- mutate(reprod, reprod_pf = reprod / F )
reprod    <- mutate(reprod, reprod_pc = reprod /TotDensity)
vs_3d_3   <- form_3d_surf(reprod,reprod_pc)

s_d       <- 0.31 #per capita seedling death
reprod    <- merge(fert,pred_till)
reprod    <- mutate(reprod, reprod    = (viable_seeds * s_d) + pred_new_t)
reprod    <- mutate(reprod, reprod_pf = reprod / F )
reprod    <- mutate(reprod, reprod_pc = reprod /TotDensity)
vs_3d_4   <- form_3d_surf(reprod,reprod_pc)

# NOTE: you need to export this file MANUALLY using R studio.
par(opar)
dev.off()

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


tiff("Results/VitalRates_4_Cade2015/figure2_seedling_surv_estimates.tiff",
     unit="in",width=6.3,height=6.3,res=600,compression="lzw")

par(new = "TRUE", plt = c(x_rng_1,y_rng_1),
    las = 1, cex.axis = 0.8, tck = -0.02, mgp = mgps) 

gray.colors_me <- function(n, start = 0, end =1, gamma = 2.2, alpha = NULL){
  gray.colors(n = n, start = start, end = end, gamma = gamma, alpha = alpha)
}

# Reproductive potential --------------------------------------------
xcoords <- unique(vs_3d_1$x)
ycoords <- unique(vs_3d_1$y)
surface.matrix <- vs_3d_1$z


# Top left plot:
filled.contour3(xcoords,
                ycoords,
                surface.matrix,
                color=gray.colors_me,
                xlab = "",        # suppress x-axis annotation
                ylab = "",        # suppress y-axis annotation
                xlim = c(min(xcoords),max(xcoords)),
                ylim = c(min(ycoords),max(ycoords)),
                zlim = c(min(surface.matrix),max(surface.matrix)) )

text_coord <- par("usr")
text_y_d   <- (text_coord[4]-text_coord[3])*0.05
text(x=text_coord[1],y=text_coord[4]+text_y_d,xpd = NA,
     expression("a) Recruits: seedling survival 100%"),
     cex = 1,font = 2, pos = 4, offset = off_m)

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
xcoords <- unique(vs_3d_2$x)
ycoords <- unique(vs_3d_2$y)
surface.matrix <- vs_3d_2$z

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
     expression("b) Recruits: seedling survival 70%"),
     cex = 1,font = 2, pos = 4, offset = off_m)


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
xcoords <- unique(vs_3d_3$x)
ycoords <- unique(vs_3d_3$y)
surface.matrix <- vs_3d_3$z

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
     expression("c) Recruits: seedling survival 20%"),
     cex = 1,font = 2, pos = 4, offset = off_m)


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
xcoords <- unique(vs_3d_4$x)
ycoords <- unique(vs_3d_4$y)
surface.matrix <- vs_3d_4$z

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
     expression("d) Recruits: seedling survival 5%"),cex = 1,pos = 4, offset = off_m)


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

mtext("Sex ratio (proportion female)",side=1, 
      line=-1, outer=T, adj = c(0.43) ,cex=1.5)

dev.off()
