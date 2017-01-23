# Code to produce Figure 2 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(dplyr)
library(boot)
library(testthat) # to test "order/correspondence of predictors" 

# Read best models --------------------------------------------------------------------
f_pan_avg   <- read.csv("Results/VitalRates_3/f_flowers_plot_best.csv", stringsAsFactors = F)
m_pan_avg   <- read.csv("Results/VitalRates_3/m_flowers_plot_best.csv", stringsAsFactors = F)
fec_beta    <- read.csv("Results/VitalRates_3/fecuntity_best.csv", stringsAsFactors = F)
germ_beta   <- read.csv("Results/VitalRates_3/germination_best.csv", stringsAsFactors = F) 

# Format beta tables -------------------------------------------------------------
format_beta = function(x){
  if(any(x$predictor == "totFlow") ){
    r           <- match( c("(Intercept)", "totFlow", "sr_f", "sr_f:totFlow"), x$predictor )
  } else{
    x$predictor <- gsub("TotDensity","N",x$predictor)
    x$predictor <- gsub("sr:N","N:sr",x$predictor)
    r           <- match( c("(Intercept)", "N", "sr", "N:sr"), x$predictor )
  }
  return(x[r,])
}
fec_beta  <- format_beta(fec_beta)
germ_beta <- format_beta(germ_beta)


# MODEL PREDICTIONS ---------------------------------------------------------------

# Design matrix
design      <- expand.grid("(Intercept)" = 1,
                           N = seq(1,48,1), sr = seq(0,1,0.01) )
design      <- mutate(design, "N:sr" = sr * N)

# Prediction (Per capita number of flowers)
expect_equal( all(names(design) == f_pan_avg$predictor), TRUE )
expect_equal( all(names(design) == m_pan_avg$predictor), TRUE )
pred_pan    <- mutate(design, 
                      n_f_flow = as.vector( exp(as.matrix(design) %*% f_pan_avg$avg )),
                      n_m_flow = as.vector( exp(as.matrix(design) %*% m_pan_avg$avg )))

# Seeds per flower (fecundity model, Figure 1b) -----------------------------------
expect_equal( all(names(design) == fec_beta$predictor), TRUE )
pred_fec <- mutate(design, seeds_per_f = 
                           as.vector(exp( as.matrix(design) %*% fec_beta$avg )) )


# Seed viability (viability model, Figure 1d) -------------------------------------
pred_v      <- pred_pan %>% mutate(totFlow     = n_f_flow + n_m_flow,
                                   sr_f        = n_f_flow / totFlow)
pred_v      <- mutate(pred_v,   "sr_f:totFlow" = totFlow * sr_f)
viab_design <- as.matrix( pred_v[,c("(Intercept)","totFlow","sr_f","sr_f:totFlow")] )
expect_equal( all(colnames(viab_design) == germ_beta$predictor), TRUE )
pred_viab   <- mutate(pred_v, 
                      pred_viab  = as.vector( inv.logit( viab_design %*% germ_beta$avg)) )


# Final data frame ----------------------------------------------------------------

# Seeds per plot
tmp           <- merge(pred_viab,
                       dplyr::select(pred_fec,   sr, N, seeds_per_f))
seed_x_plot   <- mutate(tmp, n_seeds = n_f_flow * seeds_per_f)

# fertility (viable seeds)
fert          <- mutate(seed_x_plot, viable_seeds = pred_viab * n_seeds)
fert          <- fert[order(fert$N,fert$sr),] #order requires by persp & contour


fert          <- mutate(fert, pc_viab_seeds = viable_seeds/n_f_flow )
fert$pc_viab_seeds[fert$pc_viab_seeds == Inf] = NA

# FIGURE 2 ------------------------------------------------------------------------

# Prepare data
x<- unique(fert$sr)
y<- unique(fert$N)
z<-matrix(fert$viable_seeds, nrow=length(unique(fert$sr)),
          ncol=length(unique(fert$N)))
z1<-matrix(fert$pc_viab_seeds, nrow=length(unique(fert$sr)),
          ncol=length(unique(fert$N)))



tiff("Results/VitalRates_3/figure2_pois.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

par(mfrow=c(1,1))
persp(x,y,z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",
      zlab = "Per capita viable seeds",
      main = "Fertility")

dev.off()

# Alternatives to persp() function
tiff("Results/VitalRates_3/figure2_contour_pois.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

par(mar=c(4,4,1,0.5),mgp=c(2,0.7,0))
filled.contour(x,y,z, color.palette=heat.colors, cex.lab = 1.4,
               xlab = "Proportion of female individuals", 
               ylab = "Density", 
               main = "Viable seeds")

dev.off()



# graph shows separate repronse variables ---------------------------------------------
par(mfrow=c(2,2))

# number of flowers
one <- subset(pred_fec, sr == 0.2)
two <- subset(pred_fec, sr == 1)

plot(one$N, one$seeds_per_f, ylim=c(0,300), type="l") #; plot(new=T)
par(new=T)
plot(two$N, two$seeds_per_f, ylim=c(0,300), type="l",lty=2)

# viability
one <- subset(pred_viab, sr == 0.2)
two <- subset(pred_viab, sr == 1)

plot(one$N, one$pred_viab,ylim=c(0,1),type="l") #; plot(new=T)
par(new=T)
plot(two$N, two$pred_viab,ylim=c(0,1),type="l",lty=2)

# number of (female) flowers
one <- subset(pred_pan, sr == 0.2)
two <- subset(pred_pan, sr == 1)

plot(one$N, one$n_f_flow,ylim=c(0,100),type="l") #; plot(new=T)
par(new=T)
plot(two$N, two$n_f_flow,ylim=c(0,100),type="l",lty=2)


# number of viable seeds
one <- subset(fert, sr == 0.2)
two <- subset(fert, sr == 1)

plot(one$N, one$viable_seeds,ylim=c(0,3000),type="l") #; plot(new=T)
par(new=T)
plot(two$N, two$viable_seeds,ylim=c(0,3000),type="l",lty=2)
