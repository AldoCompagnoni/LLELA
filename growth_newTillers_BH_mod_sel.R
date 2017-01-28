##Response variable: total number of tillers 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/unload_mass.R")
unload_mass()
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/ac79/Documents/CODE/LLELA/bh_util_fnc.R")


# load and format data -----------------------------------------------------
x       <- read.csv("Data/vr.csv")
d       <- format_growth(x)
d14     <- subset(d, year == 2014)
d15     <- subset(d, year == 2015)
d15     <- subset(d15, plot != 149) # Remove outlier
# omit NAs
d14 <- na.omit(select(d14, TotDensity, sr, new_t1))
d15 <- na.omit(select(d15, TotDensity, sr, new_t1))


# Model fits --------------------------------------------------------------------------------------------

# 2014
BH.14         <-optim(par=setNames(c(2,0.1,5),c("lam","b","size")),
                      fn=fit.BH, gr=NULL, d14, control=list(maxit=5000))
BH.sex.14     <-optim(par=setNames(c(10,10,0.5,0.5,5),c("lam_f","lam_m","b_f","b_m","size")),
                      fn=fit.BH.sex,gr=NULL, d14,control=list(maxit=5000))
BH.sex.lam.14 <-optim(par=setNames(c(10,10,0.5,5),c("lam_f","lam_m","b","size")),
                      fn=fit.BH.sex.lambda,gr=NULL, d14, control=list(maxit=5000))
BH.sex.b.14   <-optim(par=setNames(c(10,0.5,0.5,5),c("lam","b_f","b_m","size")),
                      fn=fit.BH.sex.b,gr=NULL, d14, control=list(maxit=5000))

# 2015
BH.15         <-optim(par=setNames(c(2,0.1,5),c("lam","b","size")),
                      fn=fit.BH, gr=NULL, d15, control=list(maxit=5000))
BH.sex.15     <-optim(par=setNames(c(10,10,0.5,0.5,5),c("lam_f","lam_m","b_f","b_m","size")),
                      fn=fit.BH.sex,gr=NULL, d15,control=list(maxit=5000))
BH.sex.lam.15 <-optim(par=setNames(c(10,10,0.5,5),c("lam_f","lam_m","b","size")),
                      fn=fit.BH.sex.lambda,gr=NULL, d15, control=list(maxit=5000))
BH.sex.b.15   <-optim(par=setNames(c(10,0.5,0.5,5),c("lam","b_f","b_m","size")),
                      fn=fit.BH.sex.b,gr=NULL, d15, control=list(maxit=5000))


# Compare models --------------------------------------------------------------------------------------

# 2014
m14 <- list()
m14[[1]]  <-2*BH.14$value+2*length(BH.14$par)
m14[[2]]  <-2*BH.sex.lam.14$value+2*length(BH.sex.lam.14$par)
m14[[3]]  <-2*BH.sex.b.14$value+2*length(BH.sex.b.14$par)
m14[[4]]  <-2*BH.sex.14$value+2*length(BH.sex.14$par)
m14       <- setNames(m14, c("AIC.BH.14", "AIC.BH.sex.lam.14", "AIC.BH.sex.b.14", "AIC.BH.sex.14"))

# 2015
m15 <- list()
m15[[1]]  <-2*BH.15$value+2*length(BH.15$par)
m15[[2]]  <-2*BH.sex.lam.15$value+2*length(BH.sex.lam.15$par)
m15[[3]]  <-2*BH.sex.b.15$value+2*length(BH.sex.b.15$par)
m15[[4]]  <-2*BH.sex.15$value+2*length(BH.sex.15$par)
m15       <- setNames(m15, c("AIC.BH.15", "AIC.BH.sex.lam.15", "AIC.BH.sex.b.15", "AIC.BH.sex.15"))

mod_weights(m14)
mod_weights(m15)


# Model selection ----------------------------------------------------------------------------------------

# Only 14 ----------------------------------------------------------------------------------
sr_seq    <- seq(0,1,0.05)
des       <- expand.grid(TotDensity = seq(1,48,1), sr = sr_seq)

# best model 1
beta_1    <- BH.sex.lam.14$par
beta_2    <- BH.14$par
beta_3    <- BH.sex.14$par

lambdas_1 <- beta_1["lam_f"]*sr_seq + beta_1["lam_m"]*(1-sr_seq)
lambdas_3 <- beta_3["lam_f"]*sr_seq + beta_3["lam_m"]*(1-sr_seq)
b_3       <- beta_3["b_f"]*sr_seq   + beta_3["b_m"]*(1-sr_seq)

lambda_df <- data.frame(sr = sr_seq, lambda_1 = lambdas,       b_1 = beta_1["b"],
                                     lambda_2 = beta_2["lam"], b_2 = beta_2["b"],
                                     lambda_3 = lambdas_3,     b_3 = b_3)



des       <- merge(des, lambda_df)
des       <- mutate(des, yhat = (lambda*TotDensity) / (1+beta["b"]*TotDensity))
des       <- mutate(des, yhat_pc = yhat / TotDensity)





till_3d <- form_3d_surf(des,yhat_pc)

persp(till_3d$x,till_3d$y,till_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",
      zlab = "New Tillers_PC",
      main = "Production of new tillers_PC")



# plot
x <- c(0:max(d14$TotDensity))
plot(d14$TotDensity,d14$new_t1, pch = 16)
lines(x,yhat_h, col = "green", lwd = 2)
lines(x,yhat_l, col = "blue", lwd = 2)


# model average ---------------------------------------------------------------------------------------
# best mod
