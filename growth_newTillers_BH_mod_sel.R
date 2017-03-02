##Response variable: total number of tillers 
setwd("C:/Users/Aldo/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)
library(dplyr)
source("C:/Users/Aldo/Documents/CODE/LLELA/unload_mass.R")
unload_mass()
source("C:/Users/Aldo/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/Aldo/Desktop/POAR/bh_util_fnc.R")


# load and format data -----------------------------------------------------
x       <- read.csv("Data/vr.csv")
d       <- format_growth(x)
d14     <- subset(d, year == 2014)
<<<<<<< HEAD
# omit NAs
d14 <- na.omit(select(d14, TotDensity, F, M, new_t1))


# Model fits ---------------------------------------------------------------
=======
d15     <- subset(d, year == 2015)
d15     <- subset(d15, plot != 149) # Remove outlier

# omit NAs
d14 <- na.omit(select(d14, TotDensity, F, M, sr, new_t1))
d15 <- na.omit(select(d15, TotDensity, F, M, sr, new_t1))

# Model fits 2014 ---------------------------------------------------------------
>>>>>>> c71ad217545fb8340e79c46e6eca994ebe95e557
res         <- list()
res[[1]]    <-optim(par=c(2,0.1,5),         fn=fit_null, gr=NULL, d14, control=list(maxit=5000))
res[[2]]    <-optim(par=c(5,5,0.5,5),       fn=fit_lam, gr=NULL, d14, control=list(maxit=5000))
res[[3]]    <-optim(par=c(5,0.5,0.5,5),     fn=fit_b, gr=NULL, d14, control=list(maxit=5000))
res[[4]]    <-optim(par=c(5,0.5,0.1,0.1,5), fn=fit_a, gr=NULL, d14, control=list(maxit=5000))
res[[5]]    <-optim(par=c(5,5,0.5,0.5,5),       fn=fit_lam_b, gr=NULL, d14, control=list(maxit=5000))
res[[6]]    <-optim(par=c(5,5,0.5,0.1,0.1,5),   fn=fit_lam_a, gr=NULL, d14, control=list(maxit=5000))
res[[7]]    <-optim(par=c(5,0.5,0.5,0.1,0.1,5), fn=fit_b_a, gr=NULL, d14, control=list(maxit=5000))
res[[8]]    <-optim(par=c(5,5,0.5,0.5,0.1,0.1,5), fn=fit_full, gr=NULL, d14, control=list(maxit=5000))
res         <-setNames(res, c("null","lam","b","a","lam_b","lam_a","b_a","full"))

# parameter names
res$null    <- par_names(res$null, null_par)
res$lam     <- par_names(res$lam, lam_par)
res$b       <- par_names(res$b, b_par)
res$a       <- par_names(res$a, a_par)
res$lam_b   <- par_names(res$lam_b, lam_b_par)
res$lam_a   <- par_names(res$lam_a, lam_a_par)
res$b_a     <- par_names(res$b_a, b_a_par)
res$full    <- par_names(res$full, full_par)

# aic weights
m14         <- lapply(res,aic_calc)
mod_w       <- mod_weights(m14)

<<<<<<< HEAD
=======
# Model selection table
write.csv(mod_w,"new_tillers_BH14_mod_sel.csv",row.names=F)

# Average AIC weight models --------------------------------------------------------------------------------------
mat <- matrix(data = 1, nrow = length(res), ncol = 6) 
for(i in 1:length(res)){
  
  params  <- c("lam.f","lam.m","b.f","b.m","a.f","a.m")
  
  lam_r   <- grep("lam.",names(res[[i]]$par) )
  b_r     <- grep("b.",names(res[[i]]$par) )
  a_r     <- grep("^a$|a.f|a.m",names(res[[i]]$par) )
  
  if( length(lam_r) > 0 ) mat[i,grep("lam.",params)] <- res[[i]]$par[lam_r]
  if( length(b_r) > 0 )   mat[i,grep("b.",params)]   <- res[[i]]$par[b_r]
  if( length(a_r) > 0 )   mat[i,grep("a.f|a.m",params)]   <- res[[i]]$par[a_r]
  
}
results   <- data.frame(mat)
results   <- setNames(results, params)
results   <- mutate(results, model = names(res) )
results   <- merge(results, mod_w[,c("model","weights")])
results   <- arrange(results, desc(weights) ) # best model first, then the rest
results   <- mutate(results, cum_weight = cumsum(weights))
# only the the models making up to 95% (included) weights
min_cum_weight <- min(results$cum_weight[results$cum_weight > 0.95])
weight_i  <- which(results$cum_weight == min_cum_weight)
weighted  <- (select(results,lam.f:a.m) * results$weights)[1:weight_i,]
summed    <- apply(weighted,2,sum) 
avg       <- summed / results$cum_weight[weight_i]
avg14     <- as.data.frame(t(avg))

write.csv(avg14, "Results/VitalRates_3/new_t_bh14_best.csv",row.names=F)


# Model fits 2015 ---------------------------------------------------------------
res         <- list()
res[[1]]    <-optim(par=c(2,0.1,5),         fn=fit_null, gr=NULL, d15, control=list(maxit=5000))
res[[2]]    <-optim(par=c(5,5,0.5,5),       fn=fit_lam, gr=NULL, d15, control=list(maxit=5000))
res[[3]]    <-optim(par=c(5,0.5,0.5,5),     fn=fit_b, gr=NULL, d15, control=list(maxit=5000))
res[[4]]    <-optim(par=c(5,0.5,0.1,0.1,5), fn=fit_a, gr=NULL, d15, control=list(maxit=5000))
res[[5]]    <-optim(par=c(5,5,0.5,0.5,5),       fn=fit_lam_b, gr=NULL, d15, control=list(maxit=5000))
res[[6]]    <-optim(par=c(5,5,0.5,0.1,0.1,5),   fn=fit_lam_a, gr=NULL, d15, control=list(maxit=5000))
res[[7]]    <-optim(par=c(5,0.5,0.5,0.1,0.1,5), fn=fit_b_a, gr=NULL, d15, control=list(maxit=5000))
res[[8]]    <-optim(par=c(5,5,0.5,0.5,0.1,0.1,5), fn=fit_full, gr=NULL, d15, control=list(maxit=5000))
res         <-setNames(res, c("null","lam","b","a","lam_b","lam_a","b_a","full"))

# parameter names
res$null    <- par_names(res$null, null_par)
res$lam     <- par_names(res$lam, lam_par)
res$b       <- par_names(res$b, b_par)
res$a       <- par_names(res$a, a_par)
res$lam_b   <- par_names(res$lam_b, lam_b_par)
res$lam_a   <- par_names(res$lam_a, lam_a_par)
res$b_a     <- par_names(res$b_a, b_a_par)
res$full    <- par_names(res$full, full_par)

# aic weights
m15         <- lapply(res,aic_calc)
mod_w       <- mod_weights(m15)

# Model selection table
write.csv(mod_w,"new_tillers_BH15_mod_sel.csv",row.names=F)


>>>>>>> c71ad217545fb8340e79c46e6eca994ebe95e557

# Average AIC weight models --------------------------------------------------------------------------------------
mat <- matrix(data = 1, nrow = length(res), ncol = 6) 
for(i in 1:length(res)){
  
  params  <- c("lam.f","lam.m","b.f","b.m","a.f","a.m")
  
  lam_r   <- grep("lam.",names(res[[i]]$par) )
  b_r     <- grep("b.",names(res[[i]]$par) )
  a_r     <- grep("^a$|a.f|a.m",names(res[[i]]$par) )
  
  if( length(lam_r) > 0 ) mat[i,grep("lam.",params)] <- res[[i]]$par[lam_r]
  if( length(b_r) > 0 )   mat[i,grep("b.",params)]   <- res[[i]]$par[b_r]
  if( length(a_r) > 0 )   mat[i,grep("a.f|a.m",params)]   <- res[[i]]$par[a_r]
  
}
results   <- data.frame(mat)
results   <- setNames(results, params)
results   <- mutate(results, model = names(res) )
results   <- merge(results, mod_w[,c("model","weights")])
results   <- arrange(results, desc(weights) ) # best model first, then the rest
results   <- mutate(results, cum_weight = cumsum(weights))
# only the the models making up to 95% (included) weights
min_cum_weight <- min(results$cum_weight[results$cum_weight > 0.95])
weight_i  <- which(results$cum_weight == min_cum_weight)
weighted  <- (select(results,lam.f:a.m) * results$weights)[1:weight_i,]
summed    <- apply(weighted,2,sum) 
avg       <- summed / results$cum_weight[weight_i]
<<<<<<< HEAD
avg       <- as.data.frame(t(avg))

write.csv(avg, "Results/VitalRates_3/new_t_bh_best.csv",row.names=F)
=======
avg15     <- as.data.frame(t(avg))

write.csv(avg15, "Results/VitalRates_3/new_t_bh15_best.csv",row.names=F)
>>>>>>> c71ad217545fb8340e79c46e6eca994ebe95e557


# Model selection ----------------------------------------------------------------------------------------

# Only 14 ----------------------------------------------------------------------------------
sr_seq    <- seq(0,1,0.05)
des       <- expand.grid(TotDensity = seq(1,48,1), sr = seq(0,1,0.05))
des       <- mutate(des, F = sr*TotDensity, M = (1-sr)*TotDensity)

# best model 1
fem_n     <- (avg["lam.f"]*des$F) / (1 + avg["b.f"]*(avg["a.m"]*des$M +            des$F))
mal_n     <- (avg["lam.m"]*des$M) / (1 + avg["b.m"]*(           des$M + avg["a.f"]*des$F))
t_pred    <- fem_n + mal_n
des       <- mutate(des, t_pred = t_pred)
des       <- mutate(des, t_pred_pc = t_pred / TotDensity)

# Graph --------------------------------------
till_3d <- form_3d_surf(des,t_pred_pc)
persp(till_3d$x,till_3d$y,till_3d$z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      xlab = "Proportion of female individuals", 
      ylab = "Density",ticktype = "detailed",
      zlab = "New Tillers_PC",
      main = "Production of new tillers_PC")
<<<<<<< HEAD

=======
>>>>>>> c71ad217545fb8340e79c46e6eca994ebe95e557
