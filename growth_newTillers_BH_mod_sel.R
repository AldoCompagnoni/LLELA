##Response variable: total number of tillers 
setwd("C:/cloud/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)
library(dplyr)
source("C:/CODE/LLELA/unload_mass.R")
unload_mass()
source("C:/CODE/LLELA/model_avg.R")
source("C:/CODE/LLELA/bh_util_fnc.R")

# load and format data -----------------------------------------------------
x   <- read.csv("Data/vr.csv")
d   <- format_new_tillers(x)
d14 <- subset(d, year == 2014)
d15 <- subset(d, year == 2015)


# model averaged prediction ------------------------------------------------

# predicting 
predict_nl <- function(x,des_mat){
  
  pars    <- x$par
  
  lam_par <- pars[grepl("lam.", names(pars))]
  b_par   <- pars[grepl("b.", names(pars))]
  a_par   <- pars[grepl("^a.$", names(pars))]
  
  # numerator --------------------------------------------------------------
  
  if(length(lam_par) == 1 ){
    num_f     <- lam_par * des_mat$F
    num_m     <- lam_par * des_mat$M
  }else{
    num_f     <- lam_par["lam.f"] * des_mat$F
    num_m     <- lam_par["lam.m"] * des_mat$M
  }
  
  # denominator --------------------------------------------------------------
  
  # No "a." parameter
  if( any(names(pars) != "a.") ){
    if( length(b_par) == 1 ){
      
      den_f     <- (1 + b_par * (des_mat$F + des_mat$M) )
      den_m     <- (1 + b_par * (des_mat$F + des_mat$M) )
      
    }else{
      
      den_f     <- (1 + b_par["b.f"] * (des_mat$F + des_mat$M) )
      den_m     <- (1 + b_par["b.m"] * (des_mat$F + des_mat$M) )
      
    }
    # with "a." parameter
  }else{
    if( length(b_par) == 1 ){
      
      den_f     <- (1 + b_par * (des_mat$F + (a_par*des_mat$M) ) )
      den_m     <- (1 + b_par * (des_mat$F + (a_par*des_mat$M) ) )
      
    }else{
      
      den_f     <- (1 + b_par["b.f"] * (des_mat$F + (a_par*des_mat$M) ) )
      den_m     <- (1 + b_par["b.m"] * (des_mat$F + (a_par*des_mat$M) ) )
      
    }
  }
  
  y_pred <- (num_f / den_f) + (num_m / den_m) 
  return(y_pred)
  
}


# Model fits 2014 ---------------------------------------------------------------
res         <- list()
res[[1]]    <- optim(par=c(2,0.1,5),         fn=fit_null, gr=NULL, d14, control=list(maxit=5000))
res[[2]]    <- optim(par=c(5,5,0.5,5),       fn=fit_lam, gr=NULL, d14, control=list(maxit=5000))
res[[3]]    <- optim(par=c(5,0.5,0.5,5),     fn=fit_b, gr=NULL, d14, control=list(maxit=5000))
res[[4]]    <- optim(par=c(5,0.5,0.5,5),     fn=fit_a, gr=NULL, d14, control=list(maxit=5000))
res[[5]]    <- optim(par=c(5,5,0.5,0.5,5),   fn=fit_lam_b, gr=NULL, d14, control=list(maxit=5000))
res[[6]]    <- optim(par=c(5,5,0.5,0.5,5),   fn=fit_lam_a, gr=NULL, d14, control=list(maxit=5000))
res[[7]]    <- optim(par=c(1,0.1,0.1,0.1,5), fn=fit_b_a, gr=NULL, d14, control=list(maxit=5000))
res[[8]]    <- optim(par=c(5,5,0.5,0.5,0.5,5), fn=fit_full, gr=NULL, d14, control=list(maxit=5000))
res         <- setNames(res, c("null","lam","b","a","lam_b","lam_a","b_a","full"))

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


# Average model output 

# design matrix
des_mat <- expand.grid( TotN = c(1:49),
                        sr = round(seq(0,1,by=0.05),2) ) %>%
  mutate( F = sr * TotN ) %>%
  mutate( M = TotN - F )

# store parameters, store predictions
par_res     <- list(null = res$null, lam = res$lam, b = res$b, a = res$a, 
                    lam_b = res$lam_b, lam_a = res$lam_a, b_a = res$b_a, full = res$full)
predict_l   <- lapply(par_res, predict_nl, des_mat) 
predict_df  <- Reduce(function(...) cbind(...), predict_l) %>% 
  as.data.frame %>%
  setNames( names(par_res) ) %>%
  bind_cols( des_mat )
id          <- max( which(cumsum(mod_w$weights) < 0.95) ) + 1
weight_v    <- mod_w$weights[1:id]
model_n     <- as.character(mod_w$model[1:id])
pred_df14   <- predict_df %>%
                  mutate( pred = (as.matrix(predict_df[,model_n]) %*% weight_v) / sum(weight_v) )

# Model outputs
write.csv(mod_w, "Results/VitalRates_4_Cade2015/new_tillers_BH14_mod_sel.csv", row.names=F)
write.csv(pred_df14, "Results/VitalRates_4_Cade2015/new_t_bh14_best.csv", row.names=F)



# Model fits 2015 ---------------------------------------------------------------
res         <- list()
res[[1]]    <-optim(par=c(2,0.1,5),         fn=fit_null, gr=NULL, d15, control=list(maxit=5000))
res[[2]]    <-optim(par=c(5,5,0.5,5),       fn=fit_lam, gr=NULL, d15, control=list(maxit=5000))
res[[3]]    <-optim(par=c(5,0.5,0.5,5),     fn=fit_b, gr=NULL, d15, control=list(maxit=5000))
res[[4]]    <-optim(par=c(5,0.5,0.5,5),     fn=fit_a, gr=NULL, d15, control=list(maxit=5000))
res[[5]]    <-optim(par=c(5,5,0.5,0.5,5),   fn=fit_lam_b, gr=NULL, d15, control=list(maxit=5000))
res[[6]]    <-optim(par=c(5,5,0.5,0.5,5),   fn=fit_lam_a, gr=NULL, d15, control=list(maxit=5000))
res[[7]]    <-optim(par=c(1,0.1,0.1,0.1,5), fn=fit_b_a, gr=NULL, d15, control=list(maxit=5000))
res[[8]]    <-optim(par=c(5,5,0.5,0.5,0.5,5), fn=fit_full, gr=NULL, d15, control=list(maxit=5000))
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

# store parameters, store predictions
par_res     <- list(null = res$null, lam = res$lam, b = res$b, a = res$a, 
                    lam_b = res$lam_b, lam_a = res$lam_a, b_a = res$b_a, full = res$full)
predict_l   <- lapply(par_res, predict_nl, des_mat) 
predict_df  <- Reduce(function(...) cbind(...), predict_l) %>% 
                    as.data.frame %>%
                    setNames( names(par_res) ) %>%
                    bind_cols( des_mat )
id          <- max( which(cumsum(mod_w$weights) < 0.95) ) + 1
weight_v    <- mod_w$weights[1:id]
model_n     <- as.character(mod_w$model[1:id])
pred_df15   <- predict_df %>%
                    mutate( pred = (as.matrix(predict_df[,model_n]) %*% weight_v) / sum(weight_v) )

# Model outputs
write.csv(mod_w, "Results/VitalRates_4_Cade2015/new_tillers_BH15_mod_sel.csv", row.names=F)
write.csv(pred_df15, "Results/VitalRates_4_Cade2015/new_t_bh15_best.csv",row.names=F)



# graph ----------------------------------------------------------------------------------------

# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(gray.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  


# Graph
tiff("Results/VitalRates_4_Cade2015/figure_newtiller_prod.tiff",unit="in",width=6.3,height=4,res=600,compression="lzw")

par(mfrow=c(1,2),mar=c(3,2.5,1,0.1),mgp=c(1.4,0.35,0),cex.lab=1.1,cex.axis=0.8,
    cex.main=1.2, oma=c(0,0,0.2,0),cex=0.9)

#  2014 --------------------------------------------------------
plot(d14$TotDensity,d14$new_t1,pch=21,ylab="Number of new tillers",
     xlab="Planting density", bg = cRamp(d14$sr), cex = 1.5,
     main = "2014")

tmp <- subset(pred_df14, sr == 0.05)
lines(tmp$TotN, tmp$pred)

tmp <- subset(pred_df14, sr == 0.95)
lines(tmp$TotN, tmp$pred)

legend(-1,56,
       c("95% female plot","  5% female plot"),
       lty = c(1,2), lwd = 2, bty="n")


#  2015 --------------------------------------------------------
plot(d15$TotDensity,d15$new_t1,pch=21,ylab="Number of new tillers",
     xlab="Planting density", bg = cRamp(d15$sr), cex = 1.5,
     main = "2015")


tmp <- subset(pred_df15, sr == 0.05)
lines(tmp$TotN, tmp$pred)

tmp <- subset(pred_df15, sr == 0.95)
lines(tmp$TotN, tmp$pred)


colfunc = colorRampPalette(cRamp(unique(arrange(d14,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(19), ncol=1))
text(x=7, y = seq(140,220,l=3), labels = seq(1,0,l=3))
rasterImage(legend_image, 0.5, 140, 5, 220)
text(8, 220, "Percent of", pos = 4)
text(8, 180, "females in", pos = 4)
text(8, 140, "plot", pos = 4)

dev.off()
