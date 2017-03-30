##Response variable: total number of tillers 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/unload_mass.R")
unload_mass()
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/ac79/Documents/CODE/LLELA/bh_util_fnc.R")

# load and format data -----------------------------------------------------
x   <- read.csv("Data/vr.csv")
d   <- format_new_tillers(x)
d14 <- subset(d, year == 2014)
d15 <- subset(d, year == 2015)


# Model fits 2014 ---------------------------------------------------------------
res         <- list()
res[[1]]    <-optim(par=c(2,0.1,5),         fn=fit_null, gr=NULL, d14, control=list(maxit=5000))
res[[2]]    <-optim(par=c(5,5,0.5,5),       fn=fit_lam, gr=NULL, d14, control=list(maxit=5000))
res[[3]]    <-optim(par=c(5,0.5,0.5,5),     fn=fit_b, gr=NULL, d14, control=list(maxit=5000))
res[[4]]    <-optim(par=c(5,0.5,0.5,5), fn=fit_a, gr=NULL, d14, control=list(maxit=5000))
res[[5]]    <-optim(par=c(5,5,0.5,0.5,5),       fn=fit_lam_b, gr=NULL, d14, control=list(maxit=5000))
res[[6]]    <-optim(par=c(5,5,0.5,0.5,5),   fn=fit_lam_a, gr=NULL, d14, control=list(maxit=5000))
res[[7]]    <-optim(par=c(1,0.1,0.1,0.1,5), fn=fit_b_a, gr=NULL, d14, control=list(maxit=5000))
res[[8]]    <-optim(par=c(5,5,0.5,0.5,0.5,5), fn=fit_full, gr=NULL, d14, control=list(maxit=5000))
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

# Model selection table
write.csv(mod_w,"Results/VitalRates_3/new_tillers_BH14_mod_sel.csv",row.names=F)

# Average AIC weight models --------------------------------------------------------------------------------------
mat <- matrix(data = 1, nrow = length(res), ncol = 5) 
for(i in 1:length(res)){
  
  params  <- c("lam.f","lam.m","b.f","b.m","a.")
  
  lam_r   <- grep("lam.",names(res[[i]]$par) )
  b_r     <- grep("b.",names(res[[i]]$par) )
  a_r     <- grep("^a.$",names(res[[i]]$par) )
  
  if( length(lam_r) > 0 ) mat[i,grep("lam.",params)] <- res[[i]]$par[lam_r]
  if( length(b_r) > 0 )   mat[i,grep("b.",params)]   <- res[[i]]$par[b_r]
  if( length(a_r) > 0 )   mat[i,grep("^a.$",params)]   <- res[[i]]$par[a_r]
  
}
results   <- mat %>%
              data.frame() %>%
              setNames(params) %>%
              mutate(model = names(res) ) %>%
              merge(mod_w[,c("model","weights")]) %>%
              arrange(desc(weights) ) %>%
              mutate(cum_weight = cumsum(weights))
# only the the models making up to 95% (included) weights
min_cum_weight <- min(results$cum_weight[results$cum_weight > 0.95])
weight_i  <- which(results$cum_weight == min_cum_weight)
weighted  <- (select(results,lam.f:a.) * results$weights)[1:weight_i,]
summed    <- apply(weighted,2,sum) 
avg       <- summed / results$cum_weight[weight_i]
avg14     <- as.data.frame(t(avg))

write.csv(avg14, "Results/VitalRates_3/new_t_bh14_best.csv",row.names=F)


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

# Model selection table
write.csv(mod_w,"Results/VitalRates_3/new_tillers_BH15_mod_sel.csv",row.names=F)

# Average AIC weight models --------------------------------------------------------------------------------------
mat <- matrix(data = 1, nrow = length(res), ncol = 5) 
for(i in 1:length(res)){
  
  #params  <- c("lam.f","lam.m","b.f","b.m","a.f","a.m")
  params  <- c("lam.f","lam.m","b.f","b.m","a.")
  
  lam_r   <- grep("lam.",names(res[[i]]$par) )
  b_r     <- grep("b.",names(res[[i]]$par) )
  a_r     <- grep("^a.$",names(res[[i]]$par) )
  
  if( length(lam_r) > 0 ) mat[i,grep("lam.",params)] <- res[[i]]$par[lam_r]
  if( length(b_r) > 0 )   mat[i,grep("b.",params)]   <- res[[i]]$par[b_r]
  if( length(a_r) > 0 )   mat[i,grep("^a.$",params)]   <- res[[i]]$par[a_r]
  
}
results   <- mat %>%
              data.frame() %>%
              setNames(params) %>%
              mutate(model = names(res) ) %>%
              merge(mod_w[,c("model","weights")]) %>%
              arrange(desc(weights) ) %>%
              mutate(cum_weight = cumsum(weights))
# only the the models making up to 95% (included) weights
min_cum_weight <- min(results$cum_weight[results$cum_weight > 0.95])
weight_i  <- which(results$cum_weight == min_cum_weight)
weighted  <- (select(results,lam.f:a.) * results$weights)[1:weight_i,]
summed    <- apply(weighted,2,sum) 
avg       <- summed / results$cum_weight[weight_i]
avg15     <- as.data.frame(t(avg))

write.csv(avg15, "Results/VitalRates_3/new_t_bh15_best.csv",row.names=F)


# graph ----------------------------------------------------------------------------------------

# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(gray.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  


# Graph
tiff("Results/VitalRates_3/figure_newtiller_prod.tiff",unit="in",width=6.3,height=4,res=600,compression="lzw")

par(mfrow=c(1,2),mar=c(3,2.5,1,0.1),mgp=c(1.4,0.35,0),cex.lab=1.1,cex.axis=0.8,
    cex.main=1.2, oma=c(0,0,0.2,0),cex=0.9)

#  2014 --------------------------------------------------------
plot(d14$TotDensity,d14$new_t1,pch=21,ylab="Number of new tillers",
     xlab="Planting density", bg = cRamp(d14$sr), cex = 1.5,
     main = "2014")
N    <- seq(0,48,1)
beta <- as.data.frame(avg14)
fem  <- N*0.95
mal  <- N*0.05
y_h  <- ((beta[,"lam.f"]*fem) / (1 + beta[,"b.f"] * (fem + beta[,"a."]*mal))) + 
        ((beta[,"lam.m"]*mal) / (1 + beta[,"b.m"] * (fem + beta[,"a."]*mal)))
fem  <- N*0.05
mal  <- N*0.95
y_l  <- ((beta[,"lam.f"]*fem) / (1 + beta[,"b.f"] * (fem + beta[,"a."]*mal))) + 
        ((beta[,"lam.m"]*mal) / (1 + beta[,"b.m"] * (fem + beta[,"a."]*mal)))
lines(N,y_h,lwd=2, lty = 1) #col="#DCDCDC",
lines(N,y_l,lwd=2, lty = 2) #col="#636363",


legend(-1,56,
       c("95% female plot","  5% female plot"),
       lty = c(1,2), lwd = 2, bty="n")


#  2015 --------------------------------------------------------
plot(d15$TotDensity,d15$new_t1,pch=21,ylab="Number of new tillers",
     xlab="Planting density", bg = cRamp(d15$sr), cex = 1.5,
     main = "2015")
N    <- seq(0,48,1)
beta <- as.data.frame(avg15)
fem  <- N*0.95
mal  <- N*0.05
y_h  <- (beta[,"lam.f"]*fem) / (1 + beta[,"b.f"] * (fem + beta[,"a."]*mal)) + 
        (beta[,"lam.m"]*mal) / (1 + beta[,"b.m"] * (fem + beta[,"a."]*mal))
fem  <- N*0.05
mal  <- N*0.95
y_l  <- (beta[,"lam.f"]*fem) / (1 + beta[,"b.f"] * (fem + beta[,"a."]*mal)) + 
  (beta[,"lam.m"]*mal) / (1 + beta[,"b.m"] * (fem + beta[,"a."]*mal))
lines(N,y_h,lwd=2, lty = 1) #col="#DCDCDC",
lines(N,y_l,lwd=2, lty = 2) #col="#636363",

colfunc = colorRampPalette(cRamp(unique(arrange(d14,sr)$sr)))
legend_image <- as.raster(matrix(colfunc(19), ncol=1))
text(x=7, y = seq(140,220,l=3), labels = seq(1,0,l=3))
rasterImage(legend_image, 0.5, 140, 5, 220)
text(8, 220, "Percent of", pos = 4)
text(8, 180, "females in", pos = 4)
text(8, 140, "plot", pos = 4)

dev.off()
