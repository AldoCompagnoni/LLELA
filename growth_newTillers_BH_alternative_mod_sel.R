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
d14 <- na.omit(select(d14, TotDensity, F, M, sr, new_t1))
d15 <- na.omit(select(d15, TotDensity, sr, new_t1))


# sex differences in lambda and 
fit_null <- function(params, N){
  
  lambda  <- params[1]
  b       <- params[2]
  yhat    <- (lambda*N$TotDensity) / (1 + b*N$TotDensity)
  NLL     <- -sum(dnbinom(N$new_t1,mu=yhat,size=params[3],log=T))
  return(NLL)
  
}

# sex diffs in lambda
fit_lam <- function(params, N){
  
  lambda.F  <-params[1]
  lambda.M  <-params[2]
  b         <-params[3]
  
  yhat      <- (lambda.F*N$F + lambda.M*N$M)/(1+ b*N$TotDensity)
  NLL       <- -sum(dnbinom(N$new_t1,mu=yhat,size=params[4],log=T))
  return(NLL)
  
}

# sex differences in lambda and 
fit_c_resp <- function(params, N){
  
  lambda.F<-params[1]
  lambda.M<-params[2]
  b.F     <-params[3]
  b.M     <-params[4]
  
  y_f     <- (lambda.F*N$F)/ (1 + b.F*(N$M + N$F))
  y_m     <- (lambda.M*N$M)/ (1 + b.M*(N$M + N$F))
  
  NLL     <- -sum(dnbinom(N$new_t1,mu=y_m+y_f,size=params[5],log=T))
  return(NLL)
  
}


# full (effect + response)
fit_full <- function(params, N){
  
  lambda.F<-params[1]
  lambda.M<-params[2]
  b.F     <-params[3]
  b.M     <-params[4]
  a.F     <-params[5]
  a.M     <-params[6]   
  
  y_f     <- (lambda.F*N$F)/ (1 + b.F*(a.M*N$M + N$F))
  y_m     <- (lambda.M*N$M)/ (1 + b.M*(N$M + a.F*N$F))
    
  NLL     <- -sum(dnbinom(N$new_t1,mu=(y_m + y_f),size=params[7],log=T))
  return(NLL)
  
}


# model fits
null    <-optim(par=setNames(c(2,0.1,5),c("lam","b","size")),
                fn=fit_null, gr=NULL, d14, control=list(maxit=50000))
lam     <-optim(par=setNames(c(10,10,0.5,5),c("lam_f","lam_m","b","size")),
                fn=fit_lam,  gr=NULL, d14, control=list(maxit=50000))
c_resp  <-optim(par=setNames(c(10,10,0.5,0.5,5),c("lam_f","lam_m","b_f","b_m","size")),
                fn=fit_c_resp,gr=NULL, d14,control=list(maxit=50000))
full    <-optim(par=setNames(c(10,10,0.5,0.5,0.5,0.5,5),c("lam_f","lam_m","b_f","b_m","a_f","a_m","size")),
                fn=fit_full, gr=NULL, d14, control=list(maxit=50000))

# AIC weights
m14 <- list()
m14[[1]]  <-2*null$value+2*length(null$par)
m14[[2]]  <-2*lam$value+2*length(lam$par)
m14[[3]]  <-2*c_resp$value+2*length(c_resp$par)
m14[[4]]  <-2*full$value+2*length(full$par)
m14       <- setNames(m14, c("null", "lam", "c_resp","full"))

mod_weights(m14)


# Model averages -----------------------------------------------------------------
sr_seq    <- seq(0,1,0.05)
des       <- expand.grid(TotDensity = seq(1,48,1), sr = sr_seq)
des       <- mutate(des, F = TotDensity*sr,
                         M = TotDensity*(1-sr)
                    )

beta   <- lam$par
beta["lam_f"]* 

des     <- mutate(des, lams = (beta["lam_f"]*F + beta["lam_m"]*M) / ( 1 + beta["b"]*TotDensity )  )
  
# model 1
beta1   <- null$par

des     <- mutate(des, lam_1 = beta1["lam"] * TotDensity,
                       comp_1 = (1 + beta1["b"]*TotDensity),
                       lam_2 = beta2["lam_m"] * M + beta2["lam_f"] * F,
                  )



des       <- merge(des, lambda_df)
des       <- mutate(des, yhat = (lambda*TotDensity) / (1+beta["b"]*TotDensity))
des       <- mutate(des, yhat_pc = yhat / TotDensity)


