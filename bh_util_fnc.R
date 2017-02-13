
# Null model 
fit_null <- function(params, N){
  
  lambda    <- params[1]
  b         <- params[2]
  yhat      <- (lambda*N$TotDensity) / (1 + b*N$TotDensity)
  NLL       <- -sum(dnbinom(N$new_t1,mu=yhat,size=params[3],log=T))
  
  null_par  <<- c("lam.", "b.", "size") 
  return(NLL)
  
}

# lambda
fit_lam <- function(params, N){
  
  lambda.F  <- params[1]
  lambda.M  <- params[2]
  b         <- params[3]
  
  yhat      <- (lambda.F*N$F + lambda.M*N$M) / (1 + b*N$TotDensity)
  NLL       <- -sum(dnbinom(N$new_t1,mu=yhat,size=params[4],log=T))
  
  lam_par   <<- c("lam.f", "lam.m", "b.", "size")
  return(NLL)
  
}

# b
fit_b <- function(params, N){
  
  lambda  <-params[1]
  b.F     <-params[2]
  b.M     <-params[3]
  
  y_f     <- (lambda*N$F) / (1 + b.F*(N$M + N$F))
  y_m     <- (lambda*N$M) / (1 + b.M*(N$M + N$F))
  NLL     <- -sum(dnbinom(N$new_t1,mu=y_f + y_m,size=params[4],log=T))
  
  b_par   <<- c("lam.", "b.f", "b.m", "size")
  return(NLL)
  
}

# alpha
fit_a <- function(params, N){
  
  lambda    <-params[1]
  b         <-params[2]
  a.F       <-params[3]
  a.M       <-params[4]
  
  y_f       <- (lambda*N$F)/ (1 + b*(a.M*N$M + N$F))
  y_m       <- (lambda*N$M)/ (1 + b*(N$M + a.F*N$F))
  NLL       <- -sum(dnbinom(N$new_t1,mu=y_f + y_m,size=params[5],log=T))
  
  a_par     <<- c("lam.", "b.", "a.f", "a.m", "size")
  return(NLL)
  
}

# lambda + b 
fit_lam_b <- function(params, N){
  
  lambda.F  <-params[1]
  lambda.M  <-params[2]
  b.F       <-params[3]
  b.M       <-params[4]
  
  y_f       <- (lambda.F*N$F) / (1 + b.F*(N$M + N$F))
  y_m       <- (lambda.M*N$M) / (1 + b.M*(N$M + N$F))
  NLL       <- -sum(dnbinom(N$new_t1,mu=y_m+y_f,size=params[5],log=T))
  
  lam_b_par <<- c("lam.f", "lam.m", "b.f", "b.m", "size")
  return(NLL)
  
}


# lambda + alpha 
fit_lam_a <- function(params, N){
  
  lambda.F  <-params[1]
  lambda.M  <-params[2]
  b         <-params[3]
  a.F       <-params[4]
  a.M       <-params[5]
  
  y_f       <- (lambda.F*N$F) / (1 + b*(a.M*N$M + N$F))
  y_m       <- (lambda.M*N$M) / (1 + b*(N$M + a.F*N$F))
  NLL       <- -sum(dnbinom(N$new_t1,mu=y_m+y_f,size=params[6],log=T))
  
  lam_a_par <<- c("lam.f", "lam.m", "b.", "a.f", "a.m", "size")
  return(NLL)
  
}


# b + alpha 
fit_b_a <- function(params, N){
  
  lambda    <- params[1]
  b.F       <- params[2]
  b.M       <- params[3]
  a.F       <- params[4]
  a.M       <- params[5]
  
  y_f       <- (lambda * N$F) / (1 + b.F * (a.M*N$M + N$F) )
  y_m       <- (lambda * N$M) / (1 + b.M * (N$M + a.F*N$F) )
  NLL       <- -sum(dnbinom(N$new_t1,mu=y_m+y_f,size=params[6],log=T))
  
  b_a_par   <<- c("lam.", "b.f", "b.m", "a.f", "a.m", "size")
  return(NLL)
  
}

# lambda + b + alpha
fit_full <- function(params, N){
  
  lambda.F  <-params[1]
  lambda.M  <-params[2]
  b.F       <-params[3]
  b.M       <-params[4]
  a.F       <-params[5]
  a.M       <-params[6]   
  
  y_f       <- (lambda.F*N$F)/ (1 + b.F*(a.M*N$M + N$F))
  y_m       <- (lambda.M*N$M)/ (1 + b.M*(N$M + a.F*N$F))
  NLL       <- -sum(dnbinom(N$new_t1,mu=(y_m + y_f),size=params[7],log=T))
  
  full_par  <<- c("lam.f", "lam.m", "b.f", "b.m", "a.f", "a.m", "size")
  return(NLL)
  
}


# model selection ------------------------------------------------------------
par_names <- function(results, par_names){
  
  results$par <- setNames(results$par, par_names)
  return(results)
  
}

# calculate aic value from skratch
aic_calc <- function(x){
  
  out <- 2*x$value + 2*length(x$par)
  return(out)
  
}

# model weights
mod_weights <- function(x){
  out <- data.frame(model=names(x),aic=unlist(x),row.names=NULL)
  out <- out[order(out$aic),]
  out <- out %>% mutate(deltaAIC= aic - aic[1])
  out <- out %>% mutate(relLik  = exp(-0.5* deltaAIC) )
  out <- out %>% mutate(weights  = round(relLik / sum(relLik),3) )
  return(out)
}
