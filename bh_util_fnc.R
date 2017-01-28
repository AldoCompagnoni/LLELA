# ML Candidate models -------------------------------------------------------

# sex independent model 
fit.BH<-function(params, N){
  yhat<-(params[1]*N$TotDensity)/(1+params[2]*N$TotDensity)
  NLL<- -sum(dnbinom(N$new_t1,mu=yhat,size=params[3],log=T))
  return(NLL)
}

## option 1 is that there may be male-specific and female-specific parameters
## and population mean is a weighted mean of the two (but this ignores their interactions)
fit.BH.sex<-function(params, N){
  lambda.F<-params[1]
  lambda.M<-params[2]
  lambda.mean<-lambda.F*N$sr+lambda.M*(1-N$sr)
  b.F<-params[3]
  b.M<-params[4]
  b.mean<-b.F*N$sr+b.M*(1-N$sr)
  
  yhat<-(lambda.mean*N$TotDensity)/(1+b.mean*N$TotDensity)
  NLL<- -sum(dnbinom(N$new_t1,mu=yhat,size=params[5],log=T))
  return(NLL)
}

# sex effect in plambda only
fit.BH.sex.lambda<-function(params, N){
  lambda.F<-params[1]
  lambda.M<-params[2]
  lambda.mean<-lambda.F*N$sr+lambda.M*(1-N$sr)
  
  yhat<-(lambda.mean*N$TotDensity)/(1+params[3]*N$TotDensity)
  NLL<- -sum(dnbinom(N$new_t1,mu=yhat,size=params[4],log=T))
  return(NLL)
}

# sex effect in competition coefficients only
fit.BH.sex.b<-function(params, N){
  lambda  <-params[1]
  b.F     <-params[2]
  b.M     <-params[3]
  b.mean<-b.F*N$sr+b.M*(1-N$sr)
  
  yhat<-(lambda*N$TotDensity)/(1+b.mean*N$TotDensity)
  NLL<- -sum(dnbinom(N$new_t1,mu=yhat,size=params[4],log=T))
  return(NLL)
}

# model selection ------------------------------------------------------------


# model weights
mod_weights <- function(x){
  out <- data.frame(model=names(x),aic=unlist(x),row.names=NULL)
  out <- out[order(out$aic),]
  out <- out %>% mutate(deltaAIC= aic - aic[1])
  out <- out %>% mutate(relLik  = exp(-0.5* deltaAIC) )
  out <- out %>% mutate(weights  = round(relLik / sum(relLik),3) )
  return(out)
}
