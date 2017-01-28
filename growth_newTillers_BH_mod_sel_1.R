##Response variable: total number of tillers 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/unload_mass.R")
unload_mass()
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")

# load and format data -----------------------------------------------------
x       <- read.csv("Data/vr.csv")
d       <- format_growth(x)
d14     <- subset(d, year == 2014)
d15     <- subset(d, year == 2015)
d15     <- subset(d15, plot != 149) # Remove outlier
# omit NAs
d14 <- na.omit(select(d14, TotDensity, sr, new_t1))
d15 <- na.omit(select(d15, TotDensity, sr, new_t1))

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
  b.F<-params[1]
  b.M<-params[2]
  b.mean<-b.F*N$sr+b.M*(1-N$sr)
  
  yhat<-(params[3]*N$TotDensity)/(1+b.mean*N$TotDensity)
  NLL<- -sum(dnbinom(N$new_t1,mu=yhat,size=params[4],log=T))
  return(NLL)
}


# Model fits --------------------------------------------------------------------------------------------

# 2014
BH.14<-optim(par=c(2,0.1,5),fn=fit.BH, gr=NULL, d15, control=list(maxit=5000))
BH.sex.14<-optim(par=c(10,10,0.5,0.5,5),fn=fit.BH.sex,gr=NULL, d14,control=list(maxit=5000))
BH.sex.lam.14<-optim(par=c(10,10,0.5,5),fn=fit.BH.sex.lambda,gr=NULL, d14, control=list(maxit=5000))
BH.sex.b.14<-optim(par=c(10,10,0.5,5),fn=fit.BH.sex.lambda,gr=NULL, d14, control=list(maxit=5000))

# 2015
BH.15<-optim(par=c(2,0.1,5),fn=fit.BH, gr=NULL, d15, control=list(maxit=5000))
BH.sex.15<-optim(par=c(10,10,0.5,0.5,5),fn=fit.BH.sex,gr=NULL, d14,control=list(maxit=5000))
BH.sex.lam.15<-optim(par=c(10,10,0.5,5),fn=fit.BH.sex.lambda,gr=NULL, d14, control=list(maxit=5000))
BH.sex.b.15<-optim(par=c(10,10,0.5,5),fn=fit.BH.sex.lambda,gr=NULL, d14, control=list(maxit=5000))


# Compare models --------------------------------------------------------------------------------------

# 2014
AIC.BH.14             <-2*BH.14$value+2*length(BH.14$par)
AIC.BH.sex.lam.14     <-2*BH.sex.lam.14$value+2*length(BH.sex.lam.14$par)
AIC.BH.sex.b.14       <-2*BH.sex.b.14$value+2*length(BH.sex.b.14$par)
AIC.BH.sex.14         <-2*BH.sex.14$value+2*length(BH.sex.14$par)

# 2015
AIC.BH.15             <-2*BH.15$value+2*length(BH.15$par)
AIC.BH.sex.lam.15     <-2*BH.sex.lam.15$value+2*length(BH.sex.lam.15$par)
AIC.BH.sex.b.15       <-2*BH.sex.b.15$value+2*length(BH.sex.b.15$par)
AIC.BH.sex.15         <-2*BH.sex.15$value+2*length(BH.sex.15$par)



# Model selection ----------------------------------------------------------------------------------------


x<-0:max(d15$TotDensity)
plot(d15$TotDensity,d15$new_t1, pch = 16)
lines(x,(MLE.BH$par[1]*x)/(1+MLE.BH$par[2]*x), lwd = 2)

# Only 14 ----------------------------------------------------------------------------------
MLE.BH_14<-optim(par=c(2,0.1,5),fn=fit.BH, gr=NULL, d14, control=list(maxit=5000))
plot(d14$TotDensity,d14$new_t1, pch = 16)
lines(x,(MLE.BH_14$par[1]*x)/(1+MLE.BH_14$par[2]*x), lwd = 2)
betas_14 <- as.data.frame(t(MLE.BH_14$par))
betas_14 <- setNames(betas_14,c("lam","b","size"))
#write.csv(betas_14, 
#          "C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/Results/VitalRates_3/new_t_best.csv",
#          row.names = F)
# ----------------------------------------------------------------------------------





x<-0:max(d15$TotDensity)
plot(d15$TotDensity,d15$new_t1,ylab="New tillers",xlab="Total density", pch = 16)
lines(x,(MLE.BH$par[1]*x)/(1+MLE.BH$par[2]*x),lwd=4)
