##Response variable: total number of tillers 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")

# load and format data -----------------------------------------------------
x       <- read.csv("Data/vr.csv")
d       <- format_growth(x)
d14     <- subset(d, year == 2014)
d15     <- subset(d, year == 2015)
d15     <- subset(d15, plot != 149) # Remove outlier
# omit NAs
d14 <- na.omit(select(d14, TotDensity, new_t1))
d15 <- na.omit(select(d15, TotDensity, new_t1))

# ML Model selection -------------------------------------------------------

# sex independent model 
fit.BH<-function(params, N){
  yhat<-(params[1]*N$TotDensity)/(1+params[2]*N$TotDensity)
  NLL<- -sum(dnbinom(N$new_t1,mu=yhat,size=params[3],log=T))
  return(NLL)
}
MLE.BH<-optim(par=c(2,0.1,5),fn=fit.BH, gr=NULL, d15, control=list(maxit=5000))

x<-0:max(d15$TotDensity)
plot(d15$TotDensity,d15$new_t1, pch = 16)
lines(x,(MLE.BH$par[1]*x)/(1+MLE.BH$par[2]*x), lwd = 2)

# Only 14 ----------------------------------------------------------------------------------
MLE.BH_14<-optim(par=c(2,0.1,5),fn=fit.BH, gr=NULL, d14, control=list(maxit=5000))
plot(d14$TotDensity,d14$new_t1, pch = 16)
lines(x,(MLE.BH_14$par[1]*x)/(1+MLE.BH_14$par[2]*x), lwd = 2)
betas_14 <- as.data.frame(t(MLE.BH_14$par))
betas_14 <- setNames(betas_14,c("lam","b","size"))
write.csv(betas_14, 
          "C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/Results/VitalRates_3/new_t_best.csv",
          row.names = F)
# ----------------------------------------------------------------------------------

#########################################################################
## option 1 is that there may be male-specific and female-specific parameters
## and population mean is a weighted mean of the two (but this ignores their interactions)
fit.BH.sex<-function(params){
  lambda.F<-params[1]
  lambda.M<-params[2]
  lambda.mean<-lambda.F*d15$sr+lambda.M*(1-d15$sr)
  b.F<-params[3]
  b.M<-params[4]
  b.mean<-b.F*d15$sr+b.M*(1-d15$sr)
  
  yhat<-(lambda.mean*d15$TotDensity)/(1+b.mean*d15$TotDensity)
  NLL<- -sum(dnbinom(d15$new_t1,mu=yhat,size=params[5],log=T))
  return(NLL)
}
MLE.BH.sex<-optim(par=c(10,10,0.5,0.5,5),fn=fit.BH.sex,control=list(maxit=50000))

# sex effect in plambda only
fit.BH.sex.lambda<-function(params){
  lambda.F<-params[1]
  lambda.M<-params[2]
  lambda.mean<-lambda.F*d15$sr+lambda.M*(1-d15$sr)
  
  yhat<-(lambda.mean*d15$TotDensity)/(1+params[3]*d15$TotDensity)
  NLL<- -sum(dnbinom(d15$new_t1,mu=yhat,size=params[4],log=T))
  return(NLL)
}
MLE.BH.sex.lambda<-optim(par=c(10,10,0.5,5),fn=fit.BH.sex.lambda,control=list(maxit=5000))

fit.BH.sex.b<-function(params){
  b.F<-params[1]
  b.M<-params[2]
  b.mean<-b.F*d15$sr+b.M*(1-d15$sr)
  
  yhat<-(params[3]*d15$TotDensity)/(1+b.mean*d15$TotDensity)
  NLL<- -sum(dnbinom(d15$new_t1,mu=yhat,size=params[4],log=T))
  return(NLL)
}
MLE.BH.sex.b<-optim(par=c(10,0.5,0.5,5),fn=fit.BH.sex.b,control=list(maxit=5000))

#####################################################
## compare models

AIC.BH<-2*MLE.BH$value+2*length(MLE.BH$par)
AIC.BH.sex.lambda<-2*MLE.BH.sex.lambda$value+2*length(MLE.BH.sex.lambda$par)
AIC.BH.sex.b<-2*MLE.BH.sex.b$value+2*length(MLE.BH.sex.b$par)
AIC.BH.sex<-2*MLE.BH.sex$value+2*length(MLE.BH.sex$par)


x<-0:max(d15$TotDensity)
plot(d15$TotDensity,d15$new_t1,ylab="New tillers",xlab="Total density", pch = 16)
lines(x,(MLE.BH$par[1]*x)/(1+MLE.BH$par[2]*x),lwd=4)
