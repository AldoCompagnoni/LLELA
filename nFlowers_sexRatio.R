##Sex predictor of flowering probability: SEX RATIO (female individuals/total individuals)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) #For AIC weights, because I'm lazy!
library(lme4) #; library(glmmADMB)


#read in data
d <- read.csv("Data/vr.csv")

#Remove resuscitated individuals (dead in spring 2014, alive in 2015)
#NOTE: possibly these are NOT DEAD in 2014, because they all have
#few leaves: resprouted from base?
d <- d[-which(d$plot == 18 & d$focalI =="m2"),]
d <- d[-which(d$plot == 38 & d$focalI =="f1"),]
d <- d[-which(d$plot == 46 & d$focalI =="f1"),]
d <- d[-which(d$plot == 83 & d$focalI =="f5"),]
d <- d[-which(d$plot == 36 & d$focalI =="m3"),]

#logtransform leaf numbers
d$log_l_t0 <- log(d$l_t0)
d$log_l_t1 <- log(d$l_t1)
#make "plot" a factor to fit models 
d$plot <- as.factor(d$plot) 

#Transform densities to SEX RATIO
d$sr <- d$F / d$TotDensity


# Select YEAR ONE data ------------------------------------------------------------------------------------
tmp                               <- subset(d,year==2014)
tmp[which(tmp$log_l_t1==-Inf),]   <- NA
#prepare the 'new_t1' column
tmp$new_t1[tmp$new_t1=="SKIPPED"] <- NA
tmp$new_t1[tmp$new_t1=="cnf"]     <- NA
tmp$new_t1[tmp$new_t1==""]        <- NA
tmp$new_t1                        <- as.numeric(as.character(tmp$new_t1))
f14                               <- na.omit(tmp[,c("plot","flow_t1","log_l_t1","log_l_t0","flowN_t1",
                                                    "sex","sr","new_t1","TotDensity")])

# Model selection ----------------------------------------------------------------------------------

nfMod=list()
#Target fitness
nfMod[[1]]=glmmadmb(flowN_t1 ~ log_l_t0 + (1 | plot),data=f14,family="poisson")
nfMod[[2]]=glmmadmb(flowN_t1 ~ log_l_t0 + sex + (1 | plot),data=f14,family="poisson")
nfMod[[3]]=glmmadmb(flowN_t1 ~ log_l_t0 * sex + (1 | plot),data=f14,family="poisson")
#Target fitness + effect = tot density 
nfMod[[4]]=glmmadmb(flowN_t1 ~ log_l_t0 + TotDensity + (1 | plot),data=f14,family="poisson")
nfMod[[5]]=glmmadmb(flowN_t1 ~ log_l_t0 + sex + TotDensity + (1 | plot),data=f14,family="poisson")
nfMod[[6]]=glmmadmb(flowN_t1 ~ log_l_t0 * sex + TotDensity + (1 | plot),data=f14,family="poisson")
#Target fitness + effect = tot density + response=sex 
nfMod[[7]]=glmmadmb(flowN_t1 ~ log_l_t0 + TotDensity + TotDensity:sex + (1 | plot),data=f14,family="poisson")
nfMod[[8]]=glmmadmb(flowN_t1 ~ log_l_t0 + sex*TotDensity + (1 | plot),data=f14,family="poisson")
nfMod[[9]]=glmmadmb(flowN_t1 ~ log_l_t0 * sex + TotDensity + TotDensity:sex + (1 | plot),data=f14,family="poisson")
#Target fitness + effect = sex 
nfMod[[10]]=glmmadmb(flowN_t1 ~ log_l_t0 + sr + TotDensity + (1 | plot),data=f14,family="poisson")
nfMod[[11]]=glmmadmb(flowN_t1 ~ log_l_t0 + sex + sr + TotDensity + (1 | plot),data=f14,family="poisson")
nfMod[[12]]=glmmadmb(flowN_t1 ~ log_l_t0 * sex + sr + TotDensity + (1 | plot),data=f14,family="poisson")
#Target fitness + effect = sex + response=sex 
nfMod[[13]]=glmmadmb(flowN_t1 ~ log_l_t0 + sr + TotDensity + sr:sex + TotDensity:sex +(1 | plot),data=f14,family="poisson")
nfMod[[14]]=glmmadmb(flowN_t1 ~ log_l_t0 + sex + sr + TotDensity + sr:sex + TotDensity:sex + (1 | plot),data=f14,family="poisson")
nfMod[[15]]=glmmadmb(flowN_t1 ~ log_l_t0 * sex + sr + TotDensity + sr:sex + TotDensity:sex  + (1 | plot),data=f14,family="poisson")
mod_select = AICtab(nfMod,weights=T)


# Model average (Models that make up more than 95% of weight)
betaList      <- list()
mod_rank <- do.call(rbind,strsplit(attributes(mod_select)$row.names , "model"))
mod_rank <- as.numeric(mod_rank[,2])

# First 11 models make up 95% of weight (sum(mod_select$weight[1:11]))
for(i in 1:11 ){

  coefficients = coef(nfMod[[mod_rank[i]]])
  betaList[[i]] <-  data.frame(predictor = names(coefficients),
                               parameter = coefficients)
  names(betaList[[i]])[2] <- paste0("parameter_",i)
  
}

# Model averages
beta_avg      <- Reduce(function(...) merge(...,all=T), betaList)
beta_avg[is.na(beta_avg)] <- 0
weights       <- mod_select$weight[1:11]
beta_avg$avg  <- as.matrix(beta_avg[,-1]) %*% weights / sum(mod_select$weight[1:11]) 

write.csv(beta_avg[,c("predictor","avg")], 
          "Results/VitalRates_3/nFlowers_best.csv", row.names = F)


# GRAPHS ---------------------------------------------------------------------------------------------------------------

beta = beta_avg[,c("predictor","avg")]$avg

plot( f14$sr, f14$flowN_t1 , pch = 16,
      ylab = "Number of flowers")

low   <- 5
high  <- 42
size  <- mean(f14$log_l_t0)
xSeq  <- seq(0,1,by = 0.1)

y_m_l <- exp( beta[1] + beta[2]*size + beta[3]*low + 
              beta[4]*low + beta[5]*xSeq + beta[6] +
              beta[7]*low + beta[8]*size + beta[9]*xSeq)
y_m_h <- exp( beta[1] + beta[2]*size + beta[3]*high + 
              beta[4]*high + beta[5]*xSeq + beta[6] +
              beta[7]*high + beta[8]*size + beta[9]*xSeq)

y_f_l <- exp( beta[1] + beta[2]*size + beta[3]*low + 
                beta[5]*xSeq)
y_f_h <- exp( beta[1] + beta[2]*size + beta[3]*high + 
                beta[5]*xSeq)

lines(xSeq,y_m_l,lty=2,lwd=2,col="red")
lines(xSeq,y_m_h,lty=1,lwd=2,col="red")

lines(xSeq,y_f_l,lty=2,lwd=2,col="blue")
lines(xSeq,y_f_h,lty=1,lwd=2,col="blue")
