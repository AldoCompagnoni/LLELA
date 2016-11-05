##Response variable: total number of tillers 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) 
library(glmmADMB) # Fit models with a Negative Binomial
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")

#read in data
d=read.csv("Data/vr.csv")
#remove dead individuals (this is a GROWTH model!)
d=subset(d, surv_t1 != 0)
#logtransform leaf numbers
d$plot=as.factor(d$plot) #glmmadmb wants plot as a factor


# Year two
tmp15=subset(d,year==2015)
#remove 
tmp15$new_t1[tmp15$new_t1=="SKIPPED"]=NA
tmp15$new_t1[tmp15$new_t1=="cnf"]=NA
tmp15$new_t1[tmp15$new_t1==""]=NA
tmp15$new_t1=as.numeric(as.character(tmp15$new_t1))
tmp15$TotDensity2 = tmp15$TotDensity^2
d15=na.omit(tmp15[,c("l_t1","log_l_t0","plot","sex","new_t1",
                     "TotDensity","TotDensity2","sr")])
d15$new_t1_pc <- d15$new_t1 / d15$TotDensity
d15 = subset(d15, plot != 149) # Remove outlier

#Response variable: total number of tillers per capita -------------------------------------------

#INITIAL STRUCTURES
#Effect of plot treatment on new tillers#################################
nt=list()
# Effect of density
nt[[1]] <- lm(new_t1_pc ~ TotDensity, data=d15)
nt[[2]] <- lm(new_t1_pc ~ TotDensity + TotDensity2, data=d15)
#Effect of density + sex ratio
nt[[3]] <- lm(new_t1_pc ~ TotDensity + sr, data=d15)
nt[[4]] <- lm(new_t1_pc ~ TotDensity + TotDensity2 + sr, data=d15)
#Interactions b/w density and sex ratio
nt[[5]] <- lm(new_t1_pc ~ sr * TotDensity, data=d15)
nt[[6]] <- lm(new_t1_pc ~ sr * TotDensity + TotDensity2, data=d15)
nt[[7]] <- lm(new_t1_pc ~ sr * TotDensity + TotDensity2*sr, data=d15)
AICtab(nt,weights=T)

# Model average
grow_new_tpc      <- AICtab(nt,weights=T)
grow_new_tpc_avg  <- model_avg(grow_new_tpc, nt)
write.csv(grow_new_tpc_avg, "Results/VitalRates_3/growth_new_tpc_best.csv", row.names = F)


#Graph: total number of tillers -----------------------------------------------------
tiff("Results/VitalRates_3/growth_newTillers_pc.tiff",unit="in",width=6.3,height=3.15,res=600,compression="lzw")

par(mfcol=c(1,2),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0.1))

plot(d15$TotDensity,d15$new_t1_pc,pch=16,ylab="Per capita new tillers (2015)",
     xlab="Planting density")

#Calculate lines
beta  <-  grow_new_tpc_avg[,c("predictor","avg")]$avg
xSeq  <-  seq(0,48,by=1)
#yL    <-  beta[1] + beta[2]*xSeq + beta[3]*xSeq^2 + beta[4]*0.1 + 
#          beta[5]*xSeq*0.1 + beta[6]*xSeq^2*0.1
#yH    <-  beta[1] + beta[2]*xSeq + beta[3]*xSeq^2 + beta[4]*0.9 + 
#          beta[5]*xSeq*0.9 + beta[6]*xSeq^2*0.9
yL    <-  beta[1] + beta[2]*0.1 + beta[3]*xSeq*0.1 + beta[4]*xSeq + 
          beta[5]*xSeq^2 + beta[6]*xSeq^2*0.1
yH    <-  beta[1] + beta[2]*0.9 + beta[3]*xSeq*0.9 + beta[4]*xSeq + 
          beta[5]*xSeq^2 + beta[6]*xSeq^2*0.9

lines(xSeq,yL,col="red",lwd=2, lty = 2)
lines(xSeq,yH,col="blue",lwd=2)

legend(10, 19, c("10% female plots", "90% female plots"),
       lty = c(1,2), lwd=2, col=c("red","blue"), bty = "n")

#dev.off()


#Response variable: total number of tillers ------------------------------------------------

#INITIAL STRUCTURES
#Effect of plot treatment on new tillers#################################
nt=list()
# Effect of density
nt[[1]] <- lm(new_t1 ~ TotDensity, data=d15)
nt[[2]] <- lm(new_t1 ~ TotDensity + TotDensity2, data=d15)
#Effect of density + sex ratio
nt[[3]] <- lm(new_t1 ~ TotDensity + sr, data=d15)
nt[[4]] <- lm(new_t1 ~ TotDensity + TotDensity2 + sr, data=d15)
#Interactions b/w density and sex ratio
nt[[5]] <- lm(new_t1 ~ sr * TotDensity, data=d15)
nt[[6]] <- lm(new_t1 ~ sr * TotDensity + TotDensity2, data=d15)
nt[[7]] <- lm(new_t1 ~ sr * TotDensity + TotDensity2*sr, data=d15)

# Model average
grow_new_t      <- AICtab(nt,weights=T)
grow_new_t_avg  <- model_avg(grow_new_t, nt)
write.csv(grow_new_t_avg, "Results/VitalRates_3/growth_new_t_best.csv", row.names = F)


#Graph: total number of tillers -----------------------------------------------------
plot(d15$TotDensity,d15$new_t1,pch=16,ylab="Number of new tillers (2015)",
     xlab="Planting density")

#Calculate lines
beta  <-  grow_new_t_avg[,c("predictor","avg")]$avg
xSeq  <-  seq(0,48,by=1)
yL    <-  beta[1] + beta[2]*0.1 + beta[3]*xSeq*0.1 + 
          beta[4]*xSeq^2*0.1 + beta[5]*xSeq + beta[6]*xSeq^2
yH    <-  beta[1] + beta[2]*0.9 + beta[3]*xSeq*0.1 + 
          beta[4]*xSeq^2*0.9 + beta[5]*xSeq + beta[6]*xSeq^2
  
lines(xSeq,yL,col="red",lwd=2, lty = 2)
lines(xSeq,yH,col="blue",lwd=2)

dev.off()

