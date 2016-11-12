##Sex predictor of flowering probability: SEX RATIO (female individuals/total individuals)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)
library(boot)
library(glmmADMB)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")

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


# YEAR ONE ----------------------------------------------------------------------------------------
tmp                               <- subset(d,year==2014)
tmp[which(tmp$log_l_t1==-Inf),]   <- NA
#prepare the 'new_t1' column
tmp$new_t1[tmp$new_t1=="SKIPPED"] <- NA
tmp$new_t1[tmp$new_t1=="cnf"]     <- NA
tmp$new_t1[tmp$new_t1==""]        <- NA
tmp$new_t1                        <- as.numeric(as.character(tmp$new_t1))
f14                               <- na.omit(tmp[,c("plot","flow_t1","log_l_t1","log_l_t0",
                                                    "sex","sr","new_t1","TotDensity")])

# Model fitting and selection --------------------------------------------------------------------------

fMod=list()
#Target fitness
fMod[[1]]=glmmadmb(flow_t1 ~ log_l_t0 + (1 | plot),data=f14,family="binomial")
fMod[[2]]=glmmadmb(flow_t1 ~ log_l_t0 + sex + (1 | plot),data=f14,family="binomial")
fMod[[3]]=glmmadmb(flow_t1 ~ log_l_t0 * sex + (1 | plot),data=f14,family="binomial")
#Target fitness + effect = tot density 
fMod[[4]]=glmmadmb(flow_t1 ~ log_l_t0 + TotDensity + (1 | plot),data=f14,family="binomial")
fMod[[5]]=glmmadmb(flow_t1 ~ log_l_t0 + sex + TotDensity + (1 | plot),data=f14,family="binomial")
fMod[[6]]=glmmadmb(flow_t1 ~ log_l_t0 * sex + TotDensity + (1 | plot),data=f14,family="binomial")
#Target fitness + effect = tot density + response=sex 
fMod[[7]]=glmmadmb(flow_t1 ~ log_l_t0 + TotDensity + sex:TotDensity + (1 | plot),data=f14,family="binomial")
fMod[[8]]=glmmadmb(flow_t1 ~ log_l_t0 + sex + TotDensity + TotDensity:sex + (1 | plot),data=f14,family="binomial")
fMod[[9]]=glmmadmb(flow_t1 ~ log_l_t0 * sex + TotDensity + TotDensity:sex + (1 | plot),data=f14,family="binomial")
#Target fitness + effect = sex 
fMod[[10]]=glmmadmb(flow_t1 ~ log_l_t0 + sr + TotDensity + (1 | plot),data=f14,family="binomial")
fMod[[11]]=glmmadmb(flow_t1 ~ log_l_t0 + sex + sr + TotDensity + (1 | plot),data=f14,family="binomial")
fMod[[12]]=glmmadmb(flow_t1 ~ log_l_t0 * sex + sr + TotDensity + (1 | plot),data=f14,family="binomial")
#Target fitness + effect = sex + response=sex 
fMod[[13]]=glmmadmb(flow_t1 ~ log_l_t0 + sr + TotDensity + sr:sex + TotDensity:sex +(1 | plot),data=f14,family="binomial")
fMod[[14]]=glmmadmb(flow_t1 ~ log_l_t0 + sex + sr + TotDensity + sr:sex + TotDensity:sex + (1 | plot),data=f14,family="binomial")
fMod[[15]]=glmmadmb(flow_t1 ~ log_l_t0 * sex + sr + TotDensity + sr:sex + TotDensity:sex  + (1 | plot),data=f14,family="binomial")

# Model average
flow_select <- AICtab(fMod,weights=T)
flow_avg    <- model_avg(flow_select, fMod)
write.csv(flow_avg, "Results/VitalRates_3/flowering_best.csv", row.names = F)


# GRAPHS -------------------------------------------------------------------------------------

tiff("Results/VitalRates_3/flowering.tiff",unit="in",width=3.5,height=3.5,res=600,compression="lzw")

par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.5),mgp=c(1.4,0.35,0),cex.lab=0.8,cex.axis=0.8,
    cex.main=0.9)

# Set up colors for plots
f14$col <- as.character(factor(as.integer(f14$sex),labels=c("blue","red")))
mal14   <- subset(f14,sex=="m")
fem14   <- subset(f14,sex=="f")

# Plot data
plot(mal14$flow_t1+0.015 ~ mal14$sr,
     xlab=expression("Sex Ratio"),ylim=c(0,1),
     ylab="Flowering probability",pch=16,col=mal14$col)
par(new=T) ; plot(fem14$flow_t1-0.015 ~ fem14$sr,pch=17,xlab="",ylab="",col=fem14$col,xaxt="n",ylim=c(0,1))

# Compute model predictions
beta  <- flow_avg[,c("predictor","avg")]$avg
low   <- 5
high  <- 42
size  <- mean(f14$log_l_t0)
xSeq  <- seq(0,1,by = 0.1)

y_m_l <- inv.logit( beta[1] + beta[2]*size + beta[3]*low + 
                    beta[4] + beta[5]*xSeq + beta[6]*low +
                    beta[7]*size + beta[8]*xSeq)
y_m_h <- inv.logit( beta[1] + beta[2]*size + beta[3]*high + 
                    beta[4] + beta[5]*xSeq + beta[6]*high +
                    beta[7]*size + beta[8]*xSeq)
y_f_l <- inv.logit( beta[1] + beta[2]*size + beta[3]*low + 
                    beta[5]*xSeq)
y_f_h <- inv.logit( beta[1] + beta[2]*size + beta[3]*high + 
                    beta[5]*xSeq)

lines(xSeq,y_m_h,lty=1,lwd=2,col="red")
lines(xSeq,y_m_l,lty=2,lwd=2,col="red")

lines(xSeq,y_f_h,lty=1,lwd=2,col="blue")
lines(xSeq,y_f_l,lty=2,lwd=2,col="blue")

legend(0,0.5,c("high density","low density","male","female"),
       lty=c(1,2,1,1),lwd=2,col=c("black","black","red","blue"),bty="n")

dev.off()
