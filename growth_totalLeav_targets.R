# Sex competition predictor: SEX RATIO (female individuals/total individuals)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle)    
library(glmmADMB) #Fit models with a Negative Binomial
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")


# load and format data ------------------------------------------------------------------
d=read.csv("Data/vr.csv")
# remove dead individuals (this is a GROWTH model!)
d=subset(d, surv_t1 != 0)
# logtransform leaf numbers
d$plot=as.factor(d$plot) # glmmadmb wants plot as a factor

# Year one
d14=subset(d, year==2014)

# Data frame for analysis
d_tot_14 <- d14 %>% group_by(plot, sr, TotDensity) %>%
  summarise(tot_l = sum(l_t1),
            tot_fI = n())
d_tot_14 <- mutate(d_tot_14, TotDensity2 = TotDensity^2,
                   l_pc = tot_l / tot_fI)

# Year two
tmp15=subset(d,year==2015)
tmp15$new_t1[tmp15$new_t1=="SKIPPED"]=NA
tmp15$new_t1[tmp15$new_t1=="cnf"]=NA
tmp15$new_t1[tmp15$new_t1==""]=NA
tmp15$new_t1=as.numeric(as.character(tmp15$new_t1))
d15=na.omit(tmp15[,c("l_t1","log_l_t0","plot","focalI","sex","new_t1","sr","TotDensity")])

# Data frame for analysis
d_tot_15 <- d15 %>% group_by(plot, sr, TotDensity) %>%
              summarise(tot_l = sum(l_t1),
                        tot_fI = n())
d_tot_15 <- mutate(d_tot_15, TotDensity2 = TotDensity^2,
                             l_pc = tot_l / tot_fI)


# Total per capita number of tillers ------------------------------------------------

# 2014 ---------------------------------------------------------------
nt=list()
# Effect of density
nt[[1]] <- lm(l_pc ~ TotDensity, data=d_tot_14)
nt[[2]] <- lm(l_pc ~ TotDensity + TotDensity2, data=d_tot_14)
#Effect of density + sex ratio
nt[[3]] <- lm(l_pc ~ TotDensity + sr, data=d_tot_14)
nt[[4]] <- lm(l_pc ~ TotDensity + TotDensity2 + sr, data=d_tot_14)
#Interactions b/w density and sex ratio
nt[[5]] <- lm(l_pc ~ sr * TotDensity, data=d_tot_14)
nt[[6]] <- lm(l_pc ~ sr * TotDensity + TotDensity2, data=d_tot_14)
nt[[7]] <- lm(l_pc ~ sr * TotDensity + TotDensity2*sr, data=d_tot_14)

# Model average
grow_tot_14      <- AICtab(nt,weights=T)
grow_tot_14_avg  <- model_avg(grow_tot_14, nt)

# 2015 ---------------------------------------------------------------
nt=list()
# Effect of density
nt[[1]] <- lm(l_pc ~ TotDensity, data=d_tot_15)
nt[[2]] <- lm(l_pc ~ TotDensity + TotDensity2, data=d_tot_15)
#Effect of density + sex ratio
nt[[3]] <- lm(l_pc ~ TotDensity + sr, data=d_tot_15)
nt[[4]] <- lm(l_pc ~ TotDensity + TotDensity2 + sr, data=d_tot_15)
#Interactions b/w density and sex ratio
nt[[5]] <- lm(l_pc ~ sr * TotDensity, data=d_tot_15)
nt[[6]] <- lm(l_pc ~ sr * TotDensity + TotDensity2, data=d_tot_15)
nt[[7]] <- lm(l_pc ~ sr * TotDensity + TotDensity2*sr, data=d_tot_15)

# Model average
grow_tot_15      <- AICtab(nt,weights=T)
grow_tot_15_avg  <- model_avg(grow_tot_15, nt)




# Graph --------------------------------------------------------------------------

# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(topo.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

# Graph
tiff("Results/VitalRates_3/growth_totalLeav_targets.tiff",unit="in",width=6.3,height=3.15,res=600,compression="lzw")

par(mfrow=c(1,2),mar=c(2.5,2.5,1,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0.1))

# total number of tillers
plot(d_tot_14$TotDensity,d_tot_14$l_pc,pch=21,ylab="Number of new leaves in plot (2015)",
     xlab="Planting density", bg = cRamp(d_tot_15$sr), cex = 1.5 )
beta  <-  grow_tot_14_avg[,c("predictor","avg")]$avg
xSeq  <-  seq(0,48,by=1)
yL    <-  beta[1] + beta[2]*0.1 + beta[3]*xSeq*0.1 + 
          beta[4]*xSeq + beta[5]*xSeq^2
yH    <-  beta[1] + beta[2]*0.9 + beta[3]*xSeq*0.9 + 
          beta[4]*xSeq + beta[5]*xSeq^2
lines(xSeq,yL,col="red",lwd=2, lty = 2)
lines(xSeq,yH,col="blue",lwd=2)
title("2014")
legend(10, 60, c("10% female plots", "90% female plots"), cex = 1,
       lty = c(1,2), lwd=2, col=c("red","blue"), bty = "n")

# total number of tillers
plot(d_tot_15$TotDensity,d_tot_15$l_pc,pch=21,ylab="Number of new leaves in plot (2015)",
     xlab="Planting density", bg = cRamp(d_tot_15$sr), cex = 1.5 )
beta  <-  grow_tot_avg[,c("predictor","avg")]$avg
xSeq  <-  seq(0,48,by=1)
yL    <-  beta[1] + beta[2]*0.1 + beta[3]*xSeq*0.1 + 
          beta[4]*xSeq + beta[5]*xSeq^2
yH    <-  beta[1] + beta[2]*0.9 + beta[3]*xSeq*0.9 + 
          beta[4]*xSeq + beta[5]*xSeq^2
lines(xSeq,yL,col="red",lwd=2, lty = 2)
lines(xSeq,yH,col="blue",lwd=2)
title("2015")

dev.off()
