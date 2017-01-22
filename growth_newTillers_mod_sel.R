##Response variable: total number of tillers 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) 
library(glmmADMB) # Fit models with a Negative Binomial
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")


# load and format data -----------------------------------------------------
x       <- read.csv("Data/vr.csv")
d       <- format_growth(x)
d14     <- subset(d, year == 2014)
d15     <- subset(d, year == 2015)
d15     <- subset(d15, plot != 149) # Remove outlier

# One sex plots 
one_sex <- one_sex_format(d)

# Model selection ##############################################################
d15     <- unique(select(d15, new_t1, new_t1_pc, plot, sr, TotDensity, TotDensity2, year))
d15     <- na.omit(d15)

# 2015 --------------------------------------------------------------------------------
# Tillers per capita------------------------------------------------------------
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
grow_new_tpc15      <- AICtab(nt,weights=T)
grow_new_tpc15_avg  <- model_avg(grow_new_tpc15, nt)
write.csv(grow_new_tpc15_avg, "Results/VitalRates_3/growth_new_tpc_15.csv", row.names = F)


# Total number of tillers ------------------------------------------------
nt=list()
# Effect of density
nt[[1]] <- glmmadmb(new_t1 ~ TotDensity, data=d15, family="nbinom2")
nt[[2]] <- glmmadmb(new_t1 ~ TotDensity + TotDensity2, data=d15, family="nbinom2")
#Effect of density + sex ratio
nt[[3]] <- glmmadmb(new_t1 ~ TotDensity + sr, data=d15, family="nbinom2")
nt[[4]] <- glmmadmb(new_t1 ~ TotDensity + TotDensity2 + sr, data=d15, family="nbinom2")
#Interactions b/w density and sex ratio
nt[[5]] <- glmmadmb(new_t1 ~ sr * TotDensity, data=d15, family="nbinom2")
nt[[6]] <- glmmadmb(new_t1 ~ sr * TotDensity + TotDensity2, data=d15, family="nbinom2")
nt[[7]] <- glmmadmb(new_t1 ~ sr * TotDensity + TotDensity2*sr, data=d15, family="nbinom2")

# Model average
grow_new_t15      <- AICtab(nt,weights=T)
grow_new_t15_avg  <- model_avg(grow_new_t15, nt)
write.csv(grow_new_t15_avg, "Results/VitalRates_3/growth_new_t_15.csv", row.names = F)


# 2014 --------------------------------------------------------------------------------
d14     <- na.omit(unique(select(d14, new_t1, new_t1_pc, plot, sr, TotDensity, TotDensity2, year)))

# Tillers per capita------------------------------------------------------------
nt=list()
# Effect of density
nt[[1]] <- lm(new_t1_pc ~ TotDensity, data=d14)
nt[[2]] <- lm(new_t1_pc ~ TotDensity + TotDensity2, data=d14)
#Effect of density + sex ratio
nt[[3]] <- lm(new_t1_pc ~ TotDensity + sr, data=d14)
nt[[4]] <- lm(new_t1_pc ~ TotDensity + TotDensity2 + sr, data=d14)
#Interactions b/w density and sex ratio
nt[[5]] <- lm(new_t1_pc ~ sr * TotDensity, data=d14)
nt[[6]] <- lm(new_t1_pc ~ sr * TotDensity + TotDensity2, data=d14)
nt[[7]] <- lm(new_t1_pc ~ sr * TotDensity + TotDensity2*sr, data=d14)

# Model average
grow_new_tpc14      <- AICtab(nt,weights=T)
grow_new_tpc14_avg  <- model_avg(grow_new_tpc14, nt)
write.csv(grow_new_tpc14_avg, "Results/VitalRates_3/growth_new_tpc_14.csv", row.names = F)


# New tillers ------------------------------------------------------------
nt=list()
# Effect of density
nt[[1]] <- glmmadmb(new_t1 ~ TotDensity, data=d14, family = "nbinom2")
nt[[2]] <- glmmadmb(new_t1 ~ TotDensity + TotDensity2, data=d14, family="nbinom2")
#Effect of density + sex ratio
nt[[3]] <- glmmadmb(new_t1 ~ TotDensity + sr, data=d14, family="nbinom2")
nt[[4]] <- glmmadmb(new_t1 ~ TotDensity + TotDensity2 + sr, data=d14, family="nbinom2")
#Interactions b/w density and sex ratio
nt[[5]] <- glmmadmb(new_t1 ~ sr * TotDensity, data=d14, family="nbinom2")
nt[[6]] <- glmmadmb(new_t1 ~ sr * TotDensity + TotDensity2, data=d14, family="nbinom2")
nt[[7]] <- glmmadmb(new_t1 ~ sr * TotDensity + TotDensity2*sr, data=d14, family="nbinom2")

# Model average
grow_new_t14      <- AICtab(nt,weights=T)
grow_new_t14_avg  <- model_avg(grow_new_t14, nt)
write.csv(grow_new_t14_avg, "Results/VitalRates_3/growth_new_t_14.csv", row.names = F)


# 2014 + 2015--------------------------------------------------------------------------------
tmp_14  <- select(d14,plot,new_t1,new_t1_pc,sr,TotDensity,TotDensity2,year)
tmp_15  <- select(d15,plot,new_t1,new_t1_pc,sr,TotDensity,TotDensity2,year)
d_all   <- rbind(tmp_14,tmp_15)
d_all   <- mutate(d_all, year = as.factor(year))

# tillers per capita
nt=list()
# Effect of density
nt[[1]] <- lm(new_t1_pc ~ TotDensity, data=d_all)
nt[[2]] <- lm(new_t1_pc ~ TotDensity + TotDensity2, data=d_all)
#Effect of density + sex ratio
nt[[3]] <- lm(new_t1_pc ~ TotDensity + sr, data=d_all)
nt[[4]] <- lm(new_t1_pc ~ TotDensity + TotDensity2 + sr, data=d_all)
#Interactions b/w density and sex ratio
nt[[5]] <- lm(new_t1_pc ~ sr * TotDensity, data=d_all)
nt[[6]] <- lm(new_t1_pc ~ sr * TotDensity + TotDensity2, data=d_all)
nt[[7]] <- lm(new_t1_pc ~ sr * TotDensity + TotDensity2*sr, data=d_all)

nt[[8]] <- lm(new_t1_pc ~ TotDensity + year, data=d_all)
nt[[9]] <- lm(new_t1_pc ~ TotDensity + TotDensity2 + year, data=d_all)
#Effect of density + sex ratio
nt[[10]] <- lm(new_t1_pc ~ TotDensity + sr + year, data=d_all)
nt[[11]] <- lm(new_t1_pc ~ TotDensity + TotDensity2 + sr + year, data=d_all)
#Interactions b/w density and sex ratio
nt[[12]] <- lm(new_t1_pc ~ sr * TotDensity + year, data=d_all)
nt[[13]] <- lm(new_t1_pc ~ sr * TotDensity + TotDensity2 + year, data=d_all)
nt[[14]] <- lm(new_t1_pc ~ sr * TotDensity + TotDensity2*sr + year, data=d_all)

# model average
till_years_tpc_sel  <- AICtab(nt, weights = T)
till_years_tpc_avg  <- model_avg(till_years_tpc_sel, nt)
write.csv(till_years_tpc_avg, "Results/VitalRates_3/till_years_tpc_avg.csv", row.names = F)


# Number of tillers
nt=list()
# Effect of density
nt[[1]] <- glmmadmb(new_t1 ~ TotDensity, data=d_all)
nt[[2]] <- glmmadmb(new_t1 ~ TotDensity + TotDensity2, data=d_all)
#Effect of density + sex ratio
nt[[3]] <- glmmadmb(new_t1 ~ TotDensity + sr, data=d_all)
nt[[4]] <- glmmadmb(new_t1 ~ TotDensity + TotDensity2 + sr, data=d_all)
#Interactions b/w density and sex ratio
nt[[5]] <- glmmadmb(new_t1 ~ sr * TotDensity, data=d_all)
nt[[6]] <- glmmadmb(new_t1 ~ sr * TotDensity + TotDensity2, data=d_all)
nt[[7]] <- glmmadmb(new_t1 ~ sr * TotDensity + TotDensity2*sr, data=d_all)

# model average
till_years_t_sel  <- AICtab(nt, weights = T)
till_years_t_avg  <- model_avg(till_years_t_sel, nt)
write.csv(till_years_t_avg, "Results/VitalRates_3/till_years_t_avg.csv", row.names = F)


# One sex plots ---------------------------------------------------------------------

nt=list()
#Effect of density
nt[[1]] <- lm(new_t1 ~ year, data=one_sex)
nt[[2]] <- lm(new_t1 ~ TotDensity + year, data=one_sex)
nt[[3]] <- lm(new_t1 ~ TotDensity + year + sex, data=one_sex)
nt[[4]] <- lm(new_t1 ~ TotDensity + sex, data=one_sex)
#Effect of density + sex ratio
nt[[5]] <- glmmadmb(new_t1 ~ year+ (1 | plot), data=one_sex, family="nbinom2")
nt[[6]] <- glmmadmb(new_t1 ~ TotDensity + year + (1 | plot), data=one_sex,family="nbinom2")
nt[[7]] <- glmmadmb(new_t1 ~ TotDensity + year + sex + (1 | plot), data=one_sex,family="nbinom2")
nt[[8]] <- glmmadmb(new_t1 ~ TotDensity + sex + (1 | plot), data=one_sex,family="nbinom2")
AICtab(nt,weights=T)

oneSex_sel      <- AICtab(nt,weights=T)
oneSex_new_avg  <- model_avg(oneSex_sel, nt)
write.csv(oneSex_new_avg, "Results/VitalRates_3/oneSex_new_avg.csv", row.names = F)




# Graph -----------------------------------------------------------------------------

# Sex ratio as dot color ------
# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(topo.colors(7))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

#Graph: total number of tillers -----------------------------------------------------
tiff("Results/VitalRates_3/growth_newTillers.tiff",unit="in",width=6.3,height=9.45,res=600,compression="lzw")

par(mfrow=c(3,2),mar=c(3,2.5,1,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0.1))

# 2015 ----------------------------------------------------------------------------
# Per capita amount of new tillers
plot(d15$TotDensity,d15$new_t1_pc,pch=21,ylab="(New tillers in plot 2015)/Planting density",
     xlab="Planting density", bg = cRamp(d15$sr), cex = 1.5)
beta  <-  grow_new_tpc15_avg[,c("predictor","avg")]$avg
xSeq  <-  seq(0,48,by=1)
yL    <-  beta[1] + beta[2]*xSeq + beta[3]*xSeq^2 + beta[4]*0.1 + beta[5]*xSeq*0.1
yH    <-  beta[1] + beta[2]*xSeq + beta[3]*xSeq^2 + beta[4]*0.9 + beta[5]*xSeq*0.9
lines(xSeq,yL,col="red",lwd=2, lty = 2)
lines(xSeq,yH,col="blue",lwd=2)
legend(10, 19, c("10% female plots", "90% female plots"), cex = 1.5,
       lty = c(1,2), lwd=2, col=c("red","blue"), bty = "n")
title("2015")

# total number of tillers
plot(d15$TotDensity,d15$new_t1,pch=21,ylab="Number of new tillers in plot (2015)",
     xlab="Planting density", bg = cRamp(d15$sr), cex = 1.5)
beta  <-  grow_new_t15_avg[,c("predictor","avg")]$avg
xSeq  <-  seq(0,48,by=1)
yL    <-  exp( beta[1] + beta[2]*xSeq + beta[3]*xSeq^2 + beta[4]*0.1 + beta[5]*xSeq*0.1 )
yH    <-  exp( beta[1] + beta[2]*xSeq + beta[3]*xSeq^2 + beta[4]*0.9 + beta[5]*xSeq*0.9 )
lines(xSeq,yL,col="red",lwd=2, lty = 2)
lines(xSeq,yH,col="blue",lwd=2)
title("2015")

# 2014 ----------------------------------------------------------------------------
# Per capita amount of new tillers
plot(d14$TotDensity,d14$new_t1_pc,pch=21,ylab="(New tillers in plot 2015)/Planting density",
     xlab="Planting density", bg = cRamp(d14$sr), cex = 1.5)
beta  <-  grow_new_tpc14_avg[,c("predictor","avg")]$avg
xSeq  <-  seq(0,48,by=1)
yL    <-  beta[1] + beta[2]*xSeq + beta[3]*xSeq^2 + 
          beta[4]*0.1 + beta[5]*xSeq*0.1 + beta[6]*xSeq^2*0.1 
yH    <-  beta[1] + beta[2]*xSeq + beta[3]*xSeq^2 + 
          beta[4]*0.9 + beta[5]*xSeq*0.9 + beta[6]*xSeq^2*0.9
lines(xSeq,yL,col="red",lwd=2, lty = 2)
lines(xSeq,yH,col="blue",lwd=2)
title("2014")

# total number of tillers
plot(d14$TotDensity,d14$new_t1,pch=21,ylab="Number of new tillers in plot (2015)",
     xlab="Planting density", bg = cRamp(d14$sr), cex = 1.5)
beta  <-  grow_new_t14_avg[,c("predictor","avg")]$avg
xSeq  <-  seq(0,48,by=1)
yL    <-  exp( beta[1] + beta[2]*xSeq + beta[3]*xSeq^2 + beta[4]*0.1 + beta[5]*xSeq*0.1 )
yH    <-  exp( beta[1] + beta[2]*xSeq + beta[3]*xSeq^2 + beta[4]*0.9 + beta[5]*xSeq*0.9 )
lines(xSeq,yL,col="red",lwd=2, lty = 2)
lines(xSeq,yH,col="blue",lwd=2)
title("2014")

# One tiller plots
par(mar=c(2.5,2.5,0.1,0.1), mgp=c(1.4,0.5,0))
boxplot(new_t1 ~ plot_sex + year, data = one_sex,  
        ylab = "Number of new tillers in one sex plots",
        col = c("blue", "red"), cex.names=0.5, 
        names = c("2014", "2014", "2015", "2015") )
legend(0.5,225, c("female","male"), fill=c("blue","red"), bty = "n")
text(0.55,180,"Data from one-sex plots only",pos=4)

dev.off()






