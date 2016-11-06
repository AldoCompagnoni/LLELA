##Response variable: total number of tillers 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) 
library(glmmADMB) # Fit models with a Negative Binomial
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/analysis/model_avg.R")

# load and format data -----------------------------------------------------
d=read.csv("Data/vr.csv")
#remove dead individuals (this is a GROWTH model!)
d=subset(d, surv_t1 != 0)
#logtransform leaf numbers
d$plot=as.factor(d$plot) #glmmadmb wants plot as a factor

# Year two
tmp15=subset(d,year==2015)
# remove 
tmp15$new_t1[tmp15$new_t1=="SKIPPED"]=NA
tmp15$new_t1[tmp15$new_t1=="cnf"]=NA
tmp15$new_t1[tmp15$new_t1==""]=NA
tmp15$new_t1=as.numeric(as.character(tmp15$new_t1))
tmp15$TotDensity2 = tmp15$TotDensity^2
d15=na.omit(tmp15[,c("l_t1","log_l_t0","plot","sex","new_t1",
                     "TotDensity","TotDensity2","sr")])
d15$new_t1_pc <- d15$new_t1 / d15$TotDensity
d15 = subset(d15, plot != 149) # Remove outlier

# One sex plots -----------------------------------------------------------
#remove 
tmp = d
tmp$new_t1[tmp$new_t1=="SKIPPED"]=NA
tmp$new_t1[tmp$new_t1=="cnf"]=NA
tmp$new_t1[tmp$new_t1==""]=NA
tmp$new_t1=as.numeric(as.character(tmp$new_t1))
d2=na.omit(tmp[,c("l_t1","log_l_t0","plot","sex","new_t1","F","M",
                  "TotDensity","sr","year")])
d2$oneSex             <- 0
d2$oneSex[d2$F==0 | 
            d2$M==0]  <- 1 # Flag treatments with only one sex 
one_sex               <- subset(d2,oneSex == 1)
# flag what is male, and what is female?
one_sex$plot_sex      <- "m"
one_sex$plot_sex[
  one_sex$M==0]       <- "f"
one_sex$plot_sex      <- factor(one_sex$plot_sex, levels = c("f", "m") )


# Model selection ##############################################################
# PER CAPITA total number of tillers -------------------------------------------

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


# Total number of tillers ------------------------------------------------

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

#Graph: total number of tillers -----------------------------------------------------
tiff("Results/VitalRates_3/growth_newTillers.tiff",unit="in",width=6.3,height=6.3,res=600,compression="lzw")

par(mfrow=c(2,2),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.4,0.5,0),oma=c(0,0,0,0.1))

# Per capita amount of new tillers
plot(d15$TotDensity,d15$new_t1_pc,pch=16,ylab="Per capita new tillers (2015)",
     xlab="Planting density")
beta  <-  grow_new_tpc_avg[,c("predictor","avg")]$avg
xSeq  <-  seq(0,48,by=1)
yL    <-  beta[1] + beta[2]*0.1 + beta[3]*xSeq*0.1 + beta[4]*xSeq + 
          beta[5]*xSeq^2 + beta[6]*xSeq^2*0.1
yH    <-  beta[1] + beta[2]*0.9 + beta[3]*xSeq*0.9 + beta[4]*xSeq + 
          beta[5]*xSeq^2 + beta[6]*xSeq^2*0.9
lines(xSeq,yL,col="red",lwd=2, lty = 2)
lines(xSeq,yH,col="blue",lwd=2)
legend(10, 19, c("10% female plots", "90% female plots"),
       lty = c(1,2), lwd=2, col=c("red","blue"), bty = "n")

# total number of tillers
plot(d15$TotDensity,d15$new_t1,pch=16,ylab="Number of new tillers (2015)",
     xlab="Planting density")
beta  <-  grow_new_t_avg[,c("predictor","avg")]$avg
xSeq  <-  seq(0,48,by=1)
yL    <-  beta[1] + beta[2]*0.1 + beta[3]*xSeq*0.1 + 
          beta[4]*xSeq^2*0.1 + beta[5]*xSeq + beta[6]*xSeq^2
yH    <-  beta[1] + beta[2]*0.9 + beta[3]*xSeq*0.1 + 
          beta[4]*xSeq^2*0.9 + beta[5]*xSeq + beta[6]*xSeq^2
lines(xSeq,yL,col="red",lwd=2, lty = 2)
lines(xSeq,yH,col="blue",lwd=2)

# One tiller plots
par(mar=c(2.5,2.5,0.1,0.1), mgp=c(1.4,0.5,0))
boxplot(new_t1 ~ plot_sex + year, data = one_sex,  ylab = "Number of new tillers",
        col = c("blue", "red"), cex.names=0.5, 
        names = c("2014", "2014", "2015", "2015") )
legend(0.5,225, c("female","male"), fill=c("blue","red"), bty = "n")
text(0.55,180,"Data from one-sex plots only",pos=4)

dev.off()

