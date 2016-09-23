# Compare growth between male-only and female-only plots
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(bbmle) #For AIC weights, because I'm lazy!
library(glmmADMB) #Fit models with a Negative Binomial
library(dplyr) ; library(testthat)

# Read and format data ======================================================================
d <- read.csv("Data/vr.csv")
# Remove dead individuals (this is a growth model)
d <- subset(d, surv_t1 != 0)
# Remove a clear mistake: missing data in 2013 and 2014, but 19 leaves in 2015
d <- d[-which(d$plot ==36 & d$focalI =="m3"),]

# Logtransform leaf numbers
d$log_l_t0 <- log(d$l_t0)
d$log_l_t1 <- log(d$l_t1)
# glmmadmb wants plot as a factor
d$plot     <- as.factor(d$plot) 
# Transform densities to SEX RATIO
d$sr = d$F / d$TotDensity


# Compare growth between female only and male only plots ==========================================
d$oneSex  <- 0
d$oneSex[d$F==0 | d$M==0] <- 1 # Flag treatments with only one sex 
oneSexD   <- subset(d,oneSex == 1)

# I start from tillers in year 2014
tmp14            <- subset(oneSexD,year == 2014)
focal_tillers_14 <- tmp14 %>% 
                    group_by(plot,sex) %>% 
                    summarise(till_t0 = sum(till_t1,na.rm=T))
new_tillers_14        <- distinct(select(tmp14,plot,sex,new_t1))
new_tillers_14$new_t1 <- as.numeric(as.character(new_tillers_14$new_t1))
# calculate total number of tillers
tillers_14        <- merge(focal_tillers_14,new_tillers_14)
tillers_14        <- mutate(tillers_14,tot_till_t0 = till_t0 + new_t1)
tillers_14        <- select(tillers_14,plot,sex,tot_till_t0)

# Total tillers per plot in 2015
tmp15     <- subset(oneSexD,year == 2015)
# Exclude plots for which we do not have "new tillers"
tmp15$new_t1[tmp15$new_t1=="SKIPPED"] <- NA
tmp15$new_t1[tmp15$new_t1=="cnf"]     <- NA
tmp15$new_t1[tmp15$new_t1==""]        <- NA
tmp15$new_t1                          <- as.numeric(as.character(tmp15$new_t1))
# Total number of tillers, only focal individuals 
focal_tillers_15  <- tmp15 %>% 
                     group_by(plot,sex) %>% 
                     summarise(till_t1 = sum(till_t1),
                               TotDensity = mean(TotDensity))
new_tillers_15    <- distinct(select(tmp15,plot,sex,new_t1))
# check for mistakes (left columns should only contain value 1)
expect_equal(aggregate(new_t1 ~ plot,length,data=new_tillers_15)$new_t1,
             rep(1,nrow(new_tillers_15)))
# calculate total number of tillers
tillers_15        <- merge(focal_tillers_15,new_tillers_15)
tillers_15        <- mutate(tillers_15,tot_till_t1 = till_t1 + new_t1)
tillers_15        <- select(tillers_15,plot,sex,tot_till_t1,TotDensity)
# finally, the tiller data set for 2015
tillers        <- merge(tillers_14,tillers_15)

 
# Model selections ==============================================================================

# single sex models 2015 -----------------------------------------------------
ssm15=list()

# Effect of sex. Analyses use "tiller numbers"
ssm15[[1]]=glmmadmb(tot_till_t1 ~ tot_till_t0, data=tillers,family="nbinom2")
ssm15[[2]]=glmmadmb(tot_till_t1 ~ tot_till_t0 + sex, data=tillers,family="nbinom2")
ssm15[[3]]=glmmadmb(tot_till_t1 ~ tot_till_t0 * sex, data=tillers,family="nbinom2")
ssm15[[4]]=glmmadmb(tot_till_t1 ~ tot_till_t0 + sex + TotDensity, data=tillers,family="nbinom2")
ssm15[[5]]=glmmadmb(tot_till_t1 ~ tot_till_t0 * sex + TotDensity, data=tillers,family="nbinom2")
ssm15[[6]]=glmmadmb(tot_till_t1 ~ sex, data=tillers,family="nbinom2")
ssm15[[7]]=glmmadmb(tot_till_t1 ~ sex + TotDensity, data=tillers,family="nbinom2")
ssm15[[8]]=glmmadmb(tot_till_t1 ~ sex * TotDensity, data=tillers,family="nbinom2")

AICtab(ssm15,weights=T)



# single sex simplified models 2015 -----------------------------------------------------
sssm15=list()

# Effect of sex. Analyses use "tiller numbers"
sssm15[[1]]=glmmadmb(tot_till_t1 ~ sex, data=tillers,family="nbinom2")
sssm15[[2]]=glmmadmb(tot_till_t1 ~ sex + TotDensity, data=tillers,family="nbinom2")
sssm15[[3]]=glmmadmb(tot_till_t1 ~ sex * TotDensity, data=tillers,family="nbinom2")

AICtab(sssm15,weights=T)



#GRAPH ==================================================================================================

# Best two single "sex" models in 2015 ----------------------------------------------------------------------------
tiff("Results/VitalRates_Sept16/growth_OneSexPlots_2015.tiff",
     unit="in",width=6.3,height=3.15,res=600,compression="lzw")

# 2015
tillers$col=as.integer(tillers$sex)
tillers$col=as.character(factor(tillers$col,labels=c("blue","red")))

# Start plotting
par(mfcol=c(1,3),mar=c(3,2.8,1,0.62),mgp=c(1.4,0.5,0))

# best model
plot(tillers$tot_till_t0,tillers$tot_till_t1, pch=16, 
     ylab="Tot. number of tillers at time t",
     xlab="Tot. number of tillers at time t-1")
title(main = "Best model (48% weight)", line=0.2,cex=0.9)

betas=coef(ssm15[[1]])
xSeq <- seq(0,max(c(tillers$tot_till_t0,tillers$tot_till_t0)),by=0.1)
y    <- exp( betas[1] + betas[2]*xSeq )
lines(xSeq,y,lwd=3,lty=1)

# 2nd best
plot(tillers$tot_till_t0,tillers$tot_till_t1, pch=16, 
     ylab="Tot. number of tillers at time t",
     xlab="Tot. number of tillers at time t-1",col=tillers$col)
title(main = "2nd best model (18% weight)", line=0.2,cex=0.9)
legend(0,235,c("Male individuals","Female individuals"),
       lty=1,lwd=3,col=c("red","blue"),bty="n")

betas=coef(ssm15[[2]])
y_f    <- exp( betas[1] + betas[2] * xSeq )
y_m    <- exp( betas[1] + betas[2] * xSeq + betas[3])
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")

# 3rd best
plot(tillers$tot_till_t0,tillers$tot_till_t1, pch=16, 
     ylab="Tot. number of tillers at time t",
     xlab="Tot. number of tillers at time t-1",col=tillers$col)
title(main = "3rd best model (18% weight)", line=0.2,cex=0.9)

betas     <- coef(ssm15[[4]])
meanDens  <- mean(tillers$TotDensity)
y_f       <- exp( betas[1] + betas[2] * xSeq + betas[4] * meanDens)
y_m       <- exp( betas[1] + betas[2] * xSeq + betas[3] + betas[4] * meanDens)
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")

dev.off()


# Best two single "sex" models in 2015 ----------------------------------------------------------------------------
tiff("Results/VitalRates_Sept16/growth_OneSexPlots_simple.tiff",
     unit="in",width=6.3,height=3.15,res=600,compression="lzw")

# Start plotting
par(mfcol=c(1,2),mar=c(3,2.5,1,0.62),mgp=c(1.4,0.5,0))

# best model
plot(tillers$TotDensity,tillers$tot_till_t1, pch=16, 
     ylab="Tot. number of tillers at time t",xlab="Total plot density",
     col=tillers$col)
title(main = "Best model (70% weight)", line=0.2,cex=0.9)

betas=coef(sssm15[[2]])
xSeq <- seq(0,max(tillers$TotDensity),by=1)
y_f  <- exp( betas[1] + betas[3]*xSeq )
y_m  <- exp( betas[1] + betas[2] + betas[3]*xSeq )
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")

legend(-2.5,190,c("Male individuals","Female individuals"),
       lty=1,lwd=3,col=c("red","blue"),bty="n")


# 2nd best model
plot(tillers$TotDensity,tillers$tot_till_t1, pch=16, 
     ylab="Tot. number of tillers at time t",xlab="Total plot density",
     col=tillers$col)
title(main = "2nd best model (27% weight)", line=0.2,cex=0.9)

betas=coef(sssm15[[3]])
xSeq <- seq(0,max(tillers$TotDensity),by=1)
y_f  <- exp( betas[1] + betas[3]*xSeq )
y_m  <- exp( betas[1] + betas[2] + betas[3]*xSeq + betas[4]*xSeq)
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")

dev.off()
