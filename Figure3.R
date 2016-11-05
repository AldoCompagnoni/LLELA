# Code to produce figure1 (panels a,b,c)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(MASS)
library(dplyr)
library(boot)

#Read in data--------------------------------------------------------------------
# raw data
d                 <- read.csv("Data/vr.csv")
tillers           <- read.csv("Data/one_sex_plots.csv", stringsAsFactors = F)

# best models
grow_avg          <- read.csv("Results/VitalRates_3/growth_best.csv")
one_sex_grow_avg  <- read.csv("Results/VitalRates_3/one_sex_growth_best.csv")


# FORMAT DATA -------------------------------------------------------------------

# remove dead individuals (this is a GROWTH model!)
d=subset(d, surv_t1 != 0)

# only use year two
tmp15=subset(d, year==2015)
tmp15$new_t1[tmp15$new_t1=="SKIPPED"]=NA
tmp15$new_t1[tmp15$new_t1=="cnf"]=NA
tmp15$new_t1[tmp15$new_t1==""]=NA
tmp15$new_t1=as.numeric(as.character(tmp15$new_t1))
d15=na.omit(tmp15[,c("l_t1","log_l_t0","plot","sex","new_t1","sr")])


# FIGURE 3 ----------------------------------------------------------------------

tiff("Results/VitalRates_3/Figure3.tiff",unit="in",width=6.3,height=3.5,res=600,compression="lzw")

par(mfrow=c(1,2),mar=c(2.5,2.5,1.3,0.5),mgp=c(1.4,0.35,0),cex.lab=1.1,cex.axis=0.8,
    cex.main=0.9)

# Panel a -  Color code sexes --------------------------------------------------------------
tillers$col=as.integer(as.factor(tillers$sex))
tillers$col=as.character(factor(tillers$col,labels=c("blue","red")))

plot(tillers$tot_till_t0,tillers$tot_till_t1, pch=16, 
     ylab="Tillers per plot (spring 2015)",
     xlab="Tillers per plot (spring 2014)",col=tillers$col)
legend(0,235,c("Male individuals","Female individuals"),
       lty=1,pch=16,lwd=3,col=c("red","blue"),bty="n")

xSeq  <- seq(0,max(c(tillers$tot_till_t0,tillers$tot_till_t0)),by=0.1)
betas <- one_sex_grow_avg[,c("predictor","avg")]$avg
y_m   <- exp( betas[1] + betas[2] * xSeq + betas[3] + betas[4]*xSeq)
y_f   <- exp( betas[1] + betas[2] * xSeq )

lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.03,
     par("usr")[4]*1.05,"a) One sex plots", cex = 1.2, xpd = T, pos = 4)

# Panel b - Set up colors for plots -------------------------------------------------------
d15$col=as.integer(d15$sex)
d15$col=as.character(factor(d15$col,labels=c("blue","red")))
d15$symb=as.integer(as.character(factor(d15$col,labels=c("17","16"))))

plot(d15$sr,d15$l_t1, pch=16, ylab="Number of leaves per individual",
     xlab="Proportion of female individuals",col=d15$col, ylim = c(0,82))
beta = grow_avg[,c("predictor","avg")]$avg
xSeq <- seq(0,1,by=0.1)
size <- mean(d15$log_l_t0)
dens <- 1
y_m <- exp( beta[1] + beta[2]*size + beta[3]*dens + beta[4] + 
              beta[5]*xSeq + beta[6]*size + beta[7]*dens + beta[8]*xSeq)
y_f <- exp( beta[1] + beta[2]*size + beta[3]*dens + beta[5]*xSeq )
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")
dens <- 48
y_m <- exp( beta[1] + beta[2]*size + beta[3]*dens + beta[4] + 
              beta[5]*xSeq + beta[6]*size + beta[7]*dens + beta[8]*xSeq)
y_f <- exp( beta[1] + beta[2]*size + beta[3]*dens + beta[5]*xSeq )
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.03,
     par("usr")[4]*1.05,"b) Intersexual competition", cex = 1.2, xpd = T, pos = 4)

dev.off()



tiff("Results/VitalRates_3/Figure3_a.tiff",unit="in",width=6.3,height=3.5,res=600,compression="lzw")

par(mfrow=c(1,2),mar=c(2.5,2.5,1.3,0.5),mgp=c(1.4,0.35,0),cex.lab=1.1,cex.axis=0.8,
    cex.main=0.9)

# Panel a -  Color code sexes --------------------------------------------------------------
tillers$col=as.integer(as.factor(tillers$sex))
tillers$col=as.character(factor(tillers$col,labels=c("blue","red")))

boxplot(tot_till_t1 ~ sex, data = tillers, col = c("blue","red"),
        names = c("Females","Males"),
        ylab = "Total number of tillers (2015)")

#legend(0.5,235,c("Male individuals","Female individuals"),fill=c("red","blue"),bty="n")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.03,
     par("usr")[4]*1.05,"a) One sex plots", cex = 1.2, xpd = T, pos = 4)

# Panel b - Set up colors for plots -------------------------------------------------------
d15$col=as.integer(d15$sex)
d15$col=as.character(factor(d15$col,labels=c("blue","red")))
d15$symb=as.integer(as.character(factor(d15$col,labels=c("17","16"))))

plot(d15$sr,d15$l_t1, pch=16, ylab="Number of leaves per individual",
     xlab="Proportion of female individuals",col=d15$col, ylim = c(0,82))
beta = grow_avg[,c("predictor","avg")]$avg
xSeq <- seq(0,1,by=0.1)
size <- mean(d15$log_l_t0)
dens <- 1
y_m <- exp( beta[1] + beta[2]*size + beta[3]*dens + beta[4] + 
              beta[5]*xSeq + beta[6]*size + beta[7]*dens + beta[8]*xSeq)
y_f <- exp( beta[1] + beta[2]*size + beta[3]*dens + beta[5]*xSeq )
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")
dens <- 48
y_m <- exp( beta[1] + beta[2]*size + beta[3]*dens + beta[4] + 
              beta[5]*xSeq + beta[6]*size + beta[7]*dens + beta[8]*xSeq)
y_f <- exp( beta[1] + beta[2]*size + beta[3]*dens + beta[5]*xSeq )
lines(xSeq,y_f,lwd=3,lty=1,col="blue")
lines(xSeq,y_m,lwd=3,lty=1,col="red")

legend(0,88.5,c("Male individuals","Female individuals"),
       lty=1,pch=16,lwd=3,col=c("red","blue"),bty="n")

text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.03,
     par("usr")[4]*1.05,"b) Intersexual competition", cex = 1.2, xpd = T, pos = 4)

dev.off()
