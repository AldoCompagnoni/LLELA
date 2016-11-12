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

# Format -----------------------------------------------------------

# Year two
tmp15=subset(d,year==2015)
#remove 
tmp15$new_t1[tmp15$new_t1=="SKIPPED"]=NA
tmp15$new_t1[tmp15$new_t1=="cnf"]=NA
tmp15$new_t1[tmp15$new_t1==""]=NA
tmp15$new_t1=as.numeric(as.character(tmp15$new_t1))
tmp15$TotDensity2 = tmp15$TotDensity^2
d15=na.omit(tmp15[,c("l_t1","log_l_t0","plot","sex","new_t1","F","M",
                     "TotDensity","TotDensity2","sr")])


#remove 
tmp = d
tmp$new_t1[tmp$new_t1=="SKIPPED"]=NA
tmp$new_t1[tmp$new_t1=="cnf"]=NA
tmp$new_t1[tmp$new_t1==""]=NA
tmp$new_t1=as.numeric(as.character(tmp$new_t1))
d2=na.omit(tmp[,c("l_t1","log_l_t0","plot","sex","new_t1","F","M",
                     "TotDensity","sr","year")])


# One sex plots -----------------------------------------------------------
d2$oneSex          <- 0
d2$oneSex[d2$F==0 | 
            d2$M==0]   <- 1 # Flag treatments with only one sex 
one_sex           <- subset(d2,oneSex == 1)
# flag what is male, and what is female?
one_sex$plot_sex  <- "m"
one_sex$plot_sex[
  one_sex$M==0]   <- "f"
one_sex$plot_sex  <- factor(one_sex$plot_sex, levels = c("f", "m") )

# Model selection --------------------------------------------------------------

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


# Graph ----------------------------------------------------------------

tiff("Results/VitalRates_3/oneSexPlots_oneSexPlots.tiff",unit="in",
     width=3.5,height=3.5,res=600,compression="lzw")

par(mar=c(2.5,2.5,0.1,0.1), mgp=c(1.4,0.5,0))
boxplot(new_t1 ~ plot_sex + year, data = one_sex,  ylab = "Number of new tillers",
        col = c("blue", "red"), cex.names=0.5, 
        names = c("2014", "2014", "2015", "2015") )
legend(0.5,230, c("female","male"), fill=c("blue","red"), bty = "n")

dev.off()
