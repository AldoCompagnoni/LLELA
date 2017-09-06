# model selection for fecundity (seeds per flower) 
setwd("C:/cloud/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle) 
library(MASS)
library(dplyr)
library(glmmADMB)
source("C:/CODE/LLELA/model_avg.R")
source("C:/CODE/LLELA/model_sel_results.R")

#read in data--------------------------------------------------------------
d         <- read.csv("Data/vr.csv")
fem_seeds <- read.csv("Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")

# plot level data ---------------------------------------------------------

# Only data from 2014
d14 <- subset(d, year == 2014)
d14 <- subset(d14, surv_t1 != 0)

# fecundity data ---------------------------------------------------------------------- 
fem_seeds$focalI    <- paste("f",fem_seeds$IndividualN,sep="")
fem_seeds           <- fem_seeds[,c("Plot","focalI","SeedN")] 
fem_seeds           <- fem_seeds[!is.na(fem_seeds$SeedN),]
names(fem_seeds)[1] <- "plot"

# merge data sets 
fecund_data      <- merge(d14,fem_seeds) 
fecund_data$plot <- as.factor(fecund_data$plot)

# Preliminary analyses ---------------------------------------------------------------------

# should I use individual size or not? 
isMod = list()
isMod[[1]]=glmmadmb(SeedN ~ (1 | plot), family = "nbinom2", data=fecund_data)
isMod[[2]]=glmmadmb(SeedN ~ log_l_t0 + (1 | plot), family = "nbinom2", data=fecund_data)

# answer: no.
AICtab(isMod, weights=T)


# MODEL SELECTION --------------------------------------------------------

# formulas for candidate models 
candidate_mods <- list(
  SeedN ~ ( 1 | plot),
  SeedN ~ TotDensity + ( 1 | plot),
  SeedN ~ TotDensity + TotDensity:sr + ( 1 | plot)
)

# fit models
fec_mod     <- lapply(candidate_mods, function(x) glmmadmb(x, family="nbinom2", data=fecund_data) )
fec_select  <- AICtab(fec_mod, weights=T)
fec_avg     <- model_avg(fec_select, fec_mod)
write.csv(fec_avg, "Results/VitalRates_4_Cade2015/fecuntity_best.csv", row.names = F)

# Model selection table
x <- fec_mod # trick to feed AICtab the models within a function
write.csv(mod_sel_res(fec_mod),
          "Results/VitalRates_4_Cade2015/mod_sel_fecuntity.csv",row.names=F)
rm(x)


# Sex ratio as dot color -------------------------------------------------

# service functions
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x){
  cols <- colorRamp(gray.colors(7,start=0,end=1))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}  

# Symbols for the sexes 
sex_symbol <- function(x){
  out <- x %>% 
    mutate(symb = factor(sex, labels=c("21","24")) ) %>%
    mutate(symb = as.character(symb) ) %>%
    mutate(symb = as.integer(symb) )
  return(out)
}

# Graph -----------------------------------------------------------------------------------

# parameters
cex_dots <- 1.5
cex_pan  <- 1.3
cex_leg  <- 1.3

# fecundity 
plot(jitter(fecund_data$TotDensity), fecund_data$SeedN, pch = 21, xlim = c(0,48.5),
     bg = cRamp(fecund_data$sr) , ylim = c(0, 1000), cex = cex_dots,
     xlab = "Planting density", ylab = "Seeds per female panicle")

fec_h <- subset(fec_avg, sr == 0.95)
fec_l <- subset(fec_avg, sr == 0.05)

lines(fec_h$TotDensity,fec_h$pred,lty=1,lwd=2,col="#ABABAB")
lines(fec_l$TotDensity,fec_l$pred,lty=2,lwd=2,col="#141414")

