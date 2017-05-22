setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(bbmle)
library(lme4)
library(dplyr)
library(glmmADMB)
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")
source("C:/Users/ac79/Documents/CODE/LLELA/model_sel_results.R")

# load and format data --------------------------------------------------------------
d             <- read.csv("Data/vr.csv")
malPanicules  <- read.csv("Data/Spring 2014/maleCounts/malePaniculesSpring2014.csv")


# MALE panicule lengths --------------------------------------------------------------
d14                     <- subset(d, year==2014)
d14                     <- subset(d14, surv_t1!=0)
m_alloc                 <- merge(d14, malPanicules)
m_alloc$pan_area        <- m_alloc$panicule_Length_cm *  m_alloc$panicule_Width_cm * pi
m_alloc$plot            <- as.factor(m_alloc$plot)


# preliminary analysis ----------------------------------------------------------------------
spike_data              <- na.omit(dplyr::select(m_alloc, CountSpikelets, 
                                                 log_l_t0, sr, TotDensity, plot))

# should I use individual size or not? 
isMod = list()
isMod[[1]]=glmmadmb(CountSpikelets ~ (1 | plot), family = "nbinom2", data=spike_data)
isMod[[2]]=glmmadmb(CountSpikelets ~ log_l_t0 + (1 | plot), family = "nbinom2", data=spike_data)

# answer: no
AICtab(isMod, weights=T)


# Model selection -------------------------------------------------

# formulas for candidate models 
candidate_mods <- list(
  CountSpikelets ~ ( 1 | plot),
  CountSpikelets ~ TotDensity + ( 1 | plot),
  CountSpikelets ~ TotDensity + TotDensity:sr + ( 1 | plot)
)

# model specification
fit_func <- function(x){
  return( glmmadmb(x, family="nbinom2", data=spike_data) )
}

# fit the models 
plMod <- lapply(candidate_mods,fit_func)

# Model average
m_alloc_select    <- AICtab(plMod,weights=T)
m_alloc_avg       <- model_avg(m_alloc_select, plMod)
write.csv(m_alloc_avg, "Results/VitalRates_3/male_spikelets.csv", row.names = F)

# Model selection table
x <- plMod
write.csv(mod_sel_res(plMod),
          "Results/VitalRates_3/mod_sel_male_spikelets.csv",row.names=F)
rm(x)