setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")

#Fit vital rates models

#Three appraches to growth models
source("Analysis/growth_oneSexPlots.R") ; rm(list=ls())
source("Analysis/growth_sexRatio.R") ; rm(list=ls())
source("Analysis/growth_newTillers.R") ; rm(list=ls())

#Survival models
