# Seed viability in female-only plots
# Data: tetrazolium scoring and germination data 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(lme4) ; library(bbmle) ; library(boot) ; library(testthat)
library(dplyr) ; library(glmmADMB)
#source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")
#source("C:/Users/ac79/Documents/CODE/LLELA/model_sel_results.R")


# Read in data and format -----------------------------------------------------------
d       <- read.csv("Data/vr.csv")
viabVr  <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                    stringsAsFactors = F)

# Format data -----------------------------------------------------------------------
viabVr  <- viabVr %>%
            mutate(totN = F + M) %>%
            mutate(sr = F / totN) %>%
            mutate(plot = as.factor(plot) )


# viability: female only vs. more plots with at least 1 male -----------------
one_f   <- subset(viabVr, F == 1 ) %>%
            mutate( male_presence = 1 ) %>%
            mutate( male_presence = replace(male_presence, M > 0, 0) )

viab_m  <- glm(cbind(germTot,germFail) ~ male_presence, 
               data = one_f, family = "binomial") 

summary(viab_m)$coefficients
