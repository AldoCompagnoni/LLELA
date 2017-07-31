# Seed viability in female-only plots
# Data: tetrazolium scoring and germination data 
setwd("C:/Cloud/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
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
            # male_presence == 0 is the effect of male ABSENCE, which is what we are estimating here
            mutate( male_presence = replace(male_presence, M > 0, 0) )

viab_of <- glm(cbind(germTot,germFail) ~ male_presence, 
               data = one_f, family = "binomial") 

summary(viab_of)$coefficients

# proportions of germination
beta_of <- coef(viab_of)
df_one_f<- data.frame( viab_one_fem = inv.logit( sum(beta_of) ),
                       viab_with_m  = inv.logit( beta_of[1] ) )


# viability: ALL female-only plots -------------------------------------------
only_f  <- subset(viabVr, F > 0 ) %>%
              mutate( male_presence = 1 ) %>%
              mutate( male_presence = replace(male_presence, M > 0, 0) )
# "all females"
viab_af <- glm(cbind(germTot,germFail) ~ male_presence, 
               data = only_f, family = "binomial") 

summary(viab_af)$coefficients

# proportions of germination
beta_af <- coef(viab_af)
df_o_f  <- data.frame( viab_only_fem = inv.logit( sum(beta_af) ),
                       viab_with_m  = inv.logit( beta_af[1] ) )
