# Code to produce figure1 (panels a,b,c)
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(dplyr)
source("C:/Users/ac79/Documents/CODE/LLELA/model_avg.R")

#Read data--------------------------------------------------------------------
d       <- read.csv("Data/vr.csv")

# FIGURE ON EXPERIMENTAL DESIGN ----------------------------------------------------------------------
exp_des <- select(d,plot,F,M) %>%
              distinct() %>%
              group_by(F,M) %>% 
              summarise(rep = n())

# plot
pt_exp  <- 1.5
tiff("Results/VitalRates_3/exp_design.tiff",unit="in",
     width=4,height=4,res=600,compression="lzw")

par(mar = c(3.5,3.8,0.1,0.1), mgp = c(2.3,1,0))

plot(exp_des$F, exp_des$M, pch = 21, lwd = 2,
     xlab = "Female density", ylab = "Male density",
     cex = exp_des$rep/pt_exp, cex.lab = 1.5)

legend(x = "topright", bty="n", pch = 21, pt.cex = c(1,2,4)/pt_exp, 
       cex = 1.2, lwd = 2, lty= NA,
       legend = c("1 replicate","2 replicates","4 replicates") )

dev.off()
