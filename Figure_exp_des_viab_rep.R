# Plots on study replication
setwd("C:/Cloud/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation")
library(lme4) ; library(bbmle) ; library(boot) ; library(testthat)
library(tidyverse) ; library(glmmADMB)


# Read in data and format -----------------------------------------------------------

# select only information on experimental design 
d       <- read.csv("Data/vr.csv") %>%
              dplyr::select(plot, F, M) %>%
              unique

# select only germination information
viabVr  <- read.csv("Data/Spring 2014/viability/tetra_germ_plot_data.csv",
                    stringsAsFactors = F) %>%
              dplyr::select(plot, germTot, germFail) %>%
              subset( !is.na(germTot) & !is.na(germFail) )
            

# calculate replications -----------------------------------------------------------

# overall experimental replication
exp_des <- dplyr::select(d,plot,F,M) %>%
              unique %>%
              group_by(F,M) %>% 
              summarise(rep = n())


# total replication 
tot_rep  <- inner_join(viabVr, d) %>%
              group_by( F, M ) %>%
              summarise( rep = n() ) %>%
              ungroup 

# plot replication
plot_rep <- inner_join(viabVr, d) %>%
              dplyr::select(plot, F, M) %>%
              unique %>%
              group_by( F, M ) %>%
              summarise( rep = n() ) %>%
              ungroup     


# plots ------------------------------------------------------------------------------

tiff("Results/VitalRates_4_Cade2015/exp_design_viab_rep.tiff",unit="in",
     width=6.3,height=6.3,res=600,compression="lzw")

# total replicates
pt_exp  <- 2
cex_lab <- 1
cex_pan <- 1.1

par(mfrow = c(2,2), mar = c(3,3,0.3,0.1), mgp = c(1.8,0.5,0) )

# experimental design
plot(M ~ F, pch = 21, lwd = 2, data = exp_des,
     xlab = "Female density", ylab = "Male density",
     cex = rep/pt_exp, cex.lab = cex_lab)
text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(a)", cex = cex_pan, xpd = T)

# plot replicates
plot(M ~ F, pch = 21, lwd = 2, data = plot_rep,
     xlab = "Female density", ylab = "Male density",
     cex = rep/pt_exp, cex.lab = cex_lab)
text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(b)", cex = cex_pan, xpd = T)

legend(x = "topright", bty="n", pch = 21, pt.cex = c(1,4,8)/pt_exp, 
       cex = cex_lab, lwd = 2, lty= NA,
       legend = c("1 replicate","2 replicates","8 replicates") )

plot(M ~ F, pch = 21, lwd = 2, data = tot_rep,
     xlab = "Female density", ylab = "Male density",
     cex = rep/pt_exp, cex.lab = cex_lab)
text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
     par("usr")[4]*0.97,"(c)", cex = cex_pan, xpd = T)

dev.off()


# 
# par(mfrow = c(1,2), mar = c(3.5,3.8,0.3,0.1), mgp = c(1.8,0.5,0) )
# 
# # experimental design
# plot(M ~ F, pch = 21, lwd = 2, data = exp_des,
#      xlab = "Female density", ylab = "Male density",
#      cex = rep/pt_exp, cex.lab = cex_lab)
# text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
#      par("usr")[4]*0.97,"(a)", cex = cex_pan, xpd = T)
# 
# # plot replicates
# plot(M ~ F, pch = 21, lwd = 2, data = plot_rep,
#      xlab = "Female density", ylab = "Male density",
#      cex = rep/pt_exp, cex.lab = cex_lab)
# text(par("usr")[1] - (par("usr")[2] - par("usr")[1])*0.11,
#      par("usr")[4]*0.97,"(b)", cex = cex_pan, xpd = T)
# 
# legend(x = "topright", bty="n", pch = 21, pt.cex = c(1,2,4)/pt_exp, 
#        cex = cex_lab, lwd = 2, lty= NA,
#        legend = c("1 replicate","2 replicates","4 replicates") )
