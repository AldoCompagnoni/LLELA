# Codes produces summaries for the manuscript 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/")
library(dplyr)

# Summaries on replication
d           <- read.csv("Data/vr.csv")
trt_rep     <-  select(d,plot,F,M) %>% 
                  distinct() %>% 
                  group_by(F,M) %>% 
                  summarise(rep = n())

# information on replicates
rep_summ    <- data.frame(total_reps      = sum(trt_rep$rep),
                          one_male_reps   = subset(trt_rep, M == 1 & F == 0)$rep,
                          one_female_reps = subset(trt_rep, M == 0 & F == 1)$rep)

# total n. of planted individuals
plant_idv   <- sum( (trt_rep$F + trt_rep$M) * trt_rep$rep )


# summaries on seed ~ weight/length relationship
seeds       <- read.csv("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")
mod_weight  <- lm(SeedN ~ Seed.weight..mg., data = seeds)
mod_length  <- lm(SeedN ~ Panicule.Length_.cm., data = seeds)

plot(log(seeds$Seed.weight..mg.), log(seeds$SeedN) )


# summary from COMADRE database -------------------------------------

# load latest version of COMADRE
load("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/Data/COMADRE_v.2.0.1.RData")

# 
#mono_mat    <- sum( is.na(comadre$metadata$StudiedSex) )
#dioec_mat   <- length( comadre$metadata$StudiedSex ) - mono_mat
#malFem_mat  <- sum( comadre$metadata$StudiedSex == "M/F", na.rm=T)

# tot n. of spp
n_spp       <- length(unique(comadre$metadata$SpeciesAccepted))

# total n. of studies
n_studies   <- length(unique(comadre$metadata$DOI.ISBN))

# number of studies wich studied sex 
stu_sex     <- unique(select(comadre$metadata, StudiedSex, DOI.ISBN))
stu_sex     <- nrow(subset(stu_sex, StudiedSex == "M/F"))



# pH buffer solution ------------------------------------------

# potassium phosphate concentration
((9.078*2) / 5000) * 100
# sodium phosphate  concentration
((9.472*3) / 5000) * 100
