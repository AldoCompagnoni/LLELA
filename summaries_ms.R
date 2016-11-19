#summaries_ms.R

# Summaries on replication
d           <- read.csv("Data/vr.csv")
trt_rep     <- distinct(select(d,plot,F,M)) %>% 
                group_by(F,M) %>% 
                  summarise(rep = n())

# information on replicates
rep_summ    <- data.frame(total_reps      = sum(trt_rep$rep),
                          one_male_reps   = subset(trt_rep, M == 1 & F == 0)$rep,
                          one_female_reps = subset(trt_rep, M == 0 & F == 1)$rep)

# summaries on seed ~ weight/length relationship
seeds       <- read.csv("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/Data/Spring 2014/SeedCount/Poa Arachnifera_seedCount.csv")
mod_weight  <- lm(SeedN ~ Seed.weight..mg., data = seeds)
mod_length  <- lm(SeedN ~ Panicule.Length_.cm., data = seeds)

plot(log(seeds$Seed.weight..mg.), log(seeds$SeedN) )


