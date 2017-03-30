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

# Mean seed weight
mean(seeds$Seed.weight..mg. / seeds$SeedN, na.rm=T)


# table 1 main text -------------------------------------------------------
nflow <- read.csv("Results/VitalRates_3/n_flowers_mod_sel.csv")
viab  <- read.csv("Results/VitalRates_3/germination_mod_sel.csv")
fec   <- read.csv("Results/VitalRates_3/fecuntity_mod_sel.csv")
male  <- read.csv("Results/VitalRates_3/male_spikelets_mod_sel.csv")

mods  <- list(nflow,viab,fec,male)

# tab 1
select_cols <- function(x){ return(select(x,Equation, dAIC, weight))} 
all_mods    <- lapply(mods, select_cols )

# 5 alternative models 
no_eq_sel   <- function(x){ return(select(x, dAIC, weight)) }
mods_5      <- lapply(all_mods[2:4], no_eq_sel)
mods_5_df   <- Reduce(function(...) cbind(...), mods_5)
new_names   <- paste0( names(mods_5_df), c(rep("_viab",2),rep("_fec",2),rep("_spk",2)))
mods_5_df   <- setNames(mods_5_df, new_names)

# order 15 alternative models
nflow         <- nflow[c(1:5,8,6,9,7,10,14,13,12,11,15),]
all_sel       <- select_cols(nflow) %>%  
                  cbind(NA,NA,NA,NA,NA,NA) %>%
                  setNames(c("Equation","dAIC_pnc","weight_pnc",
                           names(mods_5_df))
                           )
all_sel[1,4:9] <- mods_5_df[1,]
all_sel[3,4:9] <- mods_5_df[2,]
all_sel[4,4:9] <- mods_5_df[3,]
all_sel[7,4:9] <- mods_5_df[4,]
all_sel[8,4:9] <- mods_5_df[5,]

# order columns as in tab 1 
all_sel_out    <- select(all_sel, Equation, 
                         dAIC_viab, weight_viab,
                         dAIC_pnc, weight_pnc,
                         dAIC_fec, weight_fec,
                         dAIC_spk, weight_spk)
write.csv(all_sel_out, "Results/VitalRates_3/table1.csv", row.names = F)
  

# table 2 main text ------------------------------------------------------
bh_14   <- read.csv("Results/VitalRates_3/new_tillers_BH14_mod_sel.csv")
bh_15   <- read.csv("Results/VitalRates_3/new_tillers_BH15_mod_sel.csv")

sel_bh  <- function(x,year) {  namez <- paste0(c("dAIC_","weights_"),year)
                               out <- x %>%
                                        select(model,deltaAIC,weights) %>%
                                        setNames(c("model",namez))
}

tab_bh  <- sel_bh(bh_14,"14") %>%
            merge(sel_bh(bh_15,"15"), sort = F)

write.csv(tab_bh, "Results/VitalRates_3/table2.csv", row.names=F)


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
study_df    <- unique(select(comadre$metadata, StudiedSex, DOI.ISBN, Journal))
sex_df      <- subset(study_df, StudiedSex == "M/F")
stu_sex     <- nrow(sex_df)

# studies DOI
sex_doi     <- sex_df %>%
                select(Journal, DOI.ISBN) %>%
                unique()
  
cbind(c(1:31),sex_doi)

no, NA, yes/no, no, yes, no, no, yes/no, yes, no
no, NA, NA, no, no, no, no, no, no, no


# pH buffer solution ------------------------------------------

# potassium phosphate concentration
((9.078*2) / 5000) * 100
# sodium phosphate  concentration
((9.472*3) / 5000) * 100


# Poa pratensis Seedling survival estimate
chai  <- read.csv("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/Data/literature_seedling_surv/poa_seedl_surv_chai_et_al.csv")
edwd  <- read.csv("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Response-Surface experiment/Experiment/Implementation/Data/literature_seedling_surv/poa_seedl_surv_edwards_et_al.csv")

quantile(as.matrix(chai), prob = c(0.1,0.9))
min(as.matrix(chai))
