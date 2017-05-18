# model selection result table
sel_results <- function(x, candidate_num, response){

  if(candidate_num == 5)  { 
    mod_str <- read.csv("Results/VitalRates_3/5_model_structures.csv")
    mod_str <- mutate(mod_str, Equation = gsub("Response",response,Equation) )
  }
  if(candidate_num == 15) { 
    mod_str <- read.csv("Results/VitalRates_3/15_model_structures.csv")
    mod_str <- mutate(mod_str, Equation = gsub("Response",response,Equation) )
  }  

  # model ranks plus mode weights  
  mod_rank  <- do.call(rbind,strsplit(attributes(x)$row.names , "model"))
  mod_rank  <- data.frame( Model = as.numeric(mod_rank[,2]),
                           dAIC  = x$dAIC,
                           AIC_weight = x$weight )
  
  # Combine model structures and model weights
  sel_res     <- merge(mod_rank, mod_str, sort = F)
  if(candidate_num == 15){
    sel_res   <- dplyr::select(sel_res, Model, Target, Competitive.effect, 
                               Competitive.response,Equation, dAIC, AIC_weight)
  } 
  if(candidate_num == 5){
    sel_res   <- dplyr::select(sel_res, Model, Competitive.effect, Equation, dAIC, AIC_weight)
  }
  
  sel_res     <- mutate(sel_res, weight = round(AIC_weight,3))
  
  return(sel_res)
  
}

# Final model selection function - implemented with RGR on 4.26.2017
rgr_mod_sel <- function(x){
  
  form_extract <- function(x){
    
    tmp <- formula(x) %>%
      deparse() %>%
      strsplit("~") 
    # only left part of formula
    tmp <- tmp[[1]][[2]]
    tmp <- gsub("\\+ \\(1 \\| plot\\)", "", tmp)
    tmp <- trimws(tmp)
  }
  
  mod_structures  <- lapply(x, form_extract) %>%
    unlist()
  
  # model ranks plus mode weights  
  mod_sel_tab <- AICtab(x, weights = T, sort = F)
  mod_rank    <- do.call(rbind, strsplit(attributes(mod_sel_tab)$row.names , "model"))
  mod_rank    <- data.frame( Model = as.numeric(mod_rank[,2]),
                             dAIC  = mod_sel_tab$dAIC,
                             AIC_weight = mod_sel_tab$weight )
  mod_rank    <- mutate(mod_rank, Equation = mod_structures)
  
  return(mod_rank)
  
}
