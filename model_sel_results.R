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
                           AIC_weight = x$weight )
  
  # Combine model structures and model weights
  sel_res     <- merge(mod_rank, mod_str, sort = F)
  if(candidate_num == 15){
    sel_res   <- dplyr::select(sel_res, Model, Target, Effect, Equation, AIC_weight)
  } 
  if(candidate_num == 5){
    sel_res   <- dplyr::select(sel_res, Model, Effect, Equation, AIC_weight)
  }
  
  sel_res     <- mutate(sel_res, weight = round(AIC_weight,5))
  
  return(sel_res)
  
}
