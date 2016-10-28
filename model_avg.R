# Function that averages model coefficients of models that make up 95% of AIC weight
model_avg = function(model_sel,model_list){ #,

  mod_list1    <- model_list
  
  betaList <- list()
  mod_rank <- do.call(rbind,strsplit(attributes(model_sel)$row.names , "model"))
  mod_rank <- as.numeric(mod_rank[,2])

  # Models that make up more than 95% of weight
  k <- 0 ; sumW <- 0 
  while(sumW < 0.95) { 
    k <- k + 1
    sumW <- sum(model_sel$weight[1:k]) 
  }

  # Store values
  for(i in 1:k){
    
    estimates     <- coef(model_list[[mod_rank[i]]])
    
    # ANNOYING: Translate "TotDensity:sexm" into "sexm:TotDensity"
    if(any(names(estimates) == "TotDensity:sexm")){
      fixI      <- grep("TotDensity:sexm", names(estimates))
      names(estimates)[fixI] <- "sexm:TotDensity"
    }
    
    betaList[[i]] <- data.frame(predictor = names(estimates),
                                 parameter = estimates)
    names(betaList[[i]])[2] <- paste0("parameter_",i)
    
  }

  # Model averages
  beta_avg                  <- Reduce(function(...) merge(...,all=T), betaList)
  beta_avg[is.na(beta_avg)] <- 0
  mod_weights               <- model_sel$weight[1:k]
  beta_avg$avg              <- as.matrix(beta_avg[,-1]) %*% mod_weights / sum(mod_weights)

  return(beta_avg)
  
}
