# create design for prediction
new_design <- function(mod_avg){
  
  # ACTUAL PARAMETERS
  params <- mod_avg$predictor
  
  # design
  des <- expand.grid("(Intercept)" = 1, 
                     TotDensity = seq(1,48,1),
                     sr = seq(0,1,0.01),
                     sexm = c(1,0))
  des <- mutate(des, "sr:TotDensity" = sr * TotDensity)
  des <- mutate(des, "sexm:sr" = sr * sexm)
  des <- mutate(des, "sexm:TotDensity" = sexm * TotDensity)
  
  if( any(params == "sexm") ){
    return(des[,params])
  } else {
    out <- subset(des, sexm == 0) 
    return(out[,params])
  }
  
}
  
# predict 
pred <- function(design, mod_avg, pred_name){

  # give it a name
  pred_name <- deparse( substitute(pred_name) )
  #test correspondence
  expect_equal( all(names(design) == mod_avg$predictor), TRUE )
  
  # Prediction
  pred_val <- as.vector( as.matrix(design) %*% mod_avg$avg )
  design   <- cbind(design,NA)
  design[ncol(design)] <- pred_val
  names(design)[ncol(design)] <- pred_name

  return(design)
    
}

  
