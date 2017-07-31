# create design for prediction
new_design <- function(mod_avg){
  
  # ACTUAL PARAMETERS
  params <- mod_avg$predictor
  
  # design
  des <- expand.grid("(Intercept)" = 1, 
                     TotDensity = seq(1,48,1),
                     sr = seq(0,1,0.05),
                     sexm = c(1,0))
  des <- mutate(des, "TotDensity:sr" = sr * TotDensity)
  # des <- mutate(des, "sexm:sr" = sr * sexm)
  des <- mutate(des, "sexm:TotDensity" = sexm * TotDensity)
  
  if( any(params == "sexm") ){
    return(des)
  } else {
    out <- subset(des, sexm == 0) 
    return(out)
  }
  
}


# predict 
pred <- function(design, mod_avg, pred_name, link){
  
  # give it a name
  pred_name <- deparse( substitute(pred_name) )
  link_f    <- deparse( substitute(link) )
  
  #test correspondence
  # if( all(names(design) %in% mod_avg$predictor) != TRUE ) {
  #   stop("Predictors in Design matrix do not correspond to sequence of predictors")
  # } else {
  #   id_des <- match(as.character(mod_avg$predictor), names(design))
  #   design <- design[,id_des]
  # }
  
  # Prediction
  pred_val <- design %>%
                .[c(mod_avg$predictor)] %>%
                as.matrix %*% mod_avg$avg %>%
                as.vector
  
  pred_val <- eval(parse(text = paste0("sapply(pred_val,",link_f,")")) )
  design   <- design %>%
                mutate( tmp = pred_val) %>%
                rename_( .dots = setNames("tmp", pred_name ) )
  
  # pred_val <- as.vector( as.matrix(design[,mod_avg$predictor]) %*% mod_avg$avg )
  # design   <- cbind(design,NA)
  # design[ncol(design)] <- pred_val
  # names(design)[ncol(design)] <- pred_name
  
  return(design)
  
}


