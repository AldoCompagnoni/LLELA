model_avg_format <- function(x){
  
  out     <- dplyr::select(x, predictor, avg)
  out$avg <- as.vector(out$avg)
  
  out <- dplyr::mutate(out, predictor = gsub("TotDensity","N", predictor))
  out <- dplyr::mutate(out, predictor = gsub("log_l_t0","sizet", predictor))
  out <- dplyr::mutate(out, predictor = gsub("(Intercept)","intercept", predictor))
  out <- dplyr::mutate(out, predictor = gsub("sexm","male", predictor))
  
  return(out)
  
}
