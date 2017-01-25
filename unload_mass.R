unload_mass <- function(){
  
  unload_pckg <- c("lme4","glmmADMB","MASS")
  if( any(unload_pckg %in% loadedNamespaces()) ){
  
    r <- which(unload_pckg %in% loadedNamespaces()) 
    lapply(unload_pckg[r],unloadNamespace)
    return(NULL)
  } else{
  }
}
