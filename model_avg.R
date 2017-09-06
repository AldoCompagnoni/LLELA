model_avg <- function(model_sel,model_list, dat = NULL){
  
  # construct "design matrices"
  vars          <- lapply(model_list, function(x) names(fixef(x)) ) %>% unlist
  
  if( any(vars == "totFlow") ){
    # model design for panicles
    mod_des   <- expand.grid(totFlow = seq(1,49,1),
                             sr_f = round(seq(0,1,0.05),2) )
  }else{
    if( any(grepl("sex",vars)) ){
      mod_des <- expand.grid(TotDensity = seq(1,49,1),
                             sr = round(seq(0,1,0.05),2),
                             sex = as.factor( c("f","m") ) )
    }else{
      mod_des <- expand.grid(TotDensity = seq(1,49,1),
                             sr = round(seq(0,1,0.05),2) )
    }
  }
  
  # test that mod sel table is sorted.
  n <- seq_along(model_sel$weight)
  if(!identical(n, order(model_sel$weight,decreasing = T))){
    stop("Sort model selection table (sort=T)!")
  }
  
  # extract model rank
  mod_rank <- gsub("model", "", attributes(model_sel)$row.names) %>% as.numeric
    
  # Models that make up more than 95% of weight
  k <- 0 ; sumW <- 0 
  while(sumW < 0.95){ 
    k <- k + 1
    sumW <- sum(model_sel$weight[1:k])
  }
  
  # create prediction data frame
  pred_l  <- lapply(model_list[mod_rank[1:k]], 
                    function(x) predict(x, newdata = mod_des, re.form = NA, type = "response") )
  preds   <- ( do.call(cbind, pred_l) %*% model_sel$weight[1:k] ) / sum(model_sel$weight[1:k])
  pred_df <- mutate(mod_des, pred = preds)
  
  return(pred_df)
  
}


# function formats the growth data 
format_growth <- function(x){
  
  # remove dead individuals (this is a GROWTH model!)
  d       <- subset(x, surv_t1 != 0)
  # create log ratio: log(sizet1/sizet0)
  d       <- mutate(d, log_ratio = log(l_t1 / l_t0) )
  # glmmadmb wants plot as a factor
  d       <- mutate(d, plot = as.factor(plot) ) 
  
  # format the new_t1 column
  d       <- mutate(d, new_t1 = replace(new_t1, new_t1 == "SKIPPED", NA))
  d       <- mutate(d, new_t1 = replace(new_t1, new_t1 == "cnf", NA))
  d       <- mutate(d, new_t1 = replace(new_t1, new_t1 == "", NA))
  d       <- mutate(d, new_t1 = as.numeric(as.character(new_t1)))
  # new columns
  d       <- mutate(d, TotDensity2 = TotDensity^2,
                    new_t1_pc   = new_t1/TotDensity)
  
  return(d)
  
}


# function formats data on new tillers   
format_new_tillers <- function(x){
  
  # remove dead individuals (this is a GROWTH model!)
  d       <- subset(x, surv_t1 != 0)
  # create log ratio: log(sizet1/sizet0)
  d       <- mutate(d, log_ratio = l_t1 / l_t0)
  # glmmadmb wants plot as a factor
  d       <- mutate(d, plot = as.factor(plot) ) 
  
  # format the new_t1 column
  d       <- mutate(d, new_t1 = replace(new_t1, new_t1 == "SKIPPED", NA))
  d       <- mutate(d, new_t1 = replace(new_t1, new_t1 == "cnf", NA))
  d       <- mutate(d, new_t1 = replace(new_t1, new_t1 == "", NA))
  d       <- mutate(d, new_t1 = as.numeric(as.character(new_t1)))
  # new columns
  d       <- mutate(d, TotDensity2 = TotDensity^2,
                    new_t1_pc   = new_t1/TotDensity)
  
  # take unique values - new tillers refer to plots, not individuals 
  out     <- d %>%
    select(new_t1, new_t1_pc, plot, sr, M, F, TotDensity, TotDensity2, year) %>%
    unique() %>%
    na.omit()
  
  return(out)
  
}


# format data for one-sex plots only
one_sex_format <- function(form_gr_dat){
  
  # use object CREATED BY format_growth()
  tmp     <- form_gr_dat 
  tmp     <- select(tmp,new_t1,plot,sex,year,TotDensity)
  tmp     <- mutate(tmp1, oneSex = replace(oneSex, F==0 | M==0, 1) )
  
  # file with one sex-only plots
  one_sex <- subset(tmp, oneSex == 1)
  # flag what is male, and what is female
  one_sex <- mutate(one_sex, plot_sex = "m")
  one_sex <- mutate(one_sex, plot_sex = replace(plot_sex, M == 0,"f") )
  one_sex <- mutate(one_sex, plot_sex = factor(plot_sex, levels = c("f", "m")) ) 
  
  return(one_sex)
  
}


# format data for one-sex plots only
format_flower <- function(x){
  
  d       <- mutate(x, plot = as.factor(plot) )
  # Data from 2014; only live individuals  
  tmp     <- subset(d, year==2014 & surv_t1 != 0)
  f14     <- na.omit(dplyr::select(tmp,flowN_t1,plot,sex,sr,TotDensity))
  
  return(f14)
  
}


# format data for 3d graph
form_3d_surf <- function(d,response){
  
  # transform response in a string
  resp <- deparse( substitute(response) )
  
  # order data
  d  <- d[order(d$TotDensity,d$sr),] 
  
  # Prepare data
  x<-unique(d$sr)
  y<-unique(d$TotDensity)
  z<-matrix(d[,resp], nrow=length(unique(d$sr)),
            ncol=length(unique(d$TotDensity)))
  
  # substitute Inf values
  out <- list(x=x,y=y,z=z)
  
  for(i in seq_along(out)){
    if(any(out[[i]]==Inf,na.rm=T) == T){ out[[i]][out[[i]]==Inf]=NA }
    if(any(is.nan(out[[i]])) == T){ out[[i]][is.nan(out[[i]])]=NA }
  }
  
  return(out)
  
}
  