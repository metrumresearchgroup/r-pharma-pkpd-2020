getPrevalence <- function(mdl_df, mdl_pdf, data_df, data_pdf, runs=10,optim_control=list()){
  p_include <- data_pdf(mdl_df) / mdl_pdf(mdl_df) # I NEED TO CHECK THIS!
  hist_F <- purrr::partial(prevalenceFun, p_include=p_include, mdl_df=mdl_df, runs=runs, data_df=data_df)
  sf_max <-  1/max(p_include)
  sf_lower <- 0
  sf_upper <- 1000*sf_max
  # print(sf_upper)
  k= GenSA::GenSA(par=sf_max,fn=hist_F,lower=sf_lower,upper=sf_upper,control=optim_control)
  sf <- k$par
  # score <- k$value
  return(list(sf=sf,p_include=p_include))
}


prevalenceFun <- function(sf, p_include, mdl_df, runs, data_df=data_df){
  # print(sf)
  NP <- NROW(mdl_df)
  hist_score <- c()
  for(i in 1:runs){
    r = stats::runif(NP)
    # print(r)
    select <- r < (p_include * sf)
    num_vp = sum(select)
    # print(num_vp)
    # print(num_vp)
    if(num_vp > 1){
      # print('here')
      ksstat <- c()
      pval <- c()
      for(j in 1:NCOL(mdl_df)){
        ksresults <- suppressWarnings(stats::ks.test(mdl_df[select,j], data_df[,j]))
        ksstat <- append(ksstat,ksresults$statistic)
        pval <- append(pval,ksresults$p.value)
      }
      hist_score <- append(hist_score,sum(ksstat))
    }else{
      hist_score <- append(hist_score,NCOL(mdl_df))
    }
  }
  return(sum(hist_score)/runs)
}
