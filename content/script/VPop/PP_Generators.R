genPP_SA <- function(NP, model_fn, paramLims, stateLims, steadystateLims, model_args, scoreThreshold,optim_control=list()){
  df_list <- vector(mode = "list", length = NP)

  # print(df_list)
  i <- 0
  k <- 0

  RVS_f <- purrr::partial(runVariantAndScore, model_fn=model_fn, stateLims=stateLims,
                          steadystateLims=steadystateLims, model_args=model_args)
  out_df <- data.frame()
  stateLims <- split(stateLims,stateLims$Name)
  steadystateLims <- split(steadystateLims,steadystateLims$Name)
  while((i < NP) & (k < NP*100)){
    out_tmp <- list("value"=1e256)
    result <- tryCatch({
      parameters <- paramLims$Lower + (paramLims$Upper - paramLims$Lower) * stats::runif(NROW(paramLims),0,1)
      out_tmp <- GenSA::GenSA(par = parameters, fn = RVS_f, lower=paramLims$Lower,upper=paramLims$Upper,control=optim_control)

      # print(out_tmp$count)
      k <- k + 1
      # print(k)
    },error=function(e){
        print(e)
        # k <- k+1
        out_tmp$value <- 1e256
    },warning=function(w){
      print(w)
    })
    if(out_tmp$value <= scoreThreshold){
      i <- i + 1
      df_tmp <- data.frame(ID=i)
      for(ip in seq_along(paramLims$Name)){
        df_tmp[as.character(paramLims$Name[ip])] <- out_tmp$par[ip]
      }
      df_tmp["Score"] <- out_tmp["value"]
      simresults <- do.call(model_fn,lappend(model_args,"parameters"=out_tmp$par,"simulate"=1))
      out_df <- rbind(out_df,df_tmp)
      simresults$ID=i
      df_list[[i]] = as.data.frame(simresults)
    }
  }
  simresults <- suppressWarnings(dplyr::bind_rows(df_list))
  return(list(parameters=out_df,simulation_results=simresults))
}
