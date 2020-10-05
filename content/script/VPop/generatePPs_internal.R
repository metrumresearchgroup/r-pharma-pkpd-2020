generatePPs_internal <- function(model_fn=model_fn, initGuess, NP, paramLims, stateLims, method="SA", scoreThreshold=0.0, model_args=list(),optim_control=list()){
    # paramLims and stateLims are both dataframes with a column for the parameter/state name, upper and lower bounds
  paramLims <- as.data.frame(paramLims)
  stateLims <- as.data.frame(stateLims)
  ssFLAG = 0
  nssFLAG = 0
  nssTimes = c()
  if(!("Time" %in% names(stateLims))){
    stateLims$Time <-  "SS"
    ssFLAG <-  1
  }else if("SS" %in% stateLims$Time){
    # print('here')
    ssFLAG <-  1
  }
  if(NROW((stateLims %>% dplyr::filter(Time!="SS")))>0){
    nssFLAG <-  1
    nssTimes <- as.numeric(as.character((stateLims%>%dplyr::filter(Time!="SS"))$Time))
    nssTimes <- unique(nssTimes)
  }

  if(method!="SA"){
    stop("Only Simulated Annealing supported at the moment")
  }
  for(state in unique(stateLims$Name)){
    UNIQUEDF <- stateLims %>% dplyr::filter(Name==state)
    if(NROW(unique(as.character(UNIQUEDF$Time)))<NROW(UNIQUEDF)){
      stop(paste0(state," has simultaneous values"))
    }
  }


  # for(name in unique(paramLims$Name)){
  #   if(NROW(paramLims %>% dplyr::filter(Name==name))>1){
  #     stop(paste0("Parameter ",name," has multiply defined limits!"))
  #   }else if((!name %in% names(param(model)))){
  #     stop(paste0(name," is not a model parameter!"))
  #   }
  # }

  row.names(paramLims) <- paramLims$Name
  for(p in paramLims$Name){
    if(is.na(paramLims[p,]$Lower)){
      paramLims[p,]$Lower <- -Inf
    }
    if(is.na(paramLims[p,]$Upper)){
      paramLims[p,]$Upper <- Inf
    }
  }

  stateLims <- transform(stateLims, Lower = ifelse(is.na(Lower), -Inf, Lower))
  stateLims <- transform(stateLims, Upper = ifelse(is.na(Upper), -Inf, Upper))

  row.names(paramLims) <- NULL
  # Get Time varying and Steady-State Limits
  steadystateLims <- stateLims %>% dplyr::filter(Time=="SS")
  steadystateLims$Time <- as.character(steadystateLims$Time)
  stateLims <- stateLims %>% dplyr::filter(Time!="SS")
  stateLims$Time <- as.numeric(as.character(stateLims$Time))
  # GENERATE TEST OUTPUT
  test_output <- do.call(model_fn,model_args)
  if(class(test_output)!="list"){
    stop("Model output must be a list containing some combination of steady state (SS) and non steady state (NSS) values")
  }

  # Check for SS output
  if(!("SS" %in% names(test_output)) & ssFLAG==1){
    stop("Steady steady state outputs must be provided for steady state limits!")
  }else if(ssFLAG==1){
    for(state in steadystateLims$Name){
      if(!(state %in% names(test_output$SS))){
        stop(paste0("State ", state, " must be in model steady state output!"))
      }
    }
  }
  if(nssFLAG==1){
    if(!("NSS" %in% names(test_output))){
      stop("Model output must contain non-steady state values!")
    }else{
      test_output$NSS <- mrgsolve::lctran(test_output$NSS)
      for(state in stateLims$Name){
        if(!(state %in% names(test_output$NSS))){
          stop(paste0("State ", state, " must be present in non-steady-state output"))
        }
      }
      for(t in unique(stateLims$Time)){
        if(!(t %in% test_output$NSS$time)){
          stop(paste0("Time ",t," found in state limits but not in non-steady-state model outputs"))
        }
        tunique <- unique(test_output$NSS$time)
        if(length(tunique)!=NROW(test_output$NSS)){
          stop("Duplicated times found in model output! Please ensure model output times are unique!")
        }
      }
    }
  }
 # Drop unused levels
  stateLims <- droplevels(stateLims)
  steadystateLims <- droplevels(steadystateLims)
  if(method=="SA"){
    gen_pp_f <- purrr::partial(genPP_SA, model_fn=model_fn,paramLims=paramLims,stateLims=stateLims,steadystateLims=steadystateLims, model_args=model_args, scoreThreshold=scoreThreshold,optim_control=optim_control)
  }
  out_df <- gen_pp_f(NP)
  return(out_df)
}
