
scorePP <- function(mdl_output,stateLims, steadystateLims, SSatol, SSrtol){
  SCORE <- 0.0
  SSout <- data.frame()
  NSSout <- data.frame()
  if("SS" %in% names(mdl_output)){
    SSout <- as.data.frame(mrgsolve::lctran(mdl_output$SS))
    row.names(SSout) <- SSout$time
  }
  if("NSS" %in% names(mdl_output)){
    NSSout <- as.data.frame(mrgsolve::lctran(mdl_output$NSS))
    row.names(NSSout) <- as.numeric(NSSout$time)
  }
  for(state in names(stateLims)){
    SLtmp <- stateLims[state][[1]]
    mdl_state <- NSSout[state]
    mdl_state <- mdl_state[as.character(SLtmp$Time),]
    sMP <- (SLtmp$Lower + SLtmp$Upper)/2
    P1 = (sMP - mdl_state)^2
    P2 = (sMP - SLtmp$Lower)^2
    if(is.infinite(P1)){
      P1 = sign(P1) * .Machine$double.xmax
    }
    if(is.infinite(P2)){
      P2 = sign(P2) * .Machine$double.xmax
    }
    stmp <- P1 - P2
    stmp <- stmp / sMP^2
    stmp <- stmp[stmp>=0]
    SCORE <- SCORE + sum(stmp,na.rm=FALSE)
  }
  for(state in names(steadystateLims)){
    mdl_state <- SSout[state][[1]]
    SLtmp <- steadystateLims[state][[1]]
    # mdl_state <- mdl_state[length(mdl_state)]
    sMP <- (SLtmp$Lower + SLtmp$Upper)/2
    stmp <- (sMP - mdl_state)^2 - (sMP - SLtmp$Lower)^2
    stmp <- stmp / sMP^2
    stmp <- stmp[stmp>=0]
    SCORE <- SCORE + sum(stmp, na.rm=FALSE)
  }
  if(is.na(SCORE)){
    SCORE <- 1e256
  }
  return(SCORE)
}
