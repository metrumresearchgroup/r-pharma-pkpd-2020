#' Resample NP plausible patients to generate a virtual population which matches distribution (uni- or multivariate) of observed data
#'
#' @param plausiblePatients Named list of plausible patient simulations and plausible patient parameters generated from \code{\link{generatePPs}}
#' \describe{
#' \item{\code{parameters}}{dataframe with \code{NP} rows giving the values of all plausible patient parameters}
#' \item{\code{simulation_results}}{dataframe with \code{NP x nTimes} rows giving the values of plausible patient simulation results at specified times}
#' }
#' @param data dataframe containing \code{time} and observed values at times in \code{time}. Column names must be a subset of column names in \code{plausiblePatients}
#' @param runs Number of runs to average virtual population inclusion statistic over.
#' @param plausible_pdf Either \code{"auto"} or a function returning a vector probability density estimates for a given multivariate dataframe
#' \describe{
#' \item{\code{"auto"}}{Automatically computes probability density function estimates for the plausible patient simulations using a Kernel Density Estimate (see \code{\link{ks}})}
#' \item{\code{...}}{Function which must take a dataframe as the input argument and return probability densities (generated from the plausible patient simulations) for each row in the data frame}
#' }
#' @param data_pdf Either "auto" or a function returning a vector probability density estimates for a given multivariate dataframe
#' \describe{
#' \item{"auto"}{Automatically computes probability density function estimates for observed data using a Kernel Density Estimate (see \code{\link{ks}})}
#' \item{\code{...}}{Function which must take a dataframe as the input argument and return probability densities (generated from the observed data) for each row in the data frame}
#' }
#' @param plausible_timeout Numeric. Only used for automatic plausible patient kernel density estimation (\code{\link{ks}}). Number of seconds after which explicit
#' estimation of KDE bandwidth times out. If bandwidth estimation times out a binned KDE is used instead.
#' @param data_timeout Numeric. Only used for automatic data kernel density estimation (\code{\link{ks}}). Number of seconds after which explicit
#' estimation of KDE bandwidth times out. If bandwidth estimation times out a binned KDE is used instead.
#' @param alpha_algorithm String. Algorithm used to optimize the Virtual Population Size/inclusion rate into the virtual population
#' \describe{
#' \item{"SA"}{Simulated annealing function (same as used in \code{\link{generatePPs}})}
#' \item{"PSO"}{Particle swarm optimization \code{\link{hydroPSO}}. Parallel and more robust (recommended)}
#' }
#' @param optim_control Named List. Arguments to the optimization algorithm. (See \code{\link{GenSA}} or \code{\link{hydroPSO}})
#'
#'@return Named list containing the parameters and simulation outputs of all the patients within the virtual population
#'\describe{
#'\item{\code{VPs}}{Simulation results for all virtual patients include within the virtual population}
#'\item{\code{VP_params}}{Parameter sets for all virtual patients included within the virtual population}
#'}
#'@export
getVPs <- function(plausiblePatients, data, runs=20, plausible_pdf = "auto", data_pdf = "auto", plausible_timeout = 300, data_timeout = 300,
                   alpha_algorithm = "PSO",optim_control = list()){

  # Restructure Data
  data <- mrgsolve::lctran(data)
  dnames <- names(dplyr::select(data,-ID,-time))
  pp_sims <- plausiblePatients$simulation_results
  pp_sims <- dplyr::select(pp_sims,c(dnames,"time"))
  mdl_df <- data.frame()
  data_df <- data.frame()
  for(t in unique(data$time)){
    for(nm in unique(dnames)){
      tmp_df <- data.frame((dplyr::filter(pp_sims,time==t))[[nm]])
      names(tmp_df) <- paste0(nm,"_",t)
      mdl_df <- rbind(mdl_df,tmp_df)
      tmp_df_dat <- data.frame((dplyr::filter(data,time==t))[[nm]])
      names(tmp_df_dat) <- paste0(nm,"_",t)
      data_df <- rbind(data_df,tmp_df_dat)
    }
  }
  if(plausible_pdf == "auto"){
    if(length(dnames)>1){
      Hpp <- R.utils::withTimeout(expr = ks::Hpi(x=as.vector(mdl_df)), timeout = plausible_timeout)
    }else{
      Hpp <- R.utils::withTimeout(expr = ks::hpi(x=as.matrix(mdl_df)), timeout = plausible_timeout)
    }
    if(is.null(Hpp)){
      fhat_pp <- ks::kde(x=as.matrix(mdl_df),binned=TRUE)
    }else{
      fhat_pp <- ks::kde(x=as.matrix(mdl_df),H=Hpp)
      }
    pp_pdf <- function(data_in){
      return(stats::predict(fhat_pp,x=data_in))
      }
    }else{
      pp_pdf <- plausible_pdf
    }

  if(data_pdf == "auto"){
    if(length(dnames)>1){
      Hdat <- R.utils::withTimeout(expr=ks::Hpi(x=as.matrix(data_df)),timeout=data_timeout)
    }else{
      Hdat <- R.utils::withTimeout(expr=ks::hpi(x=as.matrix(data_df)),timeout=data_timeout)
    }
    if(is.null(Hdat)){
      fhat_dat <- ks::kde(x=as.matrix(data_df),binned=TRUE)
    }else{
      fhat_dat <- ks::kde(x=as.matrix(data_df),H=Hdat)
    }
    data_pdf <- function(data_in){
      return(stats::predict(fhat_dat,x=data_in))
    }
  }else{
    data_pdf <- dat_pdf
  }




  if(alpha_algorithm == "SA"){
    if(length(optim_control)==0){
      optim_control = list(threshold.stop=1e-15,maxit=150,temperature=0.5,verbose=TRUE)
    }
    sf <- getPrevalence(mdl_df=plausiblePatients, mdl_pdf=pp_pdf, data_df = data, data_pdf = data_pdf, runs=runs,optim_control)
  }else if(alpha_algorithm == "PSO"){
    if(length(optim_control)==0){
      optim_control <- list(abstol=1e-8,parallel="parallel",par.nnodes=parallel::detectCores(),verbose=TRUE)
    }
    if(!("ncores" %in% names(optim_control))){
      optim_control[["par.nnodes"]] <- parallel::detectCores()
    }
    sf <- getPrevalencePSO(mdl_df=mdl_df, mdl_pdf=pp_pdf, data_df=data_df, data_pdf=data_pdf,runs=runs,optim_control = optim_control)
  }
  p_include <- sf$p_include
  sf <- sf$sf
  NP <- NROW(p_include)
  hist_score <- c()
  r = stats::runif(NP)
  select <- r < (p_include * sf)
  mdl_output <- plausiblePatients$simulation_results
  VPs <- split(mdl_output, mdl_output$ID)
  VPs <- VPs[select]
  VPs <- do.call(rbind,VPs)
  VP_params <- (plausiblePatients$parameters)[select,]
  return(list(VPs=VPs, VP_params = VP_params))
}
