#' Generate a Set of NP Plausible Patients
#' @importFrom dplyr %>%
#' @importFrom dplyr bind_rows
#' @param model_fn User defined function to simulate the model and return a data.frame of outcomes. See \code{\link{model_fn}}.
#' @param NP  Number of plausible patients to generate.
#' @param paramLims data.frame specifying \code{upper} and \code{lower} bounds for plausible patient parameters.
#'  The dataframe must contain the following variables:
#' \describe{
#'   \item{\code{Name}}{Parameter name}
#'   \item{\code{Lower}}{Parameter lower bound}
#'   \item{\code{Upper}}{Parameter upper bound}
#' }
#' \code{Lower} and \code{Upper} must be numeric between \code{+Inf} and \code{-Inf}. Missing \code{Lower} values will be replaced with \code{-Inf}
#' and missing \code{Upper} values will be replaced with \code{Inf}.
#' @param stateLims data.frame specifying \code{upper} and \code{lower} bounds for plausible patient simulations at specified \code{time}. \code{time} may also be specified as the string \code{"SS"} to take final simulation values at steady-state which may not be at a consistent \code{time}.
#' \describe{
#'   \item{\code{Name}}{State/Output name}
#'   \item{\code{Lower}}{Upper bound of \code{Name}}
#'   \item{\code{Upper}}{Lower bound of \code{Name}}
#'   \item{\code{Time}}{Time at which \code{Upper} or \code{Lower} bound of \code{Name} is applied}
#' }
#' @param method String specifying the method for plausible patient optimization. Currently only simulated annealing (\code{"SA"}) is supported.
#' @param model_args Optional list specifying named arguments to \code{model_fn} beyond \code{parameters}
#' @param scoreThreshold Threshold for plausible patient acceptance. Must be \code{>=0}.  \code{scoreThreshold > 0} will allow for acceptance patients which violate parameter and outcome constraints to some degree in exchange for faster computation
#' @param control Optional list of named arguments which controls the behavior, including parallelization of the plausible patient generation:
#'    \describe{
#'      \item{\code{runParallel}}{
#'        String. \code{"parallel"} when local parallelization (requires \code{doParallel}) is to be used. \code{"qapply"} when parallelization across computational nodes is to be used (requires \code{qapply})}
#'      \item{\code{nCores}}{Numeric. Used to specify the number of cores to parallelize across. If \code{"qapply"} is used for parallelization \code{nCores} is used as the number of codes per computation node.
#'      }
#'      \item{\code{nNodes}}{Numeric. Used only if \code{"qapply"} is used for \code{runParallel}. Specifies the number of computational nodes to use. Note: If 48 plausible patients are desired and 4 nodes with 16 cores each are specified, then 64 plausible patients will be generated.}
#'      \item{\code{source_parallel}}{Vector of strings. Used only if \code{runParallel} is specified. Provides a vector of files which are sourced in \code{model_fn} to the child threads or nodes for parallelization}
#'      \item{\code{parallel_libs}}{Vector of strings. Used only if \code{runParallel} is specified. Provides a vector of libraries which are loaded in \code{model_fn} to the child threads or nodes for parallelization}
#'    }
#' @param optim_control Optional named list of argmuents to be passed to patient optimization routine (eg. \code{\link{GenSA}}).
#' @return Named list containing the parameters and outputs of a plausible patient:
#' \describe{
#' \item{\code{parameters}}{dataframe with \code{NP} rows giving the values of all plausible patient parameters}
#' \item{\code{simulation_results}}{dataframe with \code{NP x nTimes} rows giving the values of plausible patient simulation results at specified times}
#' }
#' @export
# utils::suppressForeignCheck(c("time","job","ID","iterator_unique","Time","Name","Lower","Upper"))
generatePPs <- function(model_fn=model_fn, NP, paramLims, stateLims, method="SA", model_args=list(), scoreThreshold=0.0, control=list(),optim_control=list()){
  if(!("runParallel" %in% names(control))){
    control[["runParallel"]] = FALSE
  }
  if((control[["runParallel"]] == "qapply")){
    if (!requireNamespace("qapply", quietly=TRUE)){
      stop("Package qapply needed for parallel argument qapply", call.=FALSE)
    }else if(!("nNodes" %in% names(control))){
      res_info <- .getResources()
      nNodes <- length(res_info$ncores_per_node)
      control[["nNodes"]] <- nNodes
      if(is.null(control[["nNodes"]])){
        stop("No cluster detected!")
      }
    }
    if(!("nCoresPerNode" %in% names(control))){
      res_info <- .getResources()
      nCores <- min(res_info$ncores_per_node)
      control[["nCoresPerNode"]] <- nCores
    }
    if(!("runTag" %in% names(control))){
      control["runTag"] <- "model"
    }
    for(arg in model_args){
      if(is.mrgmod(arg)){
        if(startsWith(dirname(soloc(arg)),tempdir())){
          stop("When using mrgsolve model with qapply soloc must be defined at compilation time. tempdir may not be used.")
        }
      }
    }
  }
  if(any(paramLims$Upper<=paramLims$Lower)){
    stop("Parameter lower bounds must be less than upper bounds")
  }
  if(any(stateLims$Upper<=stateLims$Lower)){
    stop("State lower bounds must be less than upper bounds")
  }

  if((control[["runParallel"]] == "parallel")){
    if(!requireNamespace("parallel",quietly=TRUE)){
      stop("Package parallel needed for parallel argument parallel", call.=FALSE)
    }else if(!requireNamespace("doParallel",quietly = TRUE)){
      stop("Package doParallel needed for parallel argument parallel", call.=FALSE)
    }else if(!("nCores" %in% names(control))){
      control[["nCores"]] <- parallel::detectCores(logical=TRUE)
    }else if(NP<=parallel::detectCores(logical=TRUE)){
      control[["nCores"]] <- NP
    }
  }
  if((control[["runParallel"]] == FALSE)){
    if(length(optim_control)==0){
      optim_control <- list(threshold.stop=0.005,maxit=150,temperature=0.5)
    }
    if(!("theshold.stop" %in% names(optim_control))){
      optim_control["threshold.stop"] <- 0.005
    }
    if(!("maxit" %in% names(optim_control))){
      optim_control["maxit"] <- 150
    }
    if(!("temperature" %in% names(optim_control))){
      optim_control["temperature"] <- 0.5
    }
    out <- generatePPs_internal(model_fn=model_fn,NP=NP,paramLims = paramLims, stateLims = stateLims, method=method, scoreThreshold = scoreThreshold, model_args=model_args,optim_control=optim_control)
    out$simulation_results <- dplyr::filter(out$simulation_results, time %in% stateLims$Time)
  }else if((control[["runParallel"]] == "qapply")){
    print("RUNNING QAPPLY!")
    if(length(optim_control)==0){
      optim_control <- list(threshold.stop=0.005,maxit=150,temperature=0.5)
    }
    if(!("theshold.stop" %in% names(optim_control))){
      optim_control["threshold.stop"] <- 0.005
    }
    if(!("maxit" %in% names(optim_control))){
      optim_control["maxit"] <- 150
    }
    if(!("temperature" %in% names(optim_control))){
      optim_control["temperature"] <- 0.5
    }
    if(!("clearWd" %in% names(control))){
      control["clearWd"] <- TRUE
    }
    fa <- list(model_fn = model_fn, NP = floor(max(NP/(control[["nNodes"]]*control[["nCoresPerNode"]]),1.0)),
               paramLims = paramLims, stateLims = stateLims, method = method,
               scoreThreshold=scoreThreshold, model_args = model_args,optim_control=optim_control)
    coresTotal <- control[["nNodes"]]*control[["nCoresPerNode"]]
    nPerPacket = rep(1,coresTotal)
    remainder <- NP-coresTotal*floor(max(NP/(control[["nNodes"]]*control[["nCoresPerNode"]]),1.0))
    if(remainder > 0){
      nPerPacket[1:remainder] <- nPerPacket[1:remainder] + 1
    }
    parfun <- function(i, ...){
      if("source_parallel" %in% names(control)){
        for(filename in control[["source_parallel"]]){
          source(filename)
        }
      }
      out <- generatePPs_internal(...)
      out$parameters$job <- i
      out$simulation_results <- dplyr::filter(out$simulation_results,time %in% stateLims$Time)
      out$simulation_results$job <- i
      out
    }
    # model=model_args$model
    out <- qapply::qapply(X=seq_len(floor(NP/floor(max(NP/(control[["nNodes"]]*control[["nCoresPerNode"]]),1.0)))),FUN=parfun,fargs=fa,global=FALSE,tag=control[["runTag"]],
                          clearWd = control[["clearWd"]],
                  commonData = lappend(c(list(fa=fa,model_fn=model_fn, NP = floor(max(NP/(control[["nNodes"]]*control[["nCoresPerNode"]]),1.0)),
                                    paramLims=paramLims,stateLims=stateLims,method=method,control=control,
                                    scoreThreshold=scoreThreshold,model_args=model_args),model_args)))
    sims <- lapply(out, function(x) {
      x$simulation_results
    })%>% bind_rows
    pars <- lapply(out, function(x){
      x$parameters
    }) %>% bind_rows

    out <- list(parameters=pars, simulation_results=sims)
    # print(out)
    numTimes <- out$simulation_results%>%dplyr::filter(job==1 & ID==1)
    numTimes <- length(unique(numTimes$time))
    ID_tmp <- 1:(NROW(out$simulation_results)/numTimes)
    ID_tmp <- rep(ID_tmp,numTimes)
    ID_tmp <- sort(ID_tmp)
    out$parameters <- out$parameters %>% dplyr::select(-job)
    out$parameters$ID <- 1:(NROW(out$simulation_results)/numTimes)
    out$simulation_results <- out$simulation_results%>%dplyr::select(-job)
    out$simulation_results$ID <- ID_tmp
    IDlast <- max(sort(out$parameters$ID))
    NP_par <- NROW(out$parameters)
    if((NP-NP_par)>0){
    out_add <- generatePPs_internal(model_fn=model_fn,NP=NP-NP_par,paramLims = paramLims, stateLims = stateLims, method=method, scoreThreshold = scoreThreshold, model_args=model_args,optim_control=optim_control)
    out_add$simulation_results <- dplyr::filter(out_add$simulation_results, time %in% stateLims$Time)
    out_add$simulation_results$ID <- out_add$simulation_results$ID + IDlast
    out_add$parameters$ID <- out_add$parameters$ID + IDlast
    out$parameters <- rbind(out$parameters,out_add$parameters)
    out$simulation_results <- rbind(out$simulation_results,out_add$simulation_results)
    # model_path <- dirname(soloc(mod))
    }

  }else if((control[["runParallel"]] == "parallel")){
    print("RUNNING IN PARALLEL")
    if(length(optim_control)==0){
      optim_control <- list(threshold.stop=0.005,maxit=150,temperature=0.5)
    }
    if(!("theshold.stop" %in% names(optim_control))){
      optim_control["threshold.stop"] <- 0.005
    }
    if(!("maxit" %in% names(optim_control))){
      optim_control["maxit"] <- 150
    }
    if(!("temperature" %in% names(optim_control))){
      optim_control["temperature"] <- 0.5
    }
    fa <- list(model_fn = model_fn, NP = floor(NP/control[["nCores"]]),
               paramLims = paramLims, stateLims = stateLims, method = method,
               scoreThreshold=scoreThreshold, model_args = model_args,optim_control=optim_control)
    cl <- parallel::makeCluster(control[["nCores"]])
    doParallel::registerDoParallel(cl)
    # for(item in names(model_args)){
    #   parallel::clusterExport(cl,c(item),envir=environment())
    # }
    # print("FOO!")
    lPath <- .libPaths()
    parallel::clusterCall(cl, function(){
      if("libPaths" %in% names(control)){
        .libPaths(control[["libPaths"]])
      }else{
        .libPaths(lPath)
      }
      if("parallel_libs" %in% names(control)){
        for(libname in control[["parallel_libs"]]){
          library(libname,character.only=TRUE)
        }
      }
      if("source_parallel" %in% names(control)){
        for(filename in control[["source_parallel"]]){
          source(filename)
        }
      }

    })
    nJobs <- as.integer(control[["nCores"]])
    # nJobs <- nJobs + as.integer((NP-(floor(NP/control[["nCores"]])*control[["nCores"]]))/(floor(NP/control[["nCores"]])))
    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`
    if("parallel_libs" %in% names(control)){
      extra_pkgs <- control[["parallel_libs"]]
    }else{
      extra_pkgs <- c()
    }
    # tmp <- do.call(generatePPs_internal,fa)
    out <-  foreach::foreach(iterator_unique=1:nJobs,.packages = c("purrr","dplyr",extra_pkgs,"magrittr"))  %dopar% {
          tmp <- do.call(generatePPs_internal,fa)
          tmp$parameters$job <- iterator_unique
          tmp$simulation_results <- dplyr::filter(tmp$simulation_results, time %in% stateLims$Time)
          tmp$simulation_results$job <- iterator_unique
          out <- tmp
    }
    print(out)
    parallel::stopCluster(cl)
    sims <- lapply(out, function(x) {
      x$simulation_results
    })%>%bind_rows
    pars <- lapply(out, function(x){
      x$parameters
    }) %>%bind_rows
    out <- list(parameters=pars, simulation_results=sims)
    numTimes <- out$simulation_results%>%dplyr::filter(job==1 & ID==1)
    numTimes <- length(unique(numTimes$time))
    ID_tmp <- 1:(NROW(out$simulation_results)/numTimes)
    ID_tmp <- rep(ID_tmp,numTimes)
    ID_tmp <- sort(ID_tmp)
    out$parameters <- out$parameters %>%dplyr::select(-job)
    out$parameters$ID <- 1:(NROW(out$simulation_results)/numTimes)
    out$simulation_results <- out$simulation_results%>%dplyr::select(-job)
    out$simulation_results$ID <- ID_tmp
    IDlast <- max(sort(out$parameters$ID))
    NP_par <- NROW(out$parameters)
    if((NP-NP_par)>0){
      out_add <- generatePPs_internal(model_fn=model_fn,NP=NP-NP_par,paramLims = paramLims, stateLims = stateLims, method=method, scoreThreshold = scoreThreshold, model_args=model_args,optim_control=optim_control)
      out_add$simulation_results <- dplyr::filter(out_add$simulation_results, time %in% stateLims$Time)
      out_add$simulation_results$ID <- out_add$simulation_results$ID + IDlast
      out_add$parameters$ID <- out_add$parameters$ID + IDlast
      out$parameters <- rbind(out$parameters,out_add$parameters)
      out$simulation_results <- rbind(out$simulation_results,out_add$simulation_results)
    }

  }
  closeAllConnections()
  return(out)
}
