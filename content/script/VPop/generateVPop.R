#' Generate a Virtual Population from a simulation function and observed data
#'
#' @param model_fn User defined function to simulate the model and return a data.frame of outcomes. See \code{\link{model_fn}}.
#' @param NP  Number of plausible patients to generate.
#' \describe{
#' A higher number of plausible patients results in a larger virtual population. For extremely low VP yields a very large number (>1M) plausible patients may be needed to generate a sufficienty large virtual population
#' }
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
#' @param data dataframe containing \code{time} and observed values at times in \code{time}. Column names must be a subset of column names in \code{plausiblePatients}
#' @param pp_method String specifying the method for plausible patient optimization. Currently only simulated annealing (\code{"SA"}) is supported.
#' @param model_args Optional list specifying named arguments to \code{model_fn} beyond \code{parameters}
#' @param scoreThreshold Threshold for plausible patient acceptance. Must be \code{>=0}.  \code{scoreThreshold > 0} will allow for acceptance patients which violate parameter and outcome constraints to some degree in exchange for faster computation
#' @param pp_control Optional list of named arguments which controls the behavior, including parallelization of the plausible patient generation:
#'    \describe{
#'      \item{\code{runParallel}}{
#'        String. \code{"parallel"} when local parallelization (requires \code{doParallel}) is to be used. \code{"qapply"} when parallelization across computational nodes is to be used (requires \code{qapply})}
#'      \item{\code{nCores}}{Numeric. Used to specify the number of cores to parallelize across. If \code{"qapply"} is used for parallelization \code{nCores} is used as the number of codes per computation node.
#'      }
#'      \item{\code{nNodes}}{Numeric. Used only if \code{"qapply"} is used for \code{runParallel}. Specifies the number of computational nodes to use. Note: If 48 plausible patients are desired and 4 nodes with 16 cores each are specified, then 64 plausible patients will be generated.}
#'      \item{\code{source_parallel}}{Vector of strings. Used only if \code{runParallel} is specified. Provides a vector of files which are sourced in \code{model_fn} to the child threads or nodes for parallelization}
#'      \item{\code{parallel_libs}}{Vector of strings. Used only if \code{runParallel} is specified. Provides a vector of libraries which are loaded in \code{model_fn} to the child threads or nodes for parallelization}
#'    }
#' @param pp_optim_control Optional named list of argmuents to be passed to patient optimization routine (eg. \code{\link{GenSA}}).
#' @param VP_runs Number of runs to average virtual population inclusion statistic over.
#' @param plausible_timeout Numeric. Number of seconds after which explicit estimation of KDE bandwidth times out. If bandwidth estimation times out a binned KDE is used instead.
#' @param data_timeout Numeric. Number of seconds after which explicit estimation of KDE bandwidth times out. If bandwidth estimation times out a binned KDE is used instead.
#' @param alpha_algorithm String. Algorithm used to optimize the Virtual Population Size/inclusion rate into the virtual population
#' \describe{
#' \item{"SA"}{Simulated annealing function (same as used in \code{\link{generatePPs}})}
#' \item{"PSO"}{Particle swarm optimization \code{\link{hydroPSO}}. Parallel and more robust (recommended)}
#' }
#' @param vp_optim_control Named List. Arguments to the optimization algorithm. (See \code{\link{GenSA}} or \code{\link{hydroPSO}})
#'
#'@return Named list containing the parameters and simulation outputs of all the patients within the virtual population
#'\describe{
#'\item{\code{VPs}}{Simulation results for all virtual patients include within the virtual population}
#'\item{\code{VP_params}}{Parameter sets for all virtual patients included within the virtual population}
#'}
#'@export
generateVPop <- function(model_fn=model_fn, NP, paramLims, stateLims, data, pp_method="SA", model_args=list(), scoreThreshold=0.0,
                         pp_control=list(),pp_optim_control=list(), VP_runs = 20, plausible_timeout = 300, data_timeout = 300,
                         alpha_algorithm = "PSO",vp_optim_control = list()){
PlausiblePatients <- generatePPs(model_fn=model_fn, NP=NP, paramLims=paramLims, stateLims=stateLims,
                              method="SA", model_args=model_args, scoreThreshold=scoreThreshold, control=pp_control,optim_control=pp_optim_control)
VPs <- getVPs(plausiblePatients=PlausiblePatients, data=data, runs=VP_runs, plausible_pdf = "auto", data_pdf = "auto",
              plausible_timeout = plausible_timeout, data_timeout = data_timeout,
              alpha_algorithm = alpha_algorithm, optim_control = vp_optim_control)
return(VPs)

}
