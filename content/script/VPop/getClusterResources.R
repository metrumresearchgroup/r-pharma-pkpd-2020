#' @name getResources
#' @title Get SGE resources
#' @details Find the number of (potentially) available slots to push to
#' @import stringr
.getResources <- function(){
  x <- system(" qstat -f", intern = T)
  if(length(x)==0) return(NULL)
  loc <- stringr::str_locate(x[1], "resv/used/tot.")
  y <- stringr::str_sub(x[-1], loc)
  y <- stringr::str_extract(stringr::str_trim(y[stringr::str_detect(y,"\\d/\\d/\\d")]),"\\d+$")
  ncores <- (sum(as.numeric(y)))
  nnodes <- length(y)
  return(list(nnodes=nnodes,ncores_per_node  = as.numeric(y),ncores=ncores))
}
