runVariantAndScore <- function(parameters, model_fn, stateLims, steadystateLims, model_args){
    mdl_output <- do.call(model_fn,lappend(model_args,"parameters"=parameters))
    # print(mdl_output)
    score <- scorePP(mdl_output=mdl_output,stateLims=stateLims,steadystateLims=steadystateLims)
    # print(score)
    # score <- 0
    return(score)
}
