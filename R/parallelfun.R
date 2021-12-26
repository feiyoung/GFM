single_parallel <- function(func,iterable,varlist=NULL,...){


  # library(parallel)
  cores <- parallel::detectCores(logical = FALSE)
  cl <- parallel::makeCluster(cores)
  funcname <- deparse(substitute(func))
  varlist <- c(funcname,varlist)
  parallel::clusterExport(cl, varlist = varlist, envir = environment())
  parallel::clusterCall(cl, function() library(GFM))
  result <- parallel::parSapply(cl=cl,X=iterable,FUN=func,...)
  parallel::stopCluster(cl)
  return(result)
}
