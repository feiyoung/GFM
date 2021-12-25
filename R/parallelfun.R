# parallel computing for the single dynamic parameter.
single_parallel <- function(func,iterable,varlist=NULL,...){
  "
  :param func: function to be parallel computing
  :param iteralbe:a dynamic parameter(vector、list) of func.
  :param ...: all static paramters of func.
  :return list whose length is equivalent to that of iterable.
  "
  #1.load package
  library(parallel)
  #2.count the number of cores
  cores <- detectCores(logical = FALSE)
  #3.open parallel computing nodes
  cl <- makeCluster(cores)
  #4.pass objects for each node.
  funcname <- deparse(substitute(func))
  varlist <- c(funcname,varlist)
  parallel::clusterExport(cl, varlist = varlist, envir = environment())
  # Put the reqiured functions in GFM package into all nodes.
  parallel::clusterCall(cl, function() library(GFM))
  #5.start to parallel computing
  result <- parallel::parSapply(cl=cl,X=iterable,FUN=func,...)
  #6.close parallel computing
  stopCluster(cl)
  return(result)
}
