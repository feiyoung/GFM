gfm <- function(X, group, type, q=NULL, parallel=TRUE,para.type='doSNOW', ncores = 10,
                dropout=0, dc_eps=1e-4, maxIter=50,
                q_set=1:10, output=TRUE, fast_version=FALSE){
  if(!is.matrix(X)) stop("X must be a matrix.")

  n <- nrow(X); p <- ncol(X)
  omega <- 1/p

  if(is.null(q)){
    if(parallel && para.type=='doSNOW'){
      #require(doSNOW)

      if(ncores >  parallel::detectCores() ) ncores <- parallel::detectCores()
      cl <- makeSOCKcluster(ncores) # 设定并行核
      doSNOW::registerDoSNOW(cl) # 注册该核
      nq <- length(q_set)
      pb <- txtProgressBar(min=1, max=nq, style=3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress=progress)
      r <- 1
      resq <- foreach(r = 1:nq,.packages="GFM" ,.options.snow=opts,
                        .combine='rbind') %dopar% {
                          res <- paraIC(r, X, group, type,
                                                        dropout=dropout, eps2=1e-4, maxIter=10, output=output, fast_version=TRUE)
                          res
                        }
      close(pb)
      parallel::stopCluster(cl)

    }else if(parallel && para.type=='parallel'){
      # varlist <- c('gfm_eval_intercept_init', 'ICriteria', "Factorm",
      #              "localupdateB2","family2func","ortheB","Diag",
      #              "localupdateH2","ortheH", "objfunc","signrevise")
      varlist <- NULL
      resq <- single_parallel(paraIC,iterable =q_set, varlist= varlist, XX=X, group=group, type=type,
                              dropout=dropout, eps2=dc_eps, maxIter=10,  output=output, fast_version=fast_version)
      q <- q_set[which.min(resq)]
    }else{
      n <- nrow(X)
      q_num <- length(q_set)
      allhH <- list(); allhB <- list()
      for(r in 1:q_num){

        gfm1 <- gfm_eval_intercept_init(X, group, type, r,
                                        dropout, eps2=1e-4, maxIter, output)
        if(!fast_version){
          # hH <- gfm1$hH; hB <- gfm1$hB; hmu <- gfm1$hmu
          gfm1 <-gfm_eval_intercept_osfinal(X, gfm1$hH, gfm1$hB,gfm1$hmu, group, type)
        }
        allhH[[r]] <- cbind(1, gfm1$hH)
        allhB[[r]] <- cbind(gfm1$hmu, gfm1$hB)
      }
      Vr <- matrix(0, 2, q_num)
      for(r in 1:q_num){
        Vr[,r] <- ICriteria(X, allhB[[r]], allhH[[r]], r, group, type)

      }
      resq <- apply(Vr, 2, sum)
    }
    q <- q_set[which.min(resq)]
    message('\nThe factor number q is estimated as ', q, '. \n')
  }
  message('Starting the alternate minimization algorithm...\n')
  gfm2 <- gfm_eval_intercept_init(X, group, type, q,
                                  dropout, dc_eps, maxIter,
                                  output, parallel = parallel)
  gfm2$q <- q
  message('Finish the iterative algorithm...\n')
  if(fast_version){
    class(gfm2) <- 'gfm'
    return(gfm2)
  }else{
    gfm2final <- gfm_eval_intercept_osfinal(X, gfm2$hH, gfm2$hB,gfm2$hmu, group, type)
    gfm2final$q <- q
    gfm2final$history <- gfm2$history
    class(gfm2final) <- 'gfm'
    return(gfm2final)
  }

}

singleIC <- function(X, group, type, q_set=1:10, dropout=0, dc_eps=1e-4, maxIter=10,output=FALSE, fast_version=TRUE){
  n <- nrow(X)
  q_num <- length(q_set)
  allhH <- list(); allhB <- list()
  for(r in 1:q_num){

    gfm1 <- gfm_eval_intercept_init(X, group, type, r,
                                      dropout, dc_eps, maxIter, output)
    if(!fast_version){
        # hH <- gfm1$hH; hB <- gfm1$hB; hmu <- gfm1$hmu
      gfm1 <-gfm_eval_intercept_osfinal(X, gfm1$hH, gfm1$hB,gfm1$hmu, group, type)
    }
    allhH[[r]] <- cbind(1, gfm1$hH)
    allhB[[r]] <- cbind(gfm1$hmu, gfm1$hB)
  }
  Vr <- matrix(0, 2, q_num)
  for(r in 1:q_num){
    Vr[,r] <- ICriteria(X, allhB[[r]], allhH[[r]], r, group, type)

  }
  IC <- apply(Vr, 2, sum)
  q <- q_set[which.min(IC)]
  return(q)
}

ICriteria <- function(X, hB, hH, r, group, type, criteria='IC'){
  # hB <- allhB[[r]]; hH <- allhH[[r]]
  n <- nrow(X); p <- ncol(X)
  omega <- 1/p
  ind_set <- unique(group)
  ng <- length(ind_set)
  gcell <- list()
  for(j in 1:ng){
    gcell[[j]] <- which(group ==j)
  }

  c1 <- objfunc(hH, hB, X, omega, gcell, type)
  Vr <- switch(criteria,
               IC=c(log(c1+1e-7), r/min(sqrt(n), sqrt(p))^2*log(min(sqrt(n), sqrt(p))^2)),
               PC=c(c1, r/min(sqrt(n), sqrt(p))^2*log(min(sqrt(n), sqrt(p))^2))
               # r * (n+p)/(n*p)*log(n*p/(n+p))
  )
  return(Vr)
}

# para
paraIC <- function(r, XX, group, type,
                   dropout=0, eps2=1e-4, maxIter=10, output=FALSE, fast_version=TRUE){

  gfm1 <- gfm_eval_intercept_init(XX, group, type, r,
                                    dropout, eps2, maxIter, output)
  if(!fast_version){
      # hH <- gfm1$hH; hB <- gfm1$hB; hmu <- gfm1$hmu
      gfm1 <-gfm_eval_intercept_osfinal(XX, gfm1$hH, gfm1$hB,gfm1$hmu, group, type)
  }
  hHm <- cbind(1, gfm1$hH)
  hBm <- cbind(gfm1$hmu, gfm1$hB)

  ICr <- sum(ICriteria(XX, hBm, hHm, r, group, type))
  return(ICr)
}

