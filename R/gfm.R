

# library(pkgdown)
# build_site()
# build_article('GFM.Simu')
# build_article("GFM.Brain")
gfm <- function(XList, types, q=10, offset=FALSE, dc_eps=1e-4, maxIter=30,
                verbose = TRUE, algorithm=c("VEM", "Alternate_Maximization")){

  #  algorithm=c("VEM"); offset=FALSE; dc_eps=1e-4; maxIter=30; verbose = TRUE


  if((!is.null(q)) && (q<1) ) stop("q must be NULL or other positive integer!")
  n <- nrow(XList[[1]]); p <- sum(sapply(XList, ncol))
  if(p <20) stop("ncol(X) must be at least no less than 20!")
  if(n <20) stop("nrow(X) must be at least no less than 20!")

  algorithm <- match.arg(algorithm)

  type_map <- 1:3;
  names(type_map) <- c("gaussian", "poisson", "binomial")
  typeID <- unname(type_map[types]) ## match IDs
  if(length(setdiff(types,names(type_map)))>1){
    stop("The component of string vector types must be contained in ('gaussian', 'poisson', 'binomial')!")
  }

  if(algorithm== "VEM"){

    message('Starting the varitional EM algorithm...\n')
    tic <- proc.time()
    reslist <- gfm_vb.fit(XList, types, q=q, offset=offset, epsELBO=dc_eps, maxIter=maxIter, verbose=verbose)
    toc <- proc.time()
    message('Finish the varitional EM algorithm...\n')
    gfm2 <- list()
    gfm2$hH <- reslist$H
    gfm2$hB <- reslist$B
    gfm2$hmu <- t(reslist$mu)
    gfm2$obj <- reslist$ELBO
    gfm2$history <- list(c=reslist$ELBO_seq, maxIter=maxIter, eplasedTime=toc-tic)

    try({
      X <- NULL
      group <- NULL
      for(i in seq_along(typeID)){
        X <- cbind(X, XList[[i]])
        group <- c(group, rep(i, ncol(XList[[i]])))
      }
      rm(XList)
      gfm2final <- gfm_eval_intercept_osfinal(X, gfm2$hH, gfm2$hB,gfm2$hmu, group, types)
      if(any(is.nan(gfm2final$hH)) || any(is.nan(gfm2final$hB))) error("there is NaNs produced!")
      gfm2final$history <- gfm2$history
      gfm2 <- gfm2final
    }, silent = TRUE)

  }else if(algorithm== 'Alternate_Maximization'){

    X <- NULL
    group <- NULL
    for(i in seq_along(typeID)){
      X <- cbind(X, XList[[i]])
      group <- c(group, rep(i, ncol(XList[[i]])))
    }
    rm(XList)
    omega <- 1/p
    message('Starting the alternate maximization algorithm...\n')
    dropout <- 0
    if(any(typeID==2) && length(typeID)>1) dropout <- which(typeID==2)
    gfm2 <- gfm_eval_intercept_init(X, group, types, q,
                                    dropout, dc_eps, maxIter,
                                    verbose)

    message('Finish the alternate maximization algorithm...\n')
    try({
      gfm2final <- gfm_eval_intercept_osfinal(X, gfm2$hH, gfm2$hB,gfm2$hmu, group, types)
      if(any(is.nan(gfm2final$hH)) || any(is.nan(gfm2final$hB))) error("there is NaNs produced!")
      gfm2final$history <- gfm2$history
      gfm2 <- gfm2final
    }, silent = TRUE)
  }
  ## Add identifiable condition
  res_idents <- add_identifiability(gfm2$hH, gfm2$hB, gfm2$hmu)
  gfm2$hH <- res_idents$H
  gfm2$hB <- res_idents$B
  gfm2$hmu <- res_idents$mu
  gfm2$q <- q
  class(gfm2) <- 'gfm'
  return(gfm2)
}

chooseFacNumber <- function(XList, types, q_set = 2: 10, select_method = c("ratio_test", "IC"),offset=FALSE, dc_eps=1e-4, maxIter=30,
                            verbose = TRUE, parallelList=NULL){

  select_method <- match.arg(select_method)
  if(select_method == 'IC'){
    if(is.null(parallelList$parallel)){
      parallelList$parallel <- FALSE
    }
    if(is.null(parallelList$ncores)){
      parallelList$ncores <- 5
    }
    parallel<- parallelList$parallel
    ncores <- parallelList$ncores


    type_map <- 1:3;
    names(type_map) <- c("gaussian", "poisson", "binomial")
    typeID <- unname(type_map[types]) ## match IDs
    X <- NULL
    group <- NULL
    for(i in seq_along(typeID)){
      X <- cbind(X, XList[[i]])
      group <- c(group, rep(i, ncol(XList[[i]])))
    }
    rm(XList)

    dropout <- 0
    if(any(typeID==2) && length(typeID)>1) dropout <- which(typeID==2)
    if(parallel){
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
                        res <- paraIC(r, X, group=group,type= types,
                                      dropout=dropout, eps2=1e-4, maxIter=10, output=verbose)
                        res
                      }
      close(pb)
      parallel::stopCluster(cl)

    }else{
      n <- nrow(X)
      q_num <- length(q_set)
      allhH <- list(); allhB <- list()
      for(r in 1:q_num){

        gfm1 <- gfm_eval_intercept_init(X, group, type=types, r,
                                        dropout, eps2=1e-4, maxIter=10, verbose)

        gfm1 <-gfm_eval_intercept_osfinal(X, gfm1$hH, gfm1$hB,gfm1$hmu, group, type=types)

        allhH[[r]] <- cbind(1, gfm1$hH)
        allhB[[r]] <- cbind(gfm1$hmu, gfm1$hB)
      }
      Vr <- matrix(0, 2, q_num)
      for(r in 1:q_num){
        Vr[,r] <- ICriteria(X, allhB[[r]], allhH[[r]], r, group, types)

      }
      resq <- apply(Vr, 2, sum)
    }
    q <- q_set[which.min(resq)]

    message('IC criterion estimates the factor number q  as ', q, '. \n')
  }else if(select_method == 'ratio_test'){
     q_max <- max(q_set)
     gfm2 <- gfm(XList, types, q=q_max, algorithm = 'VEM', offset=offset, dc_eps=dc_eps, maxIter=maxIter,verbose = verbose)


     svalues <- svd(gfm2$hB)$d
     svalues_use <- svalues[svalues>1e-2]
     q_max <- length(svalues_use)
     q <- which.max(svalues[-q_max] / svalues[-1])

     message('Eigevalue ratio test estimates the factor number q  as ', q, '. \n')
  }



    return(q)

}
singleIC <- function(X, group, type, q_set=1:10, dropout=0, dc_eps=1e-4, maxIter=10,output=FALSE){
  n <- nrow(X)
  q_num <- length(q_set)
  allhH <- list(); allhB <- list()
  for(r in 1:q_num){

    gfm1 <- gfm_eval_intercept_init(X, group, type, r,
                                      dropout, dc_eps, maxIter, output)


    gfm1 <-gfm_eval_intercept_osfinal(X, gfm1$hH, gfm1$hB,gfm1$hmu, group, type)

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

add_identifiability <- function(H, B, mu){
  # Load the irlba library
  #library(irlba)

  # Perform SVD decomposition with rank k = 10

  mu <- mu + B %*% colMeans(H)
  q <- ncol(H); n <- nrow(H)
  svdHB <- irlba((H- matrix(colMeans(H), n, q, byrow = TRUE)) %*% t(B), nv= q)
  signB1 <- sign(svdHB$v[1,])
  H <- sqrt(n) * svdHB$u %*% Diag(signB1)

  B <- svdHB$v %*% Diag(svdHB$d[1:q]*signB1) / sqrt(n)

  return(list(H=H, B=B, mu=mu))
}

