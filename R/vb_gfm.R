

get_initials <- function(X, q){

  n <- nrow(X); p <- ncol(X)
  mu <- colMeans(X)
  X <- X - matrix(mu, nrow=n, ncol=p, byrow=TRUE)
  svdX  <- irlba(A =X, nv = q)
  PCs <- sqrt(n) * svdX$u
  loadings <- svdX$v %*% diag(svdX$d[1:q]) / sqrt(n)
  dX <- PCs %*% t(loadings) - X
  Lam_vec <- colSums(dX^2)/n
  return(list(hH = PCs, hB = loadings, hmu=mu,sigma2vec = Lam_vec))

}

vem.fit <- function(XList, types, q, offset=FALSE, epsELBO=1e-5, maxIter=30, verbose=TRUE){
  # epsELBO=1e-5; maxIter=30; verbose=TRUE
  #typeID <- c(1, 2, 3);
  n <- nrow(XList[[1]]);
  pvec <- sapply(XList, ncol); p <- sum(pvec)
  A <- matrix(0, n, p);
  type_map <- 1:3;
  names(type_map) <- c("gaussian", "poisson", "binomial")
  typeID <- unname(type_map[types])

  if(offset && any(typeID==2)){

    AList <- XList
    A <- NULL
    for(i in seq_along(typeID)){
      dim1 <- dim(AList[[i]])
      AList[[i]] <- matrix(0, dim1[1], dim1[2])
      if(typeID[i] == 2){
        tmpvec <- rowSums(AList[[i]])
        AList[[i]] <- matrix(log(tmpvec+(tmpvec==0)), dim1[1], dim1[2])
      }
      A <- cbind(A, AList[[i]]) # re-obtain A.
    }

    rm(AList)
  }
  A <- as(A, "sparseMatrix")
  Mu_y_int <- NULL
  for(i in seq_along(typeID)){
    if(typeID[i]!=2){
      Mu_y_int <- cbind(Mu_y_int,XList[[i]])
    }else if(typeID[i] == 2){
      Mu_y_int <- cbind(Mu_y_int, log(1+ XList[[i]]))
    }

  }
  S_y_int = matrix(1, n, p);
  Fac <- get_initials(Mu_y_int, q= q)
  B_int <- Fac$hB
  Mu_h_int <- Fac$hH
  # B_int = matrix(rnorm(p*q), p, q)*0.1;
  # Mu_h_int <-  matrix(rnorm(n*q), n, q)#matrix(0, n, q)
  mu_int <- colMeans(Mu_y_int)
  #invLambda_int = rep(1, p);
  invLambda_int <- 1/ Fac$sigma2vec
  reslist <- VB_GFMcpp(XList, typeID, A,  Mu_y_int, S_y_int, invLambda_int,
                       B_int, mu_int, Mu_h_int, epsELBO, maxIter,
                       verbose)



  return(reslist)
}


overdispersedGFM <- function(XList, types, q, offset=FALSE, epsELBO=1e-5, maxIter=30, verbose=TRUE){

  #   dc_eps=1e-4; maxIter=30; verbose = TRUE


  if((!is.null(q)) && (q<1) ) stop("q must be NULL or other positive integer!")
  n <- nrow(XList[[1]]); p <- sum(sapply(XList, ncol))
  if(p <2) stop("ncol(X) must be at least no less than 2!")
  if(n <2) stop("nrow(X) must be at least no less than 2!")

  type_map <- 1:3;
  names(type_map) <- c("gaussian", "poisson", "binomial")
  typeID <- unname(type_map[types]) ## match IDs
  if(length(setdiff(types,names(type_map)))>1){
    stop("The component of string vector types must be contained in ('gaussian', 'poisson', 'binomial')!")
  }
  if(verbose){
    message('Starting the varitional EM algorithm for overdispersed GFM model...\n')
  }
  tic <- proc.time()
  reslist <- vem.fit(XList, types, q=q, offset=offset, epsELBO=epsELBO, maxIter=maxIter, verbose=verbose)
  toc <- proc.time()
  if(verbose){
    message('Finish the varitional EM algorithm\n')
  }

  gfm2 <- list()
  gfm2$hH <- reslist$H
  gfm2$hB <- reslist$B
  gfm2$hmu <- t(reslist$mu)
  gfm2$obj <- reslist$ELBO
  gfm2$history <- list(c=reslist$ELBO_seq, maxIter=maxIter, eplasedTime=toc-tic)


  ## Add identifiable condition
  try({
    res_idents <- add_identifiability(gfm2$hH, gfm2$hB, gfm2$hmu)
    gfm2$hH <- res_idents$H
    gfm2$hB <- res_idents$B
    gfm2$hmu <- res_idents$mu
  }, silent = TRUE)

  gfm2$q <- q
  class(gfm2) <- 'gfm'
  return(gfm2)
}
OverGFMchooseFacNumber <- function(XList, types, q_max=15,offset=FALSE, epsELBO=1e-4, maxIter=30,
                            verbose = TRUE, threshold= 1e-2){


    gfm2 <- overdispersedGFM(XList, types, q=q_max, offset=offset, epsELBO=epsELBO,
                             maxIter=maxIter,verbose = verbose)


    svalues <- svd(gfm2$hB)$d
    svalues_use <- svalues[svalues>threshold]
    q_max <- length(svalues_use)
    q <- which.max(svalues[-q_max] / svalues[-1])
    if(verbose){
      message('SVR estimates the factor number q  as ', q, ' for the overdispersed GFM model!\n')
    }

  return(q)

}
