

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

gfm_vb.fit <- function(XList, types, q, offset=FALSE, epsELBO=1e-5, maxIter=30, verbose=TRUE){
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
  library(Matrix)
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
  #S_h_int  <- diag(rep(1, q))
  S_h_int <- var(Mu_h_int)
  Sigma_h_int <- diag(rep(1, q))
  #invLambda_int = rep(1, p);
  invLambda_int <- 1/ Fac$sigma2vec
  reslist <- VB_GFMcpp(XList, typeID, A,  Mu_y_int, S_y_int, invLambda_int,
                       B_int, mu_int, Mu_h_int, S_h_int, Sigma_h_int, epsELBO, maxIter,
                       verbose)



  return(reslist)
}
