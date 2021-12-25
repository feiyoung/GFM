localonestepH2 <- function(X, Bm, Hm, gcell, type){
  # Bm <- hBm_init ; Hm <- hHm_init
  n <- nrow(X); p <- ncol(X)
  if(ncol(Bm)-1 ==1){
    B <- matrix(Bm[,-1], ncol=1)
  }else{
    B <- Bm[,-1]
  }

  q <- ncol(Hm) - 1
  ng <- length(type)
  mucell <- list() # mean matrix
  for(j in 1:ng){
    mucell[[j]] <- switch(type[j],
           gaussian = Hm %*% t(Bm[gcell[[j]],]),
           poisson = exp(Hm %*% t(Bm[gcell[[j]],])),
           binomial = 1/(1+exp(-Hm %*% t(Bm[gcell[[j]],])))
           )
  }
  # Score matrix
  df2 <- matrix(0,n, q)
  for(j in 1:ng){
    df2 <- df2 + (X[,gcell[[j]]] - mucell[[j]]) %*% B[gcell[[j]],]
  }
  H2update <- matrix(0, n, q)
  # Hessian matrix or information matrix
  # d2f <- list()
  for(i in 1:n){
    Bng <- matrix(0, q, q)
    for(j in 1:ng){
      dBng <- switch (type[j],
        gaussian = t(B[gcell[[j]],]) %*% B[gcell[[j]],],
        poisson = t(B[gcell[[j]],]) %*% diag(mucell[[j]][i,]) %*% B[gcell[[j]],],
        binomial = t(B[gcell[[j]],]) %*% diag(mucell[[j]][i,]*(1-mucell[[j]][i,])) %*% B[gcell[[j]],]
      )
      Bng <- Bng + dBng
    }
    H2update[i,] <- qr.solve(Bng + 1e-6*diag(q)) %*% df2[i,]
  }
  H2 <- Hm[,-1] + H2update
  return(H2)
}

gfm_eval_intercept_osfinal <- function(X, hH_init, hB_init,hmu_init, group, type){
  # hH_init <- hH; hB_init <- hB; hmu_init <- hmu
  n <- nrow(hH_init); p <- nrow(hB_init)
  omega <- 1/p
  ind_set <- unique(group)
  ng <- length(ind_set)
  if(length(setdiff(1:ng, ind_set))>0){
    stop("ID number of types must match type!")
  }
  if(ng != length(type)){
    stop("The number of groups must match with length of type!")
  }
  if(ng == 1){
    hH_final <- hH_init;
    hB_final <- hB_init;
    hmu_final <- hmu_init;
    gcell <- list();gcell[[1]] <- 1:p
    objvalue <- objfunc(cbind(1, hH_final), cbind(hmu_final, hB_final),
                        X, omega, gcell, type)
  }else{
    gcell <- list()
    for(j in 1:ng){
      gcell[[j]] <- which(group ==j)
    }
    q <- ncol(hH_init)
    # update H
    hBm_init <- cbind(hmu_init, hB_init)
    hHm_init <- cbind(1, hH_init)
    H5 <- localonestepH2(X, hBm_init, hHm_init, gcell, type)
    # update B
    hhB <- NULL
    for(j in 1:ng){
      B1 <- localupdateB2(X, gcell[[j]], H5, type[j])
      hhB <- cbind(hhB, B1)
    }
    hmu <- hhB[1,]
    if(q == 1){
      hB <- matrix(hhB[-1,], ncol=1)
    }else{
      hB <- t(hhB[-1,])
    }

    # hBm_init <- cbind(hmu, hB)
    # hHm_init <- cbind(1, H5)
    # H5 <- localonestepH2(X, hBm_init, hHm_init, gcell, type)
    hH_final <- H5 # ortheH(H5)
    hB_final <- hB #ortheB(hB)
    hmu_final <- hmu
    objvalue <- objfunc(cbind(1, hH_final), cbind(hmu_final, hB_final),
                        X, omega, gcell, type)
  }
  return(list(hH=hH_final, hB=hB_final, hmu=hmu_final, obj=objvalue))
}
