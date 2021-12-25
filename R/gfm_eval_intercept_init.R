signrevise <- function(A1, A2){
  nzid1 <- which(rowSums(A1^2)> 1e-5)[1]
  q <- ncol(A1)
  A <- sapply(1:q, function(k){
    if(sign(A1[nzid1,k]) != sign(A2[nzid1,k]))
      return(-A1[,k])
    return(A1[,k])
  })
  return(A)
}

Factorm <- function (X, q = NULL)
{
  X <- scale(X, scale=FALSE)
  n <- nrow(X)
  p <- ncol(X)
  if (p > n) {
    svdX <- eigen(X %*% t(X))
    evalues <- svdX$values
    eigrt <- evalues[1:(21 - 1)]/evalues[2:21]
    if (is.null(q)) {
      q <- which.max(eigrt)
    }
    hatF <- as.matrix(svdX$vector[, 1:q] * sqrt(n))
    B2 <- n^(-1) * t(X) %*% hatF
    sB <- sign(B2[1, ])
    hB <- B2 * matrix(sB, nrow = p, ncol = q, byrow = TRUE)
    hH <- sapply(1:q, function(k) hatF[, k] * sign(B2[1,
                                                      ])[k])
  }
  else {
    svdX <- eigen(t(X) %*% X)
    evalues <- svdX$values
    eigrt <- evalues[1:(21 - 1)]/evalues[2:21]
    if (is.null(q)) {
      q <- which.max(eigrt)
    }
    hB1 <- as.matrix(svdX$vector[, 1:q])
    hH1 <- n^(-1) * X %*% hB1
    svdH <- svd(hH1)
    hH2 <- signrevise(svdH$u * sqrt(n), hH1)
    if (q == 1) {
      hB1 <- hB1 %*% svdH$d[1:q] * sqrt(n)
    }
    else {
      hB1 <- hB1 %*% diag(svdH$d[1:q]) * sqrt(n)
    }
    sB <- sign(hB1[1, ])
    hB <- hB1 * matrix(sB, nrow = p, ncol = q, byrow = TRUE)
    hH <- sapply(1:q, function(j) hH2[, j] * sB[j])
  }
  sigma2vec <- colMeans((X - hH %*% t(hB))^2)
  res <- list()
  res$hH <- hH
  res$hB <- hB
  res$q <- q
  res$sigma2vec <- sigma2vec
  res$propvar <- sum(evalues[1:q])/sum(evalues)
  res$egvalues <- evalues
  attr(res, "class") <- "fac"
  return(res)
}

ortheB <- function(Z){
  B <- qr(Z)
  eigvs <- sqrt(sort(eigen(t(Z)%*%Z)$values, decreasing = TRUE))
  B1 <- qr.Q(B) %*% Diag(eigvs)
  B0 <- B1 %*% Diag(sign(B1[1,]))
  return(B0)
}
ortheH <- function(H){
  H1 <- qr.Q(qr(H)) * sqrt(nrow(H))
  hH <- H1 %*% Diag(sign(H[1,]) * sign(H1[1,]))
  return(hH)
}

family2func <- function(type){
  switch(type,
         gaussian= gaussian(link = "identity"),
         binomial = binomial(link = "logit"),
         poisson = poisson(link = "log"))
}
paraglmfit <- function(j, Xx, g1, XX, fun1){
  glm.fit(x=Xx, y=XX[, g1[j]], family = fun1, intercept = FALSE)$coefficients
}
localupdateB2 <- function(X, g1, hH, type1, parallel=FALSE){
  # g1 <- gcell[[1]]; type1 <- type[1]
  n <- nrow(X);q <- ncol(hH)
  p1 <- length(g1); B1 <- matrix(0, q+1, p1)
  jg <- 1:p1
  fun1 <- family2func(type1)
  if(parallel){

    # varlist <- c('gfm_eval_intercept_init', 'ICriteria', "Factorm",
    #              "localupdateB2","family2func","ortheB","Diag",
    #              "localupdateH2","ortheH", "objfunc","signrevise")
    varlist <- NULL
    B1 <- single_parallel(paraglmfit,iterable = jg, varlist= varlist,Xx=cbind(1,hH), g1=g1, XX=X, fun1=fun1)
  }else{
    for(j in jg){

      B1[,j] <- glm.fit(x=cbind(1,hH), y=X[, g1[j]], family = fun1, intercept = FALSE)$coefficients
    }
  }

  return(B1)
}

paraglmj1fit <- function(i,j, hBm, gcell, w, XX, funj){
  glm.fit(hBm[gcell[[j]], -1], XX[i, gcell[[j]]], weights = w,intercept = FALSE,
          family=funj,offset = hBm[gcell[[j]], 1])$coefficients

}
paraglmj2fit <- function(i,j, hBm, gcell, XX, funj){
  glm.fit(hBm[gcell[[j]], -1], XX[i, gcell[[j]]], intercept = FALSE,
          family=funj,offset = hBm[gcell[[j]], 1])$coefficients

}
localupdateH2 <- function(X, gcell, hBm, type, dropout, parallel=FALSE){
  n <- nrow(X); q <- ncol(hBm)-1; ng <- length(type)
  if(dropout !=0 && length(setdiff(dropout, 1:ng))>0){
    stop('dropout setting is wrong!')
  }
  idres <- setdiff(1:ng, dropout)
  Harray <- array(0, dim=c(n, q, length(idres)))
  w <- rep(1,n)
  for(jj in 1: length(idres)){
    j <- idres[jj]
    funj <- family2func(type[j])
    H2 <- matrix(rnorm(q*n), q, n)
    if(type[j] == 'gaussian'){
      w <- 1/ apply(X[,gcell[[j]]], 2, var)
      if(parallel){

        # varlist <- c('gfm_eval_intercept_init', 'ICriteria', "Factorm",
        #              "localupdateB2","family2func","ortheB","Diag",
        #              "localupdateH2","ortheH", "objfunc","signrevise")
        varlist <- NULL
        H2 <- single_parallel(paraglmj1fit,iterable = 1:n, varlist= varlist,j=j, hBm=hBm, gcell=gcell,w=w, XX=X, funj=funj)

      }else{
        for(i in 1:n){ #i <- 1
          H2[,i] <- glm.fit(hBm[gcell[[j]], -1], X[i, gcell[[j]]], weights = w,intercept = FALSE,
                            family=funj,offset = hBm[gcell[[j]], 1])$coefficients
        }
      }

    }else{
      if(parallel){

        # varlist <- c('gfm_eval_intercept_init', 'ICriteria', "Factorm",
        #              "localupdateB2","family2func","ortheB","Diag",
        #              "localupdateH2","ortheH", "objfunc","signrevise")
        varlist <- NULL
        H2 <- single_parallel(paraglmj2fit,iterable = 1:n, varlist= varlist,j=j, hBm=hBm, gcell=gcell,XX=X, funj=funj)
      }else{
        for(i in 1:n){
          try(
            H2[,i] <- glm.fit(hBm[gcell[[j]], -1], X[i, gcell[[j]]],intercept = FALSE,
                              family=funj,offset = hBm[gcell[[j]], 1])$coefficients
            , silent=TRUE
          )
        }
      }

    }
    Harray[,,jj] <- t(H2)
  }
  hH <- apply(Harray, c(1,2), mean)
  return(hH)
}

objfunc <- function(Hm, Bm, X, omega, gcell, type){
  n <- nrow(X); p <- ncol(X)
  eps1 <- 1e-20
  BHm <- Hm %*% t(Bm)
  ng <- length(type)
  Q <- matrix(0, n, p)
  for(j in 1:ng){
    if(type[j]== 'gaussian'){
      Q[,gcell[[j]]] <- (X[,gcell[[j]]] - BHm[,gcell[[j]]])^2
    }else if(type[j] == 'poisson'){
      me <- exp(BHm[,gcell[[j]]])
      Q[,gcell[[j]]] <- -log(dpois(X[, gcell[[j]]], me)+eps1)
    }else if(type[j] == 'binomial'){
      me3 <- 1 / (1 + exp(-BHm[,gcell[[j]]]))
      Q[,gcell[[j]]] <- -X[, gcell[[j]]] * log(me3+eps1) +
        (1-X[, gcell[[j]]]) * log(1-me3 + eps1)
    }
  }
  obj <- 1/n*omega*sum(Q)
  return(obj)
}
gfm_eval_intercept_init <- function(X, group, type, q,
                                    dropout, eps2, maxIter=10,
                                    output=0, parallel=FALSE){
  ind_set <- unique(group)
  ng <- length(ind_set)
  if(length(setdiff(1:ng, ind_set))>0){
    stop("ID number of types must match type!")
  }
  if(ng != length(type)){
    stop("The number of groups must match with length of type!")
  }

  gcell <- list()
  for(j in 1:ng){
    gcell[[j]] <- which(group==j)
  }
  n <- nrow(X); p <- ncol(X)
  if(length(group) != p){
    stop("The length of group must match with column of X!")
  }
  omega <- 1/p
  #initialize
  hH <- Factorm(scale(X,scale=FALSE), q)$hH;hB <- 0
  eps1 <- 1e-4
  dBm <- Inf; dH <- Inf; dc =Inf; dOmega <- max(dBm, dH)
  tmpBm <- matrix(0, p, q+1); tmpH <- hH; tmpc <- 1e7
  k <- 1; history <- list()
  tic <- proc.time()
  while(k <= maxIter && dOmega > eps1 && dc >eps2){
    hhB <- NULL
    for(j in 1:ng){
      B1 <- localupdateB2(X, gcell[[j]], hH, type[j], parallel)
      hhB <- cbind(hhB, B1)
    }
    hmu <- hhB[1,]
    if(q == 1){
      hB <- matrix(hhB[-1,], ncol=1)
    }else{
      hB <- t(hhB[-1,])
    }
    # ensure indentifiability.
    hB <- ortheB(hB)
    hBm <- cbind(hmu, hB)
    dB <- norm(hBm-tmpBm, "F") / norm(hBm, 'F')
    tmpBm <- hBm
    if(output){
      message('---------- B updation is finished!---------\n')
    }
    H4 <- localupdateH2(X, gcell, hBm, type, dropout, parallel)
    hH <- ortheH(H4)# %*% diag(sign(H4[1,]))
    dH <- norm(hH- tmpH, 'F')/norm(hH, 'F')
    tmpH <- hH
    if(output){
      message('---------- H updation is finished!---------\n')
    }
    hHm <- cbind(1, hH)
    dOmega <- max(dB, dH)
    c <- objfunc(hHm, hBm, X, omega, gcell, type)
    dc <- abs(c - tmpc)/abs(tmpc)
    tmpc <- c
    if(output){
      message('Iter=', k, ', dB=',round(dB,4), ', dH=', round(dH,4),
          ',dc=', round(dc,4), ', c=', round(c,4), '\n')
    }
    history$dB[k] <- dB; history$dH[k] <- dH; history$dc[k] <- dc;
    history$c[k] <- c
    k <- k + 1
  }
  stoc <- proc.time() - tic
  history$realIter <- k -1
  history$maxIter <- maxIter
  history$elapsedTime <- stoc
  return(list(hH=hH, hB=hB, hmu=hmu, history=history))
}

measurefun <- function(hH, H, type='ccor'){
  q <- min(ncol(H), ncol(hH))
  switch(type,
         ccor=cancor(hH, H)$cor[q],
         fnorm= norm(H-hH, 'F')^2/ prod(dim(H)))
}
