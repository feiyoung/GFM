Diag <- function(vec){
  q <- length(vec)
  if(q > 1){
    y <- diag(vec)
  }else{
    y <- matrix(vec, 1,1)
  }
  return(y)
}
cor.mat <- function (p, rho, type = "toeplitz"){
  if(p == 1) return(matrix(1,1,1))
  mat <- diag(p)
  if (type == "toeplitz") {
    for (i in 2:p) {
      for (j in 1:i) {
        mat[i, j] <- mat[j, i] <- rho^(abs(i - j))
      }
    }
  }
  if (type == "identity") {
    mat[mat == 0] <- rho
  }
  return(mat)
}
gendata <- function(seed=1, n=300, p=50, type='homonorm', q=6, rho=1){
  #type = {'homonorm', 'heternorm', 'pois', 'norm_pois', 'pois_bino', 'npb'}

  if(length(rho)==1) rho <- c(rho, 1.5)
  factor_term <- rho[2]

  set.seed(1)
  Z <- matrix(rnorm(p*q), p, q)
  B <- qr(Z)
  eigvs <- sqrt(sort(eigen(t(Z)%*%Z)$values, decreasing = T))
  B1 <- qr.Q(B) %*% Diag(sqrt(eigvs))
  B0 <- rho[1]* B1 %*% Diag(sign(B1[1,]))
  mu0 <- 0.4 * rnorm(p)
  Bm0 <- cbind(mu0, B0)
  set.seed(seed)
  H <- MASS::mvrnorm(n, mu=rep(0,q), cor.mat(q, 0.5))
  svdH <- svd(cov(H))
  H0 <- scale(H,scale=F) %*% svdH$u %*% Diag(1/sqrt(svdH$d)) %*% svdH$v

if(type == 'homonorm'){
  X <- H0 %*% t(B0) + matrix(mu0, n,p, byrow=T) + MASS::mvrnorm(n, rep(0,p), diag(p))
}else if(type == 'heternorm'){
  sigmas = 0.1 + 4* runif(p)
  X <- H0 %*% t(B0) + matrix(mu0, n,p, byrow=T) + MASS::mvrnorm(n, rep(0,p), diag(sigmas))
}else if(type == 'pois'){
  g1 <- 1:p
  B0[g1, ] <- B0[g1, ]/ max(B0[g1, ]) * factor_term
  mu <- exp(H0 %*% t(B0) + matrix(mu0, n,p, byrow=T))
  X <- matrix(rpois(n*p,lambda=mu),n,p)
}else if(type == 'norm_pois'){
  g1 <- 1:floor(p/2)
  g2 <- (floor(p/2)+1):p
  Bm0[g2, -1] <- Bm0[g2, -1]/ max(Bm0[g2, -1]) * factor_term
  mu1 <-  cbind(1,H0) %*% t(Bm0[g1,])
  mu2 <- exp(cbind(1, H0) %*% t(Bm0[g2,]) )
  X <- cbind(matrix(rnorm(prod(dim(mu1)), mu1,1), n, floor(p/2)),
        matrix(rpois(prod(dim(mu2)), mu2), n, ncol(mu2)))
}else if(type == 'pois_bino'){
  g1 <- 1:floor(p/2)
  g2 <- (floor(p/2)+1):p
  Bm0[g1, -1] <- Bm0[g1, -1]/ max(Bm0[g1, -1]) * factor_term
  mu1 <-  exp(cbind(1, H0) %*% t(Bm0[g1,]) )
  mu2 <- 1/(1+exp(-cbind(1, H0) %*% t(Bm0[g2,]) ))
  X <- cbind(matrix(rpois(prod(dim(mu1)), mu1), n, ncol(mu1)),
             matrix(rbinom(prod(dim(mu2)), 1, mu2), n, ncol(mu2)))
}
  return(list(X=X, B0=B0, H0=H0, mu0=mu0))
}
