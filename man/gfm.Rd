\name{gfm}
\alias{gfm}
\title{Generalized Factor Model}
\description{
  This function is used to conduct the Generalized Factor Model.
}
\usage{
  gfm(X, group, type, q = NULL, parallel = TRUE, para.type =
                 "doSNOW", ncores = 10, dropout = 0, dc_eps = 1e-04,
                 maxIter = 50, q_set = 1:10, output = TRUE,
                 fast_version = FALSE)
}
\arguments{
  \item{X}{a matrix with dimension of n*p(p=(p1+p2+..+p_d)),observational mixed data matrix, d is the types of variables, p_j is the dimension of j-th type variable.}
  \item{group}{a vector with length equal to p, specify each column of X belonging to which group.}
  \item{type}{a d-dimensional character vector, specify the type of variables in each group. For example, \code{type=c'poisson', 'binomial')}, and it is referred to the help file of \link{glm.fit} function for more details.}
  \item{q}{a positive integer or empty, specify the number of factors. If q is \code{NULL}, then IC criteria is used to dertemined $q$ automatically.}
  \item{parallel}{a logical value with TRUE or FALSE, indicates wheter to use prallel computating. Optional parameter with default as TRUE.}
  \item{para.type}{a character specifying the type of parallel including 'doSNOW' and 'parallel'.}
  \item{ncores}{a positive integer, specify the number of cores used for parallel computing.}
  \item{dropout}{a proper subset of $[1, 2, ..., d]$,  specify which group to be dropped in obtaining the initial estimate of factor matrix $H$, and the aim is to ensure the convergence of algorithm leaded by weak signal variable types. Optional parameter with default as 0, no group dropping.}
  \item{dc_eps}{ positive real number, specify the tolerance of varing quantity of objective function in the algorithm. Optional parameter with default as \code{1e-4}.}
  \item{maxIter}{a positive integer, specify the times of iteration. Optional parameter with default as 50.}
  \item{q_set}{a positive integer vector, specify the candidates of factor number q, (optional) default as \code{c(1:10)} according to Bai,2013.}
  \item{output}{a logical value with TRUE or FALSE, specify whether ouput the mediate information in iteration process, (optional) default as FALSE.}
  \item{fast_version}{logical value with TRUE or FALSE, \code{fast_version = TRUE}: use the fast algorithm which omit the one-step updating, but it cannot ensure the estimation efficieny; \code{fast_version = FALSE}: use the original algorithm; (optional) default as \code{FALSE}; }
}
\details{
  This function also has the MATLAB version at \url{https://github.com/feiyoung/MGFM/blob/master/gfm.m}, which runs faster in MATLAB environment.
}
\value{
  return a list with class name 'gfm' and including following components,
  \item{hH}{a n*q matrix, the estimated factor matrix.}
  \item{hB}{a p*q matrix, the estimated loading matrix.}
  \item{hmu}{a p-dimensional vector, the estimated intercept terms.}
  \item{obj}{a real number, the value of objective function when the convergence achieves.}
  \item{q}{an integer, the used or estimated factor number.}
  \item{history}{a list including the following 7 components: (1)dB: the varied quantity of B in each iteration; (2)dH: the varied quantity of H in each iteration; (3)dc: the varied quantity of the objective function in each iteration; (4)c: the objective value in each iteration; (5) realIter: the real iterations to converge; (6)maxIter: the tolerance of maximum iterations; (7)elapsedTime: the elapsed time.}
}
\references{
  Liu, W., Lin, H., Zheng, S., & Liu, J. (2021). Generalized factor model for ultra-high dimensional correlated variables with mixed types. Journal of the American Statistical Association, (just-accepted), 1-42.

  Bai, J. and Liao, Y. (2013). Statistical inferences using large esti-
mated covariances for panel data and factor models.
}
\author{
  Liu Wei
}
\note{
  nothing
}
\seealso{
  nothing
}
\examples{

## mix of normal and Poisson

dat <- gendata(seed=1, n=60, p=60, type='norm_pois', q=2, rho=2)
group <- c(rep(1,ncol(dat$X)/2), rep(2,ncol(dat$X)/2))
type <- c('gaussian','poisson')
## we set maxIter=2 for example.
gfm2 <- gfm(dat$X, group, type, dropout = 2, q=2, output = FALSE, maxIter=2, parallel =FALSE)
measurefun(gfm2$hH, dat$H0, type='ccor')
measurefun(gfm2$hB, dat$B0, type='ccor')
}

\keyword{GFM}
