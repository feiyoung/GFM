\name{overdispersedGFM}
\alias{overdispersedGFM}
\title{Overdispersed Generalized Factor Model}
\description{
  This function is to implement the overdispersed generalized factor model.
}
\usage{
  overdispersedGFM(XList, types, q, offset=FALSE, epsELBO=1e-5,
                  maxIter=30, verbose=TRUE)
}
\arguments{
  \item{XList}{a list consisting of matrices with the same rows n, and different columns (p1,p2, ..., p_d),observational mixed data matrix list, d is the types of variables, p_j is the dimension of varibles with the j-th   type.}
  \item{types}{a d-dimensional character vector, specify the type of variables. For example, \code{types=c('gaussian','poisson', 'binomial')}, implies the  components of \code{XList} are matrices with continuous, count and binomial values, respectively.}
  \item{q}{a positive integer or empty, specify the number of factors.}
  \item{offset}{a logical value, whether add an offset term (the total counts for each row in the count component of XList) when there are Poisson variables.}
  \item{epsELBO}{a positive real, specify the relative tolerance of ELBO function in the algorithm. Optional parameter with default as \code{1e-5}.}
  \item{maxIter}{a positive integer, specify the times of iteration. Optional parameter with default as 30.}
  \item{verbose}{a logical value with TRUE or FALSE, specify whether ouput the information in iteration process, (optional) default as TRUE.}
}
\details{
  Overdispersion is prevalent in practical applications, particularly in fields like biomedical and genomics studies. To address this practical demand,  we propose an overdispersed generalized factor model (OverGFM) for performing high-dimensional nonlinear factor analysis on overdispersed mixed-type data.
}
\value{
  return a list with class name 'overdispersedGFM' and including following components,
  \item{hH}{a n*q matrix, the estimated factor matrix.}
  \item{hB}{a p*q matrix, the estimated loading matrix.}
  \item{hmu}{a p-dimensional vector, the estimated intercept terms.}
  \item{obj}{a real number, the value of objective function when the convergence achieves.}
  \item{q}{an integer, the used or estimated factor number.}
  \item{history}{a list including the following 7 components: (1)dB: the varied quantity of B in each iteration; (2)dH: the varied quantity of H in each iteration; (3)dc: the varied quantity of the objective function in each iteration; (4)c: the objective value in each iteration; (5) realIter: the real iterations to converge; (6)maxIter: the tolerance of maximum iterations; (7)elapsedTime: the elapsed time.}
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
## we set maxIter=2 for example.
gfm2 <- overdispersedGFM(dat$XList,  dat$types,  q=2, verbose = FALSE, maxIter=2)
measurefun(gfm2$hH, dat$H0, type='ccor')
measurefun(gfm2$hB, dat$B0, type='ccor')
}

