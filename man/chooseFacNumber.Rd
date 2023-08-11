\name{chooseFacNumber}
\alias{chooseFacNumber}
\title{Choose the Number of factors for  Generalized Factor Models}
\description{
  This function is designed to chooose the number of factors for a generalized factor model.
}
\usage{
  chooseFacNumber(XList, types, q_set = 2: 10,
          select_method = c("SVR", "IC"),offset=FALSE,
          dc_eps=1e-4, maxIter=30, verbose = TRUE, parallelList=NULL)
}
\arguments{
  \item{XList}{a list consisting of matrices with the same rows n, and different columns (p1,p2, ..., p_d),observational mixed data matrix list, d is the types of variables, p_j is the dimension of varibles with the j-th   type.}
  \item{types}{a d-dimensional character vector, specify the type of variables. For example, \code{types=c('gaussian','poisson', 'binomial')}, implies the  components of \code{XList} are matrices with continuous, count and binomial values, respectively.}
  \item{q_set}{a positive integer vector, specify the candidates of factor number q, (optional) default as \code{c(2:10)} according to Bai (2013).}
  \item{select_method}{a string, specify the method to choose the  number of factors. Two methods are supported: the singular value ratio (SVR) and information criterion (IC) based methods, default as 'SVR'. Empirically, 'SVR' is much faster than 'IC', especially for high-dimensional large-scale data.}
   \item{offset}{a logical value, whether add an offset term (the total counts for each row in the count component of XList) when there are Poisson variables.}
  \item{dc_eps}{ positive real number, specify the tolerance of varing quantity of objective function in the algorithm. Optional parameter with default as \code{1e-4}.}
  \item{maxIter}{a positive integer, specify the times of iteration. Optional parameter with default as 50.}
  \item{verbose}{a logical value, specify whether ouput the information in iteration process, (optional) default as TRUE.}
  \item{parallelList}{a list with two components:

  (1) parallel: a logical value with TRUE or FALSE, indicates wheter to use prallel computating. Optional parameter with default as FALSE.

 (2)ncores: a positive integer, specify the number of cores when parallel computing is used.

 This argument  plays its role if only \code{select_method='IC'}.
}
}


\value{
  return an integer value, the estimated number of factors.
}
\references{
  \href{https://www.tandfonline.com/doi/abs/10.1080/01621459.2021.1999818?journalCode=uasa20}{Liu, W., Lin, H., Zheng, S., & Liu, J. (2021). Generalized factor model for ultra-high dimensional correlated variables with mixed types. Journal of the American Statistical Association, (just-accepted), 1-42.}
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
  hq <- chooseFacNumber(dat$XList, dat$types, verbose = FALSE, maxIter=2)
}
