\name{OverGFMchooseFacNumber}
\alias{OverGFMchooseFacNumber}
\title{Choose the Number of factors for  Overdispersed Generalized Factor Models}
\description{
  This function is designed to chooose the number of factors for the overdispersed generalized factor model by using the singular value ratio (SVR) based method.
}
\usage{
  OverGFMchooseFacNumber(XList, types, q_max=15,offset=FALSE, epsELBO=1e-4, maxIter=30,
                            verbose = TRUE, threshold= 1e-2)
}
\arguments{
  \item{XList}{a list consisting of matrices with the same rows n, and different columns (p1,p2, ..., p_d),observational mixed data matrix list, d is the types of variables, p_j is the dimension of varibles with the j-th   type.}
  \item{types}{a d-dimensional character vector, specify the type of variables. For example, \code{types=c('gaussian','poisson', 'binomial')}, implies the  components of \code{XList} are matrices with continuous, count and binomial values, respectively.}
  \item{q_max}{a positive integer, specify the upper bound of the number of factors, defualt as 15.}
   \item{offset}{a logical value, whether add an offset term (the total counts for each row in the count component of XList) when there are Poisson variables.}
  \item{epsELBO}{a positive real, specify the relative tolerance of ELBO function in the algorithm. Optional parameter with default as \code{1e-5}.}
  \item{maxIter}{a positive integer, specify the times of iteration. Optional parameter with default as 30.}
  \item{verbose}{a logical value with TRUE or FALSE, specify whether ouput the information in iteration process, (optional) default as TRUE.}
  \item{threshold}{a postive real, the threshold that is used to filter the small singular values in the SVR criterion.}
}


\value{
  return an integer value, the estimated number of factors.
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
  hq <- OverGFMchooseFacNumber(dat$XList, dat$types, verbose = FALSE, maxIter=2)
}
