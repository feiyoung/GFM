\name{gendata}
\alias{gendata}
\title{Generate simulated data}
\description{
  Generate simulated data from high dimensional genelized nonlinear factor model.
}
\usage{
  gendatafunction (seed = 1, n = 300, p = 50,
              type =  c('homonorm', 'heternorm', 'pois', 'bino', 'norm_pois', 'pois_bino', 'npb'),
                     q = 6, rho = 1, n_bin=1)

}
\arguments{
  \item{seed}{a nonnegative integer, the random seed, default as 1.}
  \item{n}{a positive integer, the sample size.}
  \item{p}{an positive integer, the variable dimension.}
  \item{type}{a character, specify the variables types for generated data, default as 'homonorm', representing the homogeneous gaussian variables.}
  \item{q}{a positive integer, the number of factors.}
  \item{rho}{a positive number, controlling the magnitude of loading matrix.}
  \item{n_bin}{a positive integer, specify the number of trails for the binomial variables when \code{type} is set to one of 'bino', 'pois_bino' and 'npb'.}
}
\details{
  This function provides a variaty of mix of different variable types, in which 'homonorm' represents the generated data with only homogenous normal variables; 'heternorm' represents the data with only heterogenous normal variables; 'pois' means the data with only poisson variables; 'bino' means the data with only binomial variables; 'norm_pois' means the mix of normal and poisson variables; 'pois_bino' represents the mix of poisson and binomial variables; and 'npb' means the most complex mix of normal, poisson and binomial variables.
}

\value{
  return a list including two components:
  \item{X}{a \code{n}-by-\code{p} matrix, the observed data matrix.}
  \item{XList}{a list consisting of the above observed data matrices with the same rows n (observations), and different columns (p1,p2, ..., p_d) and p columns in total, where d is the types of variables, pj is the dimension of varibles with the j-th   type.}
  \item{H0}{a \code{n}-by-\code{q} matrix, the true lantent factor matrix.}
  \item{B0}{a \code{p}-by-\code{q} matrix, the true loading matrix, the last pzero rows are vectors of zeros.}
  \item{mu0}{a p-dimensional vector, the true intercept terms.}
}



\author{
  Wei Liu
}
\note{
  nothing
}
\seealso{
  \code{\link{Factorm}; \link{gfm}}.
}
\examples{
  dat <- gendata(n=300, p = 500)
  str(dat)
}
