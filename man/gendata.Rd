\name{gendata}
\alias{gendata}
\title{Generate simulated data}
\description{
  Generate simulated data from high dimensional genelized nonlinear factor model.
}
\usage{
  gendata(seed=1, n=300, p=50, type='homonorm', q=6, rho=1)

}
\arguments{
  \item{seed}{a nonnegative integer, the random seed, default as 1.}
  \item{n}{a positive integer, the sample size.}
  \item{p}{an positive integer, the variable dimension.}
  \item{type}{a character, specify the variables type, including \code{ type = c('homonorm', 'heternorm', 'pois', 'norm_pois', 'pois_bino')}.}
  \item{q}{a positive integer, the number of factors.}
  \item{rho}{a positive number, controlling the magnitude of loading matrix.}
}

\value{
  return a list including two components:

    \item{X}{a \code{n}-by-\code{p} matrix, the observed data matrix.}
  \item{H0}{a \code{n}-by-\code{q} matrix, the true lantent factor matrix.}
  \item{B0}{a \code{p}-by-\code{q} matrix, the true loading matrix, the last pzero rows are vectors of zeros.}
  \item{ind_nz}{a integer vector, the index vector for which rows of \code{B0} not zeros.}
}
\references{
  Wei Liu, Huazhen Lin, Shurong Zheng & Jin Liu (2019) . Generalized factor model for ultra-high dimensional mixed data. Submitted.
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
\keyword{Factor}
\keyword{Feature}
