\name{measurefun}
\alias{measurefun}
\title{Assess the performance of  an estimator on a matrix}
\description{
  Evaluate the smallest cononical correlation (ccor) coefficients or trace statistic between two matrices, where a larger ccor or  trace statistic is better.
}
\usage{
  measurefun(hH, H, type=c('trace_statistic','ccor'))
}
\arguments{
  \item{hH}{a \code{n}-by-\code{q} matrix, the estimated matrix.}
  \item{H}{a \code{n}-by-\code{q} matrix, the true  matrix.}
  \item{type}{a character taking value within \code{c('trace_statistic','ccor')}, default as 'trace_statistic'.}
}

\value{
  return a real number.
}
\author{
  Liu Wei
}
\note{
  nothing
}

\examples{
  dat <- gendata(n = 100, p = 200, q=2, rho=3)
  res <- Factorm(dat$XList[[1]])
  measurefun(res$hB, dat$B0)
}
