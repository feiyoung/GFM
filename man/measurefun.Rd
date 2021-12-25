\name{measurefun}
\alias{measurefun}
\title{Assess the performance of  an estimator on a matrix}
\description{
  Evaluate the smallest cononical correlation (ccor) coefficients or F-norm (fnorm) between two matrices, where a larger ccor is better; a smaller fnorm  is better.
}
\usage{
  measurefun(hH, H, type='ccor')
}
\arguments{
  \item{hH}{a \code{n}-by-\code{q} matrix, the estimated matrix.}
  \item{H}{a \code{n}-by-\code{q} matrix, the true  matrix.}
  \item{type}{a character taking value within \code{c('ccor', 'fnorm')}, default as 'ccor'.}
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
  res <- Factorm(dat$X)
  measurefun(res$hB, dat$B0)
}
