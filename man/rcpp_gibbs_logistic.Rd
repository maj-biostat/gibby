\name{rcpp_gibbs_logistic}
\alias{rcpp_gibbs_logistic}
\title{Set of functions in example RcppArmadillo package}
\description{
  These four functions are created when
  \code{RcppArmadillo.package.skeleton()} is invoked to create a
  skeleton packages.
}
\usage{
rcpp_gibbs_logistic()
}
\arguments{
  \item{x}{a numeric vector}
}
\value{
  \code{rcpp_gibbs_logistic()} does not return a value, but displays a
  message to the console.

}
\details{
  These are example functions which should be largely
  self-explanatory. Their main benefit is to demonstrate how to write a
  function using the Armadillo C++ classes, and to have to such a
  function accessible from R.
}
\references{
  See the documentation for Armadillo, and RcppArmadillo, for more details.
}
\examples{
  x <- sqrt(1:4)
  rcpparma_innerproduct(x)
  rcpparma_outerproduct(x)
}
\author{MAJ}

