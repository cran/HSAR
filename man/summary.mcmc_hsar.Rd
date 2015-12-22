\name{summary.mcmc_hsar}
\alias{summary.mcmc_hsar}
\alias{summary.mcmc_hsar_rho_0}
\alias{summary.mcmc_hsar_lambda_0}

\title{summary method for class mcmc_hsar or mcmc_hsar_rho_0 or mcmc_hsar_lambda_0}
\description{
  Methods used for presenting the results of estimating Hierarchical SAR models.This is extended to include the two simple cases.
}
\usage{
\method{summary}{mcmc_hsar}(object, ...)
\method{summary}{mcmc_hsar_rho_0}(object, ...)
\method{summary}{mcmc_hsar_lambda_0}(object, ...)
}
  
\arguments{
  \item{object}{\code{mcmc_hsar} or \code{mcmc_hsar_rho_0}, \code{mcmc_hsar_lambda_0} object created from \code{hsar} }
  \item{\dots}{Arguments passed through}
}

\seealso{\code{\link{hsar}}}

\examples{
data(land_prices)
result <- hsar( X, y, W, M, Z, Unum, 50, 100)
summary(result)

}
\keyword{spatial econometrics}
