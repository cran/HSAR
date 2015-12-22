\name{summary.mcmc_sar}
\alias{summary.mcmc_sar}

\title{summary method for class mcmc_sar}
\description{
  Methods used for presenting the results of estimating SAR models.
}
\usage{
\method{summary}{mcmc_sar}(object, ...)
}

\arguments{
  \item{object}{\code{mcmc_sar} object from \code{sar} }
  \item{\dots}{Arguments passed through}
}

\seealso{\code{\link{sar}}}

\examples{
data(land_prices)
result <- sar( X, y, W, 50, 100)
summary(result)
}
\keyword{spatial econometrics}