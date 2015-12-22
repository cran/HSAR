\name{hsar}
\alias{hsar}
%\alias{hsar}
\title{Hierarchical SAR model estimation}
\description{
  A Bayesian MCMC approach for estimation of the Hierarchical SAR model of the form:

\deqn{Y = \rho W Y + X \beta + Z \gamma +\delta \theta + \epsilon, \theta= \lambda M \theta + u}

where Y is an Nx1 vector of the outcome variable, X is an NxK matrix of the independent variables, W and M are the weight matrices at the lower and upper level respectively,\eqn{\rho} and \eqn{\lambda} are the autoregressive parameters indicating the strength of spatial interaction at the lower and higher level respectively, \eqn{\beta} is a kx1 vector of regression coefficients to estimate, \eqn{\delta} is a NxJ block diagonal matrix and \eqn{\theta} is a Jx1 vector of the regional effects.
Two special cases of HSAR are included as well. In the first one, only the regional effects are considered independent (higher level) so the autoregressive parameter \eqn{\lambda} is zero. In the second only the lower level effects are considered independent and so the autoregressive parameter \eqn{\rho} is zero.

}

\usage{
hsar( X, y, W=NULL, M=NULL, Z, Unum, burnin=5000, Nsim=10000)
}

\arguments{
  \item{X}{matrix of independent variables}
  \item{y}{vector of outcome variable values}
  \item{W}{spatial weight matrix}
  \item{M}{spatial weight matrix at the higher level}
  \item{Z}{matrix of the higher level variables}
  \item{Unum}{Relation between higher level regions and lower level areas}
  \item{burnin}{Number of samples before start collecting points}
  \item{Nsim}{Total number of samples in MC}
}

\details{
In order to run the full HSAR model both W and M matrices need to be defined.In case W is set to NULL only the regional effects are considered independent (higher level) so the autoregressive parameter \eqn{\lambda} is zero. In case M is set to NULL only the lower level effects are considered independent and so the autoregressive parameter \eqn{\rho} is zero. Both W and M can not be set to NULL.
}

\value{
a list with
  \item{Mbetas}{a vector with the mean values of the vector of the regression coefficients estimated}
  \item{SDbetas}{a vector with the standard deviation of the vector of the regression coefficients estimated}
  \item{Mrho}{Mean values of the strength of the spatial interaction rho}
  \item{SDrho}{Standard deviation of rho}
\item{Mlamda}{Mean values of the strength of the spatial interaction at the higher level}
\item{SDlambda}{Standard deviation of lambda}
\item{Msigma2e}{Mean value of \eqn{\sigma^{2}_{e} }}
\item{SDsigma2e}{Standard deviation of \eqn{\sigma^{2}_{e} }}
\item{Msigma2u}{Mean value of \eqn{\sigma^{2}_{u} }}
\item{SDsigma2u}{Standard deviation of \eqn{\sigma^{2}_{u} }}
\item{Mus}{Mean values of \eqn{\theta} }
\item{SDus}{Standard deviation of \eqn{\theta} }
\item{DIC}{deviance information criterion (DIC)}
\item{pd}{effective number of parameters}
\item{Log_Likelihood}{Log likelihood}
\item{R_Squared}{pseudo R squared}
\item{impact_direct}{Direct effect}
\item{impact_idirect}{Indirect effect}
\item{impact_total}{Total effect}

}

\references{ Dong, D. and Harris, R. 2014. Spatial Autoregressive Models for Geographically Hierarchical Data Structures. \emph{Geographical Analysis}, 1-19. }
\author{Guanpeng Dong, Richard Harris, Angelos Mimis \email{mimis@panteion.gr} }

\seealso{\code{\link{sar}}}

\examples{
data(land_prices)
result <- hsar( X, y, W, M, Z, Unum, 50, 100)
result$Mrho
}
\keyword{spatial econometrics}
