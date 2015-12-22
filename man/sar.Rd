\name{sar}
\alias{sar}
%\alias{sar}
\title{SAR model estimation}
\description{
  A Bayesian MCMC approach for estimation of the SAR model (spatial lag model) of the form:

\deqn{Y = \rho W Y + X \beta + \epsilon}

where Y is an Nx1 vector of the outcome variable, X is an NxK matrix of the independent variables and W is the weight matrix.

}

\usage{
sar( X, y, W, burnin=5000, Nsim=10000 )
}

\arguments{
  \item{X}{matrix of independent variables}
  \item{y}{vector of outcome variable values}
  \item{W}{spatial weight matrix}
  \item{burnin}{Number of samples before start collecting points}
  \item{Nsim}{Total number of samples in MC}
}

\value{
a list with
  \item{Mbetas}{a vector with the mean values of the vector of the regression coefficients estimated}
  \item{SDbetas}{a vector with the standard deviation of the vector of the regression coefficients estimated}
  \item{Mrho}{Mean values of the strength of the spatial interaction rho}
  \item{SDrho}{Standard deviation of rho}
\item{Msigma2e}{Mean value of \eqn{\sigma^{2}_{e} }}
\item{SDsigma2e}{Standard deviation of \eqn{\sigma^{2}_{e} }}
\item{DIC}{deviance information criterion (DIC)}
\item{pd}{effective number of parameters}
\item{Log_Likelihood}{Log likelihood}
\item{R_Squared}{pseudo R squared}
\item{impact_direct}{Direct effect}
\item{impact_idirect}{Indirect effect}
\item{impact_total}{Total effect}

}

\references{ Anselin, L. 1988 \emph{Spatial econometrics: methods and models.}(Dordrecht: Kluwer); 
Dong, D. and Harris, R. 2014. Spatial Autoregressive Models for Geographically Hierarchical Data Structures. \emph{Geographical Analysis}, 1-19. }
\author{Guanpeng Dong, Richard Harris, Angelos Mimis \email{mimis@panteion.gr} }

\seealso{\code{\link{hsar}}}

\examples{
data(land_prices)
result <- sar( X, y, W, 50, 100)
result$Mrho
}
\keyword{spatial econometrics}

