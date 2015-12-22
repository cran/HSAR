\name{land_prices}
\alias{land_prices}
\docType{data}
\title{Leased residential land parcels, from 2003 to 2009 in Beijing, China}
\description{
  In the data 1117 land parcels in 111 districts are included, having 16 independent variables.
}
\usage{data(land_prices)}
\format{
  A number of variables ready for illustating the HSAR model.
  \describe{
    \item{X}{The matrix of the lower level variables (matrix) }
    \item{y}{Leasing price per square meter for each residential land parcel (numeric vector) }
    \item{M}{The weight matrix at the higher level (dgCMatrix) }
    \item{W}{The weight matrix at the lower level (a dgCMatrix) }
    \item{Z}{ The matrix of the higher level variables(a dgCMatrix) }
    \item{Unum}{ Frequency of lower level areas contained in higher level regions (numeric vector) }
  }
}

\references{
Harris, R., G.P. Dong, and W.Z. Zhang. (2013).Using COntextualised Geographically Weighted Regression to model the Spatial Heterogeneity of Land Prices in Beijing, China. Transaction in GIS 17(6), 901-19.
}

\examples{
data(land_prices)
}
\keyword{hierarchical spatial data}
